/*
This file is part of tuwcartesians.

tuwcartesians is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

tuwcartesians is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with tuwcartesians.  If not, see <http://www.gnu.org/licenses/>.
*/

/* (c) Jesper Larsson Traff, TUW, September 2014, April 2018, November 2018 */

/*
 A proposal for additional, isormorphic (now: Cartesian) collectives.
 All blocking versions implemented with combining schedules. Preparations
 for persistent, non-blocking versions. Additional, missing Cartesian
 helper functionality. Limited error handling: check for Cartesian
 communicator, return of error code (not abort)
 */

/* For schedule debugging, see stand-alone versions */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>
#include <assert.h>

#include "tuwcartesians.h"

#define TUWCARTTAG 777

typedef enum {
  regular, vector, wector
} collvariant;

typedef struct {
  int commtype;     // type information TUW_CART_GRAPH
  int blocking, trivial; // implementation choices, set by info-object
  
  int d;            // Cartesian dimension
  int t;            // size of neighborhood (number of targets)
  int *vectors;     // relative coordinates (outgoing)
  int *outrank;     // ranks of target processes
  int *inrank;      // ranks of source processes
  int *weight;

  // additional schedule and datatypes for message-combining algorithms
  // (should be: be schedule per operation, alltoall, alltoallv, alltoallw, ...)
  int precomputed;
  int phases; // number of phases in schedule
  int *phase; // number of neighbors per phase
  void *sendbuf;
  void *recvbuf;
  void *tempbuf;
  MPI_Datatype *sendtype;
  MPI_Datatype *recvtype;
  int *sendrank;
  int *recvrank;
} cartattr;

// forward declarations
static int findrank(const int rank, const int r, int ranks[], const int max, int found[]);
static int nonzeros(const int d, const int vector[]);
static int maxcoordinate(const int s, const int d, const int vectors[], const int j);
static void bucketsort(const int s, const int d, const int vectors[], const int mc, const int j, int order[]);
static void schedule_free(cartattr *neighborhood);
static void schedule_cart_alltoall(collvariant variant, cartattr *neighborhood, void *sendbuf, const int *sendcount,
    const void *senddisp, const MPI_Datatype *sendtype, void *recvbuf, const int *recvcount, const void *recvdisp,
    const MPI_Datatype *recvtype, MPI_Comm cartcomm);

static void schedule_cart_allgather(collvariant variant, cartattr *neighborhood, void *sendbuf, int sendcount,
    MPI_Datatype sendtype, void *recvbuf, const int *recvcount, const void *recvdisp, const MPI_Datatype *recvtype,
    MPI_Comm cartcomm);

static int cartdel(MPI_Comm comm, int keyval, void *attr, void *s) {
  cartattr *neighborhood = (cartattr*) attr;

  assert(neighborhood->commtype==TUW_CART_GRAPH);

  free(neighborhood->vectors);
  free(neighborhood->outrank);
  free(neighborhood->inrank);
  if (neighborhood->weight != MPI_UNWEIGHTED)
    free(neighborhood->weight);

  if (neighborhood->precomputed)
    schedule_free(neighborhood);

  free(neighborhood);

  return MPI_SUCCESS;
}

static int cartkey() {
  // hidden key value for type attributes
  static int cartkeyval = MPI_KEYVAL_INVALID;

  if (cartkeyval == MPI_KEYVAL_INVALID) {
    MPI_Comm_create_keyval(MPI_COMM_NULL_COPY_FN, cartdel, &cartkeyval, NULL);
  }

  return cartkeyval;
}

/* Compute rank of relative coordinate */
int TUW_Cart_relative_rank(MPI_Comm cartcomm, const int relative[], int *rank) {
  int commtype;
  int callrank;
  int d;
  int i;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Cartdim_get(cartcomm, &d);
  int coordinate[d];

  MPI_Comm_rank(cartcomm, &callrank);
  MPI_Cart_coords(cartcomm, callrank, d, coordinate);

  for (i = 0; i < d; i++)
    coordinate[i] = coordinate[i] + relative[i];
  MPI_Cart_rank(cartcomm, coordinate, rank);

  return MPI_SUCCESS;
}

/* Shift function with relative target coordinate */
int TUW_Cart_relative_shift(MPI_Comm cartcomm, const int relative[], int *inrank, int *outrank) {
  int commtype;
  int d;
  int i;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  TUW_Cart_relative_rank(cartcomm, relative, outrank);

  MPI_Cartdim_get(cartcomm, &d);
  int coordinate[d];
  for (i = 0; i < d; i++)
    coordinate[i] = -relative[i];
  TUW_Cart_relative_rank(cartcomm, coordinate, inrank);

  return MPI_SUCCESS;
}

/* Compute relative coordinate of given rank relative to calling process' 
 rank */
int TUW_Cart_relative_coord(MPI_Comm cartcomm, const int rank, int relative[]) {
  int commtype;
  int d;
  int i;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Cartdim_get(cartcomm, &d);

  int dims[d];
  int periods[d];
  int coordinate[d];
  MPI_Cart_get(cartcomm, d, dims, periods, coordinate);
  MPI_Cart_coords(cartcomm, rank, d, relative);

  for (i = 0; i < d; i++) {
    relative[i] = (relative[i] - coordinate[i] + dims[i]) % dims[i];
  }

  return MPI_SUCCESS;
}

/* Missing functions to count number of and determine neighbors in 
 MPI standard Cartesian communicators */
int TUW_Cart_neighbors_count(MPI_Comm cartcomm, int *outdegree, int *indegree) {
  int commtype;
  int d;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Cartdim_get(cartcomm, &d);

  // TODO: correct treatment for non-periodic communicators
  *outdegree = 2 * d;
  *indegree = 2 * d;

  return MPI_SUCCESS;
}

/* Cartesian neighbors in canonical order according to standard */
int TUW_Cart_neighbors(MPI_Comm cartcomm, const int maxindegree, int source[], const int maxoutdegree, int target[]) {
  int commtype;
  int d;
  int i;
  int inrank, outrank;
  int tars, srcs;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Cartdim_get(cartcomm, &d);

  srcs = 0;
  tars = 0;
  for (i = 0; i < d; i++) {
    MPI_Cart_shift(cartcomm, i, 1, &inrank, &outrank);
    if (i < maxindegree)
      source[srcs++] = inrank;
    if (i < maxindegree)
      source[srcs++] = outrank;
    if (i < maxoutdegree)
      target[tars++] = outrank;
    if (i < maxoutdegree)
      target[tars++] = inrank;
  }
  assert(srcs <= maxindegree);
  assert(tars <= maxoutdegree);

  return MPI_SUCCESS;
}

/* Takes list of relative coordinates and attaches list of inranks and
 outranks to new communicator */
/* In line with other MPI functions the list of relative coordinates is flat */
/* Weights, info and reorder are ignored */

int TUW_Cart_neighborhood_create(MPI_Comm cartcomm, int d, const int dims[], const int periods[], int t,
    const int targetrelative[], const int weight[], MPI_Info info, int reorder, MPI_Comm *collcomm) {
  //int commtype;
  //int d;
  int i, j;
  cartattr *neighborhood;
  int *coordinate;
  int flag;

  //MPI_Topo_test(cartcomm,&commtype); // Cartesian error check here
  //if (commtype!=MPI_CART) return MPI_ERR_TOPOLOGY;

  //MPI_Comm_get_attr(cartcomm,cartkey(),&neighborhood,&flag);
  //if (flag) return MPI_ERR_TOPOLOGY;

  //MPI_Cartdim_get(cartcomm,&d);

  MPI_Cart_create(cartcomm, d, dims, periods, 0, collcomm);

  neighborhood = (cartattr*) malloc(sizeof(cartattr));
  neighborhood->commtype = TUW_CART_GRAPH;

  neighborhood->vectors = (int*) malloc(t * d * sizeof(int));
  neighborhood->inrank = (int*) malloc(t * sizeof(int));
  neighborhood->outrank = (int*) malloc(t * sizeof(int));
  if (weight != MPI_UNWEIGHTED) {
    neighborhood->weight = (int*) malloc(t * sizeof(int));
  } else
    neighborhood->weight = MPI_UNWEIGHTED;

  neighborhood->d = d;
  neighborhood->t = t;

  // compute the ranks from list of relative coordinates
  coordinate = (int*) targetrelative;
  for (i = 0; i < t; i++, coordinate += d) {
    for (j = 0; j < d; j++)
      neighborhood->vectors[i * d + j] = coordinate[j];
    TUW_Cart_relative_shift(*collcomm, coordinate, &neighborhood->inrank[i], &neighborhood->outrank[i]);
    if (weight != MPI_UNWEIGHTED)
      neighborhood->weight[i] = weight[i];
  }

  // no schedule, remaining neighborhood values undefined
  neighborhood->precomputed = 0;

  // decode info for algorithm/implementation selection
  if (info!=MPI_INFO_NULL) {
    char keyvalue[6];
    MPI_Info_get(info,"blocking",5,keyvalue,&flag);
    if (flag&&strcmp(keyvalue,"true")==0)
      neighborhood->blocking = 1; else neighborhood->blocking = 0;
    MPI_Info_get(info,"trivial",5,keyvalue,&flag);
    if (flag&&strcmp(keyvalue,"true")==0)
      neighborhood->trivial = 1; else neighborhood->trivial = 0;
  } else {
    // default choices: non-blocking, non-trivial message-combining
    neighborhood->blocking = 0;
    neighborhood->trivial  = 0;
  }

  //MPI_Comm_dup(cartcomm,collcomm);
  MPI_Comm_set_attr(*collcomm, cartkey(), neighborhood);

  return MPI_SUCCESS;
}

/* Return (maximal) size of neighborhood */
int TUW_Cart_neighborhood_count(MPI_Comm cartcomm, int *t) {
  int commtype;
  int flag;
  cartattr *neighborhood;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  *t = neighborhood->t;

  return MPI_SUCCESS;
}

/* Return lists of (actual) sources and destinations and weights in a form 
 that can be used immediately for MPI_Dist_graph_create_adjacent */
int TUW_Cart_neighborhood_get(MPI_Comm cartcomm, const int maxindegree, int source[], int sourceweight[],
    const int maxoutdegree, int target[], int targetweight[]) {
  int commtype;
  int flag;
  cartattr *neighborhood;
  int i, j, t;
  //int *coordinate;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  t = neighborhood->t;
  if (maxindegree < t)
    t = maxindegree;

  j = 0;
  for (i = 0; i < t; i++) {
    // Beware of MPI_PROC_NULL issue for non-periodic Cartesians.
    /* As I understand the MPI standard, MPI_PROC_NULL is not a valid
     neighbor in distributed graphs */
    if (neighborhood->inrank[i] != MPI_PROC_NULL) {
      source[j] = neighborhood->inrank[i];
      if (sourceweight != MPI_UNWEIGHTED && neighborhood->weight != MPI_UNWEIGHTED)
        sourceweight[j] = neighborhood->weight[i];
      j++;
    }
  }

  t = neighborhood->t;
  if (maxoutdegree < t)
    t = maxoutdegree;

  j = 0;
  for (i = 0; i < t; i++) {
    if (neighborhood->outrank[i] != MPI_PROC_NULL) {
      target[j] = neighborhood->outrank[i];
      if (targetweight != MPI_UNWEIGHTED && neighborhood->weight != MPI_UNWEIGHTED)
        targetweight[j] = neighborhood->weight[i];
      j++;
    }
  }

  return MPI_SUCCESS;
}

int TUW_Cart_neighborhood_get_relative(MPI_Comm cartcomm, const int maxt,
				       int targetrelative[], int targetweight[]) {
  int commtype;
  int flag;
  cartattr *neighborhood;
  int i, t, d;
  //int *coordinate;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);
  
  t = neighborhood->t;
  if (maxt < t)
    t = maxt;

  d = neighborhood->d;
  
  for (i = 0; i < d*t; i++) {
    targetrelative[i] = neighborhood->vectors[i];
  }
  
  if (targetweight != MPI_UNWEIGHTED) {
    for (i = 0; i < t; i++) {
      targetweight[i] = neighborhood->weight[i];
    }
  }

  return MPI_SUCCESS;
}

// find functions by O(t) linear search
// can be improved to O(log t) by presorting and binary search

int TUW_Cart_neighborhood_find_rank(MPI_Comm cartcomm, const int rank, const int maxin, int sourceindex[], int *inmult,
    const int maxout, int targetindex[], int *outmult) {
  int commtype;
  cartattr *neighborhood;
  int flag;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  *inmult = findrank(rank, neighborhood->t, neighborhood->inrank, maxin, sourceindex);
  *outmult = findrank(rank, neighborhood->t, neighborhood->outrank, maxout, targetindex);

  return MPI_SUCCESS;
}

int TUW_Cart_neighborhood_find_coord(MPI_Comm cartcomm, const int relative[], const int maxin, int sourceindex[],
    int *inmult, const int maxout, int targetindex[], int *outmult) {
  //cartattr *neighborhood;
  //int *coordinate;
  //int flag;

  return MPI_SUCCESS;
}

// In bad MPI style, the collectives are of the -w variety

int TUW_Cart_neighbor_allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;
  int p;
  int i, j, k;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  if (neighborhood->trivial) {
    MPI_Aint lb, recvextent;
    int t;
    
    MPI_Type_get_extent(recvtype,&lb,&recvextent);
    
    t = neighborhood->t;
    if (neighborhood->blocking) {
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Sendrecv(sendbuf,sendcount,sendtype,neighborhood->outrank[i],TUWCARTTAG,
		     (void*)((char*)recvbuf+i*recvcount*recvextent),
		     recvcount,recvtype,neighborhood->inrank[i],TUWCARTTAG,
		     cartcomm,MPI_STATUS_IGNORE);
      }
    } else {
      MPI_Request request[2*t];
      for (i=0; i<t; i++) {
	MPI_Irecv((void*)((char*)recvbuf+i*recvcount*recvextent),
		  recvcount,recvtype,neighborhood->inrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i]);
	MPI_Isend(sendbuf,sendcount,sendtype,neighborhood->outrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i+1]);
      }
      MPI_Waitall(2*t,request,MPI_STATUSES_IGNORE);
    }
  } else {
    if (!neighborhood->precomputed) {
      schedule_cart_allgather(regular, neighborhood, sendbuf, sendcount, sendtype, recvbuf, &recvcount, NULL, &recvtype,
			      cartcomm);
      neighborhood->precomputed = 1;
    }
    
    p = neighborhood->phases;
    k = 0;
    if (neighborhood->blocking) {
      for (j = 0; j < p; j++) {
	for (i = 0; i < neighborhood->phase[j]; i++) {
	  MPI_Sendrecv(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG,
		       MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm, MPI_STATUS_IGNORE);
	  k++;
	}
      }
    } else {
      for (j = 0; j < p; j++) {
	MPI_Request request[2 * neighborhood->phase[j]];
	for (i = 0; i < neighborhood->phase[j]; i++) {
	  MPI_Irecv(MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm,
		    &request[2 * i]);
	  MPI_Isend(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG, cartcomm,
		    &request[2 * i + 1]);
	  k++;
	}
	MPI_Waitall(2 * neighborhood->phase[j], request, MPI_STATUSES_IGNORE);
      }
    }
  }

  return MPI_SUCCESS;
}

int TUW_Cart_neighbor_allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
    const int recvcount[], const int recvdisp[], MPI_Datatype recvtype, MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;
  int p;
  int i, j, k;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  if (neighborhood->trivial) {
    MPI_Aint lb, recvextent;
    int t;
    
    MPI_Type_get_extent(recvtype,&lb,&recvextent);

    t = neighborhood->t;
    if (neighborhood->blocking) {
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Sendrecv(sendbuf,sendcount,sendtype,neighborhood->outrank[i],TUWCARTTAG,
        (void*)((char*)recvbuf+recvdisp[i]*recvextent),
		     recvcount[i],recvtype,neighborhood->inrank[i],TUWCARTTAG,
		     cartcomm,MPI_STATUS_IGNORE);
      }
    } else {
      MPI_Request request[2*t];
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Irecv((void*)((char*)recvbuf+recvdisp[i]*recvextent),
		  recvcount[i],recvtype,neighborhood->inrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i]);
	MPI_Isend(sendbuf,sendcount,sendtype,neighborhood->outrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i+1]);
      }
      MPI_Waitall(2*t,request,MPI_STATUSES_IGNORE);
    }
  } else {
    if (!neighborhood->precomputed) {
      schedule_cart_allgather(regular, neighborhood, sendbuf, sendcount, sendtype, recvbuf, recvcount, recvdisp,
			      &recvtype, cartcomm);
      neighborhood->precomputed = 1;
    }
    
    p = neighborhood->phases;
    k = 0;
    if (neighborhood->blocking) {
      for (j = 0; j < p; j++) {
	for (i = 0; i < neighborhood->phase[j]; i++) {
	  MPI_Sendrecv(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG,
		       MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm, MPI_STATUS_IGNORE);
	  k++;
	}
      }
    } else {
      for (j = 0; j < p; j++) {
	MPI_Request request[2 * neighborhood->phase[j]];
	for (i = 0; i < neighborhood->phase[j]; i++) {
	  MPI_Irecv(MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm,
		    &request[2 * i]);
	  MPI_Isend(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG, cartcomm,
		    &request[2 * i + 1]);
	  k++;
	}
	MPI_Waitall(2 * neighborhood->phase[j], request, MPI_STATUSES_IGNORE);
      }
    }
  }

  return MPI_SUCCESS;
}

int TUW_Cart_neighbor_allgatherw(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
    const int recvcount[], const MPI_Aint recvdisp[], const MPI_Datatype recvtype[], MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;
  int p;
  int i, j, k;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  if (neighborhood->trivial) {
    int t = neighborhood->t;
    if (neighborhood->blocking) {
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Sendrecv(sendbuf,sendcount,sendtype,neighborhood->outrank[i],TUWCARTTAG,
		     (void*)((char*)recvbuf+recvdisp[i]),
		     recvcount[i],recvtype[i],neighborhood->inrank[i],TUWCARTTAG,
		     cartcomm,MPI_STATUS_IGNORE);
      }
    } else {
      MPI_Request request[2*t];
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Irecv((void*)((char*)recvbuf+recvdisp[i]),
		  recvcount[i],recvtype[i],neighborhood->inrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i]);
	MPI_Isend(sendbuf,sendcount,sendtype,neighborhood->outrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i+1]);
      }
      MPI_Waitall(2*t,request,MPI_STATUSES_IGNORE);
    }
  } else {
    if (!neighborhood->precomputed) {
      schedule_cart_allgather(regular, neighborhood, sendbuf, sendcount, sendtype, recvbuf, recvcount, recvdisp, recvtype,
			      cartcomm);
      neighborhood->precomputed = 1;
    }
    
    p = neighborhood->phases;
    k = 0;
    if (neighborhood->blocking) {
      for (j = 0; j < p; j++) {
	for (i = 0; i < neighborhood->phase[j]; i++) {
	  MPI_Sendrecv(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG,
		       MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm, MPI_STATUS_IGNORE);
	  k++;
	}
      }
    } else {
      for (j = 0; j < p; j++) {
	MPI_Request request[2 * neighborhood->phase[j]];
	for (i = 0; i < neighborhood->phase[j]; i++) {
	  MPI_Irecv(MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm,
		    &request[2 * i]);
	  MPI_Isend(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG, cartcomm,
		    &request[2 * i + 1]);
	  k++;
	}
	MPI_Waitall(2 * neighborhood->phase[j], request, MPI_STATUSES_IGNORE);
      }
    }
  }

  return MPI_SUCCESS;
}

int TUW_Cart_neighbor_alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;
  int p;
  int i, j, k;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  if (neighborhood->trivial) {
    MPI_Aint lb, sendextent, recvextent;
    int t;
    
    MPI_Type_get_extent(sendtype,&lb,&sendextent);
    MPI_Type_get_extent(recvtype,&lb,&recvextent);

    t = neighborhood->t;
    if (neighborhood->blocking) {
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Sendrecv((void*)((char*)sendbuf+i*sendcount*sendextent),
		     sendcount,sendtype,neighborhood->outrank[i],TUWCARTTAG,
		     (void*)((char*)recvbuf+i*recvcount*recvextent),
		     recvcount,recvtype,neighborhood->inrank[i],TUWCARTTAG,
		     cartcomm,MPI_STATUS_IGNORE);
      }
    } else {
      MPI_Request request[2*t];
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Irecv((void*)((char*)recvbuf+i*recvcount*recvextent),
		  recvcount,recvtype,neighborhood->inrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i]);
	MPI_Isend((void*)((char*)sendbuf+i*sendcount*sendextent),
		  sendcount,sendtype,neighborhood->outrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i+1]);
      }
      MPI_Waitall(2*t,request,MPI_STATUSES_IGNORE);
    }
  } else {
    if (!neighborhood->precomputed) {
      schedule_cart_alltoall(regular, neighborhood, sendbuf, &sendcount, NULL, &sendtype, recvbuf, &recvcount, NULL,
			     &recvtype, cartcomm);
      neighborhood->precomputed = 1;
    }
    
    p = neighborhood->phases;
    k = 0;
    if (neighborhood->blocking) {
      for (j = 0; j < p; j++) {
	for (i = 0; i < neighborhood->phase[j]; i++) {
	  MPI_Sendrecv(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG,
		       MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm, MPI_STATUS_IGNORE);
	  k++;
	}
      }
    } else {
      for (j = 0; j < p; j++) {
	MPI_Request request[2 * neighborhood->phase[j]];
	for (i = 0; i < neighborhood->phase[j]; i++) {
	  MPI_Irecv(MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm,
		    &request[2 * i]);
	  MPI_Isend(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG, cartcomm,
		    &request[2 * i + 1]);
	  k++;
	}
	MPI_Waitall(2 * neighborhood->phase[j], request, MPI_STATUSES_IGNORE);
      }
    }
  }

  return MPI_SUCCESS;
}

int TUW_Cart_neighbor_alltoallv(void *sendbuf, const int sendcount[], const int senddisp[], MPI_Datatype sendtype,
    void *recvbuf, const int recvcount[], const int recvdisp[], MPI_Datatype recvtype, MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;
  int p;
  int i, j, k;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  if (neighborhood->trivial) {
    MPI_Aint lb, sendextent, recvextent;
    int t;
    
    MPI_Type_get_extent(sendtype,&lb,&sendextent);
    MPI_Type_get_extent(recvtype,&lb,&recvextent);
    
    t = neighborhood->t;
    if (neighborhood->blocking) {
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Sendrecv((void*)((char*)sendbuf+senddisp[i]*sendextent),
		     sendcount[i],sendtype,neighborhood->outrank[i],TUWCARTTAG,
		     (void*)((char*)recvbuf+recvdisp[i]*recvextent),
		     recvcount[i],recvtype,neighborhood->inrank[i],TUWCARTTAG,
		     cartcomm,MPI_STATUS_IGNORE);
      }
    } else {
      MPI_Request request[2*t];
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Irecv((void*)((char*)recvbuf+recvdisp[i]*recvextent),
		  recvcount[i],recvtype,neighborhood->inrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i]);
	MPI_Isend((void*)((char*)sendbuf+senddisp[i]*sendextent),
		  sendcount[i],sendtype,neighborhood->outrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i+1]);
      }
      MPI_Waitall(2*t,request,MPI_STATUSES_IGNORE);
    }
  } else {
    if (!neighborhood->precomputed) {
      schedule_cart_alltoall(vector, neighborhood, sendbuf, sendcount, senddisp, &sendtype, recvbuf, recvcount, recvdisp,
			     &recvtype, cartcomm);
      neighborhood->precomputed = 1;
    }
    
    p = neighborhood->phases;
    k = 0;
    for (j = 0; j < p; j++) {
      MPI_Request request[2 * neighborhood->phase[j]];
      for (i = 0; i < neighborhood->phase[j]; i++) {
	MPI_Irecv(MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm,
		  &request[2 * i]);
	MPI_Isend(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG, cartcomm,
		  &request[2 * i + 1]);
	k++;
      }
      MPI_Waitall(2 * neighborhood->phase[j], request, MPI_STATUSES_IGNORE);
    }
  }

  return MPI_SUCCESS;
}

int TUW_Cart_neighbor_alltoallw(void *sendbuf, const int sendcount[], const MPI_Aint senddisp[],
    const MPI_Datatype sendtype[], void *recvbuf, const int recvcount[], const MPI_Aint recvdisp[],
    const MPI_Datatype recvtype[], MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;
  int p;
  int i, j, k;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  if (neighborhood->trivial) {
    int t = neighborhood->t;
    if (neighborhood->blocking) {
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Sendrecv((void*)((char*)sendbuf+senddisp[i]),
		     sendcount[i],sendtype[i],neighborhood->outrank[i],TUWCARTTAG,
		     (void*)((char*)recvbuf+recvdisp[i]),
		     recvcount[i],recvtype[i],neighborhood->inrank[i],TUWCARTTAG,
		     cartcomm,MPI_STATUS_IGNORE);
      }
    } else {
      MPI_Request request[2*t];
      for (i=0; i<t; i++) {
	// Deadlock free by specification
	MPI_Irecv((void*)((char*)recvbuf+recvdisp[i]),
		  recvcount[i],recvtype[i],neighborhood->inrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i]);
	MPI_Isend((void*)((char*)sendbuf+senddisp[i]),
		  sendcount[i],sendtype[i],neighborhood->outrank[i],TUWCARTTAG,
		  cartcomm,&request[2*i+1]);
      }
      MPI_Waitall(2*t,request,MPI_STATUSES_IGNORE);
    }
  } else {
    if (!neighborhood->precomputed) {
      schedule_cart_alltoall(wector, neighborhood, sendbuf, sendcount, senddisp, sendtype, recvbuf, recvcount, recvdisp,
			     recvtype, cartcomm);
      neighborhood->precomputed = 1;
    }
    
    p = neighborhood->phases;
    k = 0;
    if (neighborhood->blocking) {
      for (j = 0; j < p; j++) {
	for (i = 0; i < neighborhood->phase[j]; i++) {
	  MPI_Sendrecv(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG,
		       MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm, MPI_STATUS_IGNORE);
	  k++;
	}
      }
    } else {
      for (j = 0; j < p; j++) {
	MPI_Request request[2 * neighborhood->phase[j]];
	for (i = 0; i < neighborhood->phase[j]; i++) {
	  MPI_Irecv(MPI_BOTTOM, 1, neighborhood->recvtype[k], neighborhood->recvrank[k], TUWCARTTAG, cartcomm,
		    &request[2 * i]);
	  MPI_Isend(MPI_BOTTOM, 1, neighborhood->sendtype[k], neighborhood->sendrank[k], TUWCARTTAG, cartcomm,
		    &request[2 * i + 1]);
	  k++;
	}
	MPI_Waitall(2 * neighborhood->phase[j], request, MPI_STATUSES_IGNORE);
      }
    }
  }

  return MPI_SUCCESS;
}

/* Precompute schedule */
int TUW_Cart_neighbor_allgather_init(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  if (!neighborhood->trivial) {
    if (neighborhood->precomputed)
    schedule_free(neighborhood);
    schedule_cart_allgather(regular, neighborhood, sendbuf, sendcount, sendtype, recvbuf, &recvcount, NULL, &recvtype,
			    cartcomm);
    neighborhood->precomputed = 1;
  }
  
  return MPI_SUCCESS;
}

int TUW_Cart_neighbor_alltoall_init(void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  if (!neighborhood->trivial) {
    if (neighborhood->precomputed)
      schedule_free(neighborhood);
    schedule_cart_alltoall(regular, neighborhood, sendbuf, &sendcount, NULL, &sendtype, recvbuf, &recvcount, NULL,
			   &recvtype, cartcomm);
    neighborhood->precomputed = 1;
  }

  return MPI_SUCCESS;
}

int TUW_Cart_neighbor_alltoallv_init(void *sendbuf, int sendcount[], int senddisp[], MPI_Datatype sendtype,
    void *recvbuf, int recvcount[], int recvdisp[], MPI_Datatype recvtype, MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  if (!neighborhood->trivial) {
    if (neighborhood->precomputed)
      schedule_free(neighborhood);
    schedule_cart_alltoall(vector, neighborhood, sendbuf, sendcount, senddisp, &sendtype, recvbuf, recvcount, recvdisp,
			   &recvtype, cartcomm);
    neighborhood->precomputed = 1;
  }

  return MPI_SUCCESS;
}

int TUW_Cart_neighbor_alltoallw_init(void *sendbuf, int sendcount[], int senddisp[], MPI_Datatype sendtype[],
    void *recvbuf, int recvcount[], int recvdisp[], MPI_Datatype recvtype[], MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  if (!neighborhood->trivial) {
    if (neighborhood->precomputed)
      schedule_free(neighborhood);
    schedule_cart_alltoall(wector, neighborhood, sendbuf, sendcount, senddisp, sendtype, recvbuf, recvcount, recvdisp,
			   recvtype, cartcomm);
    neighborhood->precomputed = 1;
  }

  return MPI_SUCCESS;
}

int TUW_Cart_neighbor_alltoall_free(MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  schedule_free(neighborhood);

  return MPI_SUCCESS;
}

int TUW_Cart_neighbor_allgather_free(MPI_Comm cartcomm) {
  int commtype;
  int flag;
  cartattr *neighborhood;

  MPI_Topo_test(cartcomm, &commtype); // Cartesian error check here
  if (commtype != MPI_CART)
    return MPI_ERR_TOPOLOGY;

  MPI_Comm_get_attr(cartcomm, cartkey(), &neighborhood, &flag);
  if (!flag)
    return MPI_ERR_TOPOLOGY;
  assert(neighborhood->commtype==TUW_CART_GRAPH);

  schedule_free(neighborhood);

  return MPI_SUCCESS;
}

// Rank search function

static int findrank(const int rank, const int r, int ranks[], const int max, int found[]) {
  int i, j;
  int m;

  j = 0;
  m = 0;
  for (i = 0; i < r; i++) {
    if (rank == ranks[i]) {
      m++;
      if (j < max)
        found[j++] = i;
    }
  }

  return m;
}

// variant decoding functions

inline int variantcnt(collvariant variant, int i, const int *count) {
  int ret_cnt;

  switch (variant) {
  case regular:
    ret_cnt = *count;
    break;
  case vector:
  case wector:
    ret_cnt = count[i];
    break;
  default:
    assert(0);
  }

  return ret_cnt;
}

inline MPI_Aint variantbuf(collvariant variant, int i, void *buf, const int *count, const void *disp,
    const MPI_Aint extent) {

  MPI_Aint ret_buf;

  switch (variant) {
  case regular:
    ret_buf = (MPI_Aint) ((char*) buf + i * (*count) * extent);
    break;
  case vector:
    ret_buf = (MPI_Aint) ((char*) buf + ((int*) disp)[i] * extent);
    break;
  case wector:
    ret_buf = (MPI_Aint) ((char*) buf + ((MPI_Aint*) disp)[i]);
    break;
  default:
    assert(0);
  }

  return ret_buf;
}

inline MPI_Datatype varianttyp(collvariant variant, int i, const MPI_Datatype *type) {

  MPI_Datatype ret_type;

  switch (variant) {
  case regular:
  case vector:
    ret_type = *type;
    break;
  case wector:
    ret_type = type[i];
    break;
  default:
    assert(0);
  }

  return ret_type;
}

// Schedule computations

static void schedule_free(cartattr *neighborhood) {
  if (neighborhood->precomputed) {
    // Erase old schedule
    int p, i, j, k;

    p = neighborhood->phases;
    k = 0;
    for (j = 0; j < p; j++) {
      for (i = 0; i < neighborhood->phase[j]; i++) {
        MPI_Type_free(&neighborhood->sendtype[k]);
        MPI_Type_free(&neighborhood->recvtype[k]);
        k++;
      }
    }

    free(neighborhood->sendtype);
    free(neighborhood->recvtype);
    free(neighborhood->sendrank);
    free(neighborhood->recvrank);

    free(neighborhood->phase);

    free(neighborhood->tempbuf);

    neighborhood->precomputed = 0;
  }
}

static int nonzeros(const int d, const int vector[]) {
  int j;
  int nz;

  nz = 0;
  for (j = 0; j < d; j++)
    if (vector[j] != 0)
      nz++;

  return nz;
}

static int maxcoordinate(const int n, const int d, const int vectors[], const int j) {
  int i;
  int mc;

  mc = 0;
  for (i = 0; i < n; i++) {
    int a = abs(vectors[i * d + j]);
    if (a > mc)
      mc = a;
  }

  return mc;
}

static void bucketsort(const int n, const int d, const int vectors[], const int mc, const int j, int order[]) {
  int bucket[2 * (mc + 1)];
  int i;

  for (i = 0; i < 2 * (mc + 1); i++)
    bucket[i] = 0;
  for (i = 0; i < n; i++)
    bucket[mc + vectors[order[i] * d + j]]++;
  int bc = 0;
  for (i = 0; i < 2 * (mc + 1); i++) {
    int bb = bucket[i];
    bucket[i] = bc;
    bc += bb;
  }
  int neworder[n];
  for (i = 0; i < n; i++)
    neworder[bucket[mc + vectors[order[i] * d + j]]++] = order[i];
  for (i = 0; i < n; i++)
    order[i] = neworder[i];
}

/* basic alltoall schedule */
/* direct communication along dimensions */
/* complexity O(sd+c), c = vector max. coordinate */
static void schedule_cart_alltoall(collvariant variant, cartattr *neighborhood, void *sendbuf, const int *sendcount,
    const void *senddisp, const MPI_Datatype *sendtype, void *recvbuf, const int *recvcount, const void *recvdisp,
    const MPI_Datatype *recvtype, MPI_Comm cartcomm) {
  int t, d, b, e;
  int i, j, k;

  int sendrank, recvrank;
  MPI_Aint lb, sendextent, recvextent;

  void *tempdisp;
  int dispv = 0;
  MPI_Aint dispw;

  if (variant != wector) {
    MPI_Type_get_extent(*sendtype, &lb, &sendextent);
    MPI_Type_get_extent(*recvtype, &lb, &recvextent);
  }

  d = neighborhood->d;
  t = neighborhood->t;

  neighborhood->phases = d + 1;
  neighborhood->phase = (int*) malloc((d + 1) * sizeof(int));

  int *vectors = neighborhood->vectors;
  assert(vectors!=NULL);

  int firstround[t]; // first round for neighbor i
  int hops[t];       // number of hops for neighbor i

  int order[t];

  // compute number and sizes of phases

  k = 0;
  for (j = 0; j < d; j++) {
    // sort after dimension j
    int mc = maxcoordinate(t, d, vectors, j);
    for (i = 0; i < t; i++)
      order[i] = i;
    bucketsort(t, d, vectors, mc, j, order);
    i = 0;
    e = 0;
    int c = vectors[order[i] * d + j];
    for (i = 1; i < t; i++) {
      if (c != vectors[order[i] * d + j]) {
        if (c != 0) {
          k++;
          e++;
        }
        c = vectors[order[i] * d + j];
      }
    }
    if (c != 0) {
      k++;
      e++;
    }
    neighborhood->phase[j] = e;
    //fprintf(stderr,"Exchanges %d in phase %d\n",e,j);
  }
  int kk = k;
  k++; // for local copy

  neighborhood->sendtype = (MPI_Datatype*) malloc(k * sizeof(MPI_Datatype));
  neighborhood->recvtype = (MPI_Datatype*) malloc(k * sizeof(MPI_Datatype));
  neighborhood->sendrank = (int*) malloc(k * sizeof(int));
  neighborhood->recvrank = (int*) malloc(k * sizeof(int));

  switch (variant) {
  case regular:
    neighborhood->tempbuf = (void*) malloc(t * (*recvcount) * recvextent);
    break;
  case vector:
    tempdisp = (void*) malloc(t * sizeof(int));
    dispv = 0;
    for (i = 0; i < t; i++) {
      ((int*) tempdisp)[i] = dispv;
      dispv += recvcount[i];
    }
    neighborhood->tempbuf = (void*) malloc(dispv * recvextent);
    break;
  case wector:
    tempdisp = (void*) malloc(t * sizeof(MPI_Aint));
    dispw = 0;
    for (i = 0; i < t; i++) {
      MPI_Type_get_extent(recvtype[i], &lb, &recvextent);
      ((MPI_Aint*) tempdisp)[i] = dispw;
      dispv += recvcount[i] * recvextent;
    }
    neighborhood->tempbuf = (void*) malloc(dispw);
    break;
  default:
    assert(0);
  }

  int sc[t], rc[t];
  MPI_Aint sd[t], rd[t];
  MPI_Datatype st[t], rt[t];

  for (i = 0; i < t; i++) {
    for (j = 0; j < d; j++)
      if (vectors[i * d + j] != 0)
        break;
    firstround[i] = j;
    hops[i] = nonzeros(d, vectors + i * d);
  }

  k = 0;
  for (j = 0; j < d; j++) {
    // sort after dimension j
    int mc = maxcoordinate(t, d, vectors, j);
    for (i = 0; i < t; i++)
      order[i] = i;
    bucketsort(t, d, vectors, mc, j, order);

    b = 0;
    for (i = 0; i < t; i++) {
      int ii = order[i];
      int c = vectors[order[i] * d + j];

      if (c != 0) {
        if (firstround[ii] == j) {
          //printf("sendbuf[%d] to ",ii);
          //sc[b] = sendcount;
          //sd[b] = (MPI_Aint)((char*)sendbuf+ii*sendcount*sendextent);
          //st[b] = sendtype;
          sc[b] = variantcnt(variant, ii, sendcount);
          sd[b] = variantbuf(variant, ii, sendbuf, sendcount, senddisp, sendextent);
          st[b] = varianttyp(variant, ii, sendtype);
          if (hops[ii] % 2 == 0) {
            //printf("tempbuf[%d], offset %d\n",ii,c);
            //rc[b] = recvcount;
            //rd[b] = (MPI_Aint)
            //  ((char*)neighborhood->tempbuf+ii*recvcount*recvextent);
            //rt[b] = recvtype;
            rc[b] = variantcnt(variant, ii, recvcount);
            rd[b] = variantbuf(variant, ii, neighborhood->tempbuf, recvcount, tempdisp, recvextent);
            rt[b] = varianttyp(variant, ii, recvtype);
          } else {
            //printf("recvbuf[%d], offset %d\n",ii,c);
            //rc[b] = recvcount;
            //rd[b] = (MPI_Aint)((char*)recvbuf+ii*recvcount*recvextent);
            //rt[b] = recvtype;
            rc[b] = variantcnt(variant, ii, recvcount);
            rd[b] = variantbuf(variant, ii, recvbuf, recvcount, recvdisp, recvextent);
            rt[b] = varianttyp(variant, ii, recvtype);
          }
        } else {
          if (hops[ii] % 2 == 0) {
            //printf("recvbuf[%d] to tempbuf[%d], offset %d\n",ii,ii,c);
            //sc[b] = recvcount;
            //sd[b] = (MPI_Aint)((char*)recvbuf+ii*recvcount*recvextent);
            //st[b] = recvtype;
            //rc[b] = recvcount;
            //rd[b] = (MPI_Aint)
            //  ((char*)neighborhood->tempbuf+ii*recvcount*recvextent);
            //rt[b] = recvtype;
            sc[b] = variantcnt(variant, ii, recvcount);
            sd[b] = variantbuf(variant, ii, recvbuf, recvcount, recvdisp, recvextent);
            st[b] = varianttyp(variant, ii, recvtype);
            rc[b] = variantcnt(variant, ii, recvcount);
            rd[b] = variantbuf(variant, ii, neighborhood->tempbuf, recvcount, tempdisp, recvextent);
            rt[b] = varianttyp(variant, ii, recvtype);
          } else {
            //printf("tempbuf[%d] to recvbuf[%d], offset %d\n",ii,ii,c);
            //sc[b] = recvcount;
            //sd[b] = (MPI_Aint)
            //  ((char*)neighborhood->tempbuf+ii*recvcount*recvextent);
            //st[b] = recvtype;
            //rc[b] = recvcount;
            //rd[b] = (MPI_Aint)((char*)recvbuf+ii*recvcount*recvextent);
            //rt[b] = recvtype;
            sc[b] = variantcnt(variant, ii, recvcount);
            sd[b] = variantbuf(variant, ii, neighborhood->tempbuf, recvcount, tempdisp, recvextent);
            st[b] = varianttyp(variant, ii, recvtype);
            rc[b] = variantcnt(variant, ii, recvcount);
            rd[b] = variantbuf(variant, ii, recvbuf, recvcount, recvdisp, recvextent);
            rt[b] = varianttyp(variant, ii, recvtype);
          }
        }
        b++;
        hops[ii]--; // one hop less

        int cc;
        if (i < t - 1)
          cc = vectors[order[i + 1] * d + j];
        if (i == t - 1 || cc != c) {
          MPI_Cart_shift(cartcomm, j, c, &recvrank, &sendrank);

          MPI_Type_create_struct(b, sc, sd, st, &neighborhood->sendtype[k]);
          MPI_Type_commit(&neighborhood->sendtype[k]);
          neighborhood->sendrank[k] = sendrank;

          MPI_Type_create_struct(b, rc, rd, rt, &neighborhood->recvtype[k]);
          MPI_Type_commit(&neighborhood->recvtype[k]);
          neighborhood->recvrank[k] = recvrank;

          k++;
          b = 0;
        }
      }
    }
  }
  assert(k == kk);

  // local copies
  b = 0;
  for (i = 0; i < t; i++) {
    if (firstround[i] == d) {
      //printf("Local: sendbuf[%d] to recvbuf[%d]\n",i,i);
      //sc[b] = sendcount;
      //sd[b] = (MPI_Aint)((char*)sendbuf+i*sendcount*sendextent);
      //st[b] = sendtype;
      //rc[b] = recvcount;
      //rd[b] = (MPI_Aint)((char*)recvbuf+i*recvcount*recvextent);
      //rt[b] = recvtype;
      sc[b] = variantcnt(variant, i, sendcount);
      sd[b] = variantbuf(variant, i, sendbuf, sendcount, senddisp, sendextent);
      st[b] = varianttyp(variant, i, sendtype);
      rc[b] = variantcnt(variant, i, recvcount);
      rd[b] = variantbuf(variant, i, recvbuf, recvcount, recvdisp, recvextent);
      rt[b] = varianttyp(variant, i, recvtype);

      b++;
    }
  }
  if (b > 0) {
    neighborhood->phase[d] = 1;
    MPI_Comm_rank(cartcomm, &sendrank);
    recvrank = sendrank;

    MPI_Type_create_struct(b, sc, sd, st, &neighborhood->sendtype[k]);
    MPI_Type_commit(&neighborhood->sendtype[k]);
    neighborhood->sendrank[k] = sendrank;

    MPI_Type_create_struct(b, rc, rd, rt, &neighborhood->recvtype[k]);
    MPI_Type_commit(&neighborhood->recvtype[k]);
    neighborhood->recvrank[k] = recvrank;

    k++;
    b = 0;
  } else
    neighborhood->phase[d] = 0;

  switch (variant) {
  case regular:
    break;
  case vector:
  case wector:
    free(tempdisp);
    break;
  default:
    assert(0);
  }
}

// sort recursively down the dimensions
void treesort(const int n, const int j, const int d, int vectors[], int order[]) {
  if (j == d || n == 0)
    return;

  int i, ii;
  int c, mc;

  mc = 0;
  for (i = 0; i < n; i++) {
    int a = abs(vectors[order[i] * d + j]); // find max coordinate in order
    if (a > mc)
      mc = a;
  }

  bucketsort(n, d, vectors, mc, j, order);

  i = 0;
  c = vectors[order[i] * d + j];
  for (ii = 1; ii < n; ii++) {
    if (c != vectors[order[ii] * d + j]) {
      treesort(ii - i, j + 1, d, vectors, order + i);
      i = ii;
      c = vectors[order[i] * d + j];
    }
  }
  treesort(ii - i, j + 1, d, vectors, order + i);
}

// make dfs tree
void treecons(const int n, const int j, const int d, int vectors[], int order[], const int parent, int **pp) {
  if (j == d || n == 0)
    return;

  int i, ii;
  int c;

  i = 0;
  c = vectors[order[i] * d + j];
  for (ii = 1; ii < n; ii++) {
    if (c != vectors[order[ii] * d + j]) {
      pp[j][order[i]] = parent;
      treecons(ii - i, j + 1, d, vectors, order + i, order[i], pp);
      i = ii;
      c = vectors[order[i] * d + j];
    }
  }
  pp[j][order[i]] = parent;
  treecons(ii - i, j + 1, d, vectors, order + i, order[i], pp);
}

void schedule_cart_allgather(collvariant variant, cartattr *neighborhood, void *sendbuf, int sendcount,
    MPI_Datatype sendtype, void *recvbuf, const int *recvcount, const void *recvdisp, const MPI_Datatype *recvtype,
    MPI_Comm cartcomm) {
  int t, d, e;
  int i, j, k;
  int n;

  int sendrank, recvrank;
  MPI_Aint lb, sendextent, recvextent;

  void *tempdisp;
  int dispv;
  MPI_Aint dispw;

  MPI_Type_get_extent(sendtype, &lb, &sendextent);

  if (variant != wector) {
    MPI_Type_get_extent(*recvtype, &lb, &recvextent);
  }

  d = neighborhood->d;
  t = neighborhood->t;

  neighborhood->phases = d + 1;
  neighborhood->phase = (int*) malloc((d + 1) * sizeof(int));

  int *vectors = neighborhood->vectors;
  assert(vectors!=NULL);

  int **pp;
  pp = (int**) malloc(d * sizeof(int*));
  pp[0] = (int*) malloc(d * t * sizeof(int));
  for (j = 1; j < d; j++)
    pp[j] = pp[j - 1] + t;
  for (j = 0; j < d; j++) {
    for (i = 0; i < t; i++)
      pp[j][i] = -1;
  }

  int order[t];
  for (i = 0; i < t; i++)
    order[i] = i;

  // depth first search construction of virtual tree
  treesort(t, 0, d, vectors, order);
  treecons(t, 0, d, vectors, order, -1, pp);

  int firstround[t], nextround[t];

  for (i = 0; i < t; i++) {
    for (j = 0; j < d; j++)
      if (vectors[i * d + j] != 0)
        break;
    if (j == d)
      j--;
    firstround[i] = j;
    nextround[i] = firstround[i];
  }

  // compute number and sizes of phases

  k = 0;
  for (j = 0; j < d; j++) {
    int c;
    //int ii;

    int mc = maxcoordinate(t, d, vectors, j);

    int last[2 * (mc + 1)], first[2 * (mc + 1)], firsti[2 * (mc + 1)];
    int next[t];

    for (i = 0; i < 2 * (mc + 1); i++) {
      last[i] = -1;
      first[i] = -1;
    }

    // breadth first search in level
    for (i = 0; i < t;) {
      int ii = order[i];
      c = vectors[ii * d + j];
      if (last[mc + c] == -1) {
        first[mc + c] = ii;
        firsti[mc + c] = i;
        last[mc + c] = ii;
        next[ii] = -1;
      } else {
        next[last[mc + c]] = ii;
        last[mc + c] = ii;
        next[ii] = -1;
      }
      i++;
      while (i < t && c == vectors[order[i] * d + j] && pp[j][order[i]] == -1)
        i++;
    }

    n = 0;
    e = 0;
    for (i = 0; i < 2 * (mc + 1); i++) {
      int b = first[i];
      while (b != -1) {
        if (i - mc == 0 && j > firstround[b] && j >= nextround[b]) {
          int kk = 1;
          while (j + kk < d && vectors[b * d + j + kk] == 0)
            kk++;
          nextround[b] = j + kk - (kk % 2);
        }

        if (j >= nextround[b]) {
          n++;
        }

        b = next[b];
      }

      if (n > 0) {
        k++;
        e++;
        n = 0;
      }
    }
    neighborhood->phase[j] = e;
  }
  int kk = k;
  k++; // for local copy

  neighborhood->sendtype = (MPI_Datatype*) malloc(k * sizeof(MPI_Datatype));
  neighborhood->recvtype = (MPI_Datatype*) malloc(k * sizeof(MPI_Datatype));
  neighborhood->sendrank = (int*) malloc(k * sizeof(int));
  neighborhood->recvrank = (int*) malloc(k * sizeof(int));

  switch (variant) {
  case regular:
    neighborhood->tempbuf = (void*) malloc(t * (*recvcount) * recvextent);
    break;
  case vector:
    tempdisp = (void*) malloc(t * sizeof(int));
    dispv = 0;
    for (i = 0; i < t; i++) {
      ((int*) tempdisp)[i] = dispv;
      dispv += recvcount[i];
    }
    neighborhood->tempbuf = (void*) malloc(dispv * recvextent);
    break;
  case wector:
    tempdisp = (void*) malloc(t * sizeof(MPI_Aint));
    dispw = 0;
    for (i = 0; i < t; i++) {
      MPI_Type_get_extent(recvtype[i], &lb, &recvextent);
      ((MPI_Aint*) tempdisp)[i] = dispw;
      dispw += recvcount[i] * recvextent;
    }
    neighborhood->tempbuf = (void*) malloc(dispw);
    break;
  default:
    assert(0);
  }

  int sc[t], rc[t];
  MPI_Aint sd[t], rd[t];
  MPI_Datatype st[t], rt[t];

  for (i = 0; i < t; i++)
    nextround[i] = firstround[i];

  k = 0;
  for (j = 0; j < d; j++) {
    int c;
    //int ii;

    int mc = maxcoordinate(t, d, vectors, j);

    int last[2 * (mc + 1)], first[2 * (mc + 1)], firsti[2 * (mc + 1)];
    int next[t];

    for (i = 0; i < 2 * (mc + 1); i++) {
      last[i] = -1;
      first[i] = -1;
    }

    // breadth first search in level
    for (i = 0; i < t;) {
      int ii = order[i];
      c = vectors[ii * d + j];
      if (last[mc + c] == -1) {
        first[mc + c] = ii;
        firsti[mc + c] = i;
        last[mc + c] = ii;
        next[ii] = -1;
      } else {
        next[last[mc + c]] = ii;
        last[mc + c] = ii;
        next[ii] = -1;
      }
      i++;
      while (i < t && c == vectors[order[i] * d + j] && pp[j][order[i]] == -1)
        i++;
    }

    n = 0;
    for (i = 0; i < 2 * (mc + 1); i++) {
      int b = first[i];
      while (b != -1) {
        if (i - mc == 0 && j > firstround[b] && j >= nextround[b]) {
          int kk = 1;
          while (j + kk < d && vectors[b * d + j + kk] == 0)
            kk++;
          nextround[b] = j + kk - (kk % 2);
        }
        if (j >= nextround[b]) {
          if (firstround[b] == j) {
            //printf("sendbuf ");
            sc[n] = sendcount;
            sd[n] = (MPI_Aint) sendbuf;
            st[n] = sendtype;
          } else {
            assert(pp[j][b] != -1);
            int ii = pp[j][b];
            if ((d - j - 1) % 2 == 0) {
              //printf("interbuf[%d] ",pp[j][b]);
              //sc[n] = recvcount;
              //sd[n] = (MPI_Aint)
              //  ((char*)neighborhood->tempbuf+pp[j][b]*recvcount*recvextent);
              //st[n] = recvtype;
              sc[n] = variantcnt(variant, ii, recvcount);
              sd[n] = variantbuf(variant, ii, neighborhood->tempbuf, recvcount, tempdisp, recvextent);
              st[n] = varianttyp(variant, ii, recvtype);
            } else {
              //printf("recvbuf[%d] ",pp[j][b]);
              //sc[n] = recvcount;
              // sd[n] = (MPI_Aint)((char*)recvbuf+pp[j][b]*recvcount*recvextent);
              //st[n] = recvtype;
              sc[n] = variantcnt(variant, ii, recvcount);
              sd[n] = variantbuf(variant, ii, recvbuf, recvcount, recvdisp, recvextent);
              st[n] = varianttyp(variant, ii, recvtype);
            }
          }
          if ((d - j) % 2 == 1) {
            //printf("to recvbuf[%d], offset %d\n",b,i-mc);
            //rc[n] = recvcount;
            //rd[n] = (MPI_Aint)((char*)recvbuf+b*recvcount*recvextent);
            //rt[n] = recvtype;
            rc[n] = variantcnt(variant, b, recvcount);
            rd[n] = variantbuf(variant, b, recvbuf, recvcount, recvdisp, recvextent);
            rt[n] = varianttyp(variant, b, recvtype);
          } else {
            //printf("to interbuf[%d], offset %d\n",b,i-mc);
            //rc[n] = recvcount;
            //rd[n] = (MPI_Aint)
            //  ((char*)neighborhood->tempbuf+b*recvcount*recvextent);
            //rt[n] = recvtype;
            rc[n] = variantcnt(variant, b, recvcount);
            rd[n] = variantbuf(variant, b, neighborhood->tempbuf, recvcount, tempdisp, recvextent);
            rt[n] = varianttyp(variant, b, recvtype);
          }
          n++;
        }

        b = next[b];
      }

      if (n > 0) {
        MPI_Cart_shift(cartcomm, j, i - mc, &recvrank, &sendrank);

        MPI_Type_create_struct(n, sc, sd, st, &neighborhood->sendtype[k]);
        MPI_Type_commit(&neighborhood->sendtype[k]);
        neighborhood->sendrank[k] = sendrank;

        MPI_Type_create_struct(n, rc, rd, rt, &neighborhood->recvtype[k]);
        MPI_Type_commit(&neighborhood->recvtype[k]);
        neighborhood->recvrank[k] = recvrank;

        k++;
        n = 0;
      }
    }
  }
  assert(k == kk);

  // local copies
  n = 0;
  for (i = 0; i < t;) {
    int ii = order[i];
    int c = vectors[ii * d];
    i++;
    while (i < t && vectors[order[i] * d] == c && pp[d - 1][order[i]] == -1) {
      //nonzeros(d,&vectors[order[i]*d])!=0) {
      //printf("Local: recvbuf[%d] to recvbuf[%d]\n",ii,order[i]);
      //sc[n] = sendcount;
      //sd[n] = (MPI_Aint)((char*)recvbuf+ii*recvcount*recvextent);
      //st[n] = sendtype;
      sc[n] = variantcnt(variant, ii, recvcount);
      sd[n] = variantbuf(variant, ii, recvbuf, recvcount, recvdisp, recvextent);
      st[n] = varianttyp(variant, ii, recvtype);
      //rc[n] = recvcount;
      //rd[n] = (MPI_Aint)((char*)recvbuf+order[i]*recvcount*recvextent);
      //rt[n] = recvtype;
      rc[n] = variantcnt(variant, order[i], recvcount);
      rd[n] = variantbuf(variant, order[i], recvbuf, recvcount, recvdisp, recvextent);
      rt[n] = varianttyp(variant, order[i], recvtype);

      i++;
      n++;
    }
  }

  if (n > 0) {
    neighborhood->phase[d] = 1;
    MPI_Comm_rank(cartcomm, &sendrank);
    recvrank = sendrank;

    MPI_Type_create_struct(n, sc, sd, st, &neighborhood->sendtype[k]);
    MPI_Type_commit(&neighborhood->sendtype[k]);
    neighborhood->sendrank[k] = sendrank;

    MPI_Type_create_struct(n, rc, rd, rt, &neighborhood->recvtype[k]);
    MPI_Type_commit(&neighborhood->recvtype[k]);
    neighborhood->recvrank[k] = recvrank;

    k++;
    n = 0;
  } else
    neighborhood->phase[d] = 0;

  free(pp[0]);
  free(pp);

  switch (variant) {
  case regular:
    break;
  case vector:
  case wector:
    free(tempdisp);
    break;
  default:
    assert(0);
  }
}
