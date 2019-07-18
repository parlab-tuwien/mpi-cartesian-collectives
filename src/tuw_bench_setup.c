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

#include <stdio.h>
#include <stdlib.h>

#include "tuwcartesians.h"
#include "tuw_bench_setup.h"

#define REORDER 0     // no reorder
#define PERIODICITY 1 // only tori supported

const int F = -1; //(-(N/2)) // offset of first neighbor per dimension (F,F+1,...,F+N-1)
const int R = 1;        // repetitions of neighborhood (to get duplicate neighbors)
const int D = 1;        // distance multiplier

/**
 * store them here as we need them in several functions
 *
 * there are for the communicator creation
 */
static int *period;
static int *dimension;
static int *coords;
static int mmm; // this is only for validation

static int *target;
static int *source;

const char *func_names[] =
{
  "TUW_Cart_neighbor_alltoall",
  "TUW_Cart_neighbor_allgather",
  "TUW_Cart_neighbor_alltoallv",
  "MPI_Neighbor_alltoall",
  "MPI_Neighbor_allgather",
  "MPI_Neighbor_alltoallv",
  "MPI_Ineighbor_alltoall",
  "MPI_Ineighbor_allgather",
  "MPI_Ineighbor_alltoallv"
};


int get_t(int d, int nb_neighbor) {
  int t;
  int j;
  t = nb_neighbor;
  for (j = 1; j < d; j++)
    t *= nb_neighbor;
  t *= R;
  return t;
}

void setup_comm_buffers_alltoall(nh_exp_desc_t *desc) {
  desc->sendbuf = (type_t*) malloc(desc->t * desc->m * sizeof(type_t));
  desc->recvbuf = (type_t*) malloc(desc->t * desc->m * sizeof(type_t));
  desc->em = desc->t * desc->m;

  if( desc->init_buffers == 1 ) {
    int i;
    for (i=0; i<desc->t*desc->m; i++) {
      ((type_t*)desc->sendbuf)[i] = (type_t)(desc->np+i);
      ((type_t*)desc->recvbuf)[i] = (type_t)-1;
    }
  }

}

int check_buffers_alltoall(nh_exp_desc_t *desc, nh_exp_desc_t *desc2) {
  int same = 0;
  int i;
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (i=0; i<desc->t*desc->m; i++) {
    if (((type_t*)desc->recvbuf)[i]!=((type_t*)desc2->recvbuf)[i]) {
      fprintf(stderr, "Rank %d %d!=%d at pos %i\n", rank, ((type_t*) desc->recvbuf)[i], ((type_t*) desc2->recvbuf)[i],
          i);
      same = 1;
      break;
    }
  }
  return same;
}

void free_comm_buffers_alltoall(nh_exp_desc_t *desc) {
  free(desc->sendbuf);
  free(desc->recvbuf);
}


void setup_comm_buffers_allgather(nh_exp_desc_t *desc) {
  desc->sendbuf = (type_t*) malloc(desc->m * sizeof(type_t));
  desc->recvbuf = (type_t*) malloc(desc->t * desc->m * sizeof(type_t));
  desc->em = desc->m;

  if (desc->init_buffers == 1) {
    int i;
    for (i = 0; i < desc->m; i++) {
      ((type_t*) desc->sendbuf)[i] = (type_t) (desc->np + i);
    }
    for (i = 0; i < desc->t * desc->m; i++) {
      ((type_t*) desc->recvbuf)[i] = (type_t) -1;
    }
  }
}

int check_buffers_allgather(nh_exp_desc_t *desc, nh_exp_desc_t *desc2) {
  int same = 0;
  int i;
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for (i=0; i<desc->t*desc->m; i++) {
    if (((type_t*)desc->recvbuf)[i]!=((type_t*)desc2->recvbuf)[i]) {
      fprintf(stderr, "Rank %d %d!=%d at pos %i\n", rank, ((type_t*) desc->recvbuf)[i], ((type_t*) desc2->recvbuf)[i],
          i);
      same = 1;
      break;
    }
  }
  return same;
}

void free_comm_buffers_allgather(nh_exp_desc_t *desc) {
  free(desc->sendbuf);
  free(desc->recvbuf);
}



void setup_comm_buffers_alltoallv(nh_exp_desc_t *desc) {
  int mm;
  //mmm;
  int i, j, z;
//  int *coords;
  nh_exp_v_desc_t *descv;

  if (desc->exp_type == IRREGULAR) {
    descv = (nh_exp_v_desc_t*) (desc);
  } else {
    fprintf(stderr, "function need struct for irregular comm\n");
    exit(1);
  }

  mm = 1;
  for (i = 0; i < descv->desc.d - 1; i++)
    mm *= descv->desc.m;

  descv->sendcount = (int*) malloc(descv->desc.t * sizeof(int));
  descv->recvcount = (int*) malloc(descv->desc.t * sizeof(int));

  descv->senddisp = (int*) malloc(descv->desc.t * sizeof(int));
  descv->recvdisp = (int*) malloc(descv->desc.t * sizeof(int));

//  coords = (int*)calloc(descv->desc.t * descv->desc.d, sizeof(int));
//  TUW_Cart_neighborhood_get_relative(descv->desc.comm, descv->desc.t * descv->desc.d, coords, MPI_UNWEIGHTED);

  mmm = 0;
  for (i = 0; i < descv->desc.t; i++) {
    mm = 1;
    z = descv->desc.d;
    for (j = 0; j < descv->desc.d; j++) {
      if (coords[i * descv->desc.d + j] == 0) {
        mm *= descv->desc.m;
        z--;
      }
    }
    if (z == 0)
      mm = 0; // no data to self

    descv->sendcount[i] = mm;
    descv->recvcount[i] = mm;

    descv->senddisp[i] = mmm;
    descv->recvdisp[i] = mmm;
    mmm += mm;
  }

//  {
//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     if(rank==0)
//       printf("mmm=%d\n",mmm);
//   }

  descv->desc.sendbuf = (type_t*) malloc(mmm * sizeof(type_t));
  descv->desc.recvbuf = (type_t*) malloc(mmm * sizeof(type_t));
  descv->desc.em = mmm;

  if( descv->desc.init_buffers == 1 ) {
    for (i=0; i<mmm; i++) {
      ((type_t*)descv->desc.sendbuf)[i] = (type_t)(descv->desc.np+i);
      ((type_t*)descv->desc.recvbuf)[i] = (type_t)-1;
    }
  }

}


int check_buffers_alltoallv(nh_exp_desc_t *desc, nh_exp_desc_t *desc2) {
  int notsame = 0;
  int i;
  int rank;
  nh_exp_v_desc_t *descv, *descv2;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (desc->exp_type == IRREGULAR) {
    descv = (nh_exp_v_desc_t*) (desc);
  } else {
    fprintf(stderr, "function need struct for irregular comm\n");
    exit(1);
  }

  if (desc2->exp_type == IRREGULAR) {
    descv2 = (nh_exp_v_desc_t*) (desc);
  } else {
    fprintf(stderr, "function need struct for irregular comm\n");
    exit(1);
  }

  if( mmm == 0 ) {
    if( rank == 0 ) {
      fprintf(stderr, "mmm is 0, illegal\n");
      notsame = 1;
    }
  } else {
    for (i = 0; i < mmm; i++) {
      if (((type_t*) descv->desc.recvbuf)[i] != ((type_t*) descv2->desc.recvbuf)[i]) {
        fprintf(stderr, "Rank %d %d!=%d at pos %i\n", rank, ((type_t*) descv->desc.recvbuf)[i],
            ((type_t*) descv2->desc.recvbuf)[i], i);
        notsame = 1;
        break;
      }
    }
  }
  return notsame;
}

void free_comm_buffers_alltoallv(nh_exp_desc_t *desc) {
  nh_exp_v_desc_t *descv;

  if (desc->exp_type == IRREGULAR) {
    descv = (nh_exp_v_desc_t*) (desc);
  } else {
    fprintf(stderr, "function need struct for irregular comm\n");
    exit(1);
  }

  free(descv->sendcount);
  free(descv->recvcount);

  free(descv->senddisp);
  free(descv->recvdisp);

  free(descv->desc.sendbuf);
  free(descv->desc.recvbuf);
}


int init_tuw_alltoall(nh_exp_desc_t *desc) {
  return TUW_Cart_neighbor_alltoall_init(desc->sendbuf, desc->m, TYPE, desc->recvbuf, desc->m, TYPE, desc->comm);
}

int init_tuw_allgather(nh_exp_desc_t *desc) {
  return TUW_Cart_neighbor_allgather_init(desc->sendbuf, desc->m, TYPE, desc->recvbuf, desc->m, TYPE, desc->comm);
}

int init_tuw_alltoallv(nh_exp_desc_t *desc) {
  nh_exp_v_desc_t *descv;

  if (desc->exp_type == IRREGULAR) {
    descv = (nh_exp_v_desc_t*) (desc);
  } else {
    fprintf(stderr, "function need struct for irregular comm\n");
    exit(1);
  }

  return TUW_Cart_neighbor_alltoallv_init(descv->desc.sendbuf, descv->sendcount, descv->senddisp, TYPE,
      descv->desc.recvbuf, descv->recvcount, descv->recvdisp, TYPE, descv->desc.comm);
}

int init_mpi_alltoall(nh_exp_desc_t *desc) {
  // do nothing
  return 0;
}

int init_mpi_allgather(nh_exp_desc_t *desc) {
  // do nothing
  return 0;
}

int init_mpi_alltoallv(nh_exp_desc_t *desc) {
  // do nothing
  return 0;
}

int execute_tuw_alltoall(nh_exp_desc_t *desc) {
  return TUW_Cart_neighbor_alltoall(desc->sendbuf, desc->m, TYPE, desc->recvbuf, desc->m, TYPE, desc->comm);
}

int execute_tuw_allgather(nh_exp_desc_t *desc) {
  return TUW_Cart_neighbor_allgather(desc->sendbuf, desc->m, TYPE, desc->recvbuf, desc->m, TYPE, desc->comm);
}

int execute_tuw_alltoallv(nh_exp_desc_t *desc) {
  nh_exp_v_desc_t *descv;

  if (desc->exp_type == IRREGULAR) {
    descv = (nh_exp_v_desc_t*) (desc);
  } else {
    fprintf(stderr, "function need struct for irregular comm\n");
    exit(1);
  }

  return TUW_Cart_neighbor_alltoallv(descv->desc.sendbuf, descv->sendcount, descv->senddisp, TYPE, descv->desc.recvbuf,
      descv->recvcount, descv->recvdisp, TYPE, descv->desc.comm);
}

int execute_mpi_alltoall(nh_exp_desc_t *desc) {
  return MPI_Neighbor_alltoall(desc->sendbuf, desc->m, TYPE, desc->recvbuf, desc->m, TYPE, desc->comm);
}

int execute_mpi_allgather(nh_exp_desc_t *desc) {
  return MPI_Neighbor_allgather(desc->sendbuf, desc->m, TYPE, desc->recvbuf, desc->m, TYPE, desc->comm);
}

int execute_mpi_alltoallv(nh_exp_desc_t *desc) {
  nh_exp_v_desc_t *descv;

  if (desc->exp_type == IRREGULAR) {
    descv = (nh_exp_v_desc_t*) (desc);
  } else {
    fprintf(stderr, "function need struct for irregular comm\n");
    exit(1);
  }
  return MPI_Neighbor_alltoallv(descv->desc.sendbuf, descv->sendcount, descv->senddisp, TYPE, descv->desc.recvbuf,
      descv->recvcount, descv->recvdisp, TYPE, descv->desc.comm);
}

int execute_mpi_ialltoall(nh_exp_desc_t *desc) {
  MPI_Request req;
  MPI_Ineighbor_alltoall(desc->sendbuf, desc->m, TYPE, desc->recvbuf, desc->m, TYPE, desc->comm, &req);
  return MPI_Wait(&req, MPI_STATUS_IGNORE);
}

int execute_mpi_iallgather(nh_exp_desc_t *desc) {
  MPI_Request req;
  MPI_Ineighbor_allgather(desc->sendbuf, desc->m, TYPE, desc->recvbuf, desc->m, TYPE, desc->comm, &req);
  return MPI_Wait(&req, MPI_STATUS_IGNORE);
}

int execute_mpi_ialltoallv(nh_exp_desc_t *desc) {
  nh_exp_v_desc_t *descv;
  MPI_Request req;

  if (desc->exp_type == IRREGULAR) {
    descv = (nh_exp_v_desc_t*) (desc);
  } else {
    fprintf(stderr, "function need struct for irregular comm\n");
    exit(1);
  }
  MPI_Ineighbor_alltoallv(descv->desc.sendbuf, descv->sendcount, descv->senddisp, TYPE, descv->desc.recvbuf,
      descv->recvcount, descv->recvdisp, TYPE, descv->desc.comm, &req);
  return MPI_Wait(&req, MPI_STATUS_IGNORE);
}



int free_tuw_alltoall_comm(nh_exp_desc_t *desc) {
  TUW_Cart_neighbor_alltoall_free(desc->comm);
  return free_mpi_comm(desc);
}

int free_tuw_allgather_comm(nh_exp_desc_t *desc) {
  TUW_Cart_neighbor_allgather_free(desc->comm);
  return free_mpi_comm(desc);
}

int free_tuw_alltoallv_comm(nh_exp_desc_t *desc) {
  TUW_Cart_neighbor_alltoall_free(desc->comm);
  return free_mpi_comm(desc);
}

int free_mpi_comm(nh_exp_desc_t *desc) {
  free(coords);
  coords = NULL;  // avoid freeing it twice in testing mode
  return MPI_Comm_free(&(desc->comm));
}

void init_create_tuw_communicator(nh_exp_desc_t *desc) {
  int i, j;

  period    = (int*)calloc(desc->d, sizeof(int));
  dimension = (int*)calloc(desc->d, sizeof(int));

  // Create Cartesian communicator
  for (i = 0; i < desc->d; i++) {
    dimension[i] = 0;
    period[i]    = PERIODICITY;
  }
  MPI_Dims_create(desc->np, desc->d, dimension); // Let MPI decide

  coords = (int*) malloc(desc->t * desc->d * sizeof(int));

  int nd = 1;
  for (j = 0; j < desc->d; j++) {
    for (i = 0; i < desc->t; i++) {
      coords[i * desc->d + j] = F + (i / nd) % desc->nb_neighbor;
    }
    nd *= desc->nb_neighbor;
  }

}

void execute_create_tuw_communicator(nh_exp_desc_t *desc) {
  TUW_Cart_neighborhood_create(MPI_COMM_WORLD, desc->d, dimension, period, desc->t, coords, MPI_UNWEIGHTED,
  desc->cartinfo, REORDER, &(desc->comm));
}

void finalize_create_tuw_communicator(nh_exp_desc_t *desc) {
  //free(coords);
  free(period);
  free(dimension);
}

void init_create_mpi_communicator(nh_exp_desc_t *desc) {
  MPI_Comm dummy;

  source = (int*) malloc(desc->t * sizeof(int));
  target = (int*) malloc(desc->t * sizeof(int));

  init_create_tuw_communicator(desc);
  TUW_Cart_neighborhood_create(MPI_COMM_WORLD, desc->d, dimension, period, desc->t, coords, MPI_UNWEIGHTED,
                               MPI_INFO_NULL, REORDER, &dummy);
  finalize_create_tuw_communicator(desc);

  TUW_Cart_neighborhood_get(dummy, desc->t, source, MPI_UNWEIGHTED, desc->t, target, MPI_UNWEIGHTED);
}

void execute_create_mpi_communicator(nh_exp_desc_t *desc) {
  nh_exp_v_desc_t *descv;

  if (desc->exp_type == IRREGULAR) {
    descv = (nh_exp_v_desc_t*) (desc);

    MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD, descv->desc.t, source, MPI_UNWEIGHTED, descv->desc.t, target,
        MPI_UNWEIGHTED,
        MPI_INFO_NULL, REORDER, &(descv->desc.comm));

  } else {
    MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD, desc->t, source, MPI_UNWEIGHTED, desc->t, target, MPI_UNWEIGHTED,
    MPI_INFO_NULL, REORDER, &(desc->comm));
  }

}

void finalize_create_mpi_communicator(nh_exp_desc_t *desc) {
  free(source);
  free(target);
}

void get_comm_func_name(bench_func_t bench_func, char name[]) {
  sprintf(name, "%s", func_names[bench_func]);
}
