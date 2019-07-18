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

/* (c) Jesper Larsson Traff, TUW, July 2018 */

#ifndef TUW_CARTESIANS_LIB_H_
#define TUW_CARTESIANS_LIB_H_

#include <mpi.h>

#define TUW_CART_GRAPH (MPI_GRAPH+MPI_CART+MPI_DIST_GRAPH)

int TUW_Cart_relative_rank(MPI_Comm cartcomm, const int relative[], int *rank);
int TUW_Cart_relative_shift(MPI_Comm cartcomm, const int relative[], 
			    int *inrank, int *outrank);
int TUW_Cart_relative_coord(MPI_Comm cartcomm, const int rank, int relative[]);

int TUW_Cart_neighbors_count(MPI_Comm cartcomm, int *outdegree, int *indegree);
int TUW_Cart_neighbors(MPI_Comm cartcomm, 
		       const int maxindegree,  int source[],
		       const int maxoutdegree, int target[]);

int TUW_Cart_neighborhood_create(MPI_Comm cartcomm,
				 int d, const int dims[], const int periods[],
				 int t, const int targetrelative[],
				 const int weight[],
				 MPI_Info info, int reorder,
				 MPI_Comm *collcomm);

int TUW_Cart_neighborhood_count(MPI_Comm cartcomm, int *t);
int TUW_Cart_neighborhood_get(MPI_Comm cartcomm, 
			      const int maxindegree, 
			      int source[], int sourceweight[],
			      const int maxoutdegree,
			      int target[], int targetweight[]);
int TUW_Cart_neighborhood_get_relative(MPI_Comm cartcomm, const int maxt,
				       int targetrelative[],
				       int targetweight[]);

int TUW_Cart_neighborhood_find_rank(MPI_Comm cartcomm,
				    const int rank,
				    const int maxin, 
				    int sourceindex[], int *inmult,
				    const int maxout,
				    int targetindex[], int *outmult);
int TUW_Cart_neighborhood_find_coord(MPI_Comm cartcomm,
				     const int relative[],
				     const int maxin,
				     int sourceindex[], int *inmult,
				     const int maxout,
				     int targetindex[], int *outmult);

int TUW_Cart_neighbor_allgather(void *sendbuf, 
				int sendcount, MPI_Datatype sendtype,
				void *recvbuf, 
				int recvcount, MPI_Datatype recvtype,
				MPI_Comm cartcomm);
int TUW_Cart_neighbor_allgatherv(void *sendbuf, 
				 int sendcount, MPI_Datatype sendtype,
				 void *recvbuf, 
				 const int recvcount[], 
				 const int recvdisp[], 
				 MPI_Datatype recvtype,
				 MPI_Comm cartcomm);
int TUW_Cart_neighbor_allgatherw(void *sendbuf, 
				 int sendcount, MPI_Datatype sendtype,
				 void *recvbuf, 
				 const int recvcount[],
				 const MPI_Aint recvdisp[], 
				 const MPI_Datatype recvtype[],
				 MPI_Comm cartcomm);
int TUW_Cart_neighbor_alltoall(void *sendbuf, 
			       int sendcount, MPI_Datatype sendtype,
			       void *recvbuf, 
			       int recvcount, MPI_Datatype recvtype,
			       MPI_Comm cartcomm);
int TUW_Cart_neighbor_alltoallv(void *sendbuf, 
				const int sendcount[], const int senddisp[],
				MPI_Datatype sendtype,
				void *recvbuf, 
				const int recvcount[], const int recvdisp[],
				MPI_Datatype recvtype,
				MPI_Comm cartcomm);
int TUW_Cart_neighbor_alltoallw(void *sendbuf, 
				const int sendcount[],
				const MPI_Aint senddisp[], 
				const MPI_Datatype sendtype[],
				void *recvbuf, 
				const int recvcount[],
				const MPI_Aint recvdisp[], 
				const MPI_Datatype recvtype[],
				MPI_Comm cartcomm);

int TUW_Cart_neighbor_allgather_init(void *sendbuf, 
				     int sendcount, MPI_Datatype sendtype,
				     void *recvbuf, 
				     int recvcount, MPI_Datatype recvtype,
				     MPI_Comm cartcomm);
int TUW_Cart_neighbor_allgatherv_init(void *sendbuf, 
				      int sendcount, MPI_Datatype sendtype,
				      void *recvbuf, 
				      int recvcount[], int recvdisp[],
				      MPI_Datatype recvtype,
				      MPI_Comm cartcomm);
int TUW_Cart_neighbor_allgatherw_init(void *sendbuf, 
				      int sendcount, MPI_Datatype sendtype,
				      void *recvbuf, 
				      int recvcount[], int recvdisp[],
				      MPI_Datatype recvtype[],
				      MPI_Comm cartcomm);
int TUW_Cart_neighbor_alltoall_init(void *sendbuf, 
				    int sendcount, MPI_Datatype sendtype,
				    void *recvbuf, 
				    int recvcount, MPI_Datatype recvtype,
				    MPI_Comm cartcomm);
int TUW_Cart_neighbor_alltoallv_init(void *sendbuf, 
				     int sendcount[], int senddisp[],
				     MPI_Datatype sendtype,
				     void *recvbuf, 
				     int recvcount[], int recvdisp[],
				     MPI_Datatype recvtype,
				     MPI_Comm cartcomm);
int TUW_Cart_neighbor_alltoallw_init(void *sendbuf, 
				     int sendcount[], int senddisp[],
				     MPI_Datatype sendtype[],
				     void *recvbuf, 
				     int recvcount[], int recvdisp[],
				     MPI_Datatype recvtype[],
				     MPI_Comm cartcomm);

int TUW_Cart_neighbor_alltoall_free(MPI_Comm cartcomm);
int TUW_Cart_neighbor_allgather_free(MPI_Comm cartcomm);


#endif
