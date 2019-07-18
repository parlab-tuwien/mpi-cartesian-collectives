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

#ifndef SRC_TUW_BENCH_SETUP_H_
#define SRC_TUW_BENCH_SETUP_H_

#include <mpi.h>

#define TYPE MPI_INT
typedef int type_t;

extern const int D;
extern const int R;

typedef enum desc_type {REGULAR, IRREGULAR} exp_type_t;

typedef enum func_type {
  TUW_ALLTOALL = 0,
  TUW_ALLGATHER,
  TUW_ALLTOALLV,
  MPI_ALLTOALL,
  MPI_ALLGATHER,
  MPI_ALLTOALLV,
  MPI_IALLTOALL,
  MPI_IALLGATHER,
  MPI_IALLTOALLV,
  FUNCTIONS_END
} bench_func_t;

extern const char *func_names[];

typedef struct nh_exp_desc {
  bench_func_t func;
  int d;
  int t;
  int m;
  int nb_neighbor;
  int em; /* effective buf size */

  int trivial;
  int blocking;

  MPI_Comm comm;
  MPI_Info cartinfo;

  int np;
  type_t *sendbuf;
  type_t *recvbuf;
  int init_buffers; /* boolean, for testing purposes */

  int nrep;
  exp_type_t exp_type;

} nh_exp_desc_t;

typedef struct nh_exp_v_desc {
  nh_exp_desc_t desc;

  int *sendcount;
  int *senddisp;
  int *recvcount;
  int *recvdisp;

  //exp_type_t exp_type;
} nh_exp_v_desc_t;


//typedef struct nh_meas {
//  int nrep;
//} nh_meas_t;

int get_t(int d, int nb_neighbor);

void setup_comm_buffers_alltoall(nh_exp_desc_t *desc);
void setup_comm_buffers_allgather(nh_exp_desc_t *desc);
void setup_comm_buffers_alltoallv(nh_exp_desc_t *desc);

int check_buffers_alltoall(nh_exp_desc_t *desc, nh_exp_desc_t *desc2);
int check_buffers_allgather(nh_exp_desc_t *desc, nh_exp_desc_t *desc2);
int check_buffers_alltoallv(nh_exp_desc_t *desc, nh_exp_desc_t *desc2);

void free_comm_buffers_alltoall(nh_exp_desc_t *desc);
void free_comm_buffers_allgather(nh_exp_desc_t *desc);
void free_comm_buffers_alltoallv(nh_exp_desc_t *desc);

void init_create_tuw_communicator(nh_exp_desc_t *desc);
void execute_create_tuw_communicator(nh_exp_desc_t *desc);
void finalize_create_tuw_communicator(nh_exp_desc_t *desc);

void init_create_mpi_communicator(nh_exp_desc_t *desc);
void execute_create_mpi_communicator(nh_exp_desc_t *desc);
void finalize_create_mpi_communicator(nh_exp_desc_t *desc);

int init_mpi_alltoall(nh_exp_desc_t *desc);
int init_mpi_allgather(nh_exp_desc_t *desc);
int init_mpi_alltoallv(nh_exp_desc_t *desc);

int init_tuw_alltoall(nh_exp_desc_t *desc);
int init_tuw_allgather(nh_exp_desc_t *desc);
int init_tuw_alltoallv(nh_exp_desc_t *desc);

int execute_tuw_alltoall(nh_exp_desc_t *desc);
int execute_tuw_allgather(nh_exp_desc_t *desc);
int execute_tuw_alltoallv(nh_exp_desc_t *desc);

int execute_mpi_alltoall(nh_exp_desc_t *desc);
int execute_mpi_allgather(nh_exp_desc_t *desc);
int execute_mpi_alltoallv(nh_exp_desc_t *desc);

int execute_mpi_ialltoall(nh_exp_desc_t *desc);
int execute_mpi_iallgather(nh_exp_desc_t *desc);
int execute_mpi_ialltoallv(nh_exp_desc_t *desc);

int free_tuw_alltoall_comm(nh_exp_desc_t *desc);
int free_tuw_allgather_comm(nh_exp_desc_t *desc);
int free_tuw_alltoallv_comm(nh_exp_desc_t *desc);

int free_mpi_comm(nh_exp_desc_t *desc);

typedef struct nh_exp {
  void (*nh_exp_init_create_communicator)(nh_exp_desc_t *desc);
  void (*nh_exp_execute_create_communicator)(nh_exp_desc_t *desc);
  void (*nh_exp_finalize_create_communicator)(nh_exp_desc_t *desc);
  void (*nh_exp_setup_comm_buffers)(nh_exp_desc_t *desc);
  int (*nh_exp_check_comm_buffers)(nh_exp_desc_t *desc, nh_exp_desc_t *desc2);
  void (*nh_exp_free_comm_buffers)(nh_exp_desc_t *desc);
  int (*nh_exp_init_comm_func)(nh_exp_desc_t *desc);
  int (*nh_exp_exec_comm_func)(nh_exp_desc_t *desc);
  int (*nh_exp_free_comm)(nh_exp_desc_t *desc);
} nh_exp_t;


void get_comm_func_name(bench_func_t bench_func, char name[]);


#endif /* SRC_TUW_BENCH_SETUP_H_ */
