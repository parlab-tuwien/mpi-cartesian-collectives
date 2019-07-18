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

/* Jesper Larsson Traff, TUW, December 2018 */
/* Sascha Hunold, TUW, Jan 2019 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <mpi.h>

#include <assert.h>

#include "tuwcartesians.h"
#include "tuw_bench_setup.h"
#include "tuw_bench_arg_parse.h"


void run_validation(nh_exp_t *exp_tuw, nh_exp_desc_t *desc_tuw, nh_exp_t *exp_mpi, nh_exp_desc_t *desc_mpi) {
  int rank;
  int ok;
  int allok = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  exp_tuw->nh_exp_init_create_communicator(desc_tuw);
  exp_tuw->nh_exp_execute_create_communicator(desc_tuw);
  exp_tuw->nh_exp_finalize_create_communicator(desc_tuw);

  exp_tuw->nh_exp_setup_comm_buffers(desc_tuw);
  exp_tuw->nh_exp_init_comm_func(desc_tuw);
  exp_tuw->nh_exp_exec_comm_func(desc_tuw);

  exp_mpi->nh_exp_init_create_communicator(desc_mpi);
  exp_mpi->nh_exp_execute_create_communicator(desc_mpi);
  exp_mpi->nh_exp_finalize_create_communicator(desc_mpi);

  exp_mpi->nh_exp_setup_comm_buffers(desc_mpi);
  exp_mpi->nh_exp_init_comm_func(desc_mpi);
  exp_mpi->nh_exp_exec_comm_func(desc_mpi);

  ok = exp_tuw->nh_exp_check_comm_buffers(desc_tuw, desc_mpi);


  MPI_Reduce(&ok, &allok, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  if( rank == 0 ) {
    if( allok != 0 ) {
      printf("failed\n");
    } else {
      printf("passed\n");
    }
  }

  exp_tuw->nh_exp_free_comm_buffers(desc_tuw);
  exp_tuw->nh_exp_free_comm(desc_tuw);

  exp_mpi->nh_exp_free_comm_buffers(desc_mpi);
  exp_mpi->nh_exp_free_comm(desc_mpi);

}

int main(int argc, char *argv[]) {

  bench_arg_t bench_params;

  int np;
  int res;

  nh_exp_t exp_tuw;
  nh_exp_t exp_mpi;
  nh_exp_desc_t *desc_tuw = NULL;
  nh_exp_desc_t *desc_mpi = NULL;

  MPI_Init(&argc, &argv);
  res = parse_bench_args(argc, argv, &bench_params);

  if( res != 0 ) {
    MPI_Finalize();
    return -1;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if( bench_params.func == TUW_ALLTOALL || bench_params.func == MPI_ALLTOALL ) {

    exp_tuw.nh_exp_init_create_communicator = &init_create_tuw_communicator;
    exp_tuw.nh_exp_execute_create_communicator = &execute_create_tuw_communicator;
    exp_tuw.nh_exp_finalize_create_communicator = &finalize_create_tuw_communicator;
    exp_tuw.nh_exp_setup_comm_buffers = &setup_comm_buffers_alltoall;
    exp_tuw.nh_exp_free_comm_buffers  = &free_comm_buffers_alltoall;
    exp_tuw.nh_exp_check_comm_buffers = &check_buffers_alltoall;
    exp_tuw.nh_exp_init_comm_func = &init_tuw_alltoall;
    exp_tuw.nh_exp_exec_comm_func = &execute_tuw_alltoall;
    exp_tuw.nh_exp_free_comm = &free_tuw_alltoall_comm;

    desc_tuw = (nh_exp_desc_t*)calloc(1, sizeof(nh_exp_desc_t));

    desc_tuw->func = bench_params.func;
    desc_tuw->d = bench_params.d;
    desc_tuw->m = bench_params.m;
    desc_tuw->nb_neighbor = bench_params.nb_neighbor;
    desc_tuw->nrep = bench_params.nrep;
    desc_tuw->t = get_t(bench_params.d, bench_params.nb_neighbor);
    desc_tuw->exp_type = REGULAR;
    desc_tuw->np = np;
    desc_tuw->init_buffers = 1;
    desc_tuw->trivial  = bench_params.trivial;
    desc_tuw->blocking = bench_params.blocking;

    MPI_Info_create(&(desc_tuw->cartinfo));
    if( bench_params.trivial == 1 ) {
      MPI_Info_set(desc_tuw->cartinfo, "trivial", "false");
    }
    if( bench_params.blocking == 1 ) {
      MPI_Info_set(desc_tuw->cartinfo, "blocking", "false");
    }


    exp_mpi.nh_exp_init_create_communicator = &init_create_mpi_communicator;
    exp_mpi.nh_exp_execute_create_communicator = &execute_create_mpi_communicator;
    exp_mpi.nh_exp_finalize_create_communicator = &finalize_create_mpi_communicator;
    exp_mpi.nh_exp_setup_comm_buffers = &setup_comm_buffers_alltoall;
    exp_mpi.nh_exp_free_comm_buffers  = &free_comm_buffers_alltoall;
    exp_mpi.nh_exp_check_comm_buffers = &check_buffers_alltoall;
    exp_mpi.nh_exp_init_comm_func = &init_mpi_alltoall;
    exp_mpi.nh_exp_exec_comm_func = &execute_mpi_alltoall;
    exp_mpi.nh_exp_free_comm = &free_mpi_comm;

    desc_mpi = (nh_exp_desc_t*)calloc(1, sizeof(nh_exp_desc_t));

    desc_mpi->func = bench_params.func;
    desc_mpi->d = bench_params.d;
    desc_mpi->m = bench_params.m;
    desc_mpi->nb_neighbor = bench_params.nb_neighbor;
    desc_mpi->nrep = bench_params.nrep;
    desc_mpi->t = get_t(bench_params.d, bench_params.nb_neighbor);
    desc_mpi->exp_type = REGULAR;
    desc_mpi->np = np;
    desc_mpi->init_buffers = 1;
    desc_mpi->trivial  = bench_params.trivial;
    desc_mpi->blocking = bench_params.blocking;


//    printf("d=%d,m=%d,nb_neigh=%d,nrep=%d,t=%d", desc->d, desc->m, desc->nb_neighbor, desc->nrep, desc->t);

  } else if( bench_params.func == TUW_ALLGATHER || bench_params.func == MPI_ALLGATHER  ) {

    exp_tuw.nh_exp_init_create_communicator = &init_create_tuw_communicator;
    exp_tuw.nh_exp_execute_create_communicator = &execute_create_tuw_communicator;
    exp_tuw.nh_exp_finalize_create_communicator = &finalize_create_tuw_communicator;
    exp_tuw.nh_exp_setup_comm_buffers = &setup_comm_buffers_allgather;
    exp_tuw.nh_exp_free_comm_buffers  = &free_comm_buffers_allgather;
    exp_tuw.nh_exp_check_comm_buffers = &check_buffers_allgather;
    exp_tuw.nh_exp_init_comm_func = &init_tuw_allgather;
    exp_tuw.nh_exp_exec_comm_func = &execute_tuw_allgather;
    exp_tuw.nh_exp_free_comm = &free_tuw_allgather_comm;

    desc_tuw = (nh_exp_desc_t*)calloc(1, sizeof(nh_exp_desc_t));

    desc_tuw->func = bench_params.func;
    desc_tuw->d = bench_params.d;
    desc_tuw->m = bench_params.m;
    desc_tuw->nb_neighbor = bench_params.nb_neighbor;
    desc_tuw->nrep = bench_params.nrep;
    desc_tuw->t = get_t(bench_params.d, bench_params.nb_neighbor);
    desc_tuw->exp_type = REGULAR;
    desc_tuw->np = np;
    desc_tuw->init_buffers = 1;
    desc_tuw->trivial  = bench_params.trivial;
    desc_tuw->blocking = bench_params.blocking;

    MPI_Info_create(&(desc_tuw->cartinfo));
    if( bench_params.trivial == 1 ) {
      MPI_Info_set(desc_tuw->cartinfo, "trivial", "false");
    }
    if( bench_params.blocking == 1 ) {
      MPI_Info_set(desc_tuw->cartinfo, "blocking", "false");
    }


    exp_mpi.nh_exp_init_create_communicator = &init_create_mpi_communicator;
    exp_mpi.nh_exp_execute_create_communicator = &execute_create_mpi_communicator;
    exp_mpi.nh_exp_finalize_create_communicator = &finalize_create_mpi_communicator;
    exp_mpi.nh_exp_setup_comm_buffers = &setup_comm_buffers_allgather;
    exp_mpi.nh_exp_free_comm_buffers  = &free_comm_buffers_allgather;
    exp_mpi.nh_exp_check_comm_buffers = &check_buffers_allgather;
    exp_mpi.nh_exp_init_comm_func = &init_mpi_allgather;
    exp_mpi.nh_exp_exec_comm_func = &execute_mpi_allgather;
    exp_mpi.nh_exp_free_comm = &free_mpi_comm;

    desc_mpi = (nh_exp_desc_t*)calloc(1, sizeof(nh_exp_desc_t));

    desc_mpi->func = bench_params.func;
    desc_mpi->d = bench_params.d;
    desc_mpi->m = bench_params.m;
    desc_mpi->nb_neighbor = bench_params.nb_neighbor;
    desc_mpi->nrep = bench_params.nrep;
    desc_mpi->t = get_t(bench_params.d, bench_params.nb_neighbor);
    desc_mpi->exp_type = REGULAR;
    desc_mpi->np = np;
    desc_mpi->init_buffers = 1;
    desc_mpi->trivial  = bench_params.trivial;
    desc_mpi->blocking = bench_params.blocking;


  } else if( bench_params.func == TUW_ALLTOALLV || bench_params.func == MPI_ALLTOALLV ) {

    nh_exp_v_desc_t *descv_tuw =  (nh_exp_v_desc_t*)calloc(1, sizeof(nh_exp_v_desc_t));
    nh_exp_v_desc_t *descv_mpi =  (nh_exp_v_desc_t*)calloc(1, sizeof(nh_exp_v_desc_t));

    exp_tuw.nh_exp_init_create_communicator = &init_create_tuw_communicator;
    exp_tuw.nh_exp_execute_create_communicator = &execute_create_tuw_communicator;
    exp_tuw.nh_exp_finalize_create_communicator = &finalize_create_tuw_communicator;
    exp_tuw.nh_exp_setup_comm_buffers = &setup_comm_buffers_alltoallv;
    exp_tuw.nh_exp_free_comm_buffers  = &free_comm_buffers_alltoallv;
    exp_tuw.nh_exp_check_comm_buffers = &check_buffers_alltoallv;
    exp_tuw.nh_exp_init_comm_func = &init_tuw_alltoallv;
    exp_tuw.nh_exp_exec_comm_func = &execute_tuw_alltoallv;
    exp_tuw.nh_exp_free_comm = &free_tuw_alltoallv_comm;

    exp_mpi.nh_exp_init_create_communicator = &init_create_mpi_communicator;
    exp_mpi.nh_exp_execute_create_communicator = &execute_create_mpi_communicator;
    exp_mpi.nh_exp_finalize_create_communicator = &finalize_create_mpi_communicator;
    exp_mpi.nh_exp_setup_comm_buffers = &setup_comm_buffers_alltoallv;
    exp_mpi.nh_exp_free_comm_buffers  = &free_comm_buffers_alltoallv;
    exp_mpi.nh_exp_check_comm_buffers = &check_buffers_alltoallv;
    exp_mpi.nh_exp_init_comm_func = &init_mpi_alltoallv;
    exp_mpi.nh_exp_exec_comm_func = &execute_mpi_alltoallv;
    exp_mpi.nh_exp_free_comm = &free_mpi_comm;

    descv_tuw->desc.func = bench_params.func;
    descv_tuw->desc.d = bench_params.d;
    descv_tuw->desc.m = bench_params.m;
    descv_tuw->desc.nb_neighbor = bench_params.nb_neighbor;
    descv_tuw->desc.nrep = bench_params.nrep;
    descv_tuw->desc.t = get_t(bench_params.d, bench_params.nb_neighbor);
    descv_tuw->desc.exp_type = IRREGULAR;
    descv_tuw->desc.np = np;
    descv_tuw->desc.init_buffers = 1;
    descv_tuw->desc.trivial  = bench_params.trivial;
    descv_tuw->desc.blocking = bench_params.blocking;

    MPI_Info_create(&(descv_tuw->desc.cartinfo));
    if( bench_params.trivial == 1 ) {
      MPI_Info_set(descv_tuw->desc.cartinfo, "trivial", "false");
    }
    if( bench_params.blocking == 1 ) {
      MPI_Info_set(descv_tuw->desc.cartinfo, "blocking", "false");
    }

    descv_mpi->desc.func = bench_params.func;
    descv_mpi->desc.d = bench_params.d;
    descv_mpi->desc.m = bench_params.m;
    descv_mpi->desc.nb_neighbor = bench_params.nb_neighbor;
    descv_mpi->desc.nrep = bench_params.nrep;
    descv_mpi->desc.t = get_t(bench_params.d, bench_params.nb_neighbor);
    descv_mpi->desc.exp_type = IRREGULAR;
    descv_mpi->desc.np = np;
    descv_mpi->desc.init_buffers = 1;
    descv_mpi->desc.trivial  = bench_params.trivial;
    descv_mpi->desc.blocking = bench_params.blocking;

    desc_tuw = (nh_exp_desc_t*)descv_tuw;
    desc_mpi = (nh_exp_desc_t*)descv_mpi;

  } else {
    fprintf(stderr, "unknown function to benchmark... aborting...\n");
    exit(-1);
  }

  run_validation(&exp_tuw, desc_tuw, &exp_mpi, desc_mpi);

  free(desc_tuw);
  free(desc_mpi);

  MPI_Finalize();

}
