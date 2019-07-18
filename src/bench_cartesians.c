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
#include "tuw_bench_runner.h"

int main(int argc, char *argv[]) {

  bench_arg_t bench_params;

  int np;
  int res;

  nh_exp_t exp;
  nh_exp_desc_t *desc = NULL;

  MPI_Init(&argc, &argv);
  res = parse_bench_args(argc, argv, &bench_params);

  if( res != 0 ) {
    MPI_Finalize();
    return -1;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if( bench_params.func == TUW_ALLTOALL ) {

    exp.nh_exp_init_create_communicator = &init_create_tuw_communicator;
    exp.nh_exp_execute_create_communicator = &execute_create_tuw_communicator;
    exp.nh_exp_finalize_create_communicator = &finalize_create_tuw_communicator;
    exp.nh_exp_setup_comm_buffers = &setup_comm_buffers_alltoall;
    exp.nh_exp_free_comm_buffers  = &free_comm_buffers_alltoall;
    exp.nh_exp_init_comm_func = &init_tuw_alltoall;
    exp.nh_exp_exec_comm_func = &execute_tuw_alltoall;
    exp.nh_exp_free_comm = &free_tuw_alltoall_comm;

    desc = (nh_exp_desc_t*)calloc(1, sizeof(nh_exp_desc_t));

    desc->func = bench_params.func;
    desc->d = bench_params.d;
    desc->m = bench_params.m;
    desc->nb_neighbor = bench_params.nb_neighbor;
    desc->nrep = bench_params.nrep;
    desc->t = get_t(bench_params.d, bench_params.nb_neighbor);
    desc->exp_type = REGULAR;
    desc->np = np;
    desc->init_buffers = 0;
    desc->trivial  = bench_params.trivial;
    desc->blocking = bench_params.blocking;

    MPI_Info_create(&(desc->cartinfo));
    if( bench_params.trivial == 1 ) {
      MPI_Info_set(desc->cartinfo, "trivial", "true");
    }
    if( bench_params.blocking == 1 ) {
      MPI_Info_set(desc->cartinfo, "blocking", "true");
    }

//    printf("d=%d,m=%d,nb_neigh=%d,nrep=%d,t=%d", desc->d, desc->m, desc->nb_neighbor, desc->nrep, desc->t);

  } else if( bench_params.func == TUW_ALLGATHER ) {

    exp.nh_exp_init_create_communicator = &init_create_tuw_communicator;
    exp.nh_exp_execute_create_communicator = &execute_create_tuw_communicator;
    exp.nh_exp_finalize_create_communicator = &finalize_create_tuw_communicator;
    exp.nh_exp_setup_comm_buffers = &setup_comm_buffers_allgather;
    exp.nh_exp_free_comm_buffers  = &free_comm_buffers_allgather;
    exp.nh_exp_init_comm_func = &init_tuw_allgather;
    exp.nh_exp_exec_comm_func = &execute_tuw_allgather;
    exp.nh_exp_free_comm = &free_tuw_allgather_comm;

    desc = (nh_exp_desc_t*)calloc(1, sizeof(nh_exp_desc_t));

    desc->func = bench_params.func;
    desc->d = bench_params.d;
    desc->m = bench_params.m;
    desc->nb_neighbor = bench_params.nb_neighbor;
    desc->nrep = bench_params.nrep;
    desc->t = get_t(bench_params.d, bench_params.nb_neighbor);
    desc->exp_type = REGULAR;
    desc->np = np;
    desc->init_buffers = 0;
    desc->trivial  = bench_params.trivial;
    desc->blocking = bench_params.blocking;

    MPI_Info_create(&(desc->cartinfo));
    if( bench_params.trivial == 1 ) {
      MPI_Info_set(desc->cartinfo, "trivial", "true");
    }
    if( bench_params.blocking == 1 ) {
      MPI_Info_set(desc->cartinfo, "blocking", "true");
    }

  } else if( bench_params.func == TUW_ALLTOALLV ) {

    nh_exp_v_desc_t *descv =  (nh_exp_v_desc_t*)calloc(1, sizeof(nh_exp_v_desc_t));

    exp.nh_exp_init_create_communicator = &init_create_tuw_communicator;
    exp.nh_exp_execute_create_communicator = &execute_create_tuw_communicator;
    exp.nh_exp_finalize_create_communicator = &finalize_create_tuw_communicator;
    exp.nh_exp_setup_comm_buffers = &setup_comm_buffers_alltoallv;
    exp.nh_exp_free_comm_buffers  = &free_comm_buffers_alltoallv;
    exp.nh_exp_init_comm_func = &init_tuw_alltoallv;
    exp.nh_exp_exec_comm_func = &execute_tuw_alltoallv;
    exp.nh_exp_free_comm = &free_tuw_alltoallv_comm;

    descv->desc.func = bench_params.func;
    descv->desc.d = bench_params.d;
    descv->desc.m = bench_params.m;
    descv->desc.nb_neighbor = bench_params.nb_neighbor;
    descv->desc.nrep = bench_params.nrep;
    descv->desc.t = get_t(bench_params.d, bench_params.nb_neighbor);
    descv->desc.exp_type = IRREGULAR;
    descv->desc.np = np;
    descv->desc.init_buffers = 0;
    descv->desc.trivial = bench_params.trivial;
    descv->desc.blocking = bench_params.blocking;

    MPI_Info_create(&(descv->desc.cartinfo));
    if( bench_params.trivial == 1 ) {
      MPI_Info_set(descv->desc.cartinfo, "trivial", "true");
    }
    if( bench_params.blocking == 1 ) {
      MPI_Info_set(descv->desc.cartinfo, "blocking", "true");
    }

    desc = (nh_exp_desc_t*)descv;

  } else if( bench_params.func == MPI_ALLTOALL || bench_params.func == MPI_IALLTOALL ) {

    exp.nh_exp_init_create_communicator = &init_create_mpi_communicator;
    exp.nh_exp_execute_create_communicator = &execute_create_mpi_communicator;
    exp.nh_exp_finalize_create_communicator = &finalize_create_mpi_communicator;
    exp.nh_exp_setup_comm_buffers = &setup_comm_buffers_alltoall;
    exp.nh_exp_free_comm_buffers  = &free_comm_buffers_alltoall;
    exp.nh_exp_init_comm_func = &init_mpi_alltoall;
    exp.nh_exp_free_comm = &free_mpi_comm;

    if( bench_params.func == MPI_ALLTOALL ) {
      exp.nh_exp_exec_comm_func = &execute_mpi_alltoall;
    } else {
      exp.nh_exp_exec_comm_func = &execute_mpi_ialltoall;
    }


    desc = (nh_exp_desc_t*)calloc(1, sizeof(nh_exp_desc_t));

    desc->func = bench_params.func;
    desc->d = bench_params.d;
    desc->m = bench_params.m;
    desc->nb_neighbor = bench_params.nb_neighbor;
    desc->nrep = bench_params.nrep;
    desc->t = get_t(bench_params.d, bench_params.nb_neighbor);
    desc->exp_type = REGULAR;
    desc->np = np;
    desc->init_buffers = 0;
    desc->trivial = bench_params.trivial;
    desc->blocking = bench_params.blocking;


  } else if( bench_params.func == MPI_ALLGATHER || bench_params.func == MPI_IALLGATHER  ) {

    exp.nh_exp_init_create_communicator = &init_create_mpi_communicator;
    exp.nh_exp_execute_create_communicator = &execute_create_mpi_communicator;
    exp.nh_exp_finalize_create_communicator = &finalize_create_mpi_communicator;
    exp.nh_exp_setup_comm_buffers = &setup_comm_buffers_allgather;
    exp.nh_exp_free_comm_buffers  = &free_comm_buffers_allgather;
    exp.nh_exp_init_comm_func = &init_mpi_allgather;
    exp.nh_exp_free_comm = &free_mpi_comm;

    if( bench_params.func == MPI_ALLGATHER ) {
      exp.nh_exp_exec_comm_func = &execute_mpi_allgather;
    } else {
      exp.nh_exp_exec_comm_func = &execute_mpi_iallgather;
    }


    desc = (nh_exp_desc_t*)calloc(1, sizeof(nh_exp_desc_t));

    desc->func = bench_params.func;
    desc->d = bench_params.d;
    desc->m = bench_params.m;
    desc->nb_neighbor = bench_params.nb_neighbor;
    desc->nrep = bench_params.nrep;
    desc->t = get_t(bench_params.d, bench_params.nb_neighbor);
    desc->exp_type = REGULAR;
    desc->np = np;
    desc->init_buffers = 0;
    desc->trivial  = bench_params.trivial;
    desc->blocking = bench_params.blocking;


  } else if( bench_params.func == MPI_ALLTOALLV || bench_params.func == MPI_IALLTOALLV  ) {

    nh_exp_v_desc_t *descv =  (nh_exp_v_desc_t*)calloc(1, sizeof(nh_exp_v_desc_t));

    exp.nh_exp_init_create_communicator = &init_create_mpi_communicator;
    exp.nh_exp_execute_create_communicator = &execute_create_mpi_communicator;
    exp.nh_exp_finalize_create_communicator = &finalize_create_mpi_communicator;
    exp.nh_exp_setup_comm_buffers = &setup_comm_buffers_alltoallv;
    exp.nh_exp_free_comm_buffers  = &free_comm_buffers_alltoallv;
    exp.nh_exp_init_comm_func = &init_mpi_alltoallv;
    exp.nh_exp_exec_comm_func = &execute_mpi_alltoallv;
    exp.nh_exp_free_comm = &free_mpi_comm;


    if( bench_params.func == MPI_ALLGATHER ) {
      exp.nh_exp_exec_comm_func = &execute_mpi_alltoallv;
    } else {
      exp.nh_exp_exec_comm_func = &execute_mpi_ialltoallv;
    }

    descv->desc.func = bench_params.func;
    descv->desc.d = bench_params.d;
    descv->desc.m = bench_params.m;
    descv->desc.nb_neighbor = bench_params.nb_neighbor;
    descv->desc.nrep = bench_params.nrep;
    descv->desc.t = get_t(bench_params.d, bench_params.nb_neighbor);
    descv->desc.exp_type = IRREGULAR;
    descv->desc.np = np;
    descv->desc.init_buffers = 0;
    descv->desc.trivial = bench_params.trivial;
    descv->desc.blocking = bench_params.blocking;

    desc = (nh_exp_desc_t*)descv;

  } else {
    fprintf(stderr, "unknown function to benchmark... aborting...\n");
    exit(-1);
  }

  run_exp(&exp, desc);

  free(desc);

  MPI_Finalize();

}
