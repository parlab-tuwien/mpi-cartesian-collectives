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
#include <time.h>

#include "tuw_bench_runner.h"

void run_exp(nh_exp_t *exp, nh_exp_desc_t *desc) {
  int i;
  double start, stop;
  double *times;
  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  times = (double*) calloc(desc->nrep, sizeof(double));

  exp->nh_exp_init_create_communicator(desc);
  exp->nh_exp_execute_create_communicator(desc);
  exp->nh_exp_finalize_create_communicator(desc);

  exp->nh_exp_setup_comm_buffers(desc);
  exp->nh_exp_init_comm_func(desc);

  for (i = 0; i < desc->nrep; i++) {
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    exp->nh_exp_exec_comm_func(desc);
    stop = MPI_Wtime();
    times[i] = stop - start;
  }

  exp->nh_exp_free_comm_buffers(desc);

  if (rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, times, desc->nrep, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(times, times, desc->nrep, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  }

  if (rank == 0) {
    char name_buf[50];
    time_t t;
    time(&t);
    printf(
        "#############################################\n"
        "# current time %s\n"
        "#############################################\n"
        "# func - communcation function\n"
        "# nrep - number of repetitions\n"
        "#   np - communicator size\n"
        "#    d - dimension\n"
        "#    m - block size\n"
        "#   em - effective block size\n"
        "#    N - number of neighbors\n"
        "# triv - trivial implementation (yes/no)\n"
        "# block- blocking implementation (yes/no)\n"
        "# time - run-time (micros)\n"
        "#############################################\n",
        ctime(&t));
    printf("%30s;nrep;  np; d;     m;       em;N;triv;block;%16s\n", "func","time");

    get_comm_func_name(desc->func, name_buf);

    for (i = 0; i < desc->nrep; i++) {
      printf("%30s;%4d;%4d;%2d;%6d;%9d;%1d;%4d;%5d;%16.3f\n", name_buf, i, desc->np, desc->d, desc->m, desc->em,
          desc->nb_neighbor, desc->trivial, desc->blocking, times[i] * 1e6);
    }
  }

  free(times);

  exp->nh_exp_free_comm(desc);

}
