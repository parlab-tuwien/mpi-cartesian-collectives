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
#include <getopt.h>
#include <string.h>

#include "tuw_bench_arg_parse.h"

int parse_bench_args(int argc, char *argv[], bench_arg_t *bench_params) {
  int ret = 0;
  int c;

  int function_found = 0;

  bench_params->d = -1;
  bench_params->m = -1;
  bench_params->nb_neighbor = -1;
  bench_params->nrep = -1;
  bench_params->trivial  = 0;
  bench_params->blocking = 0;

  while (1) {
    static struct option long_options[] =
      {
        {"function",required_argument, 0, 'f'},
        {"nrep",    required_argument, 0, 'r'},
        {"nneigh",  required_argument, 0, 'n'},
        {"dim",     required_argument, 0, 'd'},
        {"msize",   required_argument, 0, 'm'},
        {"trivial", optional_argument, 0, 't'},
        {"blocking", optional_argument, 0, 'b'},
        {0, 0, 0, 0}
      };

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, "f:r:n:d:m:tb", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 'f':
    {
      int i;
      for(i=0; i<FUNCTIONS_END; i++) {
        if( strcmp(optarg, func_names[i]) == 0 ) {
          function_found = 1;
          bench_params->func = i;
          break;
        }
      }
    }
    break;

    case 'r':
      bench_params->nrep = atoi(optarg);
      break;

    case 'n':
      bench_params->nb_neighbor = atoi(optarg);
      break;

    case 'd':
      bench_params->d = atoi(optarg);
      break;

    case 'm':
      bench_params->m = atoi(optarg);
      break;

    case 't':
      bench_params->trivial = 1;
      break;

    case 'b':
      bench_params->blocking = 1;
      break;

    case '?':
      /* getopt_long already printed an error message. */
      break;

    default:
      abort();
    }
  }

  if( function_found == 0 ) {
    int i;
    fprintf(stderr, "function unknown; available are: ");
    for(i=0; i<FUNCTIONS_END; i++) {
      fprintf(stderr, "%s", func_names[i]);
      if(i < FUNCTIONS_END-1) {
        fprintf(stderr, ", ");
      }
    }
    fprintf(stderr, "\n");
    return -1;
  }

  if (bench_params->d == -1) {
    fprintf(stderr, "no dimension (-d), use default\n");
    bench_params->d = 2;
  }
  if (bench_params->nb_neighbor == -1) {
    fprintf(stderr, "no nb neighbors (-n), use default\n");
    bench_params->nb_neighbor = 4;
  }
  if (bench_params->m == -1) {
    fprintf(stderr, "no block size (-m), use default\n");
    bench_params->m = 1;
  }
  if (bench_params->nrep == -1) {
    fprintf(stderr, "no nrep (-r), use default\n");
    bench_params->nrep = 1;
  }

  return ret;
}
