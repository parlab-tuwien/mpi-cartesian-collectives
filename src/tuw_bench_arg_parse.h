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

#ifndef SRC_TUW_BENCH_ARG_PARSE_H_
#define SRC_TUW_BENCH_ARG_PARSE_H_

#include "tuw_bench_setup.h"

typedef struct args {
  bench_func_t func;

  int nrep;
  int nb_neighbor;
  int d;
  int m;
  int trivial;  /* trivial = 1 -> no message combining */
  int blocking; /* use blocking send's/recv's internally */

} bench_arg_t;

int parse_bench_args(int argc, char *argv[], bench_arg_t *bench_params);

#endif /* SRC_TUW_BENCH_ARG_PARSE_H_ */
