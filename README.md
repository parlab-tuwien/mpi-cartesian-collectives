# Cartesian Collective Communication

## Citing article

**Jesper Larsson Träff** and **Sascha Hunold**. 2019. *Cartesian Collective Communication*. In 48th International Conference on Parallel Processing (ICPP 2019), August 5–8, 2019, Kyoto, Japan. ACM, New York, NY, USA, 11 pages. https://doi.org/10.1145/3337821.3337848

## Build benchmark programs

```
cmake ./
make
```

## Running on local machine

### Message-combining Alltoall 

```
mpirun --oversubscribe -np 4 ./bin/bench_cartesians -f TUW_Cart_neighbor_alltoall -d 3 -m 100 -n 5 -r 10
#############################################
# current time Thu Jul 18 16:39:29 2019

#############################################
# func - communcation function
# nrep - number of repetitions
#   np - communicator size
#    d - dimension
#    m - block size
#   em - effective block size
#    N - number of neighbors
# triv - trivial implementation (yes/no)
# block- blocking implementation (yes/no)
# time - run-time (micros)
#############################################
                          func;nrep;  np; d;     m;       em;N;triv;block;            time
    TUW_Cart_neighbor_alltoall;   0;   4; 3;   100;    12500;5;   0;    0;         459.000
    TUW_Cart_neighbor_alltoall;   1;   4; 3;   100;    12500;5;   0;    0;          73.000
    TUW_Cart_neighbor_alltoall;   2;   4; 3;   100;    12500;5;   0;    0;         138.000
    TUW_Cart_neighbor_alltoall;   3;   4; 3;   100;    12500;5;   0;    0;          34.000
    TUW_Cart_neighbor_alltoall;   4;   4; 3;   100;    12500;5;   0;    0;          33.000
    TUW_Cart_neighbor_alltoall;   5;   4; 3;   100;    12500;5;   0;    0;          33.000
    TUW_Cart_neighbor_alltoall;   6;   4; 3;   100;    12500;5;   0;    0;          33.000
    TUW_Cart_neighbor_alltoall;   7;   4; 3;   100;    12500;5;   0;    0;          33.000
    TUW_Cart_neighbor_alltoall;   8;   4; 3;   100;    12500;5;   0;    0;          34.000
    TUW_Cart_neighbor_alltoall;   9;   4; 3;   100;    12500;5;   0;    0;          33.000
```

### Message-combining Alltoall (trivial -t)
```
mpirun --oversubscribe -np 4 ./bin/bench_cartesians -f TUW_Cart_neighbor_alltoall -d 3 -m 100 -n 5 -r 10 -t 
```
### Message-combining Alltoall (use blocking calls -b)
```
mpirun --oversubscribe -np 4 ./bin/bench_cartesians -f TUW_Cart_neighbor_alltoall -d 3 -m 100 -n 5 -r 10 -b
```


## Check message-combining collectives vs. default MPI collectives

```
mpirun --oversubscribe -np 4 ./bin/check_cartesians -f TUW_Cart_neighbor_alltoall -d 3 -m 100 -n 5 -r 10
passed
```
