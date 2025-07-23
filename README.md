**📢 Repository Moved (2025-07-23)**

This repository has been moved to: https://github.com/comphy-lab/Asymmetries-in-coalescence

Please use the new location for all future development and contributions.

---

# Asymmetries-in-coalescence
Asymmetries in coalescence: size asymmetry. Still axially symmetric. 


## running the codes

There are two ways to run the codes:

1. Using the vanilla basilisk method:

```shell
qcc -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm 
./coalescenceBubble
```

2. Using the makefile (can be interactively run using bview browser):

```shell
CFLAGS=-DDISPLAY=-1 make coalescenceBubble.tst
```
Check the localhost on coalescenceBubble/display.html. something like: [http://basilisk.fr/three.js/editor/index.html?ws://localhost:7100](http://basilisk.fr/three.js/editor/index.html?ws://localhost:7100) and run interactively.

### To run using openMP, please use the flag -fopenmp

```shell
qcc -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm -fopenmp
export OMP_NUM_THREADS=8
./coalescenceBubble
```

**Note:** The code will not directly work with openmpi (with -D_MPI flag and mpirun). To do that, please follow the procedure we use: 

1. Run the following for a few timesteps (stop using tmax=1e-2 or so)
```shell
qcc -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm -fopenmp
export OMP_NUM_THREADS=8
./coalescenceBubble
```

This will generate a "dump" file

2. Do not delete the dump file. Now you can use mpirun, like:

```shell
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm
mpirun -np 8 coalescenceBubble
```
