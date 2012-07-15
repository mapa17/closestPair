closestPair - Sequencial and Distributed closest Pair Algorithm
===============================================================
15.7.2012

Both programs will generate randomly a number of points in a two dimensional plane,
will find the closest pair of points, and will return the points as well as the distance
and the time it took to find them.

This can be used to compare the execution time of the sequencial algorithm compared
to the distributed algorithm on different problem sized and different number of
computing nodes.

A description of the distributed algorithm implemented can be found at
 http://en.wikipedia.org/wiki/Closest_pair_of_points_problem

OpenMPI is used as a framework to implemented the distributed algorithm.

Dependencies
------------

OpenMPI tools ( Debian/Ubuntu use sudo apt-get install openmpi-bin openmpi-doc libopenmpi-dev )

Compilation
------------

* MPI Version

mpicc closestPair_mpi.c closestPair_tools.c -o cp_mpi

* Sequencial Version

mpicc closestPair_sequencial.c closestPair_tools.c -o cp_seq


Execution
---------

* Run MPI version on two nodes with 100 Points

mpirun -n 2 cp_mpi 100


Blame
-----

Author: Pasieka Manuel , mapa17@posgrado.upv.es
