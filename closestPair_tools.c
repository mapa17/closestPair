/*
 * nearestPoint_tools.c
 *
 *  Created on: Jul 13, 2012
 *      Author: Pasieka Manuel , mapa17@posgrado.upv.es
 */


#include "closestPair_tools.h"

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <time.h>
#include <math.h>

//Global Variables
long np;
int mpi_size, mpi_id;
struct timeval start;
MPI_Group mpi_group;
MPI_Comm mpi_comm;
MPI_Datatype PointType;
MPI_Datatype DistanceType;


void tick(void){
	gettimeofday(&start, NULL);
}

double tack(void){
	struct timeval end;
	gettimeofday(&end, NULL);
	return (double) end.tv_sec-start.tv_sec + ( end.tv_usec - start.tv_usec ) / 1000000.0;
}

inline double dist(point a, point b){
	return ( (double) (a.x - b.x)*(a.x - b.x) + (double)(a.y - b.y)*(a.y - b.y) );
}

inline char eqal(point a, point b){
	if( (a.x == b.x) && (a.y == b.y) )
		return 1;
	else
		return 0;
}

inline int cmpf(float a, float b){
        return a < b ? -1 : a > b ? 1 : 0;
}

int cmpPoint_x(const void *a, const void *b) {
        return cmpf( ((point*)a)->x, ((point*)b)->x );
}

int cmpPoint_y(const void *a, const void *b) {
        return cmpf( ((point *)a)->y, ((point *)b)->y );
}


void mpiInit(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	//Create PointType
	MPI_Datatype type[2] = { MPI_FLOAT, MPI_FLOAT };
	int blocklen[2] = { 1, 1 };
	MPI_Aint disp[2];
	point Point;
	disp[0] = (void*)&Point.x - (void*)&Point;
	disp[1] = (void*)&Point.y - (void*)&Point;
	MPI_Type_create_struct(2, blocklen, disp, type, &PointType);
    MPI_Type_commit(&PointType);

    //Create DistanceType
    MPI_Datatype type2[3] = { MPI_DOUBLE, PointType, PointType };
    int blocklen2[3] = { 1, 1, 1 };
    MPI_Aint disp2[3];
    distance Distance;
    disp[0] = (void*)&Distance.d - (void*)&Distance;
    disp[1] = (void*)&Distance.a - (void*)&Distance;
    disp[2] = (void*)&Distance.b - (void*)&Distance;
    MPI_Type_create_struct(3, blocklen2, disp2, type2, &DistanceType);
    MPI_Type_commit(&DistanceType);
}

//Set global mpi_comm and mpi_group
//Check if the number of nodes is to the power of 2 ( 2, 4, 8, 16, ...)
//and if not create a new group that contains all nodes up to the closest lowest valid number of nodes.
//Return 1 if the node is part of the active comm world
//Return 0 if not
char prepareMPIComm(void)
{
	float tf;
	int ti;
	int idOld;
	char activeNode;

	//Get Rank and GroupSize
	MPI_Comm_rank (MPI_COMM_WORLD, &idOld);        /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);        /* get number of processes */

	if(mpi_size < 2){
		printf("Error! At least two nodes are needed!\n");
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}

	//Test for multiples of 2
	tf = log(mpi_size)/log(2);
	ti = tf;
	tf = tf - (float)ti;
	if(tf > 0)
	{
		if(idOld == 0){
			printf("Cant use all possible nodes! Will only use the first [%d] of [%d].\n",1<<ti, mpi_size); fflush(stdout);
		}
		mpi_size = 1 << ti;

		if(idOld < mpi_size)
		{
		   activeNode = 1;

		   // Split comm into two group
		   MPI_Comm_split(MPI_COMM_WORLD, 0, idOld, &mpi_comm);
		   MPI_Comm_group( mpi_comm, &mpi_group);

			//printf("[%d] Creating new group with [%d] members.\n",idOld, mpi_size);
			//printf("Created new group of size [%d] with id [%d] -> [%d]\n", mpi_size, idOld, mpi_id); fflush(stdout);

		} else {
			activeNode = 0;
			MPI_Comm_split(MPI_COMM_WORLD, 1, idOld - mpi_size, &mpi_comm);
			MPI_Comm_group( mpi_comm, &mpi_group);

			//printf("Node [%d] wont participate.\n", idOld);
		}
	} else {
		activeNode = 1;
		//If the number of nodes is alright, use the standard world comm
		MPI_Comm_group( MPI_COMM_WORLD, &mpi_group);
		MPI_Comm_create( MPI_COMM_WORLD, mpi_group, &mpi_comm);
	}

	MPI_Comm_size( mpi_comm, &mpi_size );
	MPI_Comm_rank( mpi_comm, &mpi_id );
	return activeNode;
}
