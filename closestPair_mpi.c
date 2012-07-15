/*
 * nearestPoints.c
 *
 *  Created on: Jul 5, 2012
 *      Author: Pasieka Manuel , mapa17@posgrado.upv.es
 *
 *      Good OpenMPI Reference: http://mpi.deino.net/mpi_functions/index.htm
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <time.h>
#include "mpi.h"
#include <math.h>
#include <string.h>

#include "closestPair_tools.h"

#define DEBUG_HIGH

void inputCheck(int argc, char** argv);
void mpiFinish(void);
//void mpiInit(int argc, char** argv);
void printUsage(int argc, char** argv);
char generatePoints(long np, point** points);
void trivalSearch(long np, point* points, distance* solution);
void multiSearch(long nTotalPoints, point* points, distance* solution);
void joinPointBlocks(long np1, long np2, point* block1, point* block2, distance* bestsolution);

int main(int argc, char** argv)
{
	point* points;
	distance solution;
	double elapsedTime;

	points = NULL;

	mpiInit(argc, argv);

	if( prepareMPIComm() )
	{
		inputCheck(argc, argv);

		if(mpi_id == 0)
		{
			printf("Generating [%d] points\n", np);

			if( generatePoints(np, &points) != EXIT_SUCCESS ){
				printf("Generating Points failed!\n");
				exit(EXIT_FAILURE);
			}

			printf("Starting search ...");
		}

		tick();
			multiSearch(np, points, &solution);
		elapsedTime = tack();

		printf("Found Solution a[%f,%f] , b[%f,%f] distance [%0.10f]\n", solution.a.x, solution.a.y, solution.b.x, solution.b.y, solution.d);
		printf("Operation took %f seconds \n", elapsedTime);

		if(mpi_id == 0){
			free(points);
		}

	} else {

		}

	mpiFinish();
	exit(EXIT_SUCCESS);
}


void mpiFinish(void){
	MPI_Finalize();
}


void inputCheck(int argc, char** argv)
{
	if ( argc < 2){
		printUsage(argc, argv);
		exit(EXIT_FAILURE);
	}

	errno = 0; // To distinguish success/failure after call
	np = atol(argv[1]);

	/* Check for various possible errors */

	if ((errno == ERANGE && (np == LONG_MAX || np == LONG_MIN)) || (errno != 0 && np == 0)) {
	   perror("atol");
	   exit(EXIT_FAILURE);
	}

	if(np <= 0){
		printf("Invalid number of points! Number of points must be bigger than zero.\n");
		exit(EXIT_FAILURE);
	}
}

void printUsage(int argc, char** argv){
	printf("Usage: %s NUMER_OF_POINTS\n",argv[0]);
}


//Get NP unique points
char generatePoints(long np, point** points)
{
	time_t t1;
	struct tm* t2 = NULL;
	int sec;
	long i,j;
	long equals;

	*points = malloc(sizeof(points) * np);
	if(*points == NULL){
		perror("malloc");
		return(EXIT_FAILURE);
	}

	t1 = time( NULL );
	t2 = gmtime(&t1);
	sec = t2->tm_sec;
	equals = 0;
	//printf("Seed [%d]\n", sec);

	 srand(sec);
	 for(i = 0; i < np; i++){
		 (*points)[i].x = (float) rand() / RAND_MAX;
		 (*points)[i].y = (float) rand() / RAND_MAX;

		 //Check if this point is unique
		 for(j = 0; j < i; j++){
			 if(eqal((*points)[i], (*points)[j])){
				 //printf("\ne[%d == %d]\n",i, j);
				 i--; //Do not use this point
				 equals++;
			 	 break;
			 }
		 }
	 }

	 printf("Eliminated %ld equal points.\n", equals);

#ifdef DEBUG_HIGH
	 printf("--- Generated Points ---\n");
	 for(i=0; i < np; i++){
		 printf("[%g, %g]\n", (*points)[i].x, (*points)[i].y );
	 }
	 printf("--- END ---\n");
#endif

	 return(EXIT_SUCCESS);
}


//Follow solution on
// http://en.wikipedia.org/wiki/Closest_pair_of_points_problem
void joinPointBlocks(long np1, long np2, point* block1, point* block2, distance* bestsolution)
{

	//Do nothing yet

}


//Follow solution on
// http://en.wikipedia.org/wiki/Closest_pair_of_points_problem
void multiSearch(long nTotalPoints, point* points, distance* solution)
{
	char b[2000];
	int i;
	long np;
	int height, depth;
	int parent;
	point* localPoints;
	distance localSolution;
	distance remoteSolution;
	int zeroOffset;
	long nlocalPoints;
	long localPointOffset;



	np = nTotalPoints / mpi_size;
	zeroOffset = nTotalPoints % mpi_size;
	localPointOffset = mpi_id * np + zeroOffset;

	if(mpi_id == 0)
	{
		//Sort points by their x value
		qsort(points, nTotalPoints, sizeof(point), cmpPoint_x);

		printf("Distributing Data ...\n");

		localPoints = points;
		MPI_Bcast( &localPoints[0], nTotalPoints, PointType, 0, mpi_comm);

		//MPI_Scatter(&points[zeroOffset], np, PointType, &localPoints[localPointOffset] , np, PointType, 0, mpi_comm);
		nlocalPoints = np + zeroOffset;
	}
	else
	{
		localPoints = malloc(sizeof(point) * nTotalPoints );
		if(localPoints == NULL){
			printf("Out of Memory! Cant locate buffer for local points!");
			exit(EXIT_FAILURE);
		}

		MPI_Bcast( &localPoints[0], nTotalPoints, PointType, 0, mpi_comm);
		//MPI_Scatter(NULL, np, PointType, &localPointOffset[localPointOffset] , np, PointType, 0, mpi_comm);
		nlocalPoints = np;
	}

#ifdef DEBUG_HIGH
	sprintf(b,"--- Received Points , Node [%d] ---\n", mpi_id );
		for(i=0; i < nTotalPoints; i++){
		 sprintf( &b[strlen(b)], "[%g, %g]\n", localPoints[i].x, localPoints[i].y );
		}
	sprintf(&b[strlen(b)], "--- END ---\n");
	printf(b); fflush(stdout);
#endif

	trivalSearch(nlocalPoints, &localPoints[localPointOffset], &localSolution);

	//Start binary Tree data collection
	height = log(mpi_size)/log(2);
	depth = 0;

	while(depth<height){
		parent = mpi_id - (mpi_id % (2<<depth)) ;
		if(mpi_id == parent){
			printf("Collect :[%d] [%d] <- [%d]\n",depth, mpi_id, mpi_id + (1<<depth) );

			//MPI_Recv(&remoteSolution, 1, DistanceType, mpi_id + (1<<depth), MPI_ANY_TAG, mpi_comm, MPI_STATUS_IGNORE );
			MPI_Recv(&remoteSolution, 1, DistanceType, 1, MPI_ANY_TAG, mpi_comm, MPI_STATUS_IGNORE );
			//MPI_Recv(&depth, 1, MPI_INT, mpi_id + (1<<depth), MPI_ANY_TAG, mpi_comm, MPI_STATUS_IGNORE );

			/*
			printf("Bleeeep\n");fflush(stdout);

			if(remoteSolution.d < localSolution.d){
				memcpy(&localSolution, &remoteSolution, sizeof(solution));
			}

			printf("Joining blocks"); fflush(stdout);
			joinPointBlocks(nlocalPoints, np, &localPoints[localPointOffset], &localPoints[zeroOffset + np*(mpi_id + (1<<depth)) ], &localSolution);
*/
			printf("[%d] New best Solution [%g]\n", mpi_id, localSolution.d); fflush(stdout);

		} else {
			printf("Send :[%d] [%d] -> [%d] , best Solution [%g]\n", depth, mpi_id, parent, localSolution.d); fflush(stdout);

			//MPI_Send(&localSolution, 1, DistanceType, parent, depth, mpi_comm);
			MPI_Send(&localSolution, 1, DistanceType, 0, 100, mpi_comm);

			//MPI_Send(&depth, 1, MPI_INT, parent, depth, mpi_comm);
			depth = height; //Exit after this
		}
		depth++;
	}

	if(mpi_id == 0){
		printf("Node [%d] , lowerstDistance [%g]\n", mpi_id, localSolution.d);
	} else {
		printf("Node [%d] ending ...\n", mpi_id);
	}

}

void trivalSearch(long np, point* points, distance* solution)
{
	long i,j;
	long t1;
	point* a;
	point* b;

	solution->d = 1.0;

	for(i = 0; i < np; i++){
		a = &points[i];
		for(j = i+1; j<np; j++){
			b = &(points[j]);
			if( (eqal(*a, *b) == 0) && ( dist(*a, *b) < solution->d) ){
				solution->d = dist(*a, *b);
				solution->a = *a;
				solution->b = *b;
			}
		}
	}


	//Solution should be set
}
