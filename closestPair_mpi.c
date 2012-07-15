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
#include <assert.h>

#include "closestPair_tools.h"

//#define DEBUG_HIGH
//#define DEBUG_MEDIUM

void inputCheck(int argc, char** argv);
void mpiFinish(void);
//void mpiInit(int argc, char** argv);
void printUsage(int argc, char** argv);
char generatePoints(long np, point** points);
void trivalSearch(long np, point* points, distance* solution);
void multiSearch(long nTotalPoints, point* points, distance* solution);
void searchBlocks(long np1, long np2, point* block1, point* block2, distance* bestsolution);

int randInit;

int main(int argc, char** argv)
{
	point* points;
	distance solution;
	double elapsedTime;

	points = NULL;

	inputCheck(argc, argv);

	printf("Generating [%d] points\n", np);

	if( generatePoints(np, &points) != EXIT_SUCCESS ){
		printf("Generating Points failed!\n");
		exit(EXIT_FAILURE);
	}

	tick();

	mpiInit(argc, argv);

	if( prepareMPIComm() )
	{


		if(mpi_id == 0)
		{
			printf("Starting search ...");
		}


			multiSearch(np, points, &solution);
		elapsedTime = tack();

		//printf("Found Solution a[%f,%f] , b[%f,%f] distance [%0.10f]\n", solution.a.x, solution.a.y, solution.b.x, solution.b.y, solution.d);


		if(mpi_id == 0){
			printf("Completed Search and found closest points at [%g, %g] , [%g, %g] with a distance of [%g]\n", mpi_id\
							, solution.a.x ,solution.a.y, solution.b.x, solution.b.y
							, solution.d);
			printf("Operation took %f seconds \n", elapsedTime);
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

	//Check for second argument being the random init value
	if(argc == 3)
	{
		errno = 0; // To distinguish success/failure after call
		randInit = atol(argv[2]);

		/* Check for various possible errors */

		if ((errno == ERANGE && (np == LONG_MAX || np == LONG_MIN)) || (errno != 0 && np == 0)) {
		   perror("atol");
		   exit(EXIT_FAILURE);
		}
	} else {
		randInit = -1;
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

	if(randInit < 0)
	{
		t1 = time( NULL );
		t2 = gmtime(&t1);
		sec = t2->tm_sec;

		srand(sec);
	} else {
		printf("Using random init [%d]\n", randInit);
		srand(randInit);
	}

	equals = 0;
	for(i = 0; i < np; i++)
	{
		(*points)[i].x = (float) rand() / RAND_MAX;
		(*points)[i].y = (float) rand() / RAND_MAX;

		//Check if this point is unique
		for(j = 0; j < i; j++){
			if(eqal((*points)[i], (*points)[j]))
			{
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
//
// How it works
// * Find the x separating Value
// * Create two working subset of points that reside in the region around X +- minDistance
// * Order each subset by the y value its points
// * Start a search with the lowest y value in the left set over all points in the right that are in a y distance of +- minDistance
void searchBlocks(long np1, long np2, point* block1, point* block2, distance* bestsolution)
{
	float xSep;
	double minD;
	double newMinD;
	double d;
	long nL, nR;
	point *L, *R;
	point* pt;
	point *lp, *rp;
	long i, j;

	assert( (np1 > 0) && "Block has no elements!" );
	assert( (np2 > 0) && "Block has no elements!" );
	assert( (block1 != NULL ) && "Illegal Block reference" );
	assert( (block2 != NULL ) && "Illegal Block reference" );

	xSep = ( block1[np1-1].x + block2[0].x ) / 2.0;
	minD = bestsolution->d;

	//Get number of candidates on the left and on the right side of the separation
	nL = 0;
	pt = &block1[np1-1];
	while( (pt->x > (xSep-minD)) && (pt != block1) ){
		nL++;
		pt--;
	}

	nR = 0;
	pt = &block2[0];
	while( (pt->x < (xSep+minD)) && (pt != &block2[np2])){
		nR++;
		pt++;
	}

	printf("Found xSep [%g] , minD [%g], nL [%d], nR [%d]\n", xSep, minD, nL, nR);

	//Create subset of left and right points , and sort it
	L = (point*) malloc(sizeof(point) * nL);
	R = (point*) malloc(sizeof(point) * nR);

	if( (L == NULL) || (R == NULL) ){
		printf("Out of Memory!");
		MPI_Abort(mpi_comm, 23);
		exit(EXIT_FAILURE);
	}

	memcpy(L, &block1[np1 - nL], sizeof(point) * nL);
	memcpy(R, &block2[0], sizeof(point) * nR);

	qsort(L, nL, sizeof(point), cmpPoint_y );
	qsort(R, nR, sizeof(point), cmpPoint_y );


#ifdef DEBUG_MEDIUM
	char b[2000];
	sprintf(b,"--- Candidate points on node [%d] ---\n", mpi_id );
	sprintf(&b[strlen(b)],"Left Side \n" );
		for(i=0; i < nL; i++){
		 sprintf( &b[strlen(b)], "[%g, %g]\n", L[i].x, L[i].y );
		}
	sprintf(&b[strlen(b)],"Right Side \n" );
			for(i=0; i < nR; i++){
			 sprintf( &b[strlen(b)], "[%g, %g]\n", R[i].x, R[i].y );
			}
	sprintf(&b[strlen(b)], "--- END ---\n");
	printf(b); fflush(stdout);
#endif

	//Start Search: For each Point on the Left calculate its distance to candidates on the right

	newMinD = minD;

	for(i = 0; i < nL; i++)
	{
		//printf("Test L %d\n",i);
		lp = &L[i];
		rp = &R[0];
		j = 0;

		//All points lower than test region
		//Skip all points that are lower than the test region
		while( ((lp->y - minD) > rp->y) && (j < nR) ){
			rp++;
			j++;
		}

		//Test all points that are not higher than the test region
		while( ((lp->y + minD) > rp->y) && (j < nR)){
			//Check distance

			d = dist(*lp, *rp);
			d = sqrt(d);

			//printf("Cmp [%g,%g] with [%g,%g], d[%g]\n",lp->x, lp->y, rp->x, rp->y, d);

			if(d < newMinD){
				memcpy(&(bestsolution->a), lp, sizeof(point));
				memcpy(&(bestsolution->b), rp, sizeof(point));
				bestsolution->d = d;
			}
			rp++;
			j++;
		}
	}

	printf("Found closest points to be [%g, %g] [%g, %g] distance [%g]\n"\
			, bestsolution->a.x, bestsolution->a.y \
			, bestsolution->b.x, bestsolution->b.y \
			, bestsolution->d);
}


//Follow solution on
// http://en.wikipedia.org/wiki/Closest_pair_of_points_problem
void multiSearch(long nTotalPoints, point* points, distance* solution)
{

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

		localPoints = points;
		MPI_Bcast( &localPoints[0], nTotalPoints, PointType, 0, mpi_comm);

		//MPI_Scatter(&points[zeroOffset], np, PointType, &localPoints[localPointOffset] , np, PointType, 0, mpi_comm);
		nlocalPoints = np + zeroOffset;
		localPointOffset = 0;
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
	char b[2000];
	sprintf(b,"--- Received Points , Node [%d] ---\n", mpi_id );
		for(i=0; i < nTotalPoints; i++){
		 sprintf( &b[strlen(b)], "[%g, %g]\n", localPoints[i].x, localPoints[i].y );
		}
	sprintf(&b[strlen(b)], "--- END ---\n");
	printf(b); fflush(stdout);
#endif

	trivalSearch(nlocalPoints, &localPoints[localPointOffset], &localSolution);
	printf("[%d] Closest pair in points [%d - %d] has a distance of [%g] \n", mpi_id, localPointOffset, localPointOffset + nlocalPoints - 1, localSolution.d);

	//Start binary Tree data collection
	height = log(mpi_size)/log(2);
	depth = 0;

	/*
	if(mpi_id == 0){
		printf("Collect :[%d] [%d] <- [%d]\n",depth, mpi_id, 1 );
		MPI_Recv(&remoteSolution, 1, DistanceType, 1, 100, mpi_comm, MPI_STATUS_IGNORE );
		printf("[%d] New best Solution [%g]\n", mpi_id, localSolution.d); fflush(stdout);
	} else {
		printf("Send :[%d] [%d] -> [%d] , best Solution [%g]\n", depth, mpi_id, 0, localSolution.d); fflush(stdout);
		MPI_Send(&localSolution, 1, DistanceType, 0, 100, mpi_comm);
		printf("Finished Send\n"); fflush(stdout);
	}
	*/


	while(depth<height){
		parent = mpi_id - (mpi_id % (2<<depth)) ;
		if(mpi_id == parent){

			//Get Solution of remove block
			MPI_Recv(&remoteSolution, 1, DistanceType, mpi_id + (1<<depth), MPI_ANY_TAG, mpi_comm, MPI_STATUS_IGNORE );

			if(remoteSolution.d < localSolution.d){
				memcpy(&localSolution, &remoteSolution, sizeof(distance));
			}

			printf("[%d] Joining blocks [%d] <- [%d] , points [%d-%d] [%d-%d]\n", mpi_id, mpi_id, mpi_id + (1<<depth) \
											,localPointOffset, localPointOffset + nlocalPoints - 1\
											, zeroOffset + np*(mpi_id + (1<<depth)), zeroOffset + np*(mpi_id + (1<<depth)) + (1<<depth)*np - 1\
											 ) ; fflush(stdout);

			searchBlocks(nlocalPoints, (1<<depth)*np, &localPoints[localPointOffset], &localPoints[zeroOffset + np*(mpi_id + (1<<depth)) ], &localSolution);

			printf("[%d] Joining blocks [%d] <- [%d] , points [%d-%d] [%d-%d] , best Solution [%g]\n", mpi_id, mpi_id, mpi_id + (1<<depth) \
								,localPointOffset, localPointOffset + nlocalPoints - 1\
								, zeroOffset + np*(mpi_id + (1<<depth)), zeroOffset + np*(mpi_id + (1<<depth)) + (1<<depth)*np - 1\
								, localSolution.d ) ; fflush(stdout);

			//Adapt nlocalPoints
			nlocalPoints += (1<<depth)*np;
		} else {
			//printf("Send :[%d] [%d] -> [%d] , best Solution [%g]\n", depth, mpi_id, parent, localSolution.d); fflush(stdout);
			MPI_Send(&localSolution, 1, DistanceType, parent, depth, mpi_comm);
			depth = height; //Exit after this
		}
		depth++;
	}


//	if(mpi_id == 0){
//		printf("Completed Search and found closest points at [%g, %g] , [%g, %g] with a distance of [%g]\n", mpi_id\
//				, localSolution.a.x ,localSolution.a.y, localSolution.b.x, localSolution.b.y
//				, localSolution.d);
//	} else {
//		//printf("Node [%d] Best solution found is [%g] ...\n", mpi_id, localSolution.d);
//	}
	memcpy(solution, &localSolution, sizeof(distance));
}
