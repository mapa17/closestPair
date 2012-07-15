/*
 * nearestPoints.c
 *
 *  Created on: Jul 5, 2012
 *      Author: Pasieka Manuel , mapa17@posgrado.upv.es
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <time.h>

#include "closestPair_tools.h"

void inputCheck(int argc, char** argv);
void printUsage(int argc, char** argv);
char generatePoints(long np, point** points);
void trivalSearch(long np, point* points, distance* solution);

int main(int argc, char** argv)
{
	point* points;
	distance solution;
	double elapsedTime;

	inputCheck(argc, argv);

	printf("Generating [%d] points\n", np);

	if( generatePoints(np, &points) != EXIT_SUCCESS ){
		printf("Generating Points failed!\n");
		exit(EXIT_FAILURE);
	}

	printf("Starting search ...");

	tick();
		trivalSearch(np, points, &solution);
	elapsedTime = tack();

	printf("Found Solution a[%f,%f] , b[%f,%f] distance [%0.10f]\n", solution.a.x, solution.a.y, solution.b.x, solution.b.y, solution.d);
    printf("Operation took %f seconds \n", elapsedTime);

	free(points);
	exit(EXIT_SUCCESS);
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
		 //printf("p %f, %f\n", (*points)[i].x, (*points)[i].y);

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

	 return(EXIT_SUCCESS);
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
