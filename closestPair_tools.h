/*
 * nearestPoints_tools.h
 *
 *  Created on: Jul 13, 2012
 *     Author: Pasieka Manuel , mapa17@posgrado.upv.es
 */

#ifndef NEARESTPOINTS_TOOLS_H_
#define NEARESTPOINTS_TOOLS_H_

#include "mpi.h"

//Types
typedef struct { float x; float y; } point;
typedef struct { point a; point b; double d; } distance;
extern MPI_Datatype PointType;
extern MPI_Datatype DistanceType;


//Global Variables
extern long np;
extern int mpi_size, mpi_id;
extern struct timeval start;
extern MPI_Group mpi_group;
extern MPI_Comm mpi_comm;


//Functions
double dist(point a, point b);
char eqal(point a, point b);
int cmpf(float a, float b);
int cmpPoint_x(const void *a, const void *b);
int cmpPoint_y(const void *a, const void *b);
void mpiInit(int argc, char** argv);
char prepareMPIComm(void);
void trivalSearch(long np, point* points, distance* solution);

void tick(void);
double tack(void);




#endif /* NEARESTPOINTS_TOOLS_H_ */
