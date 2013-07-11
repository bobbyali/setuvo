#ifndef __matrix_h
#define __matrix_h

#include <stdio.h>
#include <stdlib.h>

#define SQUARE(a) ( (a) * (a) )
#define EPSILON 2.2204e-16

/***********************************************************/
/****************** 3D matrix of integers ******************/
/***********************************************************/
struct TMatrix{
	long dx, dy, dz;
	int ***data;
};
typedef struct TMatrix MATRIX;

MATRIX allocateMatrix(long dx, long dy, long dz);
void freeMatrix(MATRIX M);
MATRIX readMatrix(char *fname);
void initializeMatrix(MATRIX *M, int c);
void write2DMatrixIntoFile(MATRIX M, int dim3, char *filename, int headerFlag);
void append2DMatrixIntoFile(FILE *fid, MATRIX M, int dim3);
void copyMatrix(MATRIX *dest, MATRIX src);
void copyMatrixWithLimits(MATRIX *dest, MATRIX src, int minx, int maxx, int miny, int maxy, int minz, int maxz);

/***********************************************************/
/****************** 3D matrix of doubles *******************/
/***********************************************************/
struct TMatrixD{
	long dx, dy, dz;
	double ***data;
};
typedef struct TMatrixD MATRIXD;

MATRIXD allocateMatrixD(long dx, long dy, long dz);
void freeMatrixD(MATRIXD M);
MATRIXD readMatrixD(char *fname);
void initializeMatrixD(MATRIXD *M, double c);
void write2DMatrixDIntoFile(MATRIXD M, int dim3, char *filename, int precision, int headerFlag);
void append2DMatrixDIntoFile(FILE *fid, MATRIXD M, int dim3, int precision);

void incrementMatrixD(MATRIXD *res, MATRIXD inc);
void decrementMatrixD(MATRIXD *res, MATRIXD dec);
double findMaxMatrixD(MATRIXD M);

/***********************************************************/
/***********************************************************/
/***********************************************************/
// assumes that dist has already been allocated
// there is a slight difference between the outputs of this function and the matlab bwdist function
// for example, this function gives 1.41421356237309514547 whereas 
//        the matlab function gives 1.41421353816986083984
//
void bwdist3D(MATRIX M, int background, MATRIXD *dist, int *ptsx, int *ptsy, int *ptsz, int *queuex, int *queuey, int *queuez, MATRIX labels);

#endif
