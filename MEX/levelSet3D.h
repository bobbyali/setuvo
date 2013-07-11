#ifndef __levelSet3D_h
#define __levelSet3D_h

#include "matrix.h"

/***********************************************************/
/***********************************************************/
/***********************************************************/
void initializeGlobals4bwdist(int dx, int dy, int dz, int **px, int **py, int **pz, 
				int **qx, int **qy, int **qz, MATRIX *labels);
void freeGlobals4bwdist(int *px, int *py, int *pz, int *qx, int *qy, int *qz, MATRIX labels);

void createSignedDistanceMap(MATRIX mask, MATRIXD *dist, int *px, int *py, int *pz, int *qx, int *qy, int *qz, MATRIX labels, MATRIXD tmp);
void reinitializeDistanceFunction(MATRIXD *dist, MATRIX tmpMask, int *px, int *py, int *pz, int *qx, int *qy, int *qz, MATRIX labels, MATRIXD tmp);

// We have these separate files since we do not want to use if statements 
// For example, we do not want to check whether i-1 >= 0 everytime we use phi.data[i-1][j][k]
// Thus, we write a separate function for the following boundary cases:
// i-1 >= 0, i+1 < dim1, j-1 >= 0, j+1 < dim2, k-1 >= 0, k+1 < dim3	
void calculateCurvatureWithBoundaryChecks(MATRIXD phi, MATRIXD *curv, int startx, int endx, int starty, int endy, int startz, int endz, 
				int minx, int maxx, int miny, int maxy, int minz, int maxz);
void calculateCurvature(MATRIXD phi, MATRIXD *curv, int minx, int maxx, int miny, int maxy, int minz, int maxz);
void calculateF(MATRIXD im, MATRIXD curv, MATRIXD *res, double E, double T, double alpha, int minx, int maxx, int miny, int maxy, int minz, int maxz);
void calculateGradPhi(MATRIXD phi, MATRIXD F, MATRIXD *gradPhi, int minx, int maxx, int miny, int maxy, int minz, int maxz);
void calculateGradPhiWithBoundaryChecks(MATRIXD phi, MATRIXD F, MATRIXD *gradPhi, int startx, int endx, int starty, int endy, int startz, int endz,
				int minx, int maxx, int miny, int maxy, int minz, int maxz);

void evolveCurve(MATRIXD *phi, MATRIXD F, MATRIXD gradPhi, int minx, int maxx, int miny, int maxy, int minz, int maxz);
long selectNonpositive(MATRIXD phi, MATRIX *seg, int minx, int maxx, int miny, int maxy, int minz, int maxz);
void updateMaps(MATRIX *tmap, MATRIX seg1, MATRIX seg0, int no, int minx, int maxx, int miny, int maxy, int minz, int maxz);

void levelSet3D(MATRIXD im, MATRIX init_mask, int max_its, double E, double T, double alpha, 
				MATRIX *seg0, MATRIX *tmap, MATRIXD *phi, long **ls_vols);


/************************************************************************************/
/* These are the functions to implement the fast version of the level set algorithm */
/* This version considers the regions within +/- offset of the initial mask in the  */
/* first iteration and +/- offset of the segmented part in the next iterations      */
/************************************************************************************/
void findLimits(MATRIX mask, int offset, int *minx, int *maxx, int *miny, int *maxy, int *minz, int *maxz);

void levelSet3Dfast(MATRIXD im, MATRIX init_mask, int max_its, double E, double T, double alpha, 
				MATRIX *seg0, MATRIX *tmap, MATRIXD *phi, long **ls_vols, int offset);

#endif
