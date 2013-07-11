#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

/***********************************************************/
/****************** 3D matrix of integers ******************/
/***********************************************************/
MATRIX allocateMatrix(long dx, long dy, long dz){
	MATRIX M; 
	long i, j;
	
	if (dx <= 0 || dy <= 0 || dz <= 0){
		M.dx = M.dy = M.dz = 0;
		M.data = NULL;
	}
	else {
		M.dx = dx; 
		M.dy = dy;
		M.dz = dz;
		M.data = (int ***) malloc(dx * sizeof(int**));
		for (i = 0; i < dx; i++){
			M.data[i] = (int **) malloc(dy * sizeof(int *));
			for (j = 0; j < dy; j++)
				M.data[i][j] = (int *)malloc(dz * sizeof(int));
		}
	}
	return M;
}
void freeMatrix(MATRIX M){
	long i, j;
	
	if (M.data == NULL)
		return;
	
	for (i = 0; i < M.dx; i++){
		for (j = 0; j < M.dy; j++)
			free(M.data[i][j]);
		free(M.data[i]);
	}
	free(M.data);
}
MATRIX readMatrix(char *fname){
	MATRIX M;
	FILE *id = fopen(fname,"r");
	int dx, dy, dz, i, j, k;
	
	if (id == NULL){
		printf("\nError: File %s does not exist\n\n",fname);
		exit(1);
	}
	fscanf(id,"%d%d%d",&dx,&dy,&dz);
	M = allocateMatrix(dx,dy,dz);
	
	for (i = 0; i < dx; i++)
		for (j = 0; j < dy; j++)
			for (k = 0; k < dz; k++)
				fscanf(id,"%d",&(M.data[i][j][k]));
	return M;
}
void initializeMatrix(MATRIX *M, int c){
	long i, j, k;
	
	for (i = 0; i < M->dx; i++)
		for (j = 0; j < M->dy; j++)
			for (k = 0; k < M->dz; k++)
				M->data[i][j][k] = c;
}
void write2DMatrixIntoFile(MATRIX M, int dim3, char *filename, int headerFlag){
    long i, j;
    FILE *id = fopen(filename,"w");
    
	if (headerFlag)
		fprintf(id,"%ld\t%ld\n",M.dx,M.dy);
	for (i = 0; i < M.dx; i++){
		for (j = 0; j < M.dy; j++)
			fprintf(id,"%d ",M.data[i][j][dim3]);
		fprintf(id,"\n");
	}
	fclose(id);
}
void append2DMatrixIntoFile(FILE *fid, MATRIX M, int dim3){
	long i, j;
	
	for (i = 0; i < M.dx; i++){
		for (j = 0; j < M.dy; j++)
			fprintf(fid,"%d ",M.data[i][j][dim3]);
		fprintf(fid,"\n");
	}
}
void copyMatrix(MATRIX *dest, MATRIX src){
	long i, j, k;
	
	if (dest->dx != src.dx  || dest->dy != src.dy  || dest->dz != src.dz){
		printf("\nError: Dimensions mismatch (copyMatrix)\n\n");
		exit(1);
	}
	for (i = 0; i < src.dx; i++)
		for (j = 0; j < src.dy; j++)
			for (k = 0; k < src.dz; k++)
				dest->data[i][j][k] = src.data[i][j][k];	
}
void copyMatrixWithLimits(MATRIX *dest, MATRIX src, int minx, int maxx, int miny, int maxy, int minz, int maxz){
	long i, j, k;
	
	if (dest->dx != src.dx  || dest->dy != src.dy  || dest->dz != src.dz){
		printf("\nError: Dimensions mismatch (copyMatrix)\n\n");
		exit(1);
	}
	for (i = minx; i <= maxx; i++)
		for (j = miny; j <= maxy; j++)
			for (k = minz; k <= maxz; k++)
				dest->data[i][j][k] = src.data[i][j][k];	
	
}
/***********************************************************/
/****************** 3D matrix of doubles *******************/
/***********************************************************/
MATRIXD allocateMatrixD(long dx, long dy, long dz){
	MATRIXD M; 
	long i, j;
	
	if (dx <= 0 || dy <= 0 || dz <= 0){
		M.dx = M.dy = M.dz = 0;
		M.data = NULL;
	}
	else {
		M.dx = dx; 
		M.dy = dy;
		M.dz = dz;
		M.data = (double ***) malloc(dx * sizeof(double **));
		for (i = 0; i < dx; i++){
			M.data[i] = (double **) malloc(dy * sizeof(double *));
			for (j = 0; j < dy; j++)
				M.data[i][j] = (double *)malloc(dz * sizeof(double));
		}
	}
	return M;
}
void freeMatrixD(MATRIXD M){
	long i, j;
	
	if (M.data == NULL)
		return;
	
	for (i = 0; i < M.dx; i++){
		for (j = 0; j < M.dy; j++)
			free(M.data[i][j]);
		free(M.data[i]);
	}
	free(M.data);
}
MATRIXD readMatrixD(char *fname){
	MATRIXD M;
	FILE *id = fopen(fname,"r");
	int dx, dy, dz, i, j, k;
	
	if (id == NULL){
		printf("\nError: File %s does not exist\n\n",fname);
		exit(1);
	}
	fscanf(id,"%d%d%d",&dx,&dy,&dz);
	M = allocateMatrixD(dx,dy,dz);
	
	for (i = 0; i < dx; i++)
		for (j = 0; j < dy; j++)
			for (k = 0; k < dz; k++)
				fscanf(id,"%lf",&(M.data[i][j][k]));
	return M;
}
void initializeMatrixD(MATRIXD *M, double c){
	long i, j, k;
	
	for (i = 0; i < M->dx; i++)
		for (j = 0; j < M->dy; j++)
			for (k = 0; k < M->dz; k++)
				M->data[i][j][k] = c;
}
void write2DMatrixDIntoFile(MATRIXD M, int dim3, char *filename, int precision, int headerFlag){
	long i, j;
	char temp[100];
    FILE *id = fopen(filename,"w");
    
	if (headerFlag)
		fprintf(id,"%ld\t%ld\n",M.dx,M.dy);
	sprintf(temp,"%c.%d%c ",37,precision,'f');
	for (i = 0; i < M.dx; i++){
		for (j = 0; j < M.dy; j++)
			fprintf(id,temp,M.data[i][j][dim3]);
		fprintf(id,"\n");
	}
	fclose(id);
}
void append2DMatrixDIntoFile(FILE *fid, MATRIXD M, int dim3, int precision){
	long i, j;
	char temp[100];
    
	sprintf(temp,"%c.%d%c ",37,precision,'f');
	for (i = 0; i < M.dx; i++){
		for (j = 0; j < M.dy; j++)
			fprintf(fid,temp,M.data[i][j][dim3]);
		fprintf(fid,"\n");
	}
}
void incrementMatrixD(MATRIXD *res, MATRIXD inc){
	long i, j, k;
	
	if (res->dx != inc.dx  || res->dy != inc.dy  || res->dz != inc.dz){
		printf("\nError: Dimensions mismatch (incrementMatrixD)\n\n");
		exit(1);
	}
	for (i = 0; i < inc.dx; i++)
		for (j = 0; j < inc.dy; j++)
			for (k = 0; k < inc.dz; k++)
				res->data[i][j][k] += inc.data[i][j][k];
}
void decrementMatrixD(MATRIXD *res, MATRIXD dec){
	long i, j, k;

	if (res->dx != dec.dx  || res->dy != dec.dy  || res->dz != dec.dz){
		printf("\nError: Dimensions mismatch (decrementMatrixD)\n\n");
		exit(1);
	}
	for (i = 0; i < dec.dx; i++)
		for (j = 0; j < dec.dy; j++)
			for (k = 0; k < dec.dz; k++)
				res->data[i][j][k] -= dec.data[i][j][k];	
}
double findMaxMatrixD(MATRIXD M){
	long i, j, k;
	double maxEntry = M.data[0][0][0];

	for (i = 0; i < M.dx; i++)
		for (j = 0; j < M.dy; j++)
			for (k = 0; k < M.dz; k++)
				if (M.data[i][j][k] > maxEntry)
					maxEntry = M.data[i][j][k];
	return maxEntry;
}
/***********************************************************/
/***********************************************************/
/***********************************************************/
void bwdist3D(MATRIX M, int background, MATRIXD *dist, int *ptsx, int *ptsy, int *ptsz, int *queuex, int *queuey, int *queuez, MATRIX labels){
	long dx = M.dx;
	long dy = M.dy;
	long dz = M.dz;
	long i, j, k, no, qstart, qend, currx, curry, currz, currlb;
	
	initializeMatrix(&labels,-2);
	
	// initialize the boundary points 
	no = 0;
	for (i = 0; i < dx; i++)
		for (j = 0; j < dy; j++)
			for (k = 0; k < dz; k++) {
				if (M.data[i][j][k] == background)
					continue;
			
				labels.data[i][j][k] = -1;
				if ((i-1 >= 0) && (M.data[i-1][j][k] == background))			labels.data[i][j][k] = 0;
				else if ((i+1 < dx) && (M.data[i+1][j][k] == background))		labels.data[i][j][k] = 0;
				else if ((j-1 >= 0) && (M.data[i][j-1][k] == background))		labels.data[i][j][k] = 0;
				else if ((j+1 < dy) && (M.data[i][j+1][k] == background))		labels.data[i][j][k] = 0;
				else if ((k-1 >= 0) && (M.data[i][j][k-1] == background))		labels.data[i][j][k] = 0;
				else if ((k+1 < dz) && (M.data[i][j][k+1] == background))		labels.data[i][j][k] = 0;
			
				if (labels.data[i][j][k] == 0){
					ptsx[no] = queuex[no] = i;
					ptsy[no] = queuey[no] = j;
					ptsz[no] = queuez[no] = k;
					labels.data[i][j][k] = no;
					no++;
				}
			}

	// grow the boundary points
	qstart = 0;
	qend = no - 1;
	while (qstart <= qend){
		currx = queuex[qstart];
		curry = queuey[qstart];
		currz = queuez[qstart];
		currlb = labels.data[currx][curry][currz];
		qstart++;
		
		if ((currx-1 >= 0) && (labels.data[currx-1][curry][currz] == -2)){
			qend++;
			queuex[qend] = currx-1;
			queuey[qend] = curry;
			queuez[qend] = currz;
			labels.data[currx-1][curry][currz] = currlb;
		}
		if ((currx+1 < dx) && (labels.data[currx+1][curry][currz] == -2)){
			qend++;
			queuex[qend] = currx+1;
			queuey[qend] = curry;
			queuez[qend] = currz;
			labels.data[currx+1][curry][currz] = currlb;
		}
		if ((curry-1 >= 0) && (labels.data[currx][curry-1][currz] == -2)){
			qend++;
			queuex[qend] = currx;
			queuey[qend] = curry-1;
			queuez[qend] = currz;
			labels.data[currx][curry-1][currz] = currlb;
		}
		if ((curry+1 < dy) && (labels.data[currx][curry+1][currz] == -2)){
			qend++;
			queuex[qend] = currx;
			queuey[qend] = curry+1;
			queuez[qend] = currz;
			labels.data[currx][curry+1][currz] = currlb;
		}
		if ((currz-1 >= 0) && (labels.data[currx][curry][currz-1] == -2)){
			qend++;
			queuex[qend] = currx;
			queuey[qend] = curry;
			queuez[qend] = currz-1;
			labels.data[currx][curry][currz-1] = currlb;
		}
		if ((currz+1 < dz) && (labels.data[currx][curry][currz+1] == -2)){
			qend++;
			queuex[qend] = currx;
			queuey[qend] = curry;
			queuez[qend] = currz+1;
			labels.data[currx][curry][currz+1] = currlb;
		}
	}
	
	// calculate the distances
	for (i = 0; i < dx; i++)
		for (j = 0; j < dy; j++)
			for (k = 0; k < dz; k++){
				if (labels.data[i][j][k] == -1)
					dist->data[i][j][k] = 0.0;
				else
					dist->data[i][j][k] = sqrt(SQUARE(i - ptsx[labels.data[i][j][k]]) + 
											   SQUARE(j - ptsy[labels.data[i][j][k]]) +
											   SQUARE(k - ptsz[labels.data[i][j][k]])); // Euclidean distance
			}
}
