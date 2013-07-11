#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mex.h>
#include "matrix.h"

MATRIXD convertMxArray2MatrixD(const mxArray *xData){
	double *xValues = mxGetPr(xData);
	const mwSize *a = mxGetDimensions(xData);
	long dim1 = a[0];
	long dim2 = a[1];
	long dim3 = a[2];
	MATRIXD res = allocateMatrixD(dim1,dim2,dim3);
	long i, j, k, cnt = 0;
	
	for (k = 0; k < dim3; k++)
		for (j = 0; j < dim2; j++)
			for (i = 0; i < dim1; i++)
				res.data[i][j][k] = xValues[cnt++];	
	return res;
}
MATRIX convertMxArray2Matrix(const mxArray *xData){
	double *xValues = mxGetPr(xData);
	const mwSize *a = mxGetDimensions(xData);
	long dim1 = a[0];
	long dim2 = a[1];
	long dim3 = a[2];
	MATRIX res = allocateMatrix(dim1,dim2,dim3);
	long i, j, k, cnt = 0;
	
	for (k = 0; k < dim3; k++)
		for (j = 0; j < dim2; j++)
			for (i = 0; i < dim1; i++)
				res.data[i][j][k] = (int)xValues[cnt++];	
	return res;
}
mxArray *convertMatrixD2MxArray(MATRIXD M){
	mwSize dims[3];
	mxArray *res;
	double *xValues;
	long i, j, k, cnt = 0;
	
	dims[0] = M.dx;
	dims[1] = M.dy;
	dims[2] = M.dz;

	res = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
	xValues = mxGetPr(res);
	
	for (k = 0; k < M.dz; k++)
		for (j = 0; j < M.dy; j++)
			for (i = 0; i < M.dx; i++)
				xValues[cnt++] = M.data[i][j][k];
	
	return res;
}
mxArray *convertMatrix2MxArray(MATRIX M){
	mwSize dims[3];
	mxArray *res;
	double *xValues;
	long i, j, k, cnt = 0;
	
	dims[0] = M.dx;
	dims[1] = M.dy;
	dims[2] = M.dz;
	
	res = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
	xValues = mxGetPr(res);
	
	for (k = 0; k < M.dz; k++)
		for (j = 0; j < M.dy; j++)
			for (i = 0; i < M.dx; i++)
				xValues[cnt++] = M.data[i][j][k];	
	
	return res;
}
mxArray *convertArray2MxArray(long *M, int no){
	mxArray *res;
	double *xValues;
	long i;
	
	res = mxCreateDoubleMatrix(1,no,mxREAL);
	xValues = mxGetPr(res);
	
	for (i = 0; i < no; i++)
		xValues[i] = M[i];
	
	return res;
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){	
	MATRIXD im = convertMxArray2MatrixD(prhs[0]);
	MATRIX  init_mask = convertMxArray2Matrix(prhs[1]);
	MATRIX seg, tmap;
	MATRIXD phi;
	long *ls_vols, i;
	int max_its	 = (int)(mxGetScalar(prhs[2]));
	double E	 = (double)(mxGetScalar(prhs[3]));
	double T	 = (double)(mxGetScalar(prhs[4]));
	double alpha = (double)(mxGetScalar(prhs[5]));
	int offset;
		
//	time_t   start, finish;
//	time( &start );

	if (nrhs == 7){
		offset	 = (int)(mxGetScalar(prhs[6]));
		levelSet3Dfast(im,init_mask,max_its,E,T,alpha,&seg,&tmap,&phi,&ls_vols,offset);
	}
	else
		levelSet3D(im,init_mask,max_its,E,T,alpha,&seg,&tmap,&phi,&ls_vols);
		
//	for (i = 1; i <= max_its; i++)
//		printf("Iteration %d of %d seg vol = %ld\n",i,max_its,ls_vols[i-1]);
	
//	time( &finish );
//	printf( "\nProgram takes %6.0f seconds.\n", difftime( finish, start ) );

	plhs[0] = convertMatrix2MxArray(seg);
	plhs[1] = convertMatrixD2MxArray(phi);
	plhs[2] = convertArray2MxArray(ls_vols,max_its);
	plhs[3] = convertMatrix2MxArray(tmap);
	
	free(ls_vols);
	freeMatrix(seg);
	freeMatrix(tmap);
	freeMatrix(init_mask);
	freeMatrixD(phi);
	freeMatrixD(im);
}



