#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "levelSet3D.h"

void initializeGlobals4bwdist(int dx, int dy, int dz, int **px, int **py, int **pz, 
							  int **qx, int **qy, int **qz, MATRIX *labels){
	long cnt = dx * dy * dz;
	
	*px = (int *)calloc(cnt,sizeof(int));
	*py = (int *)calloc(cnt,sizeof(int));
	*pz = (int *)calloc(cnt,sizeof(int));
	*qx = (int *)calloc(cnt,sizeof(int));
	*qy = (int *)calloc(cnt,sizeof(int));
	*qz = (int *)calloc(cnt,sizeof(int));
	*labels = allocateMatrix(dx,dy,dz);	
}
void freeGlobals4bwdist(int *px, int *py, int *pz, int *qx, int *qy, int *qz, MATRIX labels){
	free(px);
	free(py);
	free(pz);
	free(qx);
	free(qy);
	free(qz);
	freeMatrix(labels);
}
void createSignedDistanceMap(MATRIX mask, MATRIXD *dist, int *px, int *py, int *pz, int *qx, int *qy, int *qz, MATRIX labels, MATRIXD tmp){
	long i, j, k;
	
	bwdist3D(mask,0,dist,px,py,pz,qx,qy,qz,labels);
	bwdist3D(mask,1,&tmp,px,py,pz,qx,qy,qz,labels);

	for (i = 0; i < dist->dx; i++)
		for (j = 0; j < dist->dy; j++)
			for (k = 0; k < dist->dz; k++)
				dist->data[i][j][k] = dist->data[i][j][k] - tmp.data[i][j][k] - 0.5;
}
void calculateF(MATRIXD im, MATRIXD curv, MATRIXD *res, double E, double T, double alpha, int minx, int maxx, int miny, int maxy, int minz, int maxz){
	long i, j, k;
	double oneMinusAlpha = 1 - alpha, minusAlpha = -alpha;
	
	for (i = minx; i <= maxx; i++)
		for (j = miny; j <= maxy; j++)
			for (k = minz; k <= maxz; k++)
				res->data[i][j][k] = minusAlpha * (E - fabs(im.data[i][j][k] - T)) + oneMinusAlpha * curv.data[i][j][k];
}
void calculateGradPhiWithBoundaryChecks(MATRIXD phi, MATRIXD F, MATRIXD *gradPhi, int startx, int endx, int starty, int endy, int startz, int endz,
										int minx, int maxx, int miny, int maxy, int minz, int maxz){
	double ***d = phi.data;
	double u1, u3, u4, u4m, u4p, u5, u7;
	double dxplus, dxminus, dyplus, dyminus, dzplus, dzminus;
	double gradPhiMinx, gradPhiMiny, gradPhiMinz, gradPhiMaxx, gradPhiMaxy, gradPhiMaxz, gradPhiMax, gradPhiMin; 
	long i, j, k, im1, ip1, jm1, jp1, km1, kp1;
	
	for (i = startx; i <= endx; i++)
		for (j = starty; j <= endy; j++)
			for (k = startz; k <= endz; k++){
				if (i-1 >= minx)	im1 = i-1;	else	im1 = i;
				if (i+1 <= maxx)	ip1 = i+1;	else	ip1 = i;
				if (j-1 >= miny)	jm1 = j-1;	else	jm1 = j;
				if (j+1 <= maxy)	jp1 = j+1;	else	jp1 = j;
				if (k-1 >= minz)	km1 = k-1;	else	km1 = k;
				if (k+1 <= maxz)	kp1 = k+1;	else	kp1 = k;
				
				u1  = d[im1][j][k];
				u3  = d[i][jm1][k];
				u4  = d[i][j][k];
				u4m = d[i][j][km1];
				u4p = d[i][j][kp1];
				u5  = d[i][jp1][k];
				u7  = d[ip1][j][k];				
				
				dxplus  = u5 - u4;				dxminus = u4 - u3;
				dyplus  = u7 - u4;				dyminus = u4 - u1;
				dzplus  = u4p - u4;				dzminus = u4 - u4m;
				
				gradPhiMinx = gradPhiMiny = gradPhiMinz = 0.0;
				gradPhiMaxx = gradPhiMaxy = gradPhiMaxz = 0.0;
				if (dxplus > 0)			gradPhiMaxx += SQUARE(dxplus);
				else					gradPhiMinx += SQUARE(dxplus);
				
				if (-dxminus > 0)		gradPhiMaxx += SQUARE(dxminus);
				else					gradPhiMinx += SQUARE(dxminus);
				
				if (dyplus > 0)			gradPhiMaxy += SQUARE(dyplus);
				else					gradPhiMiny += SQUARE(dyplus);
				
				if (-dyminus > 0)		gradPhiMaxy += SQUARE(dyminus);
				else					gradPhiMiny += SQUARE(dyminus);
				
				if (dzplus > 0)			gradPhiMaxz += SQUARE(dzplus);
				else					gradPhiMinz += SQUARE(dzplus);
				
				if (-dzminus > 0)		gradPhiMaxz += SQUARE(dzminus);
				else					gradPhiMinz += SQUARE(dzminus);
				
				gradPhiMax = sqrt(gradPhiMaxx + gradPhiMaxy + gradPhiMaxz);
				gradPhiMin = sqrt(gradPhiMinx + gradPhiMiny + gradPhiMinz);
				
				if (F.data[i][j][k] > 0)
					gradPhi->data[i][j][k] = gradPhiMax;
				else 
					gradPhi->data[i][j][k] = 0.0;
			}
}
void calculateGradPhi(MATRIXD phi, MATRIXD F, MATRIXD *gradPhi, int minx, int maxx, int miny, int maxy, int minz, int maxz){
	double ***d = phi.data;
	double u1, u3, u4, u4m, u4p, u5, u7;
	double dxplus, dxminus, dyplus, dyminus, dzplus, dzminus;
	double gradPhiMinx, gradPhiMiny, gradPhiMinz, gradPhiMaxx, gradPhiMaxy, gradPhiMaxz, gradPhiMax, gradPhiMin; 
	long i, j, k;
	
	calculateGradPhiWithBoundaryChecks(phi,F,gradPhi,minx,minx,miny,maxy,minz,maxz,minx,maxx,miny,maxy,minz,maxz);
	calculateGradPhiWithBoundaryChecks(phi,F,gradPhi,maxx,maxx,miny,maxy,minz,maxz,minx,maxx,miny,maxy,minz,maxz);
	calculateGradPhiWithBoundaryChecks(phi,F,gradPhi,minx,maxx,miny,miny,minz,maxz,minx,maxx,miny,maxy,minz,maxz);
	calculateGradPhiWithBoundaryChecks(phi,F,gradPhi,minx,maxx,maxy,maxy,minz,maxz,minx,maxx,miny,maxy,minz,maxz);
	calculateGradPhiWithBoundaryChecks(phi,F,gradPhi,minx,maxx,miny,maxy,minz,minz,minx,maxx,miny,maxy,minz,maxz);
	calculateGradPhiWithBoundaryChecks(phi,F,gradPhi,minx,maxx,miny,maxy,maxz,maxz,minx,maxx,miny,maxy,minz,maxz);
	
	for (i = minx+1; i < maxx; i++)
		for (j = miny+1; j < maxy; j++)
			for (k = minz+1; k < maxz; k++){
				u1  = d[i-1][j][k];
				u3  = d[i][j-1][k];
				u4  = d[i][j][k];
				u4m = d[i][j][k-1];
				u4p = d[i][j][k+1];
				u5  = d[i][j+1][k];
				u7  = d[i+1][j][k];
				
				dxplus  = u5 - u4;				dxminus = u4 - u3;
				dyplus  = u7 - u4;				dyminus = u4 - u1;
				dzplus  = u4p - u4;				dzminus = u4 - u4m;
								
				gradPhiMinx = gradPhiMiny = gradPhiMinz = 0.0;
				gradPhiMaxx = gradPhiMaxy = gradPhiMaxz = 0.0;
				if (dxplus > 0)			gradPhiMaxx += SQUARE(dxplus);
				else					gradPhiMinx += SQUARE(dxplus);

				if (-dxminus > 0)		gradPhiMaxx += SQUARE(dxminus);
				else					gradPhiMinx += SQUARE(dxminus);

				if (dyplus > 0)			gradPhiMaxy += SQUARE(dyplus);
				else					gradPhiMiny += SQUARE(dyplus);

				if (-dyminus > 0)		gradPhiMaxy += SQUARE(dyminus);
				else					gradPhiMiny += SQUARE(dyminus);

				if (dzplus > 0)			gradPhiMaxz += SQUARE(dzplus);
				else					gradPhiMinz += SQUARE(dzplus);

				if (-dzminus > 0)		gradPhiMaxz += SQUARE(dzminus);
				else					gradPhiMinz += SQUARE(dzminus);
				
				gradPhiMax = sqrt(gradPhiMaxx + gradPhiMaxy + gradPhiMaxz);
				gradPhiMin = sqrt(gradPhiMinx + gradPhiMiny + gradPhiMinz);
				
				if (F.data[i][j][k] > 0)
					gradPhi->data[i][j][k] = gradPhiMax;
				else 
					gradPhi->data[i][j][k] = gradPhiMin;
			}
}
void evolveCurve(MATRIXD *phi, MATRIXD F, MATRIXD gradPhi, int minx, int maxx, int miny, int maxy, int minz, int maxz){
	long i, j, k;
	double dt = 0.0;
	
	for (i = minx; i <= maxx; i++)
		for (j = miny; j <= maxy; j++)
			for (k = minz; k <= maxz; k++)
				if (fabs(F.data[i][j][k] * gradPhi.data[i][j][k]) > dt)
					dt = F.data[i][j][k] * gradPhi.data[i][j][k];
	dt = 0.5 / dt;
	
	for (i = minx; i <= maxx; i++)
		for (j = miny; j <= maxy; j++)
			for (k = minz; k <= maxz; k++)
				phi->data[i][j][k] += dt * F.data[i][j][k] * gradPhi.data[i][j][k];
}
long selectNonpositive(MATRIXD phi, MATRIX *seg, int minx, int maxx, int miny, int maxy, int minz, int maxz){
	long cnt = 0, i, j, k;
	
	for (i = minx; i <= maxx; i++)
		for (j = miny; j <= maxy; j++)
			for (k = minz; k <= maxz; k++)
				if (phi.data[i][j][k] <= 0){
					seg->data[i][j][k] = 1;
					cnt++;
				}
				else
					seg->data[i][j][k] = 0;
	
	return cnt;
}
void updateMaps(MATRIX *tmap, MATRIX seg1, MATRIX seg0, int no, int minx, int maxx, int miny, int maxy, int minz, int maxz){
	long i, j, k;
	
	for (i = minx; i <= maxx; i++)
		for (j = miny; j <= maxy; j++)
			for (k = minz; k <= maxz; k++)
				if (seg1.data[i][j][k] != seg0.data[i][j][k])
					tmap->data[i][j][k] = no;
}
void reinitializeDistanceFunction(MATRIXD *dist, MATRIX tmpMask, int *px, int *py, int *pz, int *qx, int *qy, int *qz, MATRIX labels, MATRIXD tmp){
	long i, j, k;
	
	for (i = 0; i < dist->dx; i++)
		for (j = 0; j < dist->dy; j++)
			for (k = 0; k < dist->dz; k++)
				if (dist->data[i][j][k] < 0)
					tmpMask.data[i][j][k] = 1;
				else
					tmpMask.data[i][j][k] = 0;
	bwdist3D(tmpMask,0,&tmp,px,py,pz,qx,qy,qz,labels);
	
	for (i = 0; i < dist->dx; i++)
		for (j = 0; j < dist->dy; j++)
			for (k = 0; k < dist->dz; k++)
				if (dist->data[i][j][k] >= 0)
					tmpMask.data[i][j][k] = 1;
				else
					tmpMask.data[i][j][k] = 0;
	bwdist3D(tmpMask,0,dist,px,py,pz,qx,qy,qz,labels);
 
	for (i = 0; i < dist->dx; i++)
		for (j = 0; j < dist->dy; j++)
			for (k = 0; k < dist->dz; k++)
				dist->data[i][j][k] = tmp.data[i][j][k] - dist->data[i][j][k];
}
void calculateCurvatureWithBoundaryChecks(MATRIXD phi, MATRIXD *curv, int startx, int endx, int starty, int endy, int startz, 
										  int endz, int minx, int maxx, int miny, int maxy, int minz, int maxz){
	double ***d = phi.data;
	double u0, u1, u1m, u1p, u2, u3, u3m, u3p, u4, u4m, u4p, u5, u5m, u5p, u6, u7, u7m, u7p, u8;
	double dx, dy, dz, dxplus, dxminus, dyplus, dyminus, dzplus, dzminus;
	double dxplusy, dxminusy, dxplusz, dxminusz, dyplusx, dyminusx;
	double dyplusz, dyminusz, dzplusx, dzminusx, dzplusy, dzminusy;
	double nplusx, nplusy, nplusz, nminusx, nminusy, nminusz;
	long i, j, k, im1, ip1, jm1, jp1, km1, kp1;
	
	for (i = startx; i <= endx; i++)
		for (j = starty; j <= endy; j++)
			for (k = startz; k <= endz; k++){
				
				if (i-1 >= minx)	im1 = i-1;	else	im1 = i;
				if (i+1 <= maxx)	ip1 = i+1;	else	ip1 = i;
				if (j-1 >= miny)	jm1 = j-1;	else	jm1 = j;
				if (j+1 <= maxy)	jp1 = j+1;	else	jp1 = j;
				if (k-1 >= minz)	km1 = k-1;	else	km1 = k;
				if (k+1 <= maxz)	kp1 = k+1;	else	kp1 = k;
				
				u0  = d[im1][jm1][k];
				u1  = d[im1][j][k];				u1m = d[im1][j][km1];			u1p = d[im1][j][kp1];
				u2  = d[im1][jp1][k];
				u3  = d[i][jm1][k];				u3m = d[i][jm1][km1];			u3p = d[i][jm1][kp1];
				u4  = d[i][j][k];				u4m = d[i][j][km1];				u4p = d[i][j][kp1];
				u5  = d[i][jp1][k];				u5m = d[i][jp1][km1];			u5p = d[i][jp1][kp1];				
				u6  = d[ip1][jm1][k];
				u7  = d[ip1][j][k];				u7m = d[ip1][j][km1];			u7p = d[ip1][j][kp1];
				u8  = d[ip1][jp1][k];
				
				dx = (u5 - u3) / 2;
				dy = (u7 - u1) / 2;
				dz = (u4p - u4m) / 2;
				
				dxplus  = u5 - u4;				dxminus = u4 - u3;
				dyplus  = u7 - u4;				dyminus = u4 - u1;
				dzplus  = u4p - u4;				dzminus = u4 - u4m;
				
				dxplusy  = (u8 - u6)   / 2;		dxminusy = (u2 - u0)   / 2;
				dxplusz  = (u5p - u3p) / 2;		dxminusz = (u5m - u3m) / 2;
				dyplusx  = (u8 - u2)   / 2;		dyminusx = (u6 - u0)   / 2;
				dyplusz  = (u7p - u1p) / 2;		dyminusz = (u7m - u1m) / 2;
				dzplusx  = (u5p - u5m) / 2;		dzminusx = (u3p - u3m) / 2;
				dzplusy  = (u7p - u7m) / 2;		dzminusy = (u1p - u1m) / 2;
				
				nplusx = dxplus / sqrt(EPSILON + SQUARE(dxplus) + SQUARE((dyplusx + dy) / 2) + SQUARE((dzplusx + dz) / 2));
				nplusy = dyplus / sqrt(EPSILON + SQUARE(dyplus) + SQUARE((dxplusy + dx) / 2) + SQUARE((dzplusy + dz) / 2));
				nplusz = dzplus / sqrt(EPSILON + SQUARE(dzplus) + SQUARE((dxplusz + dx) / 2) + SQUARE((dyplusz + dy) / 2));
				
				nminusx = dxminus / sqrt(EPSILON + SQUARE(dxminus) + SQUARE((dyminusx + dy) / 2) + SQUARE((dzminusx + dz) / 2));
				nminusy = dyminus / sqrt(EPSILON + SQUARE(dyminus) + SQUARE((dxminusy + dx) / 2) + SQUARE((dzminusy + dz) / 2));
				nminusz = dzminus / sqrt(EPSILON + SQUARE(dzminus) + SQUARE((dxminusz + dx) / 2) + SQUARE((dyminusz + dy) / 2));
				
				curv->data[i][j][k] = (nplusx - nminusx + nplusy - nminusy + nplusz - nminusz) / 2;
			}
}
void calculateCurvature(MATRIXD phi, MATRIXD *curv, int minx, int maxx, int miny, int maxy, int minz, int maxz){
	double ***d = phi.data;
	double u0, u1, u1m, u1p, u2, u3, u3m, u3p, u4, u4m, u4p, u5, u5m, u5p, u6, u7, u7m, u7p, u8;
	double dx, dy, dz, dxplus, dxminus, dyplus, dyminus, dzplus, dzminus;
	double dxplusy, dxminusy, dxplusz, dxminusz, dyplusx, dyminusx;
	double dyplusz, dyminusz, dzplusx, dzminusx, dzplusy, dzminusy;
	double nplusx, nplusy, nplusz, nminusx, nminusy, nminusz;
	long i, j, k;
	
	calculateCurvatureWithBoundaryChecks(phi,curv,minx,minx,miny,maxy,minz,maxz,minx,maxx,miny,maxy,minz,maxz);
	calculateCurvatureWithBoundaryChecks(phi,curv,maxx,maxx,miny,maxy,minz,maxz,minx,maxx,miny,maxy,minz,maxz);
	calculateCurvatureWithBoundaryChecks(phi,curv,minx,maxx,miny,miny,minz,maxz,minx,maxx,miny,maxy,minz,maxz);
	calculateCurvatureWithBoundaryChecks(phi,curv,minx,maxx,maxy,maxy,minz,maxz,minx,maxx,miny,maxy,minz,maxz);
	calculateCurvatureWithBoundaryChecks(phi,curv,minx,maxx,miny,maxy,minz,minz,minx,maxx,miny,maxy,minz,maxz);
	calculateCurvatureWithBoundaryChecks(phi,curv,minx,maxx,miny,maxy,maxz,maxz,minx,maxx,miny,maxy,minz,maxz);
	
	for (i = minx+1; i < maxx; i++)
		for (j = miny+1; j < maxy; j++)
			for (k = minz+1; k < maxz; k++){
				
				u0  = d[i-1][j-1][k];
				u1  = d[i-1][j][k];				u1m = d[i-1][j][k-1];			u1p = d[i-1][j][k+1];
				u2  = d[i-1][j+1][k];
				u3  = d[i][j-1][k];				u3m = d[i][j-1][k-1];			u3p = d[i][j-1][k+1];
				u4  = d[i][j][k];				u4m = d[i][j][k-1];				u4p = d[i][j][k+1];
				u5  = d[i][j+1][k];				u5m = d[i][j+1][k-1];			u5p = d[i][j+1][k+1];				
				u6  = d[i+1][j-1][k];
				u7  = d[i+1][j][k];				u7m = d[i+1][j][k-1];			u7p = d[i+1][j][k+1];
				u8  = d[i+1][j+1][k];
				
				dx = (u5 - u3) / 2;
				dy = (u7 - u1) / 2;
				dz = (u4p - u4m) / 2;
				
				dxplus  = u5 - u4;				dxminus = u4 - u3;
				dyplus  = u7 - u4;				dyminus = u4 - u1;
				dzplus  = u4p - u4;				dzminus = u4 - u4m;
				
				dxplusy  = (u8 - u6)   / 2;		dxminusy = (u2 - u0)   / 2;
				dxplusz  = (u5p - u3p) / 2;		dxminusz = (u5m - u3m) / 2;
				dyplusx  = (u8 - u2)   / 2;		dyminusx = (u6 - u0)   / 2;
				dyplusz  = (u7p - u1p) / 2;		dyminusz = (u7m - u1m) / 2;
				dzplusx  = (u5p - u5m) / 2;		dzminusx = (u3p - u3m) / 2;
				dzplusy  = (u7p - u7m) / 2;		dzminusy = (u1p - u1m) / 2;
				
				nplusx = dxplus / sqrt(EPSILON + SQUARE(dxplus) + SQUARE((dyplusx + dy) / 2) + SQUARE((dzplusx + dz) / 2));
				nplusy = dyplus / sqrt(EPSILON + SQUARE(dyplus) + SQUARE((dxplusy + dx) / 2) + SQUARE((dzplusy + dz) / 2));
				nplusz = dzplus / sqrt(EPSILON + SQUARE(dzplus) + SQUARE((dxplusz + dx) / 2) + SQUARE((dyplusz + dy) / 2));

				nminusx = dxminus / sqrt(EPSILON + SQUARE(dxminus) + SQUARE((dyminusx + dy) / 2) + SQUARE((dzminusx + dz) / 2));
				nminusy = dyminus / sqrt(EPSILON + SQUARE(dyminus) + SQUARE((dxminusy + dx) / 2) + SQUARE((dzminusy + dz) / 2));
				nminusz = dzminus / sqrt(EPSILON + SQUARE(dzminus) + SQUARE((dxminusz + dx) / 2) + SQUARE((dyminusz + dy) / 2));
				
				curv->data[i][j][k] = (nplusx - nminusx + nplusy - nminusy + nplusz - nminusz) / 2;
			}
}
void levelSet3D(MATRIXD im, MATRIX init_mask, int max_its, double E, double T, double alpha, 
				MATRIX *seg0, MATRIX *tmap, MATRIXD *phi, long **ls_vols){
	// some of these variables will be used as the local variables for other functions
	// they are defined outside (like global variables) because we want to decrease the overhead
	// corresponding to allocation and deallocation (since we have multiple calls for the same function
	long dx = im.dx, dy = im.dy, dz = im.dz;
	int *px, *py, *pz, *qx, *qy, *qz, k, tmap_it = 1;
	MATRIX intTmp;
	MATRIXD dblTmp = allocateMatrixD(dx,dy,dz);	
	MATRIXD gradPhi = allocateMatrixD(dx,dy,dz);
	
	*phi = allocateMatrixD(dx,dy,dz);
	*seg0 = allocateMatrix(dx,dy,dz);
	*tmap = allocateMatrix(dx,dy,dz);
	*ls_vols = (long *)calloc(max_its,sizeof(long)); // initialized with 0
	
	initializeMatrix(tmap,0);
	initializeGlobals4bwdist(dx,dy,dz,&px,&py,&pz,&qx,&qy,&qz,&intTmp);

	createSignedDistanceMap(init_mask,phi,px,py,pz,qx,qy,qz,intTmp,dblTmp); // temps: intTmp, dblTmp

	for (k = 1; k <= max_its; k++){
		calculateCurvature(*phi,&dblTmp,0,dx-1,0,dy-1,0,dz-1); // dblTmp keeps the curvatures
		calculateF(im,dblTmp,&dblTmp,E,T,alpha,0,dx-1,0,dy-1,0,dz-1); // first dblTmp (input): curvatures, second dblTmp (output): F values
		calculateGradPhi(*phi,dblTmp,&gradPhi,0,dx-1,0,dy-1,0,dz-1);
		evolveCurve(phi,dblTmp,gradPhi,0,dx-1,0,dy-1,0,dz-1);
		
		if (k % 50 == 0)
			reinitializeDistanceFunction(phi,init_mask,px,py,pz,qx,qy,qz,intTmp,dblTmp); // temps: init_mask, intTmp, dblTmp
		
		(*ls_vols)[k-1] = selectNonpositive(*phi,&intTmp,0,dx-1,0,dy-1,0,dz-1);
		
		if (k == 1){
			selectNonpositive(*phi,seg0,0,dx-1,0,dy-1,0,dz-1);
			copyMatrix(tmap,*seg0);
		}
		else if (k % 10 == 0){
			selectNonpositive(*phi,&intTmp,0,dx-1,0,dy-1,0,dz-1);
			updateMaps(tmap,intTmp,*seg0,tmap_it,0,dx-1,0,dy-1,0,dz-1);
			copyMatrix(seg0,intTmp);
			tmap_it++;
		}
	}
	selectNonpositive(*phi,seg0,0,dx-1,0,dy-1,0,dz-1);

	freeMatrixD(gradPhi);
	freeMatrixD(dblTmp);
	freeGlobals4bwdist(px,py,pz,qx,qy,qz,intTmp);
}
/************************************************************************************/
/* These are the functions to implement the fast version of the level set algorithm */
/* This version considers the regions within +/- offset of the initial mask in the  */
/* first iteration and +/- offset of the segmented part in the next iterations      */
/************************************************************************************/
void findLimits(MATRIX mask, int offset, int *minx, int *maxx, int *miny, int *maxy, int *minz, int *maxz){
	long i, j, k;
	
	*minx = mask.dx;		*maxx = -1;
	*miny = mask.dy;		*maxy = -1;
	*minz = mask.dz;		*maxz = -1;
	
	for (i = 0; i < mask.dx; i++)
		for (j = 0; j < mask.dy; j++)
			for (k = 0; k < mask.dz; k++)
				if (mask.data[i][j][k]){
					if (*minx > i)		*minx = i;
					if (*maxx < i)		*maxx = i;
					if (*miny > j)		*miny = j;
					if (*maxy < j)		*maxy = j;
					if (*minz > k)		*minz = k;
					if (*maxz < k)		*maxz = k;
				}
	*minx -= offset;
	*maxx += offset;
	*miny -= offset;
	*maxy += offset;
	*minz -= offset;
	*maxz += offset;
	
	if (*minx < 0)				*minx = 0;
	if (*maxx >= mask.dx)		*maxx = mask.dx - 1;
	if (*miny < 0)				*miny = 0;
	if (*maxy >= mask.dy)		*maxy = mask.dy - 1;
	if (*minz < 0)				*minz = 0;
	if (*maxz >= mask.dz)		*maxz = mask.dz - 1;
}
void levelSet3Dfast(MATRIXD im, MATRIX init_mask, int max_its, double E, double T, double alpha, 
					MATRIX *seg0, MATRIX *tmap, MATRIXD *phi, long **ls_vols, int offset){
	long dx = im.dx, dy = im.dy, dz = im.dz;
	int *px, *py, *pz, *qx, *qy, *qz, k, tmap_it = 1;
	MATRIX intTmp;
	MATRIX intTmp2 = allocateMatrix(dx,dy,dz);	// add a new temporary matrix
	MATRIXD dblTmp = allocateMatrixD(dx,dy,dz);	
	MATRIXD gradPhi = allocateMatrixD(dx,dy,dz);
	int minx, maxx, miny, maxy, minz, maxz;
	
	*phi = allocateMatrixD(dx,dy,dz);
	*seg0 = allocateMatrix(dx,dy,dz);
	*tmap = allocateMatrix(dx,dy,dz);
	*ls_vols = (long *)calloc(max_its,sizeof(long)); // initialized with 0
	
	initializeMatrix(tmap,0);
	initializeMatrix(seg0,0);		// initialized with 0
	initializeMatrix(&intTmp2,0);	// initialized with 0
	initializeGlobals4bwdist(dx,dy,dz,&px,&py,&pz,&qx,&qy,&qz,&intTmp);
	
	createSignedDistanceMap(init_mask,phi,px,py,pz,qx,qy,qz,intTmp,dblTmp); // temps: intTmp, dblTmp

	findLimits(init_mask,offset,&minx,&maxx,&miny,&maxy,&minz,&maxz);
	
	for (k = 1; k <= max_its; k++){
		calculateCurvature(*phi,&dblTmp,minx,maxx,miny,maxy,minz,maxz); // dblTmp keeps the curvatures
		calculateF(im,dblTmp,&dblTmp,E,T,alpha,minx,maxx,miny,maxy,minz,maxz); // first dblTmp (input): curvatures, second dblTmp (output): F values
		calculateGradPhi(*phi,dblTmp,&gradPhi,minx,maxx,miny,maxy,minz,maxz);
		evolveCurve(phi,dblTmp,gradPhi,minx,maxx,miny,maxy,minz,maxz);
		
		if (k % 50 == 0)
			reinitializeDistanceFunction(phi,init_mask,px,py,pz,qx,qy,qz,intTmp,dblTmp); // temps: init_mask, intTmp, dblTmp
		
		(*ls_vols)[k-1] = selectNonpositive(*phi,&intTmp2,minx,maxx,miny,maxy,minz,maxz); // use the second temporarymatrix
		
		if (k == 1){
			selectNonpositive(*phi,seg0,minx,maxx,miny,maxy,minz,maxz);
			copyMatrixWithLimits(tmap,*seg0,minx,maxx,miny,maxy,minz,maxz);
			findLimits(*tmap,offset,&minx,&maxx,&miny,&maxy,&minz,&maxz);
		}
		else if (k % 10 == 0){
			selectNonpositive(*phi,&intTmp2,minx,maxx,miny,maxy,minz,maxz);			// use the second temporarymatrix
			updateMaps(tmap,intTmp2,*seg0,tmap_it,minx,maxx,miny,maxy,minz,maxz);	// use the second temporarymatrix
			copyMatrixWithLimits(seg0,intTmp2,minx,maxx,miny,maxy,minz,maxz);		// use the second temporarymatrix
			tmap_it++;
			
			findLimits(*tmap,offset,&minx,&maxx,&miny,&maxy,&minz,&maxz);
		}
	}
	selectNonpositive(*phi,seg0,0,dx-1,0,dy-1,0,dz-1);
	//selectNonpositive(*phi,seg0,minx,maxx,miny,maxy,minz,maxz); ???? (consider this modification if there are more problems)
	
	freeMatrixD(gradPhi);
	freeMatrixD(dblTmp);
	freeMatrix(intTmp2);
	freeGlobals4bwdist(px,py,pz,qx,qy,qz,intTmp);
}





