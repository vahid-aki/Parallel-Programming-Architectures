//ONLY MODIFY THIS FILE!
//YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "bmm.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

// TILEX and TILEY are used to set the number of threads in a CUDA block 
#define TILEX 32
#define TILEY 16

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!
#define minTile (TILEX < TILEY ? TILEX : TILEY)

dim3 getDimGrid(const int m, const int n) {
	dim3 dimGrid(n/TILEX,n/TILEY);
	return dimGrid;
}
dim3 getDimBlock(const int m, const int n) {
	dim3 dimBlock(TILEX,TILEY);
	return dimBlock;
}
__global__ void kernelFunc(float* ad, float* bd, float* cd, const int m, const int n) {
	//const int T = ((TILEX/TILEY)>4 || (TILEY/TILEX)>4) ? 8*minTile : 4*minTile;
	const int T = ((TILEX==32) && (TILEY==32)) ? 4*minTile : 8*minTile;

	// write your GPU kernel function here
	__shared__ float ads[TILEY][T];
	__shared__ float bds[T][TILEX];
	
	int Row = by * TILEY + ty;
	int Col = bx * TILEX + tx;
	
	int lx = T/TILEX;
	int ly = T/TILEY;
	
	float sum = 0.0;

	for(int p=0; p<n/T; p++){
		for(int k = 0; k < lx; k++)
			ads[ty][tx + k*TILEX] = ad[Row*n + p*T + tx + k*TILEX];		

		for(int k = 0; k < ly; k++)
			bds[ty + k*TILEY][tx] = bd[(p*T + ty + k*TILEY)*n + Col];
		__syncthreads();

		for (int k = 0; k < T; k++)
			sum += ads[ty][k] * bds[k][tx];
		__syncthreads();
	}
	cd[Row*n + Col] = sum;
}

