// ONLY MODIFY THIS FILE

#include "scan2.h"
#include "gpuerrors.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!


__global__ void kernelFunc(float* ad, float* cd, float* a1d, const int BlockSize,  int inclusive_EN) {
	__shared__ float ads[1024];
	int i = tx;
	int j = ((1<<25)/2)*by + (BlockSize/2)*bx + tx;

	ads[2*i] = ad[2*j];
	ads[2*i +1] = ad[2*j +1];
	//__syncthreads();

	for(int k=1; k<=BlockSize/2; k*=2) {
		//__syncthreads();
		if( i < BlockSize/(2*k) )
			ads[(k*2)*(i+1) -1] += ads[k*(2*i+1) -1];

		__syncthreads();
	}

	float sum;
	if(i == BlockSize/2 -1) {
		sum = ads[BlockSize -1];
		a1d[(1<<15)*by + bx] = sum;
	}
	__syncthreads();
	if(i==0)
		ads[BlockSize -1] = 0;

	float temp;
	for(int k=BlockSize/2; k>=1; k= k/2){
		//__syncthreads();
		if( i < BlockSize/(2*k) ){
			temp = ads[(k*2)*(i+1) -1];
			ads[(k*2)*(i+1) -1] += ads[k*(2*i+1) -1];
			ads[k*(2*i+1) -1] = temp;
		}
		__syncthreads();
	}

	if(inclusive_EN) {
		if(i == BlockSize/2 -1) {
			cd[2*j] = ads[2*i +1];
			cd[2*j +1] = sum;
		}
		else {
			cd[2*j] = ads[2*i +1];
			cd[2*j +1] = ads[2*i +2];
		}
	}
	else {
		cd[2*j] = ads[2*i];
		cd[2*j +1] = ads[2*i+1];
	}
}


__global__ void kernelFunc2(float* cd, float* ad, float OldSum=0.0) {
	cd[(1<<25)*by + 1024*bx + tx] += ad[(1<<15)*by + bx] + OldSum;
}


void gpuKernel(float* a, float* c,int n) {
	float* ad;
	float* a1d;
	float* a2d;
	float* cd;

    	HANDLE_ERROR(cudaMalloc((void**)&cd, sizeof(float)));

			if(n < 1<<28) {
					HANDLE_ERROR(cudaMalloc((void**)&ad, n * sizeof(float)));
		    	HANDLE_ERROR(cudaMalloc((void**)&a1d, (n/1024) * sizeof(float)));
		    	HANDLE_ERROR(cudaMalloc((void**)&a2d, (n/(1024*1024)) * sizeof(float)));

		    	HANDLE_ERROR(cudaMemcpy(ad, a, n * sizeof(float), cudaMemcpyHostToDevice));

					if(n == 1<<20) {
						kernelFunc <<< n/1024, 512 >>> (ad, ad, a1d, 1024, 1);
						kernelFunc <<< n/(1024*1024), 512 >>> (a1d, a1d, a2d, 1024, 0);
						kernelFunc2 <<< n/1024, 1024 >>> (ad, a1d);
					}
					else if(n < 1<<26) {
						kernelFunc <<< n/1024, 512 >>> (ad, ad, a1d, 1024, 1);
						kernelFunc <<< n/(1024*1024), 512 >>> (a1d, a1d, a2d, 1024, 0);
						kernelFunc <<< 1, 512 >>> (a2d, a2d, cd, 1024, 0);
						kernelFunc2 <<< n/(1024*1024), 1024 >>> (a1d, a2d);
						kernelFunc2 <<< n/1024, 1024 >>> (ad, a1d);
					}
					else {
						dim3 dimGrid(1<<15, n/(1<<25));
						kernelFunc <<< dimGrid, 512 >>> (ad, ad, a1d, 1024, 1);
						kernelFunc <<< n/(1024*1024), 512 >>> (a1d, a1d, a2d, 1024, 0);
						kernelFunc <<< 1, 512 >>> (a2d, a2d, cd, 1024, 0);
						kernelFunc2 <<< n/(1024*1024), 1024 >>> (a1d, a2d);
						kernelFunc2 <<< dimGrid, 1024 >>> (ad, a1d);
					}

					HANDLE_ERROR(cudaMemcpy(c, ad, n * sizeof(float), cudaMemcpyDeviceToHost));
			}

			else {
					int MaxMemSize = 1<<27;
					float OldSum=0;
					HANDLE_ERROR(cudaMalloc((void**)&ad, MaxMemSize * sizeof(float)));
		    	HANDLE_ERROR(cudaMalloc((void**)&a1d, (MaxMemSize/1024) * sizeof(float)));
		    	HANDLE_ERROR(cudaMalloc((void**)&a2d, (MaxMemSize/(1024*1024)) * sizeof(float)));

					int NumberofIter = n/MaxMemSize;
					for(int i=0; i<NumberofIter; i++) {
						HANDLE_ERROR(cudaMemcpy(ad, a + i*MaxMemSize, MaxMemSize * sizeof(float), cudaMemcpyHostToDevice));

						dim3 dimGrid(1<<15, MaxMemSize/(1<<25));
						kernelFunc <<< dimGrid, 512 >>> (ad, ad, a1d, 1024, 1);
						kernelFunc <<< MaxMemSize/(1024*1024), 512 >>> (a1d, a1d, a2d, 1024, 0);
						kernelFunc <<< 1, 512 >>> (a2d, a2d, cd, 1024, 0);
						kernelFunc2 <<< MaxMemSize/(1024*1024), 1024 >>> (a1d, a2d);
						kernelFunc2 <<< dimGrid, 1024 >>> (ad, a1d, OldSum);

						HANDLE_ERROR(cudaMemcpy(c + i*MaxMemSize, ad, MaxMemSize * sizeof(float), cudaMemcpyDeviceToHost));
						OldSum = c[(i+1)*MaxMemSize-1];
					}
			}

    	HANDLE_ERROR(cudaFree(ad));
			HANDLE_ERROR(cudaFree(a1d));
			HANDLE_ERROR(cudaFree(a2d));
			HANDLE_ERROR(cudaFree(cd));
}
