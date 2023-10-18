//ONLY MODIFY THIS FILE!
//YOU CAN MODIFY EVERYTHING IN THIS FILE!

#include "fft.h"

#define tx threadIdx.x
#define ty threadIdx.y
#define tz threadIdx.z

#define bx blockIdx.x
#define by blockIdx.y
#define bz blockIdx.z

// you may define other parameters here!
// you may define other macros here!
// you may define other functions here!

#define	R2 2
#define	R4 4
#define	R8 8


__global__ void KernelFunc1(float* x_r_d, float* x_i_d, float* X_r_d, float* X_i_d, const unsigned int N, const unsigned int M, const unsigned int Ns) 
{	
	int j = bx * blockDim.x + tx;
	int idxS = j;
	float v_r[R2]; float v_i[R2];
	float angle = -2*PI*(j%Ns) / (Ns*R2);
	for (int r=0; r<R2; r++)
	{	
		v_r[r] = (x_r_d[idxS+r*N/R2]*cos(r*angle) - x_i_d[idxS+r*N/R2]*sin(r*angle));
		v_i[r] = (x_r_d[idxS+r*N/R2]*sin(r*angle) + x_i_d[idxS+r*N/R2]*cos(r*angle));			
	}
	
	float v0_r = v_r[0]; float v0_i = v_i[0];
	v_r[0] = v0_r + v_r[1]; v_i[0] = v0_i + v_i[1];
	v_r[1] = v0_r - v_r[1]; v_i[1] = v0_i - v_i[1];
	
	int idxD = (j/Ns)*Ns*R2 + (j%Ns);
	
	for (int r=0; r<R2; r++)
	{
		X_r_d[idxD + r*Ns] = v_r[r];
		X_i_d[idxD + r*Ns] = v_i[r];
	}
}

__global__ void KernelFunc2(float* x_r_d, float* x_i_d, float* X_r_d, float* X_i_d, const unsigned int N, const unsigned int M, const unsigned int Ns) 
{	
	int j = bx * blockDim.x + tx;
	int idxS = j;
	float v_r[R4]; float v_i[R4];
	float angle = -2*PI*(j%Ns) / (Ns*R4);
	for (int r=0; r<R4; r++)
	{		
		v_r[r] = x_r_d[idxS+r*N/R4]*cos(r*angle) - x_i_d[idxS+r*N/R4]*sin(r*angle);
		v_i[r] = x_r_d[idxS+r*N/R4]*sin(r*angle) + x_i_d[idxS+r*N/R4]*cos(r*angle);			
	}

	float v0_r = v_r[0]; float v0_i = v_i[0];
	float v1_r = v_r[1]; float v1_i = v_i[1];
	float v2_r = v_r[2]; float v2_i = v_i[2];
	float v3_r = v_r[3]; float v3_i = v_i[3];
	v_r[0] = v0_r + v1_r + v2_r + v3_r; 	v_i[0] = v0_i + v1_i + v2_i + v3_i;
	v_r[1] = v0_r + v1_i - v2_r - v3_i;		v_i[1] = v0_i - v1_r - v2_i + v3_r;
	v_r[2] = v0_r - v1_r + v2_r - v3_r;		v_i[2] = v0_i - v1_i + v2_i - v3_i;
	v_r[3] = v0_r - v1_i - v2_r + v3_i;		v_i[3] = v0_i + v1_r - v2_i - v3_r;
	
	int idxD = (j/Ns)*Ns*R4 + (j%Ns);
	for (int r=0; r<R4 ; r++)
	{
		X_r_d[idxD + r*Ns] = v_r[r];
		X_i_d[idxD + r*Ns] = v_i[r];
	}
}

__global__ void KernelFunc3(float* x_r_d, float* x_i_d, float* X_r_d, float* X_i_d, const unsigned int N, const unsigned int M, const unsigned int Ns) 
{	
	int j = bx * blockDim.x + tx;
	int idxS = j;
	float v_r[R8]; float v_i[R8];
	float angle = -2*PI*(j%Ns) / (Ns*R8);
	for (int r=0; r<R8; r++)
	{		
		v_r[r] = x_r_d[idxS+r*N/R8]*cos(r*angle) - x_i_d[idxS+r*N/R8]*sin(r*angle);
		v_i[r] = x_r_d[idxS+r*N/R8]*sin(r*angle) + x_i_d[idxS+r*N/R8]*cos(r*angle);			
	}

	float a= 0.7071;
	float v0_r = v_r[0]; float v0_i = v_i[0];
	float v1_r = v_r[1]; float v1_i = v_i[1];
	float v2_r = v_r[2]; float v2_i = v_i[2];
	float v3_r = v_r[3]; float v3_i = v_i[3];
	float v4_r = v_r[4]; float v4_i = v_i[4];
	float v5_r = v_r[5]; float v5_i = v_i[5];
	float v6_r = v_r[6]; float v6_i = v_i[6];
	float v7_r = v_r[7]; float v7_i = v_i[7];
	v_r[0] = v0_r + v1_r + v2_r + v3_r + v4_r + v5_r + v6_r + v7_r;
	v_i[0] = v0_i + v1_i + v2_i + v3_i + v4_i + v5_i + v6_i + v7_i;
	
	v_r[1] = v0_r + a*(v1_r+v1_i) + v2_i + a*(v3_i-v3_r) + v4_r - a*(v5_r+v5_i) - v6_i + a*(v7_r-v7_i);
	v_i[1] = v0_i + a*(v1_i-v1_r) - v2_r - a*(v3_r+v3_i) + v4_i + a*(v5_r-v5_i) + v6_r + a*(v7_r+v7_i);
	
	v_r[2] = v0_r + v1_i + v2_r - v3_i + v4_r + v5_i - v6_r - v7_i;
	v_i[2] = v0_i - v1_r + v2_i + v3_r + v4_i - v5_r - v6_i + v7_r;
	
	v_r[3] = v0_r + a*(v1_i-v1_r) - v2_i + a*(v3_r+v3_i) + v4_r - a*(v5_i-v5_r) + v6_i - a*(v7_r+v7_i);
	v_i[3] = v0_i - a*(v1_r+v1_i) + v2_r + a*(v3_i-v3_r) + v4_i + a*(v5_r+v5_i) - v6_r + a*(v7_r-v7_i);
	
	v_r[4] = v0_r - v1_r + v2_r - v3_r + v4_r - v5_r + v6_r - v7_r;
	v_i[4] = v0_i - v1_i + v2_i - v3_i + v4_i - v5_i + v6_i - v7_i;
	
	v_r[5] = v0_r - a*(v1_r+v1_i) + v2_i - a*(v3_i-v3_r) + v4_r + a*(v5_r+v5_i) - v6_i - a*(v7_r-v7_i);
	v_i[5] = v0_i - a*(v1_i-v1_r) - v2_r + a*(v3_r+v3_i) + v4_i - a*(v5_r-v5_i) + v6_r - a*(v7_r+v7_i);
	
	v_r[6] = v0_r - v1_i + v2_r + v3_i + v4_r - v5_i - v6_r + v7_i;
	v_i[6] = v0_i + v1_r + v2_i - v3_r + v4_i + v5_r - v6_i - v7_r;
	
	v_r[7] = v0_r - a*(v1_i-v1_r) - v2_i - a*(v3_r+v3_i) + v4_r + a*(v5_i-v5_r) + v6_i + a*(v7_r+v7_i);
	v_i[7] = v0_i + a*(v1_r+v1_i) + v2_r - a*(v3_i-v3_r) + v4_i - a*(v5_r+v5_i) - v6_r - a*(v7_r-v7_i);
	
	int idxD = (j/Ns)*Ns*R8 + (j%Ns);
	for (int r=0; r<R8 ; r++)
	{
		X_r_d[idxD + r*Ns] = v_r[r];
		X_i_d[idxD + r*Ns] = v_i[r];
	}
}

__global__ void KernelCopy (float* x_r_d, float* x_i_d, float* X_r_d, float* X_i_d)
{	
	int t_Id = bx * blockDim.x + tx;
	x_r_d[t_Id] = X_r_d[t_Id];
	x_i_d[t_Id] = X_i_d[t_Id];	
}


//-----------------------------------------------------------------------------
void gpuKernel(float* x_r_d, float* x_i_d, /*float* X_r_d, float* X_i_d,*/ const unsigned int N, const unsigned int M)
{
	if(M>25)
	{
		gpuKernel(x_r_d,x_i_d,N/2,M-1);
		gpuKernel(&x_r_d[N/2],&x_i_d[N/2],N/2,M-1);
		
		return;
	}
	
	float* X_r_d;
	float* X_i_d;	
	HANDLE_ERROR(cudaMalloc((void**)&X_r_d, N * sizeof(float)));
	HANDLE_ERROR(cudaMalloc((void**)&X_i_d, N * sizeof(float)));
	
	int Ns;
	if ((M%4==0) && M>10)
	{
		for (Ns=1; Ns<N ; Ns*=R4)
		{	
			KernelFunc2 <<< N/(1024*R4), 1024 >>>(x_r_d, x_i_d, X_r_d, X_i_d, N, M, Ns);
			Ns = Ns*R4;
			KernelFunc2 <<< N/(1024*R4), 1024 >>>(X_r_d, X_i_d, x_r_d, x_i_d, N, M, Ns);
		}
	}
	else if(M==11)
	{
		for (Ns=1; Ns<N ; Ns*=R2)
			{							
				KernelFunc1 <<< 1, 1024 >>>(x_r_d, x_i_d, X_r_d, X_i_d, N, M, Ns);
				KernelCopy <<< 2, 1024 >>>(x_r_d, x_i_d, X_r_d, X_i_d);
			}
	}
	else if (M>10)
	{
		if (M%2)
		{
			for (Ns=1; Ns<(N/2) ; Ns*=R4)
			{							
				KernelFunc2 <<< N/(1024*R4), 1024 >>>(x_r_d, x_i_d, X_r_d, X_i_d, N, M, Ns);
				KernelCopy <<< (1<<M-10), 1024 >>>(x_r_d, x_i_d, X_r_d, X_i_d);
			}
			KernelFunc1 <<< N/(1024*R2), 1024 >>>(x_r_d, x_i_d, X_r_d, X_i_d, N, M, Ns);
			KernelCopy <<< (1<<M-10), 1024 >>>(x_r_d, x_i_d, X_r_d, X_i_d);
		}
		else
			for (Ns=1; Ns<N ; Ns*=R4)
			{							
				KernelFunc2 <<< N/(1024*R4), 1024 >>>(x_r_d, x_i_d, X_r_d, X_i_d, N, M, Ns);
				KernelCopy <<< (1<<M-10), 1024 >>>(x_r_d, x_i_d, X_r_d, X_i_d);
			}		
	}
	else 
	{
		if (M%2)
		{
			for (Ns=1; Ns<(N/2) ; Ns*=R4)
			{							
				KernelFunc2 <<< 1, N/R4 >>>(x_r_d, x_i_d, X_r_d, X_i_d, N, M, Ns);
				KernelCopy <<< 1, N >>>(x_r_d, x_i_d, X_r_d, X_i_d);
			}
			KernelFunc1 <<< 1, N/R2 >>>(x_r_d, x_i_d, X_r_d, X_i_d, N, M, Ns);
			KernelCopy <<< 1, N >>>(x_r_d, x_i_d, X_r_d, X_i_d);
		}
		else
			for (Ns=1; Ns<N ; Ns*=R4)
			{							
				KernelFunc2 <<< 1, N/R4 >>>(x_r_d, x_i_d, X_r_d, X_i_d, N, M, Ns);
				KernelCopy <<< 1, N >>>(x_r_d, x_i_d, X_r_d, X_i_d);
			}		
	}
	
    HANDLE_ERROR(cudaFree(X_r_d));
    HANDLE_ERROR(cudaFree(X_i_d));
}
