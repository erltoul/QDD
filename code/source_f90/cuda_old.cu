#include <stdio.h>
#include <math.h>
#include <cufft.h>
#include <math_functions.h> //math library for GPU
#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include "define_cuda.h"
//#include <cutil.h>

#define BATCH 1 //The number of batched ffts 

#define a0 0.458165293
#define da0 0.119086804

#define a1  2.2170586
#define da1 0.615740256

#define a2  0.740555173
#define da2 0.157420151

#define a3  0.019682278
#define da3 0.003532336

#define b1  1.000000000
#define db1 0.000000000

#define b2  4.504130959
#define db2 0.2361297

#define b3  1.1106363
#define db3 0.205200460

#define b4  0.023592917
#define db4 0.004200005

#define cr2 1.2599210498948731906665444
#define cr3 1.4422495703074083017725115
#define cr4 1.5874010519681993613971827

/*__device__ __constant__ double a0  = 0.458165293;
__device__ __constant__ double da0 = 0.119086804;

__device__ __constant__ double a1  = 2.2170586;
__device__ __constant__ double da1 = 0.615740256;

__device__ __constant__ double a2  = 0.740555173;
__device__ __constant__ double da2 = 0.157420151;

__device__ __constant__ double a3  = 0.019682278;
__device__ __constant__ double da3 = 0.003532336;

__device__ __constant__ double b1  = 1.000000000;
__device__ __constant__ double db1 = 0.000000000;

__device__ __constant__ double b2  = 4.504130959;
__device__ __constant__ double db2 = 0.2361297;

__device__ __constant__ double b3  = 1.1106363;
__device__ __constant__ double db3 = 0.205200460;

__device__ __constant__ double b4  = 0.023592917;
__device__ __constant__ double db4 = 0.004200005;

__device__ __constant__ double cr2 = 1.2599210498948731906665444;
__device__ __constant__ double cr3 = 1.4422495703074083017725115;
__device__ __constant__ double cr4 = 1.5874010519681993613971827;*/

//#define pi 3.141592653589793; //PI is already defined in math.h

//typedef float2 cuDoubleComplex;

extern "C" void cuda_plan_1d_(cufftHandle *plan, int *dim, int *batch)
{
        cufftPlan1d(plan,*dim,CUFFT_Z2Z,*batch);
}

extern "C" void cuda_plan_3d_(cufftHandle *plan, int *n1, int *n2, int *n3)
{
        cufftPlan3d(plan,*n1,*n2,*n3,CUFFT_Z2Z);
}

extern "C" void kill_plan_(cufftHandle * plan)
{
        cufftDestroy(*plan);
}

extern "C" void run_fft_for_(cufftHandle *plan, cufftDoubleComplex *a, cufftDoubleComplex *b, int *Np)
{
        int N = *Np;
        //Allocate device memory
        cufftDoubleComplex *d_signal;

        cudaMalloc((void**)&d_signal,sizeof(cufftDoubleComplex)*N*BATCH);
        cudaMemcpy(d_signal,a,sizeof(cufftDoubleComplex)*N*BATCH,cudaMemcpyHostToDevice);
        cufftExecZ2Z(*plan,(cufftDoubleComplex *)d_signal,(cufftDoubleComplex *)d_signal, CUFFT_FORWARD); //forward in-place fft
        cudaMemcpy(b, d_signal, sizeof(cufftDoubleComplex)*N*BATCH, cudaMemcpyDeviceToHost);
        cudaFree(d_signal);
}

extern "C" void run_fft_back_(cufftHandle *plan, cufftDoubleComplex *a, cufftDoubleComplex *b,int *Np)
{
        int N = *Np;
        //Allocate device memory
        cufftDoubleComplex *d_signal;

        cudaMalloc((void**)&d_signal,sizeof(cufftDoubleComplex)*N*BATCH);
        cudaMemcpy(d_signal,a,sizeof(cufftDoubleComplex)*N*BATCH,cudaMemcpyHostToDevice);
        cufftExecZ2Z(*plan, (cufftDoubleComplex *)d_signal,(cufftDoubleComplex *)d_signal, CUFFT_INVERSE);
        cudaMemcpy(b, d_signal, sizeof(cufftDoubleComplex)*N*BATCH, cudaMemcpyDeviceToHost);
        cudaFree(d_signal); //free the signal on the device.
}

extern "C" void run_fft_for3d_(cufftHandle *plan, cufftDoubleComplex *in, cufftDoubleComplex *out,int *Npx,int *Npy,int *Npz,int *ind)
{
        int Nx = *Npx;
        int Ny = *Npy;
        int Nz = *Npz;
        //Allocate device memory
        cufftDoubleComplex *d_ffta;

        cudaMalloc((void**)&d_ffta,sizeof(cufftDoubleComplex)*Nx*Ny*Nz);
        cudaMemcpy(d_ffta,in,sizeof(cufftDoubleComplex)*Nx*Ny*Nz,cudaMemcpyHostToDevice);
        cufftExecZ2Z(*plan, (cufftDoubleComplex *)d_ffta,(cufftDoubleComplex *)d_ffta, CUFFT_FORWARD); //forward in-place fft
        cudaMemcpy(out,d_ffta, sizeof(cufftDoubleComplex)*Nx*Ny*Nz,cudaMemcpyDeviceToHost);
        cudaFree(d_ffta);
}

extern "C" void run_fft_back3d_(cufftHandle *plan, cufftDoubleComplex *a, cufftDoubleComplex *b,int *Npx,int *Npy,int *Npz)
{
        int Nx = *Npx;
        int Ny = *Npy;
        int Nz = *Npz;

        //Allocate device memory
        cufftDoubleComplex *d_signal;
        
        cudaMalloc((void**)&d_signal,sizeof(cufftDoubleComplex)*Nx*Ny*Nz);
        cudaMemcpy(d_signal,a,sizeof(cufftDoubleComplex)*Nx*Ny*Nz,cudaMemcpyHostToDevice);
        cufftExecZ2Z(*plan, (cufftDoubleComplex *)d_signal,(cufftDoubleComplex *)d_signal, CUFFT_INVERSE);
        cudaMemcpy(b, d_signal, sizeof(cufftDoubleComplex)*Nx*Ny*Nz,cudaMemcpyDeviceToHost);
        cudaFree(d_signal); //free the signal on the device.
}

/*#if(directenergy)
__global__ void lda(double *d_rho,double *d_chpdft,double *d_ec,int nxyz,double e2, double enerpw)
#else
__global__ void lda(double *d_rho,double *d_chpdft,double *d_ec,unsigned int *d_ind,int nxyz,double e2)
#endif
{
       unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

       //double t1,t2,t3,t4,t5,t6,t7,t8,t10,t11,t12,t13,t15,t17,t22,t23,t24,t25,t26,t28,t29,t34,t35,t36;
       //double t37,t42,t44,t48,t53,t58,t63,t64,t65,t,t68,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;
       double t36,t48,t53,t58,t63,t64,t,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;
       //double rp,xi;
       if (ind<nxyz)
       {
         d_ind[ind]=ind;
         //rp = max(d_rho[ind],1.0e-16);
         //xi = d_rho[ind+nxyz];
         
         //t2 = max(d_rho[ind],1.0e-16);
         //t1 = d_rho[ind+nxyz]*t2;
         //t3 = 1.0/t2;
         //t4 = d_rho[ind+nxyz];
         //t6 = pow((1.0+t4),(1.0/3.0));
         //t6 = pow((1.0+d_rho[ind+nxyz]),(1.0/3.0));
         //t7 = t6*t6;
         //t8 = t7*t7;
         //t10 = pow((1.0-t4),(1.0/3.0));
         //t10 = pow((1.0-d_rho[ind+nxyz]),(1.0/3.0));
         //t11 = t10*t10;
         //t12 = t11*t11;
         //t13 = t8+t12-2.0;
         //t13 = (t6*t6*t6*t6+t10*t10*t10*t10-2.0);
         //t15 = cr2;
         //t17 = 1.0/(2.0*cr2-2.0);
         //t22 = cr3;
         //t23 = ((a1+da1*(t6*t6*t6*t6+t10*t10*t10*t10-2.0)*1.0/(2.0*cr2-2.0))*cr3);
         //t24 = cr4;
         //t25 = cr4*cr4;
         //t26 = 1.0/M_PI;
         //t28 = pow((M_1_PI/t2),(1.0/3.0));
         //t29 = cr4*cr4*pow((M_1_PI/t2),(1.0/3.0));
         //t34 = cr3*cr3;
         //t35 = (a2+da2*(t6*t6*t6*t6+t10*t10*t10*t10-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3;
         t36 = pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0));
         //t37 = cr4*t36;
         //t42 = (a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI;
         //t44 = (a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16));
         t48 = (b1+db1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3;
         t53 = (b2+db2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3;
         t58 = (b3+db3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI;
         t63 = (b4+db4*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3;
         t64 = pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0));
         t = t48*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+t53*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*t58/max(d_rho[ind],1.0e-16)+3.0/16.0*t63*cr4*cr4*t64;
         //t68 = 1.0/t;
         t70 = max(d_rho[ind],1.0e-16)*max(d_rho[ind],1.0e-16);
         t71 = 1.0/t70;
         //t72 = t1*t71;
         t72 = d_rho[ind+nxyz]*max(d_rho[ind],1.0e-16)*t71;
         t77 = 4.0/3.0*pow((1.0+d_rho[ind+nxyz]),(1.0/3.0))*(1.0/max(d_rho[ind],1.0e-16)-t72)+4.0/3.0*pow((1.0-d_rho[ind+nxyz]),(1.0/3.0))*(-1.0/max(d_rho[ind],1.0e-16)+t72);
         t82 = cr3*cr4*cr4;
         t83 = t82*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0));
         t88 = 1.0/t36*M_1_PI*t71;
         t93 = cr3*cr3*cr4*t36;
         t98 = 1.0/pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*M_1_PI*t71;
         t102 = 1.0/(2.0*cr2-2.0)*M_1_PI/max(d_rho[ind],1.0e-16);
         t109 = t*t;
         t135 = (a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16))*1.0/t+max(d_rho[ind],1.0e-16)*(da0*t77*1.0/(2.0*cr2-2.0)+da1*t77*1.0/(2.0*cr2-2.0)*t83/4.0-((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*t88/12.0+da2*t77*1.0/(2.0*cr2-2.0)*t93/4.0-(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*t98/6.0+3.0/4.0*da3*t77*t102-3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI*t71)*1.0/t-max(d_rho[ind],1.0e-16)*(a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16))/t109*(db1*t77*1.0/(2.0*cr2-2.0)*t83/4.0-t48*cr4*cr4*t88/12.0+db2*t77*1.0/(2.0*cr2-2.0)*t93/4.0-t53*cr4*t98/6.0+3.0/4.0*db3*t77*t102-3.0/4.0*t58*t71+3.0/16.0*db4*t77*1.0/(2.0*cr2-2.0)*t82*t64-t63*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*M_1_PI*t71/4.0);

         d_chpdft[ind] = -t135*e2;
  
         t77 = (4.0/3.0*pow((1.0+d_rho[ind+nxyz]),(1.0/3.0))*(-1.0/max(d_rho[ind],1.0e-16)-t72))+(4.0/3.0*pow((1.0-d_rho[ind+nxyz]),(1.0/3.0))*(1.0/max(d_rho[ind],1.0e-16)+t72));
         //t82 = cr3*cr4*cr4;
         //t83 = t82*t28;
         //t88 = 1.0/t36*M_1_PI*t71;
         //t93 = cr3*cr3*cr4*t36;
         //t98 = 1.0/t28*M_1_PI*t71;
         //t102 = t17*M_1_PI/t2;
         //t109 = t*t;
         t135 = (a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16))*1.0/t+max(d_rho[ind],1.0e-16)*(da0*t77*1.0/(2.0*cr2-2.0)+da1*t77*1.0/(2.0*cr2-2.0)*t83/4.0-((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*t88/12.0+da2*t77*1.0/(2.0*cr2-2.0)*t93/4.0-(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*t98/6.0+3.0/4.0*da3*t77*t102-3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI*t71)*1.0/t-max(d_rho[ind],1.0e-16)*(a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16))/t109*(db1*t77*1.0/(2.0*cr2-2.0)*t83/4.0-t48*cr4*cr4*t88/12.0+db2*t77*1.0/(2.0*cr2-2.0)*t93/4.0-t53*cr4*t98/6.0+3.0/4.0*db3*t77*t102-3.0/4.0*t58*t71+3.0/16.0*db4*t77*1.0/(2.0*cr2-2.0)*t82*t64-t63*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*M_1_PI*t71/4.0);

         d_chpdft[ind+nxyz] = -t135*e2;
  
         //t1= t2;
         //t4 = d_rho[ind+nxyz];
         //t6 = pow((1.0+t4),(1.0/3.0));
         //t7 = t6*t6;
         //t8 = t7*t7;
         //t10 = pow((1.0-t4),(1.0/3.0));
         //t11 = t10*t10;
         //t12 = t11*t11;
         //t13 = t8+t12-2.0;
         //t15 = pow(2.0,(1.0/3.0));
         //t17 = 1.0/(2.0*t15-2.0);
         //t22 = pow(3.0,(1.0/3.0));
         //t24 = pow(4.0,(1.0/3.0));
         //t25 = t24*t24;
         //t26 = 1.0/M_PI;
         //t28 = pow((M_1_PI*t3),(1.0/3.0));
         //t29 = t25*t28;
         //t34 = t22*t22;
         //t36 = t28*t28;
         //t37 = t24*t36;
         //t65 = t36*t36;
         t70 = max(d_rho[ind],1.0e-16)*(a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+(a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*t36/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16))/((b1+db1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(b2+db2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*t36/4.0+3.0/4.0*(b3+db3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16)+3.0/16.0*(b4+db4*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr4*cr4*t36*t36);
         //t70 = t1*(a0+da0*t13*t17+(a1+da1*t13*t17)*t22*t29/4.0+(a2+da2*t13*t17)*t34*t37/4.0+3.0/4.0*(a3+da3*t13*t17)*t26*t3)/((b1+db1*t13*t17)*t22*t29/4.0+(b2+db2*t13*t17)*t34*t37/4.0+3.0/4.0*(b3+db3*t13*t17)*t26*t3+3.0/16.0*(b4+db4*t13*t17)*t22*t25*t65);
  
         //t1 = d_rho[ind+nxyz]*t2;
         //t2 = max(d_rho[ind],1.0e-16);
         //t3=fabs((t1+t2)/2.0);  //  *e2
         //t4=fabs((t1-t2)/2.0);  //  *e2
         //t3=fabs(t2*(d_rho[ind+nxyz]+1)/2.0);  //  *e2
         //t4=fabs(t2*(d_rho[ind+nxyz]-1)/2.0);  //  *e2
         //t5= d_chpdft[ind]*t3+d_chpdft[ind+nxyz]*t4;
         //t5= (d_chpdft[ind]*fabs(max(d_rho[ind],1.0e-16)*(d_rho[ind+nxyz]+1)/2.0)+d_chpdft[ind+nxyz]*fabs(max(d_rho[ind],1.0e-16)*(d_rho[ind+nxyz]-1)/2.0));

         d_ec[ind]=-t70*e2-0.5*(d_chpdft[ind]*fabs(max(d_rho[ind],1.0e-16)*(d_rho[ind+nxyz]+1)/2.0)+d_chpdft[ind+nxyz]*fabs(max(d_rho[ind],1.0e-16)*(d_rho[ind+nxyz]-1)/2.0));
  
#if(directenergy)
         enerpw = -t70*e2 + enerpw;
#endif

      }
      __syncthreads();
}*/

//#if(directenergy)
//__global__ void lda(double *d_rho,double *d_chpdft,double *d_ec,int nxyz,double e2, double enerpw)
//#else
//__global__ void lda(double *d_rho,double *d_chpdft,double *d_ec,unsigned int *d_ind,int nxyz,double e2)
__global__ void lda(int *d_ind,int nxyz)
//#endif
{
       int ind = blockIdx.x*blockDim.x+threadIdx.x;

       /*double t1,t2,t3,t4,t5,t6,t7,t8,t10,t11,t12,t13,t17,t23,t25,t28,t29,t34,t35,t36,t37,t42,t44;
       double t48,t53,t58,t63,t64,t65,t,t68,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;
       double rp,xi;*/

       if (ind<nxyz)
       {
         d_ind[ind]=ind;
         /*rp = max(d_rho[ind],1e-16);
         xi = d_rho[ind+nxyz];
  
         t1 = xi*rp;
         t2 = rp;
         t3 = 1.0/t2;
         t4 = xi;
         t6 = pow((1.0+t4),(1.0/3.0));
         t7 = t6*t6;
         t8 = t7*t7;
         t10 = pow((1.0-t4),(1.0/3.0));
         t11 = t10*t10;
         t12 = t11*t10;
         t13 = t8+t12-2.0;
         //t15 = 2.0**(1.0/3.0);
         t17 = 1.0/(2.0*cr2-2.0);
         //t22 = 3.0**(1.0/3.0);
         t23 = (a1+da1*t13*t17)*cr3;
         //t24 = 4.0**(1.0/3.0);
         t25 = cr4*cr4;
         //t26 = 1/pi;
         t28 = pow((M_1_PI/t2),(1.0/3.0));
         t29 = t25*t28;
         t34 = cr3*cr3;
         t35 = (a2+da2*t13*t17)*t34;
         t36 = t28*t28;
         t37 = cr4*t36;
         t42 = (a3+da3*t13*t17)*M_1_PI;
         t44 = a0+da0*t13*t17+t23*t29/4.0+t35*t37/4+3.0/4.0*t42*t3;
         t48 = (b1+db1*t13*t17)*cr3;
         t53 = (b2+db2*t13*t17)*t34;
         t58 = (b3+db3*t13*t17)*M_1_PI;
         t63 = (b4+db4*t13*t17)*cr3;
         t64 = t36*t36;
         t = t48*t29/4.0+t53*t37/4.0+3.0/4.0*t58*t3+3.0/16.0*t63*t25*t64;
         t68 = 1.0/t;
         t70 = t2*t2;
         t71 = 1.0/t70;
         t72 = t1*t71;
         t77 = 4.0/3.0*t6*(t3-t72)+4.0/3.0*t10*(-t3+t72);
         t82 = cr3*t25;
         t83 = t82*t28;
         t88 = 1.0/t36*M_1_PI*t71;
         t93 = t34*cr4*t36;
         t98 = 1.0/t28*M_1_PI*t71;
         t102 = t17*M_1_PI*t3;
         t109 = t*t;
         t135 = t44*t68+t2*(da0*t77*t17+da1*t77*t17*t83/4-t23*t25*t88/12+da2*t77*t17*t93/4-t35*cr4*t98/6+3.0/4.0*da3*t77*t102-3.0/4.0*t42*t71)*t68-t2*t44/t109*(db1*t77*t17*t83/4-t48*t25*t88/12+db2*t77*t17*t93/4-t53*cr4*t98/6+3.0/4.0*db3*t77*t102-3.0/4.0*t58*t71+3.0/16.0*db4*t77*t17*t82*t64-t63*t25*t28*M_1_PI*t71/4);

         d_chpdft[ind] = -t135  * e2;

         t77 = 4.0/3.0*t6*(-t3-t72)+4.0/3.0*t10*(t3+t72);
         t82 = cr3*t25;
         t83 = t82*t28;
         t88 = 1.0/t36*M_1_PI*t71;
         t93 = t34*cr4*t36;
         t98 = 1.0/t28*M_1_PI*t71;
         t102 = t17*M_1_PI*t3;
         t109 = t*t;
         t135 = t44*t68+t2*(da0*t77*t17+da1*t77*t17*t83/4-t23*t25*t88/12+da2*t77*t17*t93/4-t35*cr4*t98/6+3.0/4.0*da3*t77*t102-3.0/4.0*t42*t71)*t68-t2*t44/t109*(db1*t77*t17*t83/4-t48*t25*t88/12+db2*t77*t17*t93/4-t53*cr4*t98/6+3.0/4.0*db3*t77*t102-3.0/4.0*t58*t71+3.0/16.0*db4*t77*t17*t82*t64-t63*t25*t28*M_1_PI*t71/4);

         d_chpdft[ind+nxyz] = -t135  * e2;

         t1=rp;
         t4 = xi;
         t6 = pow((1.0+t4),(1.0/3.0));
         t7 = t6*t6;
         t8 = t7*t7;
         t6 = pow((1.0-t4),(1.0/3.0));;
         t11 = t10*t10;
         t12 = t11*t11;
         t13 = t8+t12-2.0;
         //t15 = 2**(1.0/3.0);
         t17 = 1.0/(2.0*cr2-2.0);
         //t22 = 3**(1.0/3.0);
         //t24 = 4**(1.0/3.0);
         t25 = cr4*cr4;
         //t26 = 1/pi;
         t28 = pow((M_1_PI/t2),(1.0/3.0));;
         t29 = t25*t28;
         t34 = cr3*cr3;
         t36 = t28*t28;
         t37 = cr4*t36;
         t65 = t36*t36;
         t70 = t1*(a0+da0*t13*t17+(a1+da1*t13*t17)*cr3*t29/4+(a2+da2*t13*t17)*t34*t37/4+3.0/4.0*(a3+da3*t13*t17)*M_1_PI*t3)/((b1+db1*t13*t17)*cr3*t29/4+(b2+db2*t13*t17)*t34*t37/4+3.0/4.0*(b3+db3*t13*t17)*M_1_PI*t3+3.0/16.0*(b4+db4*t13*t17)*cr3*t25*t65);

         t1 = xi*rp;
         t2 = rp;
         t3=fabs((t1+t2)/2.0);   //  *e2;
         t4=fabs((t1-t2)/2.0);   //  *e2;
         t5= d_chpdft[ind]*t3+d_chpdft[ind+nxyz]*t4;
         
         d_ec[ind] = (-t70*e2 - 0.5*t5);*/
//#if(directenergy)
//         enerpw = -t70*e2 + enerpw;
//#endif

      }
}

/*void ldac(double *d_rho,double *d_chpdft,double *d_ec,unsigned int *d_ind,int nxyz,double e2)
{
       //double t1,t2,t3,t4,t5,t6,t7,t8,t10,t11,t12,t13,t15,t17,t22,t23,t24,t25,t26,t28,t29,t34,t35,t36;
       //double t37,t42,t44,t48,t53,t58,t63,t64,t65,t,t68,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;
       double t36,t48,t53,t58,t63,t64,t,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;
       //double rp,xi;
       
       for (int ind=0;ind<nxyz;ind++)
       {
         d_ind[ind]=ind;
         //rp = max(d_rho[ind],1.0e-16);
         //xi = d_rho[ind+nxyz];
         
         //t2 = max(d_rho[ind],1.0e-16);
         //t1 = d_rho[ind+nxyz]*t2;
         //t3 = 1.0/t2;
         //t4 = d_rho[ind+nxyz];
         //t6 = pow((1.0+t4),(1.0/3.0));
         //t6 = pow((1.0+d_rho[ind+nxyz]),(1.0/3.0));
         //t7 = t6*t6;
         //t8 = t7*t7;
         //t10 = pow((1.0-t4),(1.0/3.0));
         //t10 = pow((1.0-d_rho[ind+nxyz]),(1.0/3.0));
         //t11 = t10*t10;
         //t12 = t11*t11;
         //t13 = t8+t12-2.0;
         //t13 = (t6*t6*t6*t6+t10*t10*t10*t10-2.0);
         //t15 = cr2;
         //t17 = 1.0/(2.0*cr2-2.0);
         //t22 = cr3;
         //t23 = ((a1+da1*(t6*t6*t6*t6+t10*t10*t10*t10-2.0)*1.0/(2.0*cr2-2.0))*cr3);
         //t24 = cr4;
         //t25 = cr4*cr4;
         //t26 = 1.0/M_PI;
         //t28 = pow((M_1_PI/t2),(1.0/3.0));
         //t29 = cr4*cr4*pow((M_1_PI/t2),(1.0/3.0));
         //t34 = cr3*cr3;
         //t35 = (a2+da2*(t6*t6*t6*t6+t10*t10*t10*t10-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3;
         t36 = pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0));
         //t37 = cr4*t36;
         //t42 = (a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI;
         //t44 = (a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16));
         t48 = (b1+db1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3;
         t53 = (b2+db2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3;
         t58 = (b3+db3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI;
         t63 = (b4+db4*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3;
         t64 = pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0));
         t = t48*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+t53*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*t58/max(d_rho[ind],1.0e-16)+3.0/16.0*t63*cr4*cr4*t64;
         //t68 = 1.0/t;
         t70 = max(d_rho[ind],1.0e-16)*max(d_rho[ind],1.0e-16);
         t71 = 1.0/t70;
         //t72 = t1*t71;
         t72 = d_rho[ind+nxyz]*max(d_rho[ind],1.0e-16)*t71;
         t77 = 4.0/3.0*pow((1.0+d_rho[ind+nxyz]),(1.0/3.0))*(1.0/max(d_rho[ind],1.0e-16)-t72)+4.0/3.0*pow((1.0-d_rho[ind+nxyz]),(1.0/3.0))*(-1.0/max(d_rho[ind],1.0e-16)+t72);
         t82 = cr3*cr4*cr4;
         t83 = t82*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0));
         t88 = 1.0/t36*M_1_PI*t71;
         t93 = cr3*cr3*cr4*t36;
         t98 = 1.0/pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*M_1_PI*t71;
         t102 = 1.0/(2.0*cr2-2.0)*M_1_PI/max(d_rho[ind],1.0e-16);
         t109 = t*t;
         t135 = (a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16))*1.0/t+max(d_rho[ind],1.0e-16)*(da0*t77*1.0/(2.0*cr2-2.0)+da1*t77*1.0/(2.0*cr2-2.0)*t83/4.0-((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*t88/12.0+da2*t77*1.0/(2.0*cr2-2.0)*t93/4.0-(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*t98/6.0+3.0/4.0*da3*t77*t102-3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI*t71)*1.0/t-max(d_rho[ind],1.0e-16)*(a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16))/t109*(db1*t77*1.0/(2.0*cr2-2.0)*t83/4.0-t48*cr4*cr4*t88/12.0+db2*t77*1.0/(2.0*cr2-2.0)*t93/4.0-t53*cr4*t98/6.0+3.0/4.0*db3*t77*t102-3.0/4.0*t58*t71+3.0/16.0*db4*t77*1.0/(2.0*cr2-2.0)*t82*t64-t63*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*M_1_PI*t71/4.0);

         d_chpdft[ind] = -t135*e2;
  
         t77 = (4.0/3.0*pow((1.0+d_rho[ind+nxyz]),(1.0/3.0))*(-1.0/max(d_rho[ind],1.0e-16)-t72))+(4.0/3.0*pow((1.0-d_rho[ind+nxyz]),(1.0/3.0))*(1.0/max(d_rho[ind],1.0e-16)+t72));
         //t82 = cr3*cr4*cr4;
         //t83 = t82*t28;
         //t88 = 1.0/t36*M_1_PI*t71;
         //t93 = cr3*cr3*cr4*t36;
         //t98 = 1.0/t28*M_1_PI*t71;
         //t102 = t17*M_1_PI/t2;
         //t109 = t*t;
         t135 = (a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16))*1.0/t+max(d_rho[ind],1.0e-16)*(da0*t77*1.0/(2.0*cr2-2.0)+da1*t77*1.0/(2.0*cr2-2.0)*t83/4.0-((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*t88/12.0+da2*t77*1.0/(2.0*cr2-2.0)*t93/4.0-(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*t98/6.0+3.0/4.0*da3*t77*t102-3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI*t71)*1.0/t-max(d_rho[ind],1.0e-16)*(a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+((a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3)*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16))/t109*(db1*t77*1.0/(2.0*cr2-2.0)*t83/4.0-t48*cr4*cr4*t88/12.0+db2*t77*1.0/(2.0*cr2-2.0)*t93/4.0-t53*cr4*t98/6.0+3.0/4.0*db3*t77*t102-3.0/4.0*t58*t71+3.0/16.0*db4*t77*1.0/(2.0*cr2-2.0)*t82*t64-t63*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))*M_1_PI*t71/4.0);

         d_chpdft[ind+nxyz] = -t135*e2;
  
         //t1= t2;
         //t4 = d_rho[ind+nxyz];
         //t6 = pow((1.0+t4),(1.0/3.0));
         //t7 = t6*t6;
         //t8 = t7*t7;
         //t10 = pow((1.0-t4),(1.0/3.0));
         //t11 = t10*t10;
         //t12 = t11*t11;
         //t13 = t8+t12-2.0;
         //t15 = pow(2.0,(1.0/3.0));
         //t17 = 1.0/(2.0*t15-2.0);
         //t22 = pow(3.0,(1.0/3.0));
         //t24 = pow(4.0,(1.0/3.0));
         //t25 = t24*t24;
         //t26 = 1.0/M_PI;
         //t28 = pow((M_1_PI*t3),(1.0/3.0));
         //t29 = t25*t28;
         //t34 = t22*t22;
         //t36 = t28*t28;
         //t37 = t24*t36;
         //t65 = t36*t36;
         t70 = max(d_rho[ind],1.0e-16)*(a0+da0*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0)+(a1+da1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(a2+da2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*t36/4.0+3.0/4.0*(a3+da3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16))/((b1+db1*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr4*cr4*pow((M_1_PI/max(d_rho[ind],1.0e-16)),(1.0/3.0))/4.0+(b2+db2*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr3*cr4*t36/4.0+3.0/4.0*(b3+db3*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*M_1_PI/max(d_rho[ind],1.0e-16)+3.0/16.0*(b4+db4*(pow(pow((1.0+d_rho[ind+nxyz]),(1.0/3.0)),4.0)+pow(pow((1.0-d_rho[ind+nxyz]),(1.0/3.0)),4.0)-2.0)*1.0/(2.0*cr2-2.0))*cr3*cr4*cr4*t36*t36);
         //t70 = t1*(a0+da0*t13*t17+(a1+da1*t13*t17)*t22*t29/4.0+(a2+da2*t13*t17)*t34*t37/4.0+3.0/4.0*(a3+da3*t13*t17)*t26*t3)/((b1+db1*t13*t17)*t22*t29/4.0+(b2+db2*t13*t17)*t34*t37/4.0+3.0/4.0*(b3+db3*t13*t17)*t26*t3+3.0/16.0*(b4+db4*t13*t17)*t22*t25*t65);
  
         //t1 = d_rho[ind+nxyz]*t2;
         //t2 = max(d_rho[ind],1.0e-16);
         //t3=fabs((t1+t2)/2.0);  //  *e2
         //t4=fabs((t1-t2)/2.0);  //  *e2
         //t3=fabs(t2*(d_rho[ind+nxyz]+1)/2.0);  //  *e2
         //t4=fabs(t2*(d_rho[ind+nxyz]-1)/2.0);  //  *e2
         //t5= d_chpdft[ind]*t3+d_chpdft[ind+nxyz]*t4;
         //t5= (d_chpdft[ind]*fabs(max(d_rho[ind],1.0e-16)*(d_rho[ind+nxyz]+1)/2.0)+d_chpdft[ind+nxyz]*fabs(max(d_rho[ind],1.0e-16)*(d_rho[ind+nxyz]-1)/2.0));

         d_ec[ind]=-t70*e2-0.5*(d_chpdft[ind]*fabs(max(d_rho[ind],1.0e-16)*(d_rho[ind+nxyz]+1)/2.0)+d_chpdft[ind+nxyz]*fabs(max(d_rho[ind],1.0e-16)*(d_rho[ind+nxyz]-1)/2.0));
      }
}*/

/*__global__ void init(double *d_rho,double *d_chpdft,double *d_ec,unsigned int *d_ind,int nxyz,double e2)
{
       unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;
       
       if(ind<nxyz){
         d_ind[ind]=ind;
         d_chpdft[ind]=ind*e2;
         d_chpdft[ind+nxyz]=(ind+nxyz)*e2*2.0;
         d_ec[ind]=d_ind[ind]*d_chpdft[ind]*d_rho[ind];

       }
      __syncthreads();
}*/

/*__global__ void reduce(double *d_ec,double *ec)
{
      // each thread loads one element from global to shared mem
      unsigned int tid = threadIdx.x;
      unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
      extern __shared__ double sdata[];

      sdata[tid] = d_ec[i];
      __syncthreads();

      // do reduction in shared mem
      for (unsigned int s=1; s < blockDim.x; s *= 2) {
        int index = 2 * s * tid;
        if (index < blockDim.x) {
          sdata[index] += sdata[index + s];
        }
        __syncthreads();
      }

      // write result for this block to global mem
      if (tid == 0) ec[blockIdx.x] = sdata[0];
}*/

//#if(directenergy)
//extern "C" void calc_lda_gpu_(double *rho,double *chpdft,int *N,double *e2,double *ec,double *enpw)
//#else
//extern "C" void calc_lda_gpu_(double *rho,double *chpdft,int *N,double *e2,double *ec)
extern "C" void calc_lda_gpu_(int *N)
//#endif
{
        int nxyz =*N;

        /*double *d_rho,*d_chpdft,*d_ec;
        double e=*e2;
#if(directenergy)
        double enerpw = *enpw;
#endif*/
        int blocksize=1024;
        dim3 dimgrid(ceil(nxyz/(float)blocksize),1,1);
        dim3 dimblock(blocksize,1,1);
        //double * h_ec;
        int * d_ind;
        int * h_ind;

        //h_ec=(double*)calloc(nxyz,sizeof(double));
        h_ind=(int*)calloc(nxyz,sizeof(int));
        /*cudaMalloc((void**)&d_rho,sizeof(double)*(nxyz*2));
        cudaMalloc((void**)&d_chpdft,sizeof(double)*(nxyz*2));
        cudaMalloc((void**)&d_ec,sizeof(double)*nxyz);*/
        cudaMalloc((void**)&d_ind,sizeof(int)*nxyz);
        /*cudaMemcpy(d_rho,rho,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
        cudaMemcpy(d_chpdft,chpdft,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
        cudaMemcpy(d_ec,h_ec,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);*/
        cudaMemcpy(d_ind,h_ind,sizeof(int)*(nxyz),cudaMemcpyHostToDevice);

#if(directenergy)
        lda<<<dimgrid,dimblock>>>(d_rho,d_chpdft,ec,nxyz,e,enerpw);
#else
        //init<<<dimgrid,dimblock>>>(d_rho,d_chpdft,d_ec,d_ind,nxyz,e);
        //lda<<<dimgrid,dimblock>>>(d_rho,d_chpdft,d_ec,d_ind,nxyz,e);
        lda<<<dimgrid,dimblock>>>(d_ind,nxyz);
#endif
        //cudaMemcpy(h_ec,d_ec, sizeof(double)*(nxyz),cudaMemcpyDeviceToHost);
        cudaMemcpy(h_ind,d_ind, sizeof(unsigned int)*(nxyz),cudaMemcpyDeviceToHost);
        /*cudaMemcpy(chpdft,d_chpdft, sizeof(double)*(nxyz*2),cudaMemcpyDeviceToHost);
        cudaMemcpy(rho,d_rho, sizeof(double)*(nxyz*2),cudaMemcpyDeviceToHost);
        cudaFree(d_rho);
        cudaFree(d_chpdft);
        cudaFree(d_ec);*/

        printf("%u \n",h_ind[0]);
        /*printf("%u %1.16e %1.16e %1.16e %1.16e %1.16e\n",h_ind[0],h_ec[0],rho[0],rho[nxyz],chpdft[0],chpdft[nxyz]);
        printf("%u %1.16e %1.16e %1.16e %1.16e %1.16e\n",h_ind[1],h_ec[1],rho[1],rho[1+nxyz],chpdft[1],chpdft[1+nxyz]);
        printf("%u %1.16e %1.16e %1.16e %1.16e %1.16e\n",h_ind[2],h_ec[2],rho[2],rho[2+nxyz],chpdft[2],chpdft[2+nxyz]);

        for (int ii=0;ii<nxyz;ii++){
          *ec+=h_ec[ii];
        }*/
        /*d_rho=(double *)malloc(nxyz*sizeof(double));
        d_chpdft=(double *)malloc(nxyz*sizeof(double));
        d_ec=(double *)malloc(nxyz*sizeof(double));
        d_ind=(unsigned int *)malloc(nxyz*sizeof(unsigned int));
        d_rho=rho;
        d_chpdft=chpdft;
        ldac(d_rho,d_chpdft,d_ec,d_ind,nxyz,e);
        printf("%u %1.16e %1.16e %1.16e %1.16e %1.16e\n",d_ind[0],d_ec[0],d_rho[0],d_rho[nxyz],d_chpdft[0],d_chpdft[nxyz]);
        printf("%u %1.16e %1.16e %1.16e %1.16e %1.16e\n",d_ind[1],d_ec[1],d_rho[1],d_rho[1+nxyz],d_chpdft[1],d_chpdft[1+nxyz]);
        printf("%u %1.16e %1.16e %1.16e %1.16e %1.16e\n",d_ind[2],d_ec[2],d_rho[2],d_rho[2+nxyz],d_chpdft[2],d_chpdft[2+nxyz]);
        for (int ii=0;ii<nxyz;ii++){
          *ec+=d_ec[ii];
        }*/
//        reduce<<<512,512>>>(d_ec,ec);
}
