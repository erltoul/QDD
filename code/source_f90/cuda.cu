#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <functional>
#include <numeric>
#include <cufft.h>
#include <math_functions.h>
//#include <cutil.h>
#include "define_cuda.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/device_free.h>
#if(parayes)
#include "mpif.h"
#endif

using namespace std;

#include "params_gpu.h"
#include "gpu_compute.cu"
#include "coulex.cu"
#include "static_gpu.cu"
#define BATCH 1 //The number of batched ffts 

#if(lda_gpu)
//Work in progress for lda_gpu
__device__ __constant__ double e2  = 2.0;

__device__ __constant__ double a0  = 0.458165293;
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
__device__ __constant__ double cr4 = 1.5874010519681993613971827;
__device__ __constant__ double t17 = 1.923661051;

double *d_chpdft;
/*thrust::device_vector<double> d_ec(1);
#if(directenergy)
thrust::device_vector<double> d_enerpw(1);
#endif*/
#endif

//extern "C" void cuda_gpu_init_(int *Npx,int *Npy,int *Npz) //initialize some variables usefull for GPU computing
extern "C" void cuda_gpu_init_() //initialize some variables usefull for GPU computing
{
#if(lda_gpu)
//Work in progress for lda_gpu
	int size_lda = 2*params_mp_nxyz_*sizeof(double);
        error=cudaMalloc((void**)&d_chpdft,size_lda);
	Check_CUDA_Error(error);
	/*d_ec.resize(nxyz);
#if(directenergy)
	d_enerpw.resize(nxyz);
#endif*/
#endif
	cudaStreamCreate(&stream1);
	cudaStreamCreate(&stream2);
	cudaStreamCreate(&stream3);
}

extern "C" void cuda_end_() //Destroy variables set above
{
#if(lda_gpu)
	cudaFree(d_chpdft);
#endif
	cudaStreamDestroy(stream1);
	cudaStreamDestroy(stream2);
	cudaStreamDestroy(stream3);
}

extern "C" void cuda_plan_1d_(cufftHandle *plan, int *dim, int *batch)
{
	//Creation of a 1D FFT plan for GPU (the plan is the same for FORWARD and BACKWARD FFT)
        if(cufftPlan1d(plan,*dim,CUFFT_Z2Z,*batch) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Plan Creation failed"<<endl;
	  exit(-1);
	}
	//Associate the plan to a stream
        if(cufftSetStream(*plan,stream2) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Streamed FFT Creation failed"<<endl;
	  exit(-1);
	}
}

extern "C" void cuda_plan_3d_(cufftHandle *plan, int *n1, int *n2, int *n3)
{
	//Creation of a 3D FFT plan for GPU (the plan is the same for FORWARD and BACKWARD FFT)
        if(cufftPlan3d(plan,*n3,*n2,*n1,CUFFT_Z2Z) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Plan Creation failed"<<endl;
	  exit(-1);
	}
	//Associate the plan to a stream
        if(cufftSetStream(*plan,stream2) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Streamed FFT Creation failed"<<endl;
	  exit(-1);
	}
}

extern "C" void cuda_plan_3d_r2c_(cufftHandle *plan, int *n1, int *n2, int *n3)
{
	//Creation of a 3D FFT plan for GPU (the plan is the same for FORWARD and BACKWARD FFT)
        if(cufftPlan3d(plan,*n3,*n2,*n1,CUFFT_D2Z) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Plan Creation failed"<<endl;
	  exit(-1);
	}
	//Associate the plan to a stream
        if(cufftSetStream(*plan,stream2) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Streamed FFT Creation failed"<<endl;
	  exit(-1);
	}
}

extern "C" void cuda_plan_3d_c2r_(cufftHandle *plan, int *n1, int *n2, int *n3)
{
	//Creation of a 3D FFT plan for GPU (the plan is the same for FORWARD and BACKWARD FFT)
        if(cufftPlan3d(plan,*n3,*n2,*n1,CUFFT_Z2D) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Plan Creation failed"<<endl;
	  exit(-1);
	}
	//Associate the plan to a stream
        if(cufftSetStream(*plan,stream2) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Streamed FFT Creation failed"<<endl;
	  exit(-1);
	}
}

extern "C" void kill_plan_(cufftHandle *plan)
{
	//Destroy FFT plan
        cufftDestroy(*plan);
}

extern "C" void run_fft_for_(cufftHandle *plan, cufftDoubleComplex *in, cufftDoubleComplex *out, int *Np) //1D FFT FORWARD
{
        int N = *Np;
        cufftDoubleComplex *d_ffta_kin;

	//Allocate device memory
        error=cudaMalloc((void**)&d_ffta_kin,sizeof(cufftDoubleComplex)*N*BATCH);
	Check_CUDA_Error(error);
	//Copy data from CPU to GPU asynchonously, so we can perform copy and FFT at the same time
        error=cudaMemcpyAsync(d_ffta_kin,in,sizeof(cufftDoubleComplex)*N*BATCH,cudaMemcpyHostToDevice,stream1);
	Check_CUDA_Error(error);
	//Perform 1D FFT FORWARD
        if(cufftExecZ2Z(*plan,d_ffta_kin,d_ffta_kin, CUFFT_FORWARD) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Exec Z2Z forward failed"<<endl;
	  exit(-1);
	}
	//Copy data from GPU to CPU
        error=cudaMemcpy(out,d_ffta_kin, sizeof(cufftDoubleComplex)*N*BATCH, cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
	//Free device memory
        cudaFree(d_ffta_kin);
}

extern "C" void run_fft_back_(cufftHandle *plan, cufftDoubleComplex *in, cufftDoubleComplex *out,int *Np)
{
        int N = *Np;
        cufftDoubleComplex *d_ffta_kin;

        //Allocate device memory
        error=cudaMalloc((void**)&d_ffta_kin,sizeof(cufftDoubleComplex)*N*BATCH);
	Check_CUDA_Error(error);
	//Copy data from CPU to GPU asynchonously, so we can perform copy and FFT at the same time
        error=cudaMemcpyAsync(d_ffta_kin,in,sizeof(cufftDoubleComplex)*N*BATCH,cudaMemcpyHostToDevice,stream1);
	Check_CUDA_Error(error);
	//Perform 1D FFT BACKWARD
        if(cufftExecZ2Z(*plan,d_ffta_kin,d_ffta_kin, CUFFT_INVERSE) != CUFFT_SUCCESS)
	{
	  cout<<("CUFFT error : Exec Z2Z backward failed")<<endl;
	  exit(-1);
	}
	//Copy data from GPU to CPU
        error=cudaMemcpy(out, d_ffta_kin, sizeof(cufftDoubleComplex)*N*BATCH,cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
	//Free device memory
        cudaFree(d_ffta_kin);
}

extern "C" void gpu_to_gpu_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta_int,int *N)
//Copy a complex array from GPU to another array on GPU
{
	int nxyz = *N;
	int size_cp = nxyz*sizeof(cufftDoubleComplex);

	error=cudaMemcpyAsync(d_ffta_int,d_ffta,size_cp,cudaMemcpyDeviceToDevice,stream1);
	Check_CUDA_Error(error);
}

extern "C" void copy_on_gpu_(cufftDoubleComplex *mat,cufftDoubleComplex *d_mat,int *N)
//Copy a complex array from CPU to GPU
{

	int nxyz = *N;
	int size_cp = nxyz*sizeof(cufftDoubleComplex);

	error=cudaMemcpyAsync(d_mat,mat,size_cp,cudaMemcpyHostToDevice,stream1);
	Check_CUDA_Error(error);
}

extern "C" void copy_real_on_gpu_(cufftDoubleReal *mat,cufftDoubleReal *d_mat,int *N)
//Copy a real array from CPU to GPU
{
	int nxyz = *N;
	int size_cp = nxyz*sizeof(cufftDoubleReal);

	error=cudaMemcpyAsync(d_mat,mat,size_cp,cudaMemcpyHostToDevice,stream1);
	Check_CUDA_Error(error);
}

extern "C" void copy_from_gpu_(cufftDoubleComplex *mat,cufftDoubleComplex *d_mat,int *N)
{
//Copy a complex array from GPU to CPU
	int nxyz = *N;
	int size_cp = nxyz*sizeof(cufftDoubleComplex);

	error=cudaMemcpy(mat,d_mat,size_cp,cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
}

extern "C" void copy_real_from_gpu_(cufftDoubleReal *mat,cufftDoubleReal *d_mat,int *N)
//Copy a real array from CPU to GPU
{
	int nxyz = *N;
	int size_cp = nxyz*sizeof(cufftDoubleReal);

	error=cudaMemcpy(mat,d_mat,size_cp,cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
}

extern "C" void run_fft_for3d_(cufftHandle *plan,cufftDoubleComplex *d_ffta,int *ty)
{
	int type = *ty;
	//Perform 3D FTT FORWARD on the GPU, with different error messages to localize easily a problem
        if(cufftExecZ2Z(*plan,d_ffta,d_ffta, CUFFT_FORWARD) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Exec Z2Z forward failed ";
	  if(type==1){
	    cout<<"in fftf"<<endl;}
	  if(type==2){
	    cout<<"in rftf"<<endl;}
	  if(type==3){
	    cout<<"in coulexf"<<endl;}
	  if(type==4){
	    cout<<"in kinprop "<<endl;}
	  exit(-1);
	}
}

extern "C" void run_fft_back3d_(cufftHandle *plan,cufftDoubleComplex *d_ffta,int *ty)
{
	int type = *ty;

	//Perform 3D FTT BACKWARD on the GPU, with different error messages to localize easily a problem
        if(cufftExecZ2Z(*plan,d_ffta,d_ffta, CUFFT_INVERSE) != CUFFT_SUCCESS)
	{
	  cout<<"CUFFT error : Exec Z2Z backward failed ";
	  if(type==1){
	    cout<<"in fftback"<<endl;}
	  if(type==2){
	    cout<<"in rfftback"<<endl;}
	  if(type==3){
	    cout<<"in coulexback"<<endl;}
	  if(type==4){
	    cout<<"in kinprop"<<endl;}
	  exit(-1);
	}
}

#if(lda_gpu)
//Work in progress
__global__ void lda_enerpw(double *d_chpdft,double *d_ec,double *d_enerpw, int nxyz)
{

       double t2,t3,t4,t5,t6,t10,t13,t23,t25,t28,t29,t34,t35,t36,t37,t42,t44,t48;
       double t53,t58,t63,t64,t65,t,t68,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;

       unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

       if (ind<nxyz)
       {
         t2 = max(d_chpdft[ind],1e-16);
         t4 = d_chpdft[ind+nxyz];
         t3 = 1.0/t2;
         t6 = pow((1.0+t4),(1.0/3.0));
         t10 = pow((1.0-t4),(1.0/3.0));
         t13 = t6*t6*t6*t6+t10*t10*t10*t10-2.0;
         t23 = (a1+da1*t13*t17)*cr3;
         t25 = cr4*cr4;
         t28 = pow((M_1_PI*t3),(1.0/3.0));
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
         t72 = t2*t4*t71;
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

         t65 = t36*t36;
         t70 = t2*(a0+da0*t13*t17+(a1+da1*t13*t17)*cr3*t29/4+(a2+da2*t13*t17)*t34*t37/4+3.0/4.0*(a3+da3*t13*t17)*M_1_PI*t3)/((b1+db1*t13*t17)*cr3*t29/4+(b2+db2*t13*t17)*t34*t37/4+3.0/4.0*(b3+db3*t13*t17)*M_1_PI*t3+3.0/16.0*(b4+db4*t13*t17)*cr3*t25*t65);
         t3=fabs((t2*(t4+1.0))/2.0);   //  *e2;
         t4=fabs((t2*(t4-1.0))/2.0);   //  *e2;
         t5= d_chpdft[ind]*t3+d_chpdft[ind+nxyz]*t4;
         
         d_ec[ind] = (-t70*e2 - 0.5*t5);
         d_enerpw[ind] = -t70*e2;
      }
}

extern "C" void calc_lda_enerpw_gpu_(double *rho,double *chpdft,double *ec,double *enerpw)
{
        int blocksize=192;
	int gridx=(int)ceil(params_mp_nxyz_/(float)blocksize);
        dim3 dimgrid(gridx,1,1);
        dim3 dimblock(blocksize,1,1);
	thrust::device_vector<double> d_ec(params_mp_nxyz_);
	thrust::device_vector<double> d_enerpw(params_mp_nxyz_);

        error=cudaMemcpyAsync(d_chpdft,rho,sizeof(double)*(params_mp_nxyz_*2),cudaMemcpyHostToDevice,stream1);
	Check_CUDA_Error(error);

        lda_enerpw<<<dimgrid,dimblock,0,stream2>>>(d_chpdft,raw_pointer_cast(&d_ec[0]),raw_pointer_cast(&d_enerpw[0]),params_mp_nxyz_);

	*ec=thrust::reduce(d_ec.begin(),d_ec.end(),(double)0.0);
	*enerpw=thrust::reduce(d_enerpw.begin(),d_enerpw.end(),(double)0.0);
        error=cudaMemcpy(chpdft,d_chpdft, sizeof(double)*(params_mp_nxyz_*2),cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
}

__global__ void lda(double *d_chpdft,double *d_ec,int nxyz)
{
       double t2,t3,t4,t5,t6,t10,t13,t23,t25,t28,t29,t34,t35,t36,t37,t42,t44,t48;
       double t53,t58,t63,t64,t65,t,t68,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;

       unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

       if (ind<nxyz)
       {
         t2 = max(d_chpdft[ind],1e-16);
         t4 = d_chpdft[ind+nxyz];
         t3 = 1.0/t2;
         t6 = pow((1.0+t4),(1.0/3.0));
         t10 = pow((1.0-t4),(1.0/3.0));
         t13 = t6*t6*t6*t6+t10*t10*t10*t10-2.0;
         t23 = (a1+da1*t13*t17)*cr3;
         t25 = cr4*cr4;
         t28 = pow((M_1_PI*t3),(1.0/3.0));
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
         t72 = t2*t4*t71;
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

         t65 = t36*t36;
         t70 = t2*(a0+da0*t13*t17+(a1+da1*t13*t17)*cr3*t29/4+(a2+da2*t13*t17)*t34*t37/4+3.0/4.0*(a3+da3*t13*t17)*M_1_PI*t3)/((b1+db1*t13*t17)*cr3*t29/4+(b2+db2*t13*t17)*t34*t37/4+3.0/4.0*(b3+db3*t13*t17)*M_1_PI*t3+3.0/16.0*(b4+db4*t13*t17)*cr3*t25*t65);

         t3=fabs((t2*(t4+1.0))/2.0);   //  *e2;
         t4=fabs((t2*(t4-1.0))/2.0);   //  *e2;
         t5= d_chpdft[ind]*t3+d_chpdft[ind+nxyz]*t4;
         
         d_ec[ind] = (-t70*e2 - 0.5*t5);
      }
}

extern "C" void calc_lda_gpu_(double *rho,double *chpdft,double *ec)
{

        int blocksize=192;
	int gridx=(int)ceil(params_mp_nxyz_/(float)blocksize);
        dim3 dimgrid(gridx,1,1);
        dim3 dimblock(blocksize,1,1);
	thrust::device_vector<double> d_ec(params_mp_nxyz_);

        error=cudaMemcpyAsync(d_chpdft,rho,2*params_mp_nxyz_*sizeof(double),cudaMemcpyHostToDevice,stream1);
	Check_CUDA_Error(error);

        lda<<<dimgrid,dimblock,0,stream2>>>(d_chpdft,raw_pointer_cast(&d_ec[0]),params_mp_nxyz_);

	*ec=thrust::reduce(d_ec.begin(),d_ec.end(),(double)0.0);
        error=cudaMemcpy(chpdft,d_chpdft,2*params_mp_nxyz_*sizeof(double),cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
}
#endif
