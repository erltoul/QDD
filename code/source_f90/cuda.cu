#include <stdio.h>
#include <math.h>
#include <cufft.h>
#include <math_functions.h>
//#include <cutil.h>
#include "define_cuda.h"
#if(lda_gpu)
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/device_free.h>
#endif

#define BATCH 1 //The number of batched ffts 

cudaError_t error;
cudaStream_t stream1,stream2;

#if(lda_gpu)
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

void Check_CUDA_Error(cudaError_t error)
{
	if(error!=cudaSuccess) 
	{
		printf("Cuda error: %s\n",cudaGetErrorString(error));
		exit(-1);
	}
}

extern "C" void cuda_gpu_init_(int *Npx,int *Npy,int *Npz)
{

#if(lda_gpu)
	int size_lda= sizeof(double)*Nx*Ny*Nz*2;

        error = cudaMalloc((void**)&d_chpdft,size_lda);
	Check_CUDA_Error(error);
	/*d_ec.resize(nxyz);
#if(directenergy)
	d_enerpw.resize(nxyz);
#endif*/
#endif
	cudaStreamCreate(&stream1);
	cudaStreamCreate(&stream2);
}

extern "C" void cuda_end_()
{
#if(lda_gpu)
	cudaFree(d_chpdft);
#endif
	cudaStreamDestroy(stream1);
	cudaStreamDestroy(stream2);
}

extern "C" void cuda_plan_1d_(cufftHandle *plan, int *dim, int *batch)
{
        if(cufftPlan1d(plan,*dim,CUFFT_Z2Z,*batch) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Plan Creation failed\n");
	  exit(-1);
	}
        /*if(cufftSetStream(*plan,stream2) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Streamed FFT Creation failed\n");
	  exit(-1);
	}*/
}

extern "C" void cuda_plan_3d_(cufftHandle *plan, int *n1, int *n2, int *n3)
{
        if(cufftPlan3d(plan,*n3,*n2,*n1,CUFFT_Z2Z) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Plan Creation failed\n");
	  exit(-1);
	}
        /*if(cufftSetStream(*plan,stream2) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Streamed FFT Creation failed\n");
	  exit(-1);
	}*/
}

extern "C" void kill_plan_(cufftHandle *plan)
{
        cufftDestroy(*plan);
}

extern "C" void run_fft_for_(cufftHandle *plan, cufftDoubleComplex *in, cufftDoubleComplex *out, int *Np)
{
        int N = *Np;

        //Allocate device memory
        cufftDoubleComplex *d_ffta_kin;

        cudaMalloc((void**)&d_ffta_kin,sizeof(cufftDoubleComplex)*N*BATCH);
	Check_CUDA_Error(error);
        cudaMemcpyAsync(d_ffta_kin,in,sizeof(cufftDoubleComplex)*N*BATCH,cudaMemcpyHostToDevice,stream1);
	Check_CUDA_Error(error);
        if(cufftExecZ2Z(*plan,d_ffta_kin,d_ffta_kin, CUFFT_FORWARD) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Exec Z2Z forward failed\n");
	  exit(-1);
	}
        cudaMemcpy(out,d_ffta_kin, sizeof(cufftDoubleComplex)*N*BATCH, cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
        cudaFree(d_ffta_kin);
}

extern "C" void run_fft_back_(cufftHandle *plan, cufftDoubleComplex *in, cufftDoubleComplex *out,int *Np)
{
        int N = *Np;

        //Allocate device memory
        cufftDoubleComplex *d_ffta_kin;

        cudaMalloc((void**)&d_ffta_kin,sizeof(cufftDoubleComplex)*N*BATCH);
	Check_CUDA_Error(error);
        cudaMemcpyAsync(d_ffta_kin,in,sizeof(cufftDoubleComplex)*N*BATCH,cudaMemcpyHostToDevice,stream1);
	Check_CUDA_Error(error);
        if(cufftExecZ2Z(*plan,d_ffta_kin,d_ffta_kin, CUFFT_INVERSE) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Exec Z2Z backward failed\n");
	  exit(-1);
	}
        cudaMemcpy(out, d_ffta_kin, sizeof(cufftDoubleComplex)*N*BATCH,cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
        cudaFree(d_ffta_kin); //free the signal on the device.
}

__global__ void multiply_device(cufftDoubleComplex *d_ffta,int nxyz,double norm)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

	if (ind<nxyz)
	{
		d_ffta[ind].x=norm*d_ffta[ind].x;
		d_ffta[ind].y=norm*d_ffta[ind].y;
	}
}

__global__ void multiply_ak_gpu(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ak,int nxyz)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;
	cufftDoubleReal SAVE2;

	if (ind<nxyz)
	{
		SAVE2         = d_ak[ind].x*d_ffta[ind].x+d_ak[ind].y*d_ffta[ind].y;
		d_ffta[ind].y = d_ak[ind].x*d_ffta[ind].y+d_ak[ind].y*d_ffta[ind].x;
		d_ffta[ind].x = SAVE2;
	}
}

__global__ void multiply_ak_gpu2(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ak,int nxyz)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;
	cufftDoubleReal SAVE2;

	if (ind<nxyz)
	{
		SAVE2         = d_ak[ind].x*d_ffta[ind].x-d_ak[ind].y*d_ffta[ind].y;
		d_ffta[ind].y = d_ak[ind].x*d_ffta[ind].y+d_ak[ind].y*d_ffta[ind].x;
		d_ffta[ind].x = SAVE2;
	}
}

__global__ void multiply_ak_real_gpu(cufftDoubleComplex *d_ffta,cufftDoubleReal *d_ak,int nxyz)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

	if (ind<nxyz)
	{
		d_ffta[ind].x = d_ak[ind]*d_ffta[ind].x;
		d_ffta[ind].y = d_ak[ind]*d_ffta[ind].y;
	}
}

__global__ void multiply_shift_gpu(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_akx,cufftDoubleComplex *d_aky,cufftDoubleComplex *d_akz,cufftDoubleComplex shix,cufftDoubleComplex shiy,cufftDoubleComplex shiz,int nxyz)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;
	cufftDoubleReal SAVE2;

	if (ind<nxyz)
	{
                shix.x = shix.x*d_akx[ind].x;
		shiy.x = shiy.x*d_aky[ind].x;
		shiz.x = shiz.x*d_akz[ind].x;
                shix.y = shix.y*d_akx[ind].y;
		shiy.y = shiy.y*d_aky[ind].y;
		shiz.y = shiz.y*d_akz[ind].y;

		//multiplication by shix
		SAVE2         = exp(shix.x)*(d_ffta[ind].x*cos(shix.y)-sin(shix.y)*d_ffta[ind].y);
		d_ffta[ind].y = exp(shix.x)*(d_ffta[ind].x*sin(shix.y)+cos(shix.y)*d_ffta[ind].y);
		d_ffta[ind].x = SAVE2;
		//multiplication by shiy
		SAVE2         = exp(shiy.x)*(d_ffta[ind].x*cos(shiy.y)-sin(shiy.y)*d_ffta[ind].y);
		d_ffta[ind].y = exp(shiy.x)*(d_ffta[ind].x*sin(shiy.y)+cos(shiy.y)*d_ffta[ind].y);
		d_ffta[ind].x = SAVE2;
		//multiplication by shiz
		SAVE2         = exp(shiz.x)*(d_ffta[ind].x*cos(shiz.y)-sin(shiz.y)*d_ffta[ind].y);
		d_ffta[ind].y = exp(shiz.x)*(d_ffta[ind].x*sin(shiz.y)+cos(shiz.y)*d_ffta[ind].y);
		d_ffta[ind].x = SAVE2;
	}
}

extern "C" void multiply_gpu_(cufftDoubleComplex *d_ffta,int *N,double *tnorm)
{
	int nxyz = *N;
	double norm = *tnorm;
	int blocksize=512;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

        multiply_device<<<dimgrid,dimblock>>>(d_ffta,nxyz,norm);
	Check_CUDA_Error(error);
}

extern "C" void multiply_ak_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ak,int *N)
{
	int nxyz = *N;
	int blocksize=512;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

        multiply_ak_gpu<<<dimgrid,dimblock>>>(d_ffta,d_ak,nxyz);
	Check_CUDA_Error(error);
}

extern "C" void multiply_ak2_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ak,int *N)
{
	int nxyz = *N;
	int blocksize=512;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

        multiply_ak_gpu2<<<dimgrid,dimblock>>>(d_ffta,d_ak,nxyz);
	Check_CUDA_Error(error);
}

extern "C" void multiply_ak_real_(cufftDoubleComplex *d_ffta,cufftDoubleReal *d_ak,int *N)
{
	int nxyz = *N;
	int blocksize=512;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

        multiply_ak_real_gpu<<<dimgrid,dimblock>>>(d_ffta,d_ak,nxyz);
	Check_CUDA_Error(error);
}

extern "C" void multiply_shift_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_akx,cufftDoubleComplex *d_aky,cufftDoubleComplex *d_akz,double *sx,double *sy,double *sz,int *N)
{
	int nxyz = *N;
	cufftDoubleComplex shix,shiy,shiz;
	int blocksize=512;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

	shix.x = *sx;
	shix.y = *sx;
	shiy.x = *sy;
	shiy.y = *sy;
	shiz.x = *sz;
	shiz.y = *sz;

        multiply_shift_gpu<<<dimgrid,dimblock>>>(d_ffta,d_akx,d_aky,d_akz,shix,shiy,shiz,nxyz);
	Check_CUDA_Error(error);
}

extern "C" void copy_on_gpu_(cufftDoubleComplex *mat,cufftDoubleComplex *d_mat,int *N)
{
	int nxyz = *N;
	int size_cp = nxyz*sizeof(cufftDoubleComplex);

	cudaMemcpy(d_mat,mat,size_cp,cudaMemcpyHostToDevice);
	Check_CUDA_Error(error);
}

extern "C" void copy_real_on_gpu_(cufftDoubleReal *mat,cufftDoubleReal *d_mat,int *N)
{
	int nxyz = *N;
	int size_cp = nxyz*sizeof(cufftDoubleReal);

	cudaMemcpy(d_mat,mat,size_cp,cudaMemcpyHostToDevice);
	Check_CUDA_Error(error);
}

extern "C" void copy_from_gpu_(cufftDoubleComplex *mat,cufftDoubleComplex *d_mat,int *N)
{
	int nxyz = *N;
	int size_cp = nxyz*sizeof(cufftDoubleComplex);

	cudaMemcpy(mat,d_mat,size_cp,cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
}

extern "C" void copy_real_from_gpu_(cufftDoubleReal *mat,cufftDoubleReal *d_mat,int *N)
{
	int nxyz = *N;
	int size_cp = nxyz*sizeof(cufftDoubleReal);

	cudaMemcpy(mat,d_mat,size_cp,cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
}

extern "C" void run_fft_for3d_(cufftHandle *plan,cufftDoubleComplex *d_ffta,int *ty)
{
	int type = *ty;

        if(cufftExecZ2Z(*plan,d_ffta,d_ffta, CUFFT_FORWARD) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Exec Z2Z forward failed\n");
	  if(type==1){
	    printf("From fftf\n");}
	  if(type==2){
	    printf("From rftf\n");}
	  if(type==3){
	    printf("From coulexf\n");}
	  exit(-1);
	}
}

extern "C" void run_fft_back3d_(cufftHandle *plan,cufftDoubleComplex *d_ffta,int *ty)
{
	int type = *ty;

        if(cufftExecZ2Z(*plan,d_ffta,d_ffta, CUFFT_INVERSE) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Exec Z2Z backward failed\n");
	  if(type==1){
	    printf("From fftback\n");}
	  if(type==2){
	    printf("From rfftback\n");}
	  if(type==3){
	    printf("From coulexback\n");}
	  exit(-1);
	}
}

#if(lda_gpu)
__global__ void lda_enerpw(double *d_chpdft,double *d_ec,double *d_enerpw, int nxyz,double e2)
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

extern "C" void calc_lda_enerpw_gpu_(double *rho,double *chpdft,int *N,double *e2,double *ec,double *enerpw)
{
        const int nxyz =*N;
        //double *d_chpdft;
	//double *d_ec;
        //double *d_rho;
        double e=*e2;
        //double enerpw = *d_enerpw;
        int blocksize=512;
	int gridx=(int)ceil(nxyz/(float)blocksize);
        dim3 dimgrid(gridx,1,1);
        dim3 dimblock(blocksize,1,1);
        //double h_ec[nxyz];
        //unsigned int *d_ind;
        //unsigned int *h_ind;
	thrust::device_vector<double> d_ec(nxyz);
	thrust::device_vector<double> d_enerpw(nxyz);
	//thrust::host_vector<double> h_ec(nxyz);

        //h_ec=(double*)malloc(nxyz,sizeof(double));
        /*h_ind=(unsigned int*)calloc(nxyz,sizeof(unsigned int));

        cudaMalloc((void**)&d_ind,sizeof(unsigned int)*nxyz);
        cudaMemcpy(d_ind,h_ind,sizeof(unsigned int)*(nxyz),cudaMemcpyHostToDevice);*/
        //cudaMalloc((void**)&d_rho,sizeof(double)*(nxyz*2));
        //cudaMemcpy(d_rho,rho,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
        //cudaMalloc((void**)&d_chpdft,sizeof(double)*(nxyz*2));
        //cudaMemcpy(d_chpdft,chpdft,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
        cudaMemcpyAsync(d_chpdft,rho,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice,stream1);
	Check_CUDA_Error(error);
        //cudaMalloc((void**)&d_ec,sizeof(double)*nxyz);
        //cudaMalloc((void**)&d_enerpw,sizeof(double)*nxyz);
        //cudaMemcpy(d_ec,h_ec,sizeof(double)*(nxyz),cudaMemcpyHostToDevice);
        lda_enerpw<<<dimgrid,dimblock,0,stream2>>>(d_chpdft,raw_pointer_cast(&d_ec[0]),raw_pointer_cast(&d_enerpw[0]),nxyz,e);
	Check_CUDA_Error(error);
        //cudaMemcpy(h_ind,d_ind, sizeof(unsigned int)*(nxyz),cudaMemcpyDeviceToHost);
        //cudaMemcpy(h_ec,d_ec, sizeof(double)*(nxyz),cudaMemcpyDeviceToHost);
	//h_ec=d_ec;
	*ec=thrust::reduce(d_ec.begin(),d_ec.end(),(double)0.0);
	*enerpw=thrust::reduce(d_enerpw.begin(),d_enerpw.end(),(double)0.0);
        cudaMemcpy(chpdft,d_chpdft, sizeof(double)*(nxyz*2),cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
        //cudaMemcpy(rho,d_rho, sizeof(double)*(nxyz*2),cudaMemcpyDeviceToHost);
        /*for (int ii=0;ii<nxyz;ii++){
          *ec+=h_ec[ii];
        }*/
        //cudaFree(d_rho);
        //cudaFree(d_chpdft);
        //cudaFree(d_ec);
}

__global__ void lda(double *d_chpdft,double *d_ec,int nxyz,double e2)
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

extern "C" void calc_lda_gpu_(double *rho,double *chpdft,int *N,double *e2,double *ec)
{
        const int nxyz =*N;
        //double *d_chpdft;
	//double *d_ec;
        //double *d_rho;
        double e=*e2;
        int blocksize=512;
	int gridx=(int)ceil(nxyz/(float)blocksize);
        dim3 dimgrid(gridx,1,1);
        dim3 dimblock(blocksize,1,1);
        //double h_ec[nxyz];
        //unsigned int *d_ind;
        //unsigned int *h_ind;
	thrust::device_vector<double> d_ec(nxyz);
	//thrust::host_vector<double> h_ec(nxyz);

        //h_ec=(double*)malloc(nxyz,sizeof(double));
        /*h_ind=(unsigned int*)calloc(nxyz,sizeof(unsigned int));

        cudaMalloc((void**)&d_ind,sizeof(unsigned int)*nxyz);
        cudaMemcpy(d_ind,h_ind,sizeof(unsigned int)*(nxyz),cudaMemcpyHostToDevice);*/
        //cudaMalloc((void**)&d_rho,sizeof(double)*(nxyz*2));
        //cudaMemcpy(d_rho,rho,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
        //cudaMalloc((void**)&d_chpdft,sizeof(double)*(nxyz*2));
        //cudaMemcpy(d_chpdft,chpdft,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
        cudaMemcpyAsync(d_chpdft,rho,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice,stream1);
	Check_CUDA_Error(error);
        //cudaMalloc((void**)&d_ec,sizeof(double)*nxyz);
        //cudaMemcpy(d_ec,h_ec,sizeof(double)*(nxyz),cudaMemcpyHostToDevice);
        lda<<<dimgrid,dimblock,0,stream2>>>(d_chpdft,raw_pointer_cast(&d_ec[0]),nxyz,e);
	Check_CUDA_Error(error);
        //cudaMemcpy(h_ind,d_ind, sizeof(unsigned int)*(nxyz),cudaMemcpyDeviceToHost);
        //cudaMemcpy(h_ec,d_ec, sizeof(double)*(nxyz),cudaMemcpyDeviceToHost);
	//h_ec=d_ec;
	*ec=thrust::reduce(d_ec.begin(),d_ec.end(),(double)0.0);
        cudaMemcpy(chpdft,d_chpdft, sizeof(double)*(nxyz*2),cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
        //cudaMemcpy(rho,d_rho, sizeof(double)*(nxyz*2),cudaMemcpyDeviceToHost);
        /*for (int ii=0;ii<nxyz;ii++){
          *ec+=h_ec[ii];
        }*/
        //cudaFree(d_rho);
        //cudaFree(d_chpdft);
        //cudaFree(d_ec);
}
#endif
