#include <stdio.h>
#include <math.h>
#include <cufft.h>
#include <math_functions.h> //math library for GPU
#include "define_cuda.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/device_free.h>
//#include <cutil.h>

#define BATCH 1 //The number of batched ffts 

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

/*float time_fft_for;
float time_mem;
int count_fft;*/
#if(lda_gpu)
/*float time_lda;
count_lda;*/
#endif
//cufftDoubleComplex *d_ffta_coulex;
cufftDoubleComplex *d_ffta;
#if(lda_gpu)
double *d_chpdft;
/*thrust::device_vector<double> d_ec(1);
#if(directenergy)
thrust::device_vector<double> d_enerpw(1);
#endif*/
#endif
cudaError_t error;
//cudaStream_t stream1,stream2;

void Check_CUDA_Error(cudaError_t error)
{
	if(error!=cudaSuccess) 
	{
		printf("Cuda error: %s\n",cudaGetErrorString(error));
		exit(-1);
	}
}

/*extern "C" void cuda_fftvar_alloc_(int *Npx,int *Npy,int *Npz)
{
        int Nx = *Npx;
        int Ny = *Npy;
        int Nz = *Npz;
	int nxyz = Nx*Ny*Nz;

        cudaMalloc((void**)&d_ffta,sizeof(cufftDoubleComplex)*nxyz);
}*/

/*void cuda_fftvar_coulex_alloc(int *Npx,int *Npy,int *Npz)
{
        int Nx = *Npx;
        int Ny = *Npy;
        int Nz = *Npz;
	int nxyz = Nx*Ny*Nz;

        cudaMalloc((void**)&d_ffta_coulex,sizeof(cufftDoubleComplex)*nxyz);
}*/

extern "C" void cuda_gpu_init_(int *Npx,int *Npy,int *Npz)
{
        int Nx = *Npx;
        int Ny = *Npy;
        int Nz = *Npz;
	int nxyz = Nx*Ny*Nz;

        //error = cudaMalloc((void**)&d_ffta_coulex,sizeof(cufftDoubleComplex)*(nxyz*8)); //Does not work if allocated here (stange...)
	//Check_CUDA_Error(error);
        error = cudaMalloc((void**)&d_ffta,sizeof(cufftDoubleComplex)*nxyz);
	Check_CUDA_Error(error);
#if(lda_gpu)
        error = cudaMalloc((void**)&d_chpdft,sizeof(double)*(nxyz*2));
	Check_CUDA_Error(error);
	/*d_ec.resize(nxyz);
#if(directenergy)
	d_enerpw.resize(nxyz);
#endif*/
#endif
	//cudaStreamCreate(&stream1);
	//cudaStreamCreate(&stream2);
}

extern "C" void cuda_end_()
{
        cudaFree(d_ffta);
	//cudaFree(d_ffta_coulex);
#if(lda_gpu)
	cudaFree(d_chpdft);
#endif
	//cudaStreamDestroy(stream1);
	//cudaStreamDestroy(stream2);
}

extern "C" void cuda_plan_1d_(cufftHandle *plan, int *dim, int *batch)
{
        if(cufftPlan1d(plan,*dim,CUFFT_Z2Z,*batch) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Plan Creation failed");
	  exit(-1);
	}
}

extern "C" void cuda_plan_3d_(cufftHandle *plan, int *n1, int *n2, int *n3)
{
        if(cufftPlan3d(plan,*n3,*n2,*n1,CUFFT_Z2Z) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Plan Creation failed");
	  exit(-1);
	}
}

extern "C" void kill_plan_(cufftHandle * plan)
{
        cufftDestroy(*plan);
}

extern "C" void run_fft_for_(cufftHandle *plan, cufftDoubleComplex *a, cufftDoubleComplex *b, int *Np)
{
        int N = *Np;

        //Allocate device memory
        cufftDoubleComplex *d_ffta_kin;

        cudaMalloc((void**)&d_ffta_kin,sizeof(cufftDoubleComplex)*N*BATCH);
	Check_CUDA_Error(error);
        cudaMemcpy(d_ffta_kin,a,sizeof(cufftDoubleComplex)*N*BATCH,cudaMemcpyHostToDevice);
	Check_CUDA_Error(error);
        if(cufftExecZ2Z(*plan,(cufftDoubleComplex *)d_ffta_kin,(cufftDoubleComplex *)d_ffta_kin, CUFFT_FORWARD) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Exec Z2Z forward failed");
	  exit(-1);
	}
        cudaMemcpy(b, d_ffta_kin, sizeof(cufftDoubleComplex)*N*BATCH, cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
        cudaFree(d_ffta_kin);
}

extern "C" void run_fft_back_(cufftHandle *plan, cufftDoubleComplex *a, cufftDoubleComplex *b,int *Np)
{
        int N = *Np;
        //Allocate device memory
        cufftDoubleComplex *d_ffta_kin;

        cudaMalloc((void**)&d_ffta_kin,sizeof(cufftDoubleComplex)*N*BATCH);
	Check_CUDA_Error(error);
        cudaMemcpy(d_ffta_kin,a,sizeof(cufftDoubleComplex)*N*BATCH,cudaMemcpyHostToDevice);
	Check_CUDA_Error(error);
        if(cufftExecZ2Z(*plan, (cufftDoubleComplex *)d_ffta_kin,(cufftDoubleComplex *)d_ffta_kin, CUFFT_INVERSE) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Exec Z2Z backward failed");
	  exit(-1);
	}
        cudaMemcpy(b, d_ffta_kin, sizeof(cufftDoubleComplex)*N*BATCH, cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
        cudaFree(d_ffta_kin); //free the signal on the device.
}

extern "C" void run_fft_for3d_(cufftHandle *plan, cufftDoubleComplex *in, cufftDoubleComplex *out,int *Npx,int *Npy,int *Npz)
{
        int Nx = *Npx;
        int Ny = *Npy;
        int Nz = *Npz;
	int size_cp = Nx*Ny*Nz*sizeof(cufftDoubleComplex);
	//cudaEvent_t start,stop,start_mem,stop_mem;
	//float int_time,int_time_mem;
        //Allocate device memory
        //cufftDoubleComplex *d_ffta;

	//count_fft++;
        //cudaMalloc((void**)&d_ffta,sizeof(cufftDoubleComplex)*nxyz);
	/*cudaEventCreate(&start_mem);
	cudaEventCreate(&stop_mem);
	cudaEventRecord(start_mem,0);*/
        cudaMemcpy(d_ffta,in,size_cp,cudaMemcpyHostToDevice);
	Check_CUDA_Error(error);
	/*cudaEventRecord(stop_mem,0);
	cudaEventSynchronize(stop_mem);
	cudaEventElapsedTime(&int_time_mem,start_mem,stop_mem);
	time_mem=time_mem+int_time_mem;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);*/
        if(cufftExecZ2Z(*plan, (cufftDoubleComplex *)d_ffta,(cufftDoubleComplex *)d_ffta, CUFFT_FORWARD) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Exec Z2Z forward failed");
	  exit(-1);
	}
	/*cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&int_time,start,stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	time_fft_for=time_fft_for+int_time;
	if(count_fft%100==0) printf("Time fft_for: %f\n",time_fft_for);
	cudaEventRecord(start_mem,0);*/
        cudaMemcpy(out,d_ffta,size_cp,cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
	/*cudaEventRecord(stop_mem,0);
	cudaEventSynchronize(stop_mem);
	cudaEventElapsedTime(&int_time_mem,start_mem,stop_mem);
	time_mem=time_mem+int_time_mem;
	if(count_fft%100==0) printf("Time memcopy: %f\n",time_mem);
	cudaEventDestroy(start_mem);
	cudaEventDestroy(stop_mem);*/
        //cudaFree(d_ffta);
}

extern "C" void run_fft_for3d_coulex_(cufftHandle *plan, cufftDoubleComplex *in, cufftDoubleComplex *out,int *Npx,int *Npy,int *Npz)
{
        int Nx = *Npx;
        int Ny = *Npy;
        int Nz = *Npz;
	int size_cp = Nx*Ny*Nz*sizeof(cufftDoubleComplex);

        //Allocate device memory
        cufftDoubleComplex *d_ffta_coulex;

        cudaMalloc((void**)&d_ffta_coulex,size_cp);
	Check_CUDA_Error(error);
        cudaMemcpy(d_ffta_coulex,in,size_cp,cudaMemcpyHostToDevice);
	Check_CUDA_Error(error);
        if(cufftExecZ2Z(*plan, (cufftDoubleComplex *)d_ffta_coulex,(cufftDoubleComplex *)d_ffta_coulex, CUFFT_FORWARD) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Exec Z2Z forward failed");
	  exit(-1);
	}
        cudaMemcpy(out,d_ffta_coulex,size_cp,cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
        cudaFree(d_ffta_coulex);
}

extern "C" void run_fft_back3d_(cufftHandle *plan, cufftDoubleComplex *a, cufftDoubleComplex *b,int *Npx,int *Npy,int *Npz)
{
        int Nx = *Npx;
        int Ny = *Npy;
        int Nz = *Npz;
	int size_cp = Nx*Ny*Nz*sizeof(cufftDoubleComplex);

        //Allocate device memory
        //cufftDoubleComplex *d_ffta;
        
        //cudaMalloc((void**)&d_ffta,sizeof(cufftDoubleComplex)*Nx*Ny*Nz);
        cudaMemcpy(d_ffta,a,size_cp,cudaMemcpyHostToDevice);
	Check_CUDA_Error(error);
        if(cufftExecZ2Z(*plan, (cufftDoubleComplex *)d_ffta,(cufftDoubleComplex *)d_ffta, CUFFT_INVERSE) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Exec Z2Z backward failed");
	  exit(-1);
	}
        cudaMemcpy(b, d_ffta,size_cp,cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
        //cudaFree(d_ffta); //free the signal on the device.
}

extern "C" void run_fft_back3d_coulex_(cufftHandle *plan, cufftDoubleComplex *a, cufftDoubleComplex *b,int *Npx,int *Npy,int *Npz)
{
        int Nx = *Npx;
        int Ny = *Npy;
        int Nz = *Npz;
	int size_cp = Nx*Ny*Nz*sizeof(cufftDoubleComplex);

        //Allocate device memory
        cufftDoubleComplex *d_ffta_coulex;
        
        cudaMalloc((void**)&d_ffta_coulex,size_cp);
	Check_CUDA_Error(error);
        cudaMemcpy(d_ffta_coulex,a,size_cp,cudaMemcpyHostToDevice);
	Check_CUDA_Error(error);
        if(cufftExecZ2Z(*plan, (cufftDoubleComplex *)d_ffta_coulex,(cufftDoubleComplex *)d_ffta_coulex, CUFFT_INVERSE) != CUFFT_SUCCESS)
	{
	  printf("CUFFT error : Exec Z2Z backward failed");
	  exit(-1);
	}
        cudaMemcpy(b, d_ffta_coulex, size_cp,cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
        cudaFree(d_ffta_coulex); //free the signal on the device.
}

#if(lda_gpu)
__global__ void lda(double *d_rho,double *d_chpdft,double *d_ec,int nxyz,double e2, double *d_enerpw)
{

       //double t1,t2,t3,t4,t5,t6,t7,t8,t10,t11,t12,t13,t17,t23,t25,t28,t29,t34,t35,t36,t37,t42,t44;
       //double t48,t53,t58,t63,t64,t65,t,t68,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;
       double t2,t3,t4,t5,t6,t10,t13,t23,t25,t28,t29,t34,t35,t36,t37,t42,t44,t48;
       double t53,t58,t63,t64,t65,t,t68,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;
       //double rp,xi;

       unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

       if (ind<nxyz)
       {
         //d_ind[ind]=ind;
         /*rp = max(d_rho[ind],1e-16);
         xi = d_rho[ind+nxyz];*/
         //rp = max(d_chpdft[ind],1e-16);
         //xi = d_chpdft[ind+nxyz];
         t2 = max(d_chpdft[ind],1e-16);
         t4 = d_chpdft[ind+nxyz];
  
         //t1 = xi*rp;
         //t1 = t2*t4;
         //t2 = rp;
         t3 = 1.0/t2;
         //t4 = xi;
         t6 = pow((1.0+t4),(1.0/3.0));
         //t7 = t6*t6;
         //t8 = t7*t7;
         //t8 = t6*t6*t6*t6;
         t10 = pow((1.0-t4),(1.0/3.0));
         //t11 = t10*t10;
         //t12 = t10*t10*t10*t10;
         //t13 = t8+t12-2.0;
         t13 = t6*t6*t6*t6+t10*t10*t10*t10-2.0;
         //t15 = 2.0**(1.0/3.0);
         //t17 = 1.0/(2.0*cr2-2.0);
         //t22 = 3.0**(1.0/3.0);
         t23 = (a1+da1*t13*t17)*cr3;
         //t24 = 4.0**(1.0/3.0);
         t25 = cr4*cr4;
         //t26 = 1/pi;
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
         //t72 = t1*t71;
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

         //t1=rp;
         //t4 = xi;
         //t6 = pow((1.0+t4),(1.0/3.0));
         //t7 = t6*t6;
         //t8 = t7*t7;
         //t6 = pow((1.0-t4),(1.0/3.0));;
         //t11 = t10*t10;
         //t12 = t11*t11;
         //t13 = t8+t12-2.0;
         //t15 = 2**(1.0/3.0);
         //t17 = 1.0/(2.0*cr2-2.0);
         //t22 = 3**(1.0/3.0);
         //t24 = 4**(1.0/3.0);
         //t25 = cr4*cr4;
         //t26 = 1/pi;
         //t28 = pow((M_1_PI/t2),(1.0/3.0));;
         //t29 = t25*t28;
         //t34 = cr3*cr3;
         //t36 = t28*t28;
         //t37 = cr4*t36;
         t65 = t36*t36;
         //t70 = t1*(a0+da0*t13*t17+(a1+da1*t13*t17)*cr3*t29/4+(a2+da2*t13*t17)*t34*t37/4+3.0/4.0*(a3+da3*t13*t17)*M_1_PI*t3)/((b1+db1*t13*t17)*cr3*t29/4+(b2+db2*t13*t17)*t34*t37/4+3.0/4.0*(b3+db3*t13*t17)*M_1_PI*t3+3.0/16.0*(b4+db4*t13*t17)*cr3*t25*t65);
         t70 = t2*(a0+da0*t13*t17+(a1+da1*t13*t17)*cr3*t29/4+(a2+da2*t13*t17)*t34*t37/4+3.0/4.0*(a3+da3*t13*t17)*M_1_PI*t3)/((b1+db1*t13*t17)*cr3*t29/4+(b2+db2*t13*t17)*t34*t37/4+3.0/4.0*(b3+db3*t13*t17)*M_1_PI*t3+3.0/16.0*(b4+db4*t13*t17)*cr3*t25*t65);

         //t1 = xi*rp;
         //t2 = rp;
         //t3=fabs((t1+t2)/2.0);   //  *e2;
         //t4=fabs((t1-t2)/2.0);   //  *e2;
         t3=fabs((t2*(t4+1.0))/2.0);   //  *e2;
         t4=fabs((t2*(t4-1.0))/2.0);   //  *e2;
         t5= d_chpdft[ind]*t3+d_chpdft[ind+nxyz]*t4;
         
         d_ec[ind] = (-t70*e2 - 0.5*t5);
         d_enerpw[ind] = -t70*e2;
      }
}

extern "C" void calc_lda_gpu_(double *rho,double *chpdft,int *N,double *e2,double *ec,double *enerpw)
{
        const int nxyz =*N;
	//cudaEvent_t start,stop,start_mem,stop_mem;
	//float int_time,int_time_mem;
        //double *d_chpdft;
	//double *d_ec;
        //double *d_rho;
        double e=*e2;
        double enerpw = *d_enpw;
        int blocksize=512;
        dim3 dimgrid(ceil(nxyz/(float)blocksize),1,1);
        dim3 dimblock(blocksize,1,1);
        //double h_ec[nxyz];
        //unsigned int *d_ind;
        //unsigned int *h_ind;
	thrust::device_vector<double> d_ec(nxyz);
	thrust::device_vector<double> d_enerpw(nxyz);
	//count_lda++;
	//thrust::host_vector<double> h_ec(nxyz);

        //h_ec=(double*)malloc(nxyz,sizeof(double));
        /*h_ind=(unsigned int*)calloc(nxyz,sizeof(unsigned int));

        cudaMalloc((void**)&d_ind,sizeof(unsigned int)*nxyz);
        cudaMemcpy(d_ind,h_ind,sizeof(unsigned int)*(nxyz),cudaMemcpyHostToDevice);*/
        //cudaMalloc((void**)&d_rho,sizeof(double)*(nxyz*2));
        //cudaMemcpy(d_rho,rho,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
        //cudaMalloc((void**)&d_chpdft,sizeof(double)*(nxyz*2));
        //cudaMemcpy(d_chpdft,chpdft,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
	/*cudaEventCreate(&start_mem);
	cudaEventCreate(&stop_mem);
	cudaEventRecord(start_mem,0);*/
        cudaMemcpyAsync(d_chpdft,rho,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
	Check_CUDA_Error(error);
	/*cudaEventRecord(stop_mem,0);
	cudaEventSynchronize(stop_mem);
	cudaEventElapsedTime(&int_time_mem,start_mem,stop_mem);
	time_mem=time_mem+int_time_mem;*/
        //cudaMalloc((void**)&d_ec,sizeof(double)*nxyz);
        //cudaMalloc((void**)&d_enerpw,sizeof(double)*nxyz);
        //cudaMemcpy(d_ec,h_ec,sizeof(double)*(nxyz),cudaMemcpyHostToDevice);
	/*cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);*/
        lda<<<dimgrid,dimblock>>>(d_chpdft,raw_pointer_cast(&d_ec[0]),raw_pointer_cast(&d_enerpw[0]),nxyz,e);
	Check_CUDA_Error(error);
        //cudaMemcpy(h_ind,d_ind, sizeof(unsigned int)*(nxyz),cudaMemcpyDeviceToHost);
        //cudaMemcpy(h_ec,d_ec, sizeof(double)*(nxyz),cudaMemcpyDeviceToHost);
	//h_ec=d_ec;
	*ec=thrust::reduce(d_ec.begin(),d_ec.end(),(double)0.0);
	*enerpw=thrust::reduce(d_enerpw.begin(),d_enerpw.end(),(double)0.0);
	/*cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&int_time,start,stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	time_lda=time_lda+int_time;
	if(count_lda%100==0) printf("Time lda: %f\n",time_lda);*/
	//cudaEventRecord(start_mem,0);
        cudaMemcpyAsync(chpdft,d_chpdft, sizeof(double)*(nxyz*2),cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
	/*cudaEventRecord(stop_mem,0);
	cudaEventSynchronize(stop_mem);
	cudaEventElapsedTime(&int_time_mem,start_mem,stop_mem);
	time_mem=time_mem+int_time_mem;
	if(count_lda%100==0) printf("Time memcopy: %f\n",time_mem);
	cudaEventDestroy(start_mem);
	cudaEventDestroy(stop_mem);*/
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

       //double t1,t2,t3,t4,t5,t6,t7,t8,t10,t11,t12,t13,t17,t23,t25,t28,t29,t34,t35,t36,t37,t42,t44;
       //double t48,t53,t58,t63,t64,t65,t,t68,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;
       double t2,t3,t4,t5,t6,t10,t13,t23,t25,t28,t29,t34,t35,t36,t37,t42,t44,t48;
       double t53,t58,t63,t64,t65,t,t68,t70,t71,t72,t77,t82,t83,t88,t93,t98,t102,t109,t135;
       //double rp,xi;

       unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

       if (ind<nxyz)
       {
         //d_ind[ind]=ind;
         /*rp = max(d_rho[ind],1e-16);
         xi = d_rho[ind+nxyz];*/
         //rp = max(d_chpdft[ind],1e-16);
         //xi = d_chpdft[ind+nxyz];
         t2 = max(d_chpdft[ind],1e-16);
         t4 = d_chpdft[ind+nxyz];
  
         //t1 = xi*rp;
         //t1 = t2*t4;
         //t2 = rp;
         t3 = 1.0/t2;
         //t4 = xi;
         t6 = pow((1.0+t4),(1.0/3.0));
         //t7 = t6*t6;
         //t8 = t7*t7;
         //t8 = t6*t6*t6*t6;
         t10 = pow((1.0-t4),(1.0/3.0));
         //t11 = t10*t10;
         //t12 = t10*t10*t10*t10;
         //t13 = t8+t12-2.0;
         t13 = t6*t6*t6*t6+t10*t10*t10*t10-2.0;
         //t15 = 2.0**(1.0/3.0);
         //t17 = 1.0/(2.0*cr2-2.0);
         //t22 = 3.0**(1.0/3.0);
         t23 = (a1+da1*t13*t17)*cr3;
         //t24 = 4.0**(1.0/3.0);
         t25 = cr4*cr4;
         //t26 = 1/pi;
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
         //t72 = t1*t71;
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

         //t1=rp;
         //t4 = xi;
         //t6 = pow((1.0+t4),(1.0/3.0));
         //t7 = t6*t6;
         //t8 = t7*t7;
         //t6 = pow((1.0-t4),(1.0/3.0));;
         //t11 = t10*t10;
         //t12 = t11*t11;
         //t13 = t8+t12-2.0;
         //t15 = 2**(1.0/3.0);
         //t17 = 1.0/(2.0*cr2-2.0);
         //t22 = 3**(1.0/3.0);
         //t24 = 4**(1.0/3.0);
         //t25 = cr4*cr4;
         //t26 = 1/pi;
         //t28 = pow((M_1_PI/t2),(1.0/3.0));;
         //t29 = t25*t28;
         //t34 = cr3*cr3;
         //t36 = t28*t28;
         //t37 = cr4*t36;
         t65 = t36*t36;
         //t70 = t1*(a0+da0*t13*t17+(a1+da1*t13*t17)*cr3*t29/4+(a2+da2*t13*t17)*t34*t37/4+3.0/4.0*(a3+da3*t13*t17)*M_1_PI*t3)/((b1+db1*t13*t17)*cr3*t29/4+(b2+db2*t13*t17)*t34*t37/4+3.0/4.0*(b3+db3*t13*t17)*M_1_PI*t3+3.0/16.0*(b4+db4*t13*t17)*cr3*t25*t65);
         t70 = t2*(a0+da0*t13*t17+(a1+da1*t13*t17)*cr3*t29/4+(a2+da2*t13*t17)*t34*t37/4+3.0/4.0*(a3+da3*t13*t17)*M_1_PI*t3)/((b1+db1*t13*t17)*cr3*t29/4+(b2+db2*t13*t17)*t34*t37/4+3.0/4.0*(b3+db3*t13*t17)*M_1_PI*t3+3.0/16.0*(b4+db4*t13*t17)*cr3*t25*t65);

         //t1 = xi*rp;
         //t2 = rp;
         //t3=fabs((t1+t2)/2.0);   //  *e2;
         //t4=fabs((t1-t2)/2.0);   //  *e2;
         t3=fabs((t2*(t4+1.0))/2.0);   //  *e2;
         t4=fabs((t2*(t4-1.0))/2.0);   //  *e2;
         t5= d_chpdft[ind]*t3+d_chpdft[ind+nxyz]*t4;
         
         d_ec[ind] = (-t70*e2 - 0.5*t5);
      }
}

extern "C" void calc_lda_gpu_(double *rho,double *chpdft,int *N,double *e2,double *ec)
{
        const int nxyz =*N;
	//cudaEvent_t start,stop,start_mem,stop_mem;
	//float int_time,int_time_mem;
        //double *d_chpdft;
	//double *d_ec;
        //double *d_rho;
        double e=*e2;
        int blocksize=512;
        dim3 dimgrid(ceil(nxyz/(float)blocksize),1,1);
        dim3 dimblock(blocksize,1,1);
        //double h_ec[nxyz];
        //unsigned int *d_ind;
        //unsigned int *h_ind;
	thrust::device_vector<double> d_ec(nxyz);
	//count_lda++;
	//thrust::host_vector<double> h_ec(nxyz);

        //h_ec=(double*)malloc(nxyz,sizeof(double));
        /*h_ind=(unsigned int*)calloc(nxyz,sizeof(unsigned int));

        cudaMalloc((void**)&d_ind,sizeof(unsigned int)*nxyz);
        cudaMemcpy(d_ind,h_ind,sizeof(unsigned int)*(nxyz),cudaMemcpyHostToDevice);*/
        //cudaMalloc((void**)&d_rho,sizeof(double)*(nxyz*2));
        //cudaMemcpy(d_rho,rho,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
        //cudaMalloc((void**)&d_chpdft,sizeof(double)*(nxyz*2));
        //cudaMemcpy(d_chpdft,chpdft,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
	/*cudaEventCreate(&start_mem);
	cudaEventCreate(&stop_mem);
	cudaEventRecord(start_mem,0);*/
        cudaMemcpyAsync(d_chpdft,rho,sizeof(double)*(nxyz*2),cudaMemcpyHostToDevice);
	Check_CUDA_Error(error);
	/*cudaEventRecord(stop_mem,0);
	cudaEventSynchronize(stop_mem);
	cudaEventElapsedTime(&int_time_mem,start_mem,stop_mem);
	time_mem=time_mem+int_time_mem;*/
        //cudaMalloc((void**)&d_ec,sizeof(double)*nxyz);
        //cudaMemcpy(d_ec,h_ec,sizeof(double)*(nxyz),cudaMemcpyHostToDevice);
	/*cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);*/
        lda<<<dimgrid,dimblock>>>(d_chpdft,raw_pointer_cast(&d_ec[0]),nxyz,e);
	Check_CUDA_Error(error);
        //cudaMemcpy(h_ind,d_ind, sizeof(unsigned int)*(nxyz),cudaMemcpyDeviceToHost);
        //cudaMemcpy(h_ec,d_ec, sizeof(double)*(nxyz),cudaMemcpyDeviceToHost);
	//h_ec=d_ec;
	*ec=thrust::reduce(d_ec.begin(),d_ec.end(),(double)0.0);
	/*cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&int_time,start,stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	time_lda=time_lda+int_time;
	if(count_lda%100==0) printf("Time lda: %f\n",time_lda);*/
	//cudaEventRecord(start_mem,0);
        cudaMemcpyAsync(chpdft,d_chpdft, sizeof(double)*(nxyz*2),cudaMemcpyDeviceToHost);
	Check_CUDA_Error(error);
	/*cudaEventRecord(stop_mem,0);
	cudaEventSynchronize(stop_mem);
	cudaEventElapsedTime(&int_time_mem,start_mem,stop_mem);
	time_mem=time_mem+int_time_mem;
	if(count_lda%100==0) printf("Time memcopy: %f\n",time_mem);
	cudaEventDestroy(start_mem);
	cudaEventDestroy(stop_mem);*/
        //cudaMemcpy(rho,d_rho, sizeof(double)*(nxyz*2),cudaMemcpyDeviceToHost);
        /*for (int ii=0;ii<nxyz;ii++){
          *ec+=h_ec[ii];
        }*/
        //cudaFree(d_rho);
        //cudaFree(d_chpdft);
        //cudaFree(d_ec);
}
#endif
