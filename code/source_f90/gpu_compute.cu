/*
This file is a part of PW-TELEMAN project.
PW-TELEMAN is a Time-Dependent Electronic Dynamics in Molecules And Nanosystems library.
Copyright (C) 2011-2015  Paul-Gerhard Reinhard, Eric Suraud, Florent Calvayrac,
Phuong Mai Dinh, David Brusson, Philipp Wopperer, José María Escartín Esteban.

PW-Teleman is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PW-Teleman is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PW-Teleman.  If not, see <http://www.gnu.org/licenses/>.
*/

//All sum and multiplication routines on GPU
#include"values_gpu.h"

__global__ void multiply_device(cufftDoubleComplex *d_ffta,int nxyz,double norm)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

	if (ind<nxyz)
	{
		d_ffta[ind].x=norm*d_ffta[ind].x;
		d_ffta[ind].y=norm*d_ffta[ind].y;
	}
}

__global__ void multiply_device_real(cufftDoubleReal *d_ffta,int nxyz,double norm)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

	if (ind<nxyz)
	{
		d_ffta[ind]=norm*d_ffta[ind];
	}
}

extern "C" void multiply_gpu_(cufftDoubleComplex *d_ffta,int *N,double *tnorm)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	double norm = *tnorm;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);

	//Multiplication d_ffta*norm on the GPU
#if(asynclaunch)
        multiply_device<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta,nxyz,norm);
#else
        multiply_device<<<dimgrid,dimblock>>>(d_ffta,nxyz,norm);
#endif
}

__global__ void multiply_ak_gpu(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ak,int nxyz, double tnorm)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;
	cufftDoubleReal SAVE2;

	if (ind<nxyz)
	{
		SAVE2         = d_ak[ind].x*d_ffta[ind].x+d_ak[ind].y*d_ffta[ind].y;
		d_ffta[ind].y = (d_ak[ind].x*d_ffta[ind].y+d_ak[ind].y*d_ffta[ind].x)*tnorm;
		d_ffta[ind].x = SAVE2*tnorm;
	}
}

extern "C" void multiply_ak_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ak,int *N)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);

	//Multiplication d_ffta*d_ak on the GPU
#if(asynclaunch)
        multiply_ak_gpu<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta,d_ak,nxyz,1.0);
#else
        multiply_ak_gpu<<<dimgrid,dimblock>>>(d_ffta,d_ak,nxyz,1.0);
#endif
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

extern "C" void multiply_ak2_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ak,int *N)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);

	//Multiplication d_ffta*d_ak on the GPU
#if(asynclaunch)
        multiply_ak_gpu2<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta,d_ak,nxyz);
#else
        multiply_ak_gpu2<<<dimgrid,dimblock>>>(d_ffta,d_ak,nxyz);
#endif
}

__global__ void multiply_rak_gpu2(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ak,int nxyz)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;
	cufftDoubleReal SAVE2;

	if (ind<nxyz)
	{
		SAVE2         = d_ak[ind].y*d_ffta[ind].y-d_ak[ind].x*d_ffta[ind].x;
		d_ffta[ind].y = (-1.0)*(d_ak[ind].x*d_ffta[ind].y+d_ak[ind].y*d_ffta[ind].x);
		d_ffta[ind].x = SAVE2;
	}
}

extern "C" void multiply_rak2_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ak,int *N)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);

	//Multiplication d_ffta*d_ak on the GPU
#if(asynclaunch)
        multiply_rak_gpu2<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta,d_ak,nxyz);
#else
        multiply_rak_gpu2<<<dimgrid,dimblock>>>(d_ffta,d_ak,nxyz);
#endif
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

extern "C" void multiply_ak_real_(cufftDoubleComplex *d_ffta,cufftDoubleReal *d_ak,int *N)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);

	//Multiplication d_ffta*d_ak on the GPU
#if(asynclaunch)
        multiply_ak_real_gpu<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta,d_ak,nxyz);
#else
        multiply_ak_real_gpu<<<dimgrid,dimblock>>>(d_ffta,d_ak,nxyz);
#endif
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

extern "C" void multiply_shift_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_akx,cufftDoubleComplex *d_aky,cufftDoubleComplex *d_akz,double *sx,double *sy,double *sz,int *N)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	cufftDoubleComplex shix,shiy,shiz;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);

	shix.x = *sx;
	shix.y = *sx;
	shiy.x = *sy;
	shiy.y = *sy;
	shiz.x = *sz;
	shiz.y = *sz;

	//Multiplication on the GPU
#if(asynclaunch)
        multiply_shift_gpu<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta,d_akx,d_aky,d_akz,shix,shiy,shiz,nxyz);
#else
        multiply_shift_gpu<<<dimgrid,dimblock>>>(d_ffta,d_akx,d_aky,d_akz,shix,shiy,shiz,nxyz);
#endif
}

__global__ void hpsi_gpu(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,cufftDoubleReal *d_akv,int nxyz, int kdfull2)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;
	unsigned int nbe = ind/kdfull2;

	if (ind<nxyz)
	{
		d_ffta2[ind].x = d_akv[ind-(nbe*kdfull2)]*d_ffta[ind].x+d_ffta2[ind].x;
		d_ffta2[ind].y = d_akv[ind-(nbe*kdfull2)]*d_ffta[ind].y+d_ffta2[ind].y;
	}
}

extern "C" void hpsi_cuda_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,cufftDoubleReal *d_akv,int *N,int *kd)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	int kdfull2=*kd;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);

	//d_ffta*d_ak+d_ffta2 on the GPU
#if(asynclaunch)
        hpsi_gpu<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta,d_ffta2,d_akv,nxyz,kdfull2);
#else
        hpsi_gpu<<<dimgrid,dimblock>>>(d_ffta,d_ffta2,d_akv,nxyz,kdfull2);
#endif
	Check_CUDA_Error(error);
}

__global__ void d_grad_gpu1(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,cufftDoubleReal *d_akv,double epswf,double e0dmp,int nxyz)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

	if (ind<nxyz)
	{
		d_ffta[ind].x = d_ffta[ind].x-((epswf*d_ffta2[ind].x)/(d_akv[ind]+e0dmp));
		d_ffta[ind].y = d_ffta[ind].y-((epswf*d_ffta2[ind].y)/(d_akv[ind]+e0dmp));
	}
}

extern "C" void d_grad1_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,cufftDoubleReal *d_akv,double *ep,double *e0,int *N,int *nbe)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	double epswf = *ep, e0dmp = *e0;
	int kstate=*nbe-1;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);

	//grad on the GPU
#if(asynclaunch)
        d_grad_gpu1<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta+kstate*nxyz,d_ffta2+kstate*nxyz,d_akv,epswf,e0dmp,nxyz);
#else
        d_grad_gpu1<<<dimgrid,dimblock>>>(d_ffta+kstate*nxyz,d_ffta2+kstate*nxyz,d_akv,epswf,e0dmp,nxyz);
#endif
	Check_CUDA_Error(error);
}

__global__ void d_grad_gpu2(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,double epswf,int nxyz)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

	if (ind<nxyz)
	{
		d_ffta[ind].x = d_ffta[ind].x-epswf*d_ffta2[ind].x;
		d_ffta[ind].y = d_ffta[ind].y-epswf*d_ffta2[ind].y;
	}
}

extern "C" void d_grad2_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,double *ep,int *N,int *nbe)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	double epswf = *ep;
	int kstate=*nbe-1;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);

	//grad on the GPU
#if(asynclaunch)
        d_grad_gpu2<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta+kstate*nxyz,d_ffta2+kstate*nxyz,epswf,nxyz);
#else
        d_grad_gpu2<<<dimgrid,dimblock>>>(d_ffta+kstate*nxyz,d_ffta2+kstate*nxyz,epswf,nxyz);
#endif
	Check_CUDA_Error(error);
}

__global__ void d_sum_calc(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,cufftDoubleReal *d_akv,double *d_sum0,double *d_sumk,double *d_sume,double *d_sum2,int nxyz)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

	if (ind<nxyz)
	{
		d_sum0[ind] = d_ffta[ind].x*d_ffta[ind].x + d_ffta[ind].y*d_ffta[ind].y;
		d_sumk[ind] = d_sum0[ind]*d_akv[ind];
		d_sume[ind] = d_ffta2[ind].x*d_ffta[ind].x + d_ffta2[ind].y*d_ffta[ind].y;
		d_sum2[ind] = d_ffta2[ind].x*d_ffta2[ind].x + d_ffta2[ind].y*d_ffta2[ind].y;
	}
}

extern "C" void sum_calc_(double *s0,double *sk,double *se,double *s2,cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2, cufftDoubleReal *d_akv,int *N,int *nb)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	int nbe  = *nb-1;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);
	//Declaration of the vectors on the GPU
	thrust::device_vector<double> d_sum0(nxyz);
	thrust::device_vector<double> d_sumk(nxyz);
	thrust::device_vector<double> d_sume(nxyz);
	thrust::device_vector<double> d_sum2(nxyz);

	//Computation of the d_sum vectors on the GPU
#if(asynclaunch)
        d_sum_calc<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta+nbe*nxyz,d_ffta2+nbe*nxyz,d_akv,raw_pointer_cast(&d_sum0[0]),raw_pointer_cast(&d_sumk[0]),raw_pointer_cast(&d_sume[0]),raw_pointer_cast(&d_sum2[0]),nxyz);
#else
        d_sum_calc<<<dimgrid,dimblock>>>(d_ffta+nbe*nxyz,d_ffta2+nbe*nxyz,d_akv,raw_pointer_cast(&d_sum0[0]),raw_pointer_cast(&d_sumk[0]),raw_pointer_cast(&d_sume[0]),raw_pointer_cast(&d_sum2[0]),nxyz);
#endif

	//Reduction of the vectors
	*s0=thrust::reduce(d_sum0.begin(),d_sum0.end(),(double)0.0);
	*sk=thrust::reduce(d_sumk.begin(),d_sumk.end(),(double)0.0);
	*se=thrust::reduce(d_sume.begin(),d_sume.end(),(double)0.0);
	*s2=thrust::reduce(d_sum2.begin(),d_sum2.end(),(double)0.0);

}

__global__ void d_sum_calc2(cufftDoubleComplex *d_ffta,cufftDoubleReal *d_akv,double *d_sum0,double *d_sumk,int nxyz)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

	if (ind<nxyz)
	{
		d_sum0[ind] = d_ffta[ind].x*d_ffta[ind].x + d_ffta[ind].y*d_ffta[ind].y;
		d_sumk[ind] = d_sum0[ind]*d_akv[ind];
	}
}

extern "C" void sum_calc2_(double *s0,double *sk,cufftDoubleComplex *d_ffta,cufftDoubleReal *d_akv,int *N)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nxyz = *N;
	int blocksize=BLOCK1;int warpSize=WARPSIZE;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	;dim3 dimgrid(gridx+1,GRIDY,GRIDZ);
	dim3 dimblock(blocksize,BLOCK2,BLOCK3);
	//Declaration of the vectors on the GPU
	thrust::device_vector<double> d_sum0(nxyz);
	thrust::device_vector<double> d_sumk(nxyz);
#if(asynclaunch)
        d_sum_calc2<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(d_ffta,d_akv,raw_pointer_cast(&d_sum0[0]),raw_pointer_cast(&d_sumk[0]),nxyz);
#else
        d_sum_calc2<<<dimgrid,dimblock>>>(d_ffta,d_akv,raw_pointer_cast(&d_sum0[0]),raw_pointer_cast(&d_sumk[0]),nxyz);
#endif

	*s0=thrust::reduce(d_sum0.begin(),d_sum0.end(),(double)0.0);
	*sk=thrust::reduce(d_sumk.begin(),d_sumk.end(),(double)0.0);
}

__global__ void d_build_k(int nx,int ny,int nz,double h2m,double dt1,double dkx,double dky,double dkz,cufftDoubleComplex *d_ak,cufftDoubleReal *d_akv,cufftDoubleComplex *d_akx,cufftDoubleComplex *d_aky,cufftDoubleComplex *d_akz)
{
	int i1 = blockIdx.x*blockDim.x+threadIdx.x;
	int i2 = blockIdx.y*blockDim.y+threadIdx.y;
	int i3 = blockIdx.z*blockDim.z+threadIdx.z;
	int ind  = i3*nx*ny+i2*nx+i1;
	double zkx,zky,zkz;

	if((i1<nx) && (i2<ny) && (i3<nz))
	{
		if(i3 >= nz/2) zkz=(i3-nz)*dkz;
		else zkz=i3*dkz;

		if(i2 >= ny/2) zky=(i2-ny)*dky;
	    	else zky=i2*dky;

		if(i1 >= nx/2) zkx=(i1-nz)*dkx;
	      	else zkx=i1*dkx;
      
		d_ak[ind].x=cos(dt1*(zkx*zkx+zky*zky+zkz*zkz)*h2m);
		d_ak[ind].y=-sin(dt1*(zkx*zkx+zky*zky+zkz*zkz)*h2m);
		d_akv[ind]=(zkx*zkx+zky*zky+zkz*zkz)*h2m;
		d_akx[ind].x=0.0;
		d_akx[ind].y=-zkx;
		d_aky[ind].x=0.0;
		d_aky[ind].y=-zky;
		d_akz[ind].x=0.0;
		d_akz[ind].y=-zkz;
	}
}

extern "C" void build_kgpu_(int *NX,int *NY,int *NZ,double *hm,double *dt,double *dx,double *dy,double *dz,cufftDoubleComplex *d_ak,cufftDoubleReal *d_akv,cufftDoubleComplex *d_akx,cufftDoubleComplex *d_aky,cufftDoubleComplex *d_akz)
{
#if(parayes)
	cudaSetDevice(params_mp_mygpu_);
#endif
	int nx=*NX,ny=*NY,nz=*NZ;
	double h2m=*hm,dt1=*dt,dkx=*dx,dky=*dy,dkz=*dz;
	int blocksize=BLOCK1B;
	int gridx=(int)ceil(nx/(float)blocksize),gridy=(int)ceil(ny/(float)blocksize),gridz=(int)ceil(nz/(float)blocksize);
	dim3 dimgrid(gridx+1,gridy,gridz);
	dim3 dimblock(blocksize,blocksize,blocksize);
#if(asynclaunch)
	d_build_k<<<dimgrid,dimblock,0,stream2[params_mp_mygpu_]>>>(nx,ny,nz,h2m,dt1,dkx,dky,dkz,d_ak,d_akv,d_akx,d_aky,d_akz);
#else
	d_build_k<<<dimgrid,dimblock>>>(nx,ny,nz,h2m,dt1,dkx,dky,dkz,d_ak,d_akv,d_akx,d_aky,d_akz);
#endif
}
