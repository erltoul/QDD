//All sum and multiplication routines on GPU

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
	int nxyz = *N;
	double norm = *tnorm;
	int blocksize=192;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

	//Multiplication d_ffta*norm on the GPU
        multiply_device<<<dimgrid,dimblock,0,stream2>>>(d_ffta,nxyz,norm);
	Check_CUDA_Error(error);
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
	int nxyz = *N;
	int blocksize=192;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

	//Multiplication d_ffta*d_ak on the GPU
        multiply_ak_gpu<<<dimgrid,dimblock,0,stream2>>>(d_ffta,d_ak,nxyz,1.0);
	Check_CUDA_Error(error);
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
	int nxyz = *N;
	int blocksize=192;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

	//Multiplication d_ffta*d_ak on the GPU
        multiply_ak_gpu2<<<dimgrid,dimblock,0,stream2>>>(d_ffta,d_ak,nxyz);
	Check_CUDA_Error(error);
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
	int nxyz = *N;
	int blocksize=192;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

	//Multiplication d_ffta*d_ak on the GPU
        multiply_ak_real_gpu<<<dimgrid,dimblock,0,stream2>>>(d_ffta,d_ak,nxyz);
	Check_CUDA_Error(error);
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
	int nxyz = *N;
	cufftDoubleComplex shix,shiy,shiz;
	int blocksize=192;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

	shix.x = *sx;
	shix.y = *sx;
	shiy.x = *sy;
	shiy.y = *sy;
	shiz.x = *sz;
	shiz.y = *sz;

	//Multiplication on the GPU
        multiply_shift_gpu<<<dimgrid,dimblock,0,stream2>>>(d_ffta,d_akx,d_aky,d_akz,shix,shiy,shiz,nxyz);
	Check_CUDA_Error(error);
}

__global__ void hpsi_gpu(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,cufftDoubleReal *d_akv,int nxyz)
{
	unsigned int ind = blockIdx.x*blockDim.x+threadIdx.x;

	if (ind<nxyz)
	{
		d_ffta2[ind].x = d_akv[ind]*d_ffta[ind].x+d_ffta2[ind].x;
		d_ffta2[ind].y = d_akv[ind]*d_ffta[ind].y+d_ffta2[ind].y;
	}
}

extern "C" void hpsi_cuda_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,cufftDoubleReal *d_akv,int *N)
{
	int nxyz = *N;
	int blocksize=192;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

	//d_ffta*d_ak+d_ffta2 on the GPU
        hpsi_gpu<<<dimgrid,dimblock,0,stream2>>>(d_ffta,d_ffta2,d_akv,nxyz);
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

extern "C" void d_grad1_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,cufftDoubleReal *d_akv,double *ep,double *e0,int *N)
{
	int nxyz = *N;
	double epswf = *ep, e0dmp = *e0;
	int blocksize=192;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

	//grad on the GPU
        d_grad_gpu1<<<dimgrid,dimblock,0,stream2>>>(d_ffta,d_ffta2,d_akv,epswf,e0dmp,nxyz);
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

extern "C" void d_grad2_(cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2,double *ep,int *N)
{
	int nxyz = *N;
	double epswf = *ep;
	int blocksize=192;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);

	//grad on the GPU
        d_grad_gpu2<<<dimgrid,dimblock,0,stream2>>>(d_ffta,d_ffta2,epswf,nxyz);
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

extern "C" void sum_calc_(double *s0,double *sk,double *se,double *s2,cufftDoubleComplex *d_ffta,cufftDoubleComplex *d_ffta2, cufftDoubleReal *d_akv,int *N)
{
	int nxyz = *N;
	int blocksize=192;
	int gridx=(int)ceil(nxyz/(float)blocksize);
	dim3 dimgrid(gridx,1,1);
	dim3 dimblock(blocksize,1,1);
	//Declaration of the vectors on the GPU
	thrust::device_vector<double> d_sum0(nxyz);
	thrust::device_vector<double> d_sumk(nxyz);
	thrust::device_vector<double> d_sume(nxyz);
	thrust::device_vector<double> d_sum2(nxyz);

	//Computation of the d_sum vectors on the GPU
        d_sum_calc<<<dimgrid,dimblock,0,stream2>>>(d_ffta,d_ffta2,d_akv,raw_pointer_cast(&d_sum0[0]),raw_pointer_cast(&d_sumk[0]),raw_pointer_cast(&d_sume[0]),raw_pointer_cast(&d_sum2[0]),nxyz);
	Check_CUDA_Error(error);
	//Reduction of the vectors
	*s0=thrust::reduce(d_sum0.begin(),d_sum0.end(),(double)0.0);
	*sk=thrust::reduce(d_sumk.begin(),d_sumk.end(),(double)0.0);
	*se=thrust::reduce(d_sume.begin(),d_sume.end(),(double)0.0);
	*s2=thrust::reduce(d_sum2.begin(),d_sum2.end(),(double)0.0);
}
