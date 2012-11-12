// Conversion from fortran to C++/Cuda still in progressC

unsigned int kxmax,kymax,kzmax;
// kxmax must be the largest
unsigned int kdfull;
unsigned int kdred;
double dx,dy,dz,dxsp,grnorm,fnorm,tnorm;
unsigned int nxc,nyc,nzc;
unsigned int nxi,nyi,nzi;
unsigned int nxy1;

double dkx,dky,dkz,akmax,dksp,ecut;
int *inde,*pindex;

double zero=0.0;
double pi=3.141592653589793;

cufftHandle pfftforw,pfftback;
int batch=1;
cufftDoubleReal *rfftac,*akvr,*prfftac,*pakvr;
cufftDoubleComplex *gpu_fftac,*gpu_akvc;
cufftDoubleReal *gpu_rfftac,*gpu_akvr;

int res;

//-----fourfakv---------------------------------------------------------

void fourfakv(){

//     fourier forward transformation
//     I/O: pskr   real part of the wave-function
//          pski   imaginary part of the wave-function

//DATA  mini/0/              // flag for initialization
//----------------------------------------------------------------------

int blocksize=192;
int gridx=(int)ceil(kdred/(float)blocksize);
dim3 dimgrid(gridx,1,1);
dim3 dimblock(blocksize,1,1);

//test      sqh=sqrt(0.5)

tnorm=grnorm*fnorm;

cudaMemcpyAsync(gpu_akvr,pakvr,kdred*sizeof(cufftDoubleReal),cudaMemcpyHostToDevice,stream1);
Check_CUDA_Error(error);

if(cufftExecD2Z(pfftforw,gpu_akvr,gpu_akvc) != CUFFT_SUCCESS)
{
  cout<<"CUFFT error : Exec Z2Z forward failed in coulex (akv)"<<endl;
  exit(-1);
}

multiply_device<<<dimgrid,dimblock,0,stream3>>>(gpu_akvc,kdred,tnorm);
Check_CUDA_Error(error);

}

//-----fftinp------------------------------------------------------------

void fftinp() {

//     initializes work tables for FFT

//     grid parameters nxc,nyc,nzc,dx,dy,dz,ecut must have been read or
//     initialized before !

//-----------------------------------------------------------------------

int ikzero,ii;
double xz1,xz2,xy1,xy2,xx1,xx2,ak2;
int nxyfn,ind;

//test      sqh=sqrt(0.5)

//     initialize grid in coordinate space

nxi=nxc+nxc;
nyi=nyc+nyc;
nzi=nzc+nzc;
nxy1=nxi*nyi;
nxyfn = kxmax*kzmax;

//     grid lengths must match with parameters in incs

if(kxmax < nxi) {
  cout<< "ERROR: parameter   kxmax   too small"<<endl;
  cout<< " error in parameter: KXMAX in COULEX too small"<<endl;
  exit(-1);
}
else if(kymax < nyi) {
  cout<< " ERROR: parameter   kymax   too small"<<endl;
  cout<< " error in parameter: KYMAX in COULEX too small"<<endl;
  exit(-1);
}
else if(kzmax < nzi) {
  cout<< " ERROR: parameter   kzmax   too small"<<endl;
  cout<< " error in parameter: KZMAX in COULEX too small"<<endl;
  exit(-1);
}

//     initialize grid in fourier space

dkx=pi/(dx*double(nxc));
dky=pi/(dy*double(nyc));
dkz=pi/(dz*double(nzc));

dxsp=dx*dy*dz;
dksp=dkx*dky*dkz;
cout<< " dkx,dky,dkz,dksp= "<<dkx<<" "<<dky<<" "<<dkz<<" "<<dksp<<endl;

grnorm=sqrt(dxsp/dksp);
fnorm=1.0/sqrt(double(nxc*nyc*nzc));

//     built Greens function in Fourier space
//     by Fourier transformation from real space

ikzero = nxy1*(nzc-1)+nxi*(nyc-1)+nxc;
cout<< " nzi,nyi,nxi,nx,ny,nz,ikzero= "<<nzi<<" "<<nyi<<" "<<nxi<<" "<<nxc<<" "<<nyc<<" "<<nzc<<" "<<ikzero<<endl;
ii=0;
xz1=-(1.0)*nzc*dz;
for (int i3=1;i3<=nzi;i3++){
  xz1=xz1+dz;
  xz2=xz1*xz1;
  xy1=-(1.0)*nyc*dy;
  for (int i2=1;i2<=nyi;i2++){
    xy1=xy1+dy;
    xy2=xy1*xy1;
    xx1=-(1.0)*nxc*dx;
    for (int i1=1;i1<=nxi;i1++){
      xx1=xx1+dx;
      xx2=xx1*xx1;
      ak2=xx2+xy2+xz2;
      ii=ii+1;
//        cout<< " i1,i2,i3,ii= "<<i1<<" "<<i2<<" "<<i3<<" "ii;
      ind=((i3+nxc)%nxi)*nxyfn+((i2+nyc)%nyi)*kxmax+(i1+nyc)%nzi+1; //storage of the indices of the flatten 3D array
      if(ii != ikzero) {
        akvr[ind] =  1.0/sqrt(ak2);
      }
      else {
//              akvr[inde[ii]] = (6D0*pi/(dx*dy*dz))**(1D0/3D0)  // spherical approx
//              akvr[inde[ii]] = 1.19003868*(dx*dy*dz)**(-1D0/3D0)
        akvr[ind] = 2.34*1.19003868*pow((dx*dy*dz),(-1.0/3.0));  // empirical
      }
    }
  }
}

ii=0;
for (int i3=1;i3<=nzc;i3++){
  for (int i2=1;i2<=nyc;i2++){
    for (int i1=1;i1<=nxc;i1++){
      ii++;
      inde[ii]=((i3+nxc)%nxi)*nxyfn+((i2+nyc)%nyi)*kxmax+(i1+nzc)%nzi+1;
    }
  }
}

fourfakv();

cudaFreeHost(pakvr);
cudaFree(gpu_akvr);

}

//-------------------------------------------------------------------

extern "C" void init_coul_(double *dx0,double *dy0,double *dz0,unsigned int *nx0,unsigned int *ny0,unsigned int *nz0) {

//-----------------------------------------------------------------------

//     read grid parameters from file or simply initialize them
//     note that the Coulomb solver doubles the grid internally
nxc=*nx0;  ///2;
nyc=*ny0;  ///2;
nzc=*nz0;  ///2;
dx=*dx0;
dy=*dy0;
dz=*dz0;

kxmax=2*nxc;kymax=2*nyc;kzmax=2*nzc;
kdfull=nxc*nyc*nzc;
kdred=kxmax*kymax*kzmax;

//     check initialization

if(cufftPlan3d(&pfftforw,kxmax,kymax,kzmax,CUFFT_D2Z) != CUFFT_SUCCESS)
{
  cout<<"CUFFT error : Plan Creation failed"<<endl;
  exit(-1);
}
if(cufftSetStream(pfftforw,stream2) != CUFFT_SUCCESS)
{
  cout<<"CUFFT error : Streamed FFT Creation failed"<<endl;
  exit(-1);
}
if(cufftPlan3d(&pfftback,kxmax,kymax,kzmax,CUFFT_Z2D) != CUFFT_SUCCESS)
{
  cout<<"CUFFT error : Plan Creation failed"<<endl;
  exit(-1);
}
if(cufftSetStream(pfftback,stream1) != CUFFT_SUCCESS)
{
  cout<<"CUFFT error : Streamed FFT Creation failed"<<endl;
  exit(-1);
}
// Pinned memory allocation on the CPU to make CPU>GPU and GPU>CPU transfers faster

cudaMallocHost (&prfftac,kdred*sizeof(cufftDoubleReal));
cudaMallocHost (&pakvr,kdred*sizeof(cufftDoubleReal));

rfftac=prfftac-1; //rfftac points one location before prfftac, so rfftac[1]...rfftac[kdred] all exist (just sticks with the fortran convention)
akvr=pakvr-1;   //same trick as above

// Memory allocation on the GPU
cudaMalloc((void**)&gpu_fftac,kdred*sizeof(cufftDoubleComplex));
cudaMalloc((void**)&gpu_akvc,kdred*sizeof(cufftDoubleComplex));
cudaMalloc((void**)&gpu_rfftac,kdred*sizeof(cufftDoubleReal));
cudaMalloc((void**)&gpu_akvr,kdred*sizeof(cufftDoubleReal));

pindex = (int*)malloc(kdfull*sizeof(int));
inde=pindex-1;

//     call input routine fftinp, which initializes the grid and fft tabl

fftinp();
}

//-----cofows------------------------------------------------------------

void coufou2(){

int blocksize=192;
int gridx=(int)ceil(kdred/(float)blocksize);
dim3 dimgrid(gridx,1,1);
dim3 dimblock(blocksize,1,1);

//------------------------------------------------------------------------------

//     fourier transformation of the density

tnorm=grnorm*fnorm;

cudaMemcpyAsync(gpu_rfftac,prfftac,kdred*sizeof(cufftDoubleReal),cudaMemcpyHostToDevice,stream1);
Check_CUDA_Error(error);

if(cufftExecD2Z(pfftforw,gpu_rfftac,gpu_fftac) != CUFFT_SUCCESS)
{
  cout<<"CUFFT error : Exec Z2Z forward failed in coulex"<<endl;
  exit(-1);
}

//     calculation of the coulomb field (writing on the density field)

multiply_ak_gpu<<<dimgrid,dimblock,0,stream3>>>(gpu_fftac,gpu_akvc,kdred,tnorm);
Check_CUDA_Error(error);

//     fourier back transformation

tnorm=fnorm/(8.0*grnorm)*pow(pi,1.5);

if(cufftExecZ2D(pfftback,gpu_fftac,gpu_rfftac) != CUFFT_SUCCESS)
{
	  cout<<"CUFFT error : Exec Z2Z backward failed in coulex"<<endl;
	  exit(-1);
}

multiply_device_real<<<dimgrid,dimblock,0,stream2>>>(gpu_rfftac,kdred,tnorm);
Check_CUDA_Error(error);

cudaMemcpy(prfftac,gpu_rfftac,kdred*sizeof(cufftDoubleReal),cudaMemcpyDeviceToHost);
Check_CUDA_Error(error);

}

extern "C" void falr_(double *prhoinp,double *pchpfalr,int nxdum,int nydum,int nzdum,int kdum) {

double *rhoinp,*chpfalr;

rhoinp=prhoinp-1; //rhoinp points one location before pfftac, so rhoinp[1]...rhoinp[kdred] all exist (just sticks with the fortran convention)
chpfalr=pchpfalr-1; //same trick as above

//     Write your density field on the array rfftac.
//     remember not to send your original density array to the fcs.
//     in this case we have a homogeneously charged sphere .

for (int ii=1;ii<=kdred;ii++) rfftac[ii]=0.0;
for (int ii=1;ii<=kdfull;ii++) rfftac[inde[ii]]=rhoinp[ii];

//     call coufou, which contains the fcs procedure.
coufou2();

//     Write the result of the fcs on the array chpfalr

for (int i0=1;i0<=kdfull;i0++) chpfalr[i0] = 2.0*rfftac[inde[i0]];

}

extern "C" void coulsolv_end_() {

cudaFreeHost(prfftac);
cudaFree(gpu_fftac);
cudaFree(gpu_rfftac);
cudaFree(gpu_akvc);
cufftDestroy(pfftforw);
cufftDestroy(pfftback);
free(pindex);

}
