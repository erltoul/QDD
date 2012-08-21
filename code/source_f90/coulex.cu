// Conversion from fortran to C++/Cuda still in progressC

unsigned int kxmax,kymax,kzmax,ksmax;
// kxmax must be the largest
unsigned int kdfull;
unsigned int kdred;
unsigned int kfft2;
unsigned int kfft,kfftx,kffty,kfftz;
double dx,dy,dz,dxsp,grnorm,fnorm,tnorm;
unsigned int nxc,nyc,nzc,nx1,ny1,nz1;
unsigned int nxi,nyi,nzi;
unsigned int nxy1,nxyz,nxhigh;
unsigned int nxlow,nyhigh,nylow,nzhigh,nzlow;

double dkx,dky,dkz,akmax,dksp,ecut;
unsigned int nxk,nxklo,nxkhi,nksp,nkxyz;

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
int gridx=(int)ceil(nxyz/(float)blocksize);
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

multiply_device<<<dimgrid,dimblock,0,stream2>>>(gpu_akvc,kdred,tnorm);
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

nx1=nxc+1;
ny1=nyc+1;
nz1=nzc+1;
nxi=nxc+nxc;
nyi=nyc+nyc;
nzi=nzc+nzc;
nxy1=nxi*nyi;
nxyz=nxi*nyi*nzi;
nkxyz=nxi*nyi*nzi;

nxyfn = kfftx*kfftz;

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
nxk=nx1;

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
      ind=((i3+nxc)%nxi)*nxyfn+((i2+nyc)%nyi)*kfftx+(i1+nyc)%nzi+1; //storage in a flatten 3D complex array for FFT on GPU
      if(ii != ikzero) {
        akvr[ind] =  1.0/sqrt(ak2);
      }
      else {
//              akv2r(ii) = (6D0*pi/(dx*dy*dz))**(1D0/3D0)  // spherical approx
//              akv2r(ii) = 1.19003868*(dx*dy*dz)**(-1D0/3D0)
        akvr[ind] = 2.34*1.19003868*pow((dx*dy*dz),(-1.0/3.0));  // empirical
      }
    }
  }
}
nksp=ii;

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

kxmax=2*nxc;kymax=2*nyc;kzmax=2*nzc;ksmax=kxmax;
kdfull=nxc*nyc*nzc;
kdred=kxmax*kymax*kzmax;
kfft=ksmax;kfftx=kxmax;kffty=kymax;kfftz=kzmax;
kfft2=kfft*2;

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
if(cufftSetStream(pfftback,stream2) != CUFFT_SUCCESS)
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

//     call input routine fftinp, which initializes the grid and fft tabl

fftinp();
}

//-----fft--------------------------------------------------------------

void fft() {

int blocksize=192;
int gridx=(int)ceil(nxyz/(float)blocksize);
dim3 dimgrid(gridx,1,1);
dim3 dimblock(blocksize,1,1);

tnorm=grnorm*fnorm;

cudaMemcpyAsync(gpu_rfftac,prfftac,kdred*sizeof(cufftDoubleReal),cudaMemcpyHostToDevice,stream1);
Check_CUDA_Error(error);

if(cufftExecD2Z(pfftforw,gpu_rfftac,gpu_fftac) != CUFFT_SUCCESS)
{
  cout<<"CUFFT error : Exec Z2Z forward failed in coulex"<<endl;
  exit(-1);
}

multiply_device<<<dimgrid,dimblock,0,stream2>>>(gpu_fftac,nxyz,tnorm);
Check_CUDA_Error(error);

}

//-----ffb--------------------------------------------------------------

void ffb() {

//----------------------------------------------------------------------

int blocksize=192;
int gridx=(int)ceil(nxyz/(float)blocksize);
dim3 dimgrid(gridx,1,1);
dim3 dimblock(blocksize,1,1);

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

//-----cofows------------------------------------------------------------

void coufou2(){

int blocksize=192;
int gridx=(int)ceil(nxyz/(float)blocksize);
dim3 dimgrid(gridx,1,1);
dim3 dimblock(blocksize,1,1);

//------------------------------------------------------------------------------

//     fourier transformation of the density

fft();

//     calculation of the coulomb field (writing on the density field)

multiply_ak_gpu<<<dimgrid,dimblock,0,stream2>>>(gpu_fftac,gpu_akvc,kdred);
Check_CUDA_Error(error);

//     fourier back transformation

ffb();

}

//-----rhofld------------------------------------------------------------

void rhofld(double *rhoinp){

//     copy density on complex array of double extnesion in x,y,z

int nxyfn,nyfn,nnx2,nny2,nnz2,i0,ii;

nxyfn = kfftx*kfftz;
nyfn  = kfftx;
nnx2=nxc+nxc;
nny2=nyc+nyc;
nnz2=nzc+nzc;

i0=0;
for (int i3=1;i3<=kfftz;i3++){
  for (int i2=1;i2<=kffty;i2++){
    for (int i1=1;i1<=kfftx;i1++){
      ii=((i3+nxc)%nnx2)*nxyfn+((i2+nyc)%nny2)*nyfn+(i1+nzc)%nnz2+1;
      if(i3 <= nzc && i2 <= nyc && i1 <= nxc) {
        i0 = i0+1;
        rfftac[ii]=rhoinp[i0];
      }
      else rfftac[ii]=0.0;
    }
  }
}

}


//-----result------------------------------------------------------------

void result(double *chpfalr){

//     copy Coulomb field back to standard grid
int nxyfn,nyfn,nnx2,nny2,nnz2,ii;

int i0=0;

nxyfn = kfftx*kfftz;
nyfn  = kfftx;
nnx2=nxc+nxc;
nny2=nyc+nyc;
nnz2=nzc+nzc;

for (int i3=1;i3<=nzc;i3++){
  for (int i2=1;i2<=nyc;i2++){
    for (int i1=1;i1<=nxc;i1++){
        i0++;
        ii=((i3+nxc)%nnx2)*nxyfn+((i2+nyc)%nny2)*nyfn+(i1+nzc)%nnz2+1;
        chpfalr[i0] = 2.0*rfftac[ii];
    }
  }
}

}

//-------------------------------------------------------------------

extern "C" void falr_(double *prhoinp,double *pchpfalr,int nxdum,int nydum,int nzdum,int kdum) {

double *rhoinp,*chpfalr;

rhoinp=prhoinp-1; //rhoinp points one location before pfftac, so rhoinp[1]...rhoinp[kdred] all exist (just sticks with the fortran convention)
chpfalr=pchpfalr-1; //same trick as above

//     call a routine written by you which writes your density field
//     on the array rho.
//     remember not to send your original density array to the fcs.
//     in this case we have a homogeneously charged sphere .

rhofld(rhoinp);

//     call coufou, which contains the fcs procedure.
coufou2();

//     call a routine written by you which outputs the results of the fcs
//     and maybe some other things to an output file or the screen.

result(chpfalr);
}

extern "C" void coulsolv_end_() {

cudaFreeHost(prfftac);
cudaFree(gpu_fftac);
cudaFree(gpu_rfftac);
cudaFree(gpu_akvc);
cufftDestroy(pfftforw);
cufftDestroy(pfftback);

}
