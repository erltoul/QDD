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

// Conversion from fortran to C++/Cuda still in progress

extern "C" { void static_mfield_( double rho[], double aloc[] , double psir[], double qaux[] ); }

extern "C" { void calcrhor_( double rho[], double psir[] ); } 

extern "C" { void resume_( double psir[], char ); } 

extern "C" { void calcpseudo_(); } 

extern "C" { void pricm_( double rho[] ); } 

extern "C" { void infor_( double psir[], double rho[] ); } 

extern "C" { void sstep_( double psir[], double aloc[] ); }

extern "C" { void flush_(int *); }

extern "C" { void localizer_( double rho[], double psir[] ); }

extern "C" { void pri_pstat_( double psir[], double rho[] ); }

extern "C" { void printfield_(int , double [] , const char* ); }

extern "C" { void rsave_( double psir[] ); }

extern "C" { void transel_( double psir[] ); }

extern "C" { void mtv_fld_( double [] , int ); }

extern "C" { void printone_( double rho[], double aloc[] ); }

extern "C" { void printtwo_( ); }

extern "C" { void printthree_( double rho[], double aloc[] ); }

extern "C" { void print_densdiffc_( double rho[] ); }

extern "C" { void print_orb_( double psir[] ); }

extern "C" { void print_surf_( double rho[] ); }

#if(twostsic)
extern "C" { void infor_sic_( double psir[] ); }

extern "C" { void diag_lagr_( double psir[] ); }
#endif

#if(fullsic)
extern "C" { void sicstep_gen_( double psir[], double qaux[], double aloc[] ); }

extern "C" { void calc_locwf_( double psir[], double qaux[] ); }

extern "C" { void calc_fullsicr_( double psir[], double qaux[] ); }
#endif

#if(parayes)
double is[mpi_status_size];
#endif

bool tp_prints = false;
double accsum;
int ifsicpsav;
int iter1;
int iunit;

#if(fullsic)
double qaux[params_mp_kdfull2_][params_mp_kstate_];
#else
double qaux[1][1];       // dummy array
#endif


void prifld_cu(double *field,const char* comment){

#include<iostream>

bool tnopri=false;
int ind;

if(tnopri) return;
cout<< "field along x: "<<comment<<endl;
ind = 0;
for (int jz=params_mp_minz_;jz<params_mp_maxz_;jz++){
  for (int jy=params_mp_miny_;jy<params_mp_maxy_;jy++){
    for (int jx=params_mp_minx_;jx<=params_mp_maxx_;jx++){
      if(jz == params_mp_nzsh_-1 && jy == params_mp_nysh_-1) cout << ((double)jx-params_mp_nxsh_)*params_mp_dx_<<" "<<field[ind]<<endl;
      ind++;
    }
  }
}
}

extern "C" void statit_(double *psir,double *rho,double *aloc)
{
//     master routine for static iteration

//-----Check that later
/*#if(twostsic)
USE twostr
#endif
#if(fullsic)
USE localize_rad
#endif*/

// test Coulomb
calcrhor_(rho,psir);
accsum=0.0;
for (int ii=0;ii<params_mp_kdfull2_;ii++){
  accsum+=rho[ii];
}
cout<< "for charge="<<" "<< accsum*params_mp_dvol_<<endl<<endl;
falr_(rho,params_mp_chpcoul_,params_mp_nx2_,params_mp_ny2_,params_mp_nz2_,params_mp_kdfull2_);
prifld_cu(params_mp_chpcoul_,"coulomb pot");

/*   Number of pre-iterations for static solution with IFSICP=6.
     This parameter is to be set here "by hand" such that it can
     be communicated by 'all.inc'*/

cout<< " in STATIT"<<endl;

params_mp_itersicp6_=20;


    if(params_mp_istat_ == 1) {
        resume_(psir,params_mp_outnam_);
        calcpseudo_();                // initialize pseudo-potentials
	params_mp_int_pass_=0;
        static_mfield_(rho,aloc,psir,*qaux);
        pricm_(rho);
        infor_(psir,rho);
    }

//     initial computation of mean field, protocol prints, headers

calcpseudo_();                // initialize pseudo-potentials
params_mp_int_pass_=0;
static_mfield_(rho,aloc,psir,*qaux);
pricm_(rho);
if(params_mp_myn_ == 0) cout<< "rhoCM: "<<params_mp_rvectmp_[0]<<" "<<params_mp_rvectmp_[1]<<" "<<params_mp_rvectmp_[2]<<endl;
params_mp_int_pass_=0;
infor_(psir,rho);
#if(twostsic)
if(ifsicp >= 7) infor_sic_(psir);
#endif

if(params_mp_myn_ == 0) printone_(rho,aloc);

//     the static iteration starts here

ifsicpsav = params_mp_ifsicp_;       // save for later recycling
params_mp_sumvar_ = 100.0;           // initial value for terminator

for (iter1=1;iter1<params_mp_ismax_;iter1++)
{
  params_mp_int_pass_=iter1;
  if(ifsicpsav == 4) {     // switch safe pre-iterations for KLI
    if(iter1 < 40) {
      params_mp_ifsicp_ = 3;
    }
    else {
      params_mp_ifsicp_ = ifsicpsav;
    }
  }
  
  if(ifsicpsav == 6) {
#if(fullsic)
    if(iter1 <= params_mp_itersicp6_ && params_mp_sumvar_ > 10.0*params_mp_epsoro_) {
      params_mp_ifsicp_ = 3;   //   or 0 ?
    }
    else {
      params_mp_ifsicp_ = 6;
      if(iter1 <= 0) {
        cout<< " before LOCWF"<<endl;
        calc_locwf_(psir,*qaux);
        calc_fullsicr_(psir,*qaux);
        cout<< " after LOCWF"<<endl;
      }
    }
#else
    cout<< "IFSICP=6 requires LOCSIC or FULLSIC"<<endl;
#endif
  }
  
  if(params_mp_ifsicp_ != 6) {
    sstep_(psir,aloc);    // also for GSlat, DSIC
#if(fullsic)
  }
  else {
    sicstep_gen_(psir,*qaux,aloc);
#endif
  }

  static_mfield_(rho,aloc,psir,*qaux);
  pricm_(rho);

  if(params_mp_myn_ == 0) cout<< "rhoCM: "<<params_mp_rvectmp_[0]<<" "<<params_mp_rvectmp_[1]<<" "<<params_mp_rvectmp_[2]<<endl;
  
  if(iter1%params_mp_istinf_ == 0) {
    infor_(psir,rho);
#if(twostsic)
    if(ifsicp >= 7) infor_sic_(psir);
#endif
    
    if(params_mp_myn_ == 0) printtwo_();

    if(params_mp_sumvar2_ < params_mp_epsoro_) break;
  } //endif(iter1%params_mp_istinf_ == 0)
  
iunit=6;
flush_(&iunit);

}  // end iteration loop

if(params_mp_myn_ == 0) printthree_(rho,aloc);

//     compute and print localization

if(params_mp_iflocaliz_ == 1) {
#if(!parayes)
  localizer_(rho,psir);
#else
  cout<< " LOCALIZER (switch IFLOCALIZ) should not be invoked in parallele code"<<endl;   // cPW
  exit(-1);
#endif
}

#if(twostsic)

//     diagonalize Lagrange parameters
if(params_mp_ifsicp_ >= 7) diag_lagr_(psir);
#endif


//     final protocol on file 'pstat.<name>''

pri_pstat_(psir,rho);

//     save real wavefunctions for further applications

if(tp_prints && (params_mp_myn_ == 0 || params_mp_knode_ == 1)) {
  printfield_(491,aloc,"tp.aloc");
  printfield_(492,rho,"tp.density");
  printfield_(496,params_mp_chpcoul_,"tp.coulomb");
  printfield_(497,params_mp_potion_,"tp.potion");
}

if (params_mp_iplotorbitals_ != 0) print_orb_(psir);

calcrhor_(rho,psir);

if (params_mp_jplotdensitydiff_ != 0) print_densdiffc_(rho);

rsave_(psir);
params_mp_istat_=1;                       // prepare for reading in time step
if(params_mp_itmax_ == 0 && params_mp_isave_ > 0) {
  cout<< " CALL RSAVE  1. case"<<endl;
  params_mp_int_pass_=-1;
  infor_(psir,rho);
#if(parayes)
  mpi_finalize(params_mp_mpi_ierror)
#endif
  
  cout<<" terminate with static iteration"<<endl;
  exit(-1);
}

//     file for groundstate at the end of static iteration

//     optionally print KS potential along axes
//     use 'jforce' as switch

if(params_mp_jforce_ > 0 && params_mp_myn_ == 0) prifld_cu(aloc,"potential   \n");

if(params_mp_icooltyp_ == 1 && params_mp_itmax_ > 0) {
  cout<<" CALL RSAVE  2. case"<<endl;
  rsave_(psir);
}

//k check 1p-h transition matrix elements to unoccpied state:
//k attention: nstate must be bigger than number of electrons

#if(parano)
if(params_mp_iftransme_==1 && (params_mp_nstate_ > params_mp_nclust_)) transel_(psir);
#endif

#if(raregas)
if (params_mp_isurf_ != 0) print_surf(rho);
#endif


pricm_(rho);
if(params_mp_myn_ == 0) cout<< "rhoCM: "<<params_mp_rvectmp_[0]<<" "<<params_mp_rvectmp_[1]<<" "<<params_mp_rvectmp_[2]<<endl;

if(params_mp_idenspl_ != 0) mtv_fld_(rho,1);

if(params_mp_itmax_ <= 0) cout<< " terminate with static iteration "<<endl;

}
