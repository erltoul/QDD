// This file just allow using variable set in params.F90 module with C++/Cuda code
//
// general settings for dimensions, part.numbers etc
//
#include <cufft.h>
#define kmom 35
#define maxnang 100 // max. angular bins for PAD
#define maxmps 100  // max. nr. analyzing points fpr PES
// number of nodes (=1 for serial version)
extern int params_mp_knode_;
// max. nr. electron states per node
//fix// int,PARAMETER kstate=20
extern int params_mp_kstate_;
// max. total nr. electron  states
extern int params_mp_ksttot_;


// maximum number of ions
//fix// int,PARAMETER ng=8;
extern int params_mp_ng_;


//extern int params_mp_maxnang_;    // max. angular bins for PAD
//extern int params_mp_maxmps_;    // max. nr. analyzing points fpr PES


// physical units and constants (here Rydberg units)
extern double params_mp_e2_;
extern double params_mp_hbar_;
extern double params_mp_ame_;
extern double params_mp_h2m_;  //  h2m =  hbar**2/2m

// frequently used mathematical constants
extern double params_mp_zero_;
extern double params_mp_half_;
extern double params_mp_one_;
extern double params_mp_PI_;
extern double params_mp_fourpi_;
extern double params_mp_small_;
extern double params_mp_sq2_;
extern double params_mp_sqrtpi_;
extern cufftDoubleComplex params_mp_eye_;


//     maximum sizes of the box in x,y, and z
extern int params_mp_kxbox_,params_mp_kybox_,params_mp_kzbox_;
// deduced grid paramaters
extern int params_mp_kxmax_,params_mp_kymax_,params_mp_kzmax_;
extern int params_mp_nx2_,params_mp_ny2_,params_mp_nz2_;
//int,PARAMETER nxyzf=nx2*ny2*nz2_;                 //?
extern int params_mp_nxyz_,params_mp_nyf_,params_mp_nxyf_;
extern int params_mp_kdfull2_;
//int,PARAMETER kfbox=kdfull2_;                    // ??
extern int params_mp_nx_,params_mp_ny_,params_mp_nz_;
//int,PARAMETER nx1=nx+1,ny1=ny+1,nz1=nz+1,nzi=nz1,nzr=nz+nz_;  //?
//int,PARAMETER nxy1=nx1*ny1_;    //?
//int,PARAMETER ksmax=nx+1,kdfull=(nx+1)*(ny+1)*(nz+1)_;    //?
extern int params_mp_knodem_;     // =knode-1
#if(gridfft)
// bounds of loops
extern int params_mp_minx_,params_mp_maxx_;
extern int params_mp_miny_,params_mp_maxy_;
extern int params_mp_minz_,params_mp_maxz_;
// bounds of esc. el. ???
extern int params_mp_nbnx_,params_mp_nbxx_;
extern int params_mp_nbny_,params_mp_nbxy_;
extern int params_mp_nbnz_,params_mp_nbxz_;
// mid-point for x,y,z-values
extern int params_mp_nxsh_,params_mp_nysh_,params_mp_nzsh_;
#endif
#if(findiff|numerov)
// bounds of loops
extern int params_mp_minx_,params_mp_maxx_;
extern int params_mp_miny_,params_mp_maxy_;
extern int params_mp_minz_,params_mp_maxz_;
// bounds of esc. el.
extern int params_mp_nbnx_,params_mp_nbxx_;
extern int params_mp_nbny_,params_mp_nbxy_;
extern int params_mp_nbnz_,params_mp_nbxz_;
// offset for x,y,z-values   ???
extern int params_mp_nxsh_,params_mp_nysh,params_mp_nzsh_;
#endif


extern double params_mp_dx_,params_mp_dy_,params_mp_dz_,params_mp_dvol_;                  //  mesh spacing, volume
extern double *params_mp_xval_,*params_mp_yval_,*params_mp_zval_;        //  grid coordinates
extern double *params_mp_xt2_,*params_mp_yt2_,*params_mp_zt2_;           //  coordinates**2


extern double *params_mp_enonlo_;               //  non-local s.p.energy
extern double *params_mp_amoy_;                 //  single particle energies
extern double *params_mp_spvariance_;           //  s.p. energy variances
extern double *params_mp_spvariancep_;          //  s.p. en.varian. projected
extern double *params_mp_spvariancebi_;         //  s.p. en.varian. boost-inv.
extern double *params_mp_spenergybi_;           //  s.p. energy boost-inv.
extern double *params_mp_spnorm_;               //  norm of s.p. wf
extern double *params_mp_occup_;                //  occupation weight
extern int params_mp_nstate_,params_mp_nspdw_;                     //  Nr. of states
extern int params_mp_nstate_all_;                           //  total Nr. of states (parayes)
extern int **params_mp_nq_;                  //  nodes of init. states
extern int **params_mp_ipar_;                //  xyz-parities
extern int *params_mp_ispin_;                 //  spin of states
extern int *params_mp_nrel2abs_,*params_mp_nabs2rel_;  //  pointer to wfs
extern int **params_mp_nrel2abs_other_;      //  pointer to wfs
extern int *params_mp_nhome_;                 //  home node of wf


#if(parayes)
extern int *params_mp_nstate_node_;   // number of active states in a node
extern int *params_mp_nstart_node_;   // offset address for counting in a node
extern int **params_mp_ispin_node_;
#endif


//     basic parameters

extern bool params_mp_directenergy_;            // to compute energy directly
extern int params_mp_numspin_;                       // number of spin components
extern double params_mp_omeg_,params_mp_eferm_;                     //  initial Gaussians
extern double params_mp_dInMargin_;
extern int params_mp_iforce_;
extern int params_mp_ipseudo_;          // swich on/off subgrids and pseudo-densities
extern double params_mp_xfac_,params_mp_yfac_,params_mp_zfac_;                       //  initial. auxiliary
extern double params_mp_epswf_,params_mp_e0dmp_,params_mp_epsorc_;        //  convergence
extern double params_mp_b2occ_,params_mp_gamocc_,params_mp_deocc_,params_mp_osfac_; //  h.o.initalization
extern double params_mp_temp_,params_mp_occmix_,params_mp_epsoro_;        //  electron temperature
extern int params_mp_isurf_;
extern double params_mp_endcon_,params_mp_radjel_,params_mp_surjel_;     //  jellium params
extern int params_mp_itback_;                            //  jellium params
extern double params_mp_betatj_,params_mp_gammtj_,params_mp_bet4tj_;                 //      "      "
extern double params_mp_bbeta_,params_mp_gamma_,params_mp_beta4_;        //      "      "
extern double params_mp_falph_,params_mp_fbeta_,params_mp_fhexe_;                    //  jellium auxiliary
extern double params_mp_dpolx_,params_mp_dpoly_,params_mp_dpolz_;       //  static dipol potential
extern bool params_mp_tdipolxyz_;                     // switch to dipole fields
extern int params_mp_iexcit_,params_mp_irotat_,params_mp_ispidi_;          //  dynam. initialisation
extern double params_mp_phirot_;
extern int params_mp_nclust_,params_mp_nion_,params_mp_nion2_;
extern double params_mp_charge_;             //  Nr. el. & ions
extern double params_mp_scaleClust_,params_mp_scaleClustx_,params_mp_scaleClusty_;
extern double params_mp_scaleClustz_;
extern int params_mp_iemomsRel_;   // 1 = relative to c.m.,  0 = relative to box
extern double params_mp_shiftClustx_,params_mp_shiftClusty_,params_mp_shiftClustz_;
extern double params_mp_rotClustx_,params_mp_rotClusty_,params_mp_rotClustz_;
extern int params_mp_imob_;   // fix (=0) or unfix (=1) ions
extern int params_mp_ifixcmion_;
extern double params_mp_shiftWFx_,params_mp_shiftWFy_,params_mp_shiftWFz_;
extern int params_mp_ispinsep_;
extern double params_mp_comx_,params_mp_comy_,params_mp_comz_,params_mp_dion_[3],params_mp_qtion_[3][3];
extern double params_mp_apnum_,params_mp_rmsion_,params_mp_dmdistion_,params_mp_rmsel_,params_mp_rhopss_,params_mp_rhomix_;
extern double params_mp_codx_,params_mp_cody_,params_mp_codz_,params_mp_del_[3],params_mp_qtel_[3][3];
extern double params_mp_time_absinit_;
extern double params_mp_phangle_,params_mp_phphase_;         // angle and phase of ph rotation

extern int params_mp_nhstate_,params_mp_npstate_;
extern int params_mp_idebug_,ifreezekspot_;
extern int params_mp_itof_,params_mp_jescmask_,params_mp_ishiftcmtoorigin_,params_mp_jescmaskorb_;
extern int params_mp_ishutdown_;
// int icheckformessages=1,jcheckformessages=50_;
extern int params_mp_jplotdensitydiff_,params_mp_jplotdensity2d_,params_mp_jplotdensitydiff2d_;
extern int params_mp_nmptheta_,params_mp_nmpphi_,params_mp_nmps_,params_mp_jmp_,*params_mp_imps_[maxmps];
extern int params_mp_jovlp_,params_mp_jnorms_;
extern int params_mp_iscatterelectron_,params_mp_jcharges_;
extern int params_mp_iaddcluster_,params_mp_iswforce_,params_mp_iplotorbitals_,params_mp_ievaluate_;
extern double params_mp_ekin0pp_,params_mp_vxn0_,params_mp_vyn0_,params_mp_vzn0_;
extern double params_mp_eproj_,params_mp_vpx_,params_mp_vpy_,params_mp_vpz_,params_mp_taccel_;
extern double params_mp_trequest_,params_mp_timefrac_;
extern double params_mp_rheatclust_;
extern double params_mp_igeneratesurffile_;
//double,ALLOCATABLE rnormsinit(kstate)_;
extern double params_mp_ehom0_,params_mp_ehomx_,params_mp_ehomy_,params_mp_ehomz_;
extern int params_mp_ihome_;
extern double params_mp_scatterelectronenergy_,params_mp_scatterelectronw_;
extern double params_mp_scatterelectronvxn_,scatterelectronvyn_,scatterelectronvzn_;
extern double params_mp_scatterelectronx_,params_mp_scatterelectrony_,params_mp_scatterelectronz_;
extern double params_mp_drcharges_;


extern const char params_mp_outnam_;
extern int params_mp_iflocaliz_;                           // evaluate localization
extern int params_mp_myn_;                                 // nr. of actual node
extern int params_mp_ifls_,params_mp_ismax_,params_mp_itmax_,params_mp_istinf_,params_mp_ipasinf_;
extern int params_mp_idyniter_;        // number iterations to start dynamic E0DMP 
extern int params_mp_iffastpropag_,params_mp_ifexpevol_;
extern int params_mp_irest_,params_mp_istat_,params_mp_isave_,params_mp_idenspl_;
extern int params_mp_i3dz_,params_mp_i3dx_,params_mp_i3dstate_,params_mp_istream_,params_mp_modrho_;
extern int params_mp_jpos_,params_mp_jvel_,params_mp_jener_,params_mp_jesc_,params_mp_jforce_,params_mp_jposcm_,params_mp_jgeomion_;
extern int params_mp_jinfo_,params_mp_jdip_,params_mp_jquad_,params_mp_jang_,params_mp_jspdp_,params_mp_jenergy_;
extern int params_mp_jgeomel_,params_mp_jangabso_,params_mp_jelf_,params_mp_jstinf_,params_mp_jstboostinv_;
extern int params_mp_jstateoverlap_;
extern int params_mp_nabsorb_,params_mp_ifsicp_,params_mp_ifredmas_,params_mp_ionmdtyp_,params_mp_icooltyp_;
extern int params_mp_init_lcao_,params_mp_ipsptyp_,params_mp_ivdw_,params_mp_idenfunc_;
extern int params_mp_izforcecorr_,params_mp_icooltimes_,params_mp_ntref_;
extern int params_mp_jheatmod_;         // modulus for re-heating the system
extern int params_mp_ifrhoint_time_,params_mp_iangmo_,params_mp_ifspemoms_,params_mp_iftransme_;
extern int params_mp_iterat_,params_mp_itersicp6_;
extern bool params_mp_tstinf_;
extern double params_mp_rheattemp_;        // re-heat temperature
extern double params_mp_tempion_,params_mp_dt1_;
extern double params_mp_centfx_,params_mp_centfy_,params_mp_centfz_;
extern double params_mp_shiftinix_,params_mp_shiftiniy_,params_mp_shiftiniz_;

#if(parayes)
extern int params_mp_ifhamdiag_;
#endif
#if(parano)
extern int params_mp_ifhamdiag_;
#endif



//     spatial fields as densities and potentials
extern double *params_mp_rhojel_;                //  jellium density
extern double *params_mp_potion_;                //  pseudopotentials
extern double *params_mp_potFixedIon_;           //  potential from frozen ions
extern double *params_mp_chpcoul_;                 //  Coulomb potential



//     fields and variables for analysis of electron emission

extern double *params_mp_rhoabso_;           //  storage for absorbed density
extern double *params_mp_spherMask_;         //  mask for spherical absorbing bounds
extern double *params_mp_spherloss_;         //  loss factor for spher. abs. bounds
extern double params_mp_bcrad_,params_mp_powabso_;             // width & power of abs. bounds
extern int params_mp_ispherAbso_;               //  swicth to spherical abs. bounds
extern int params_mp_iangabso_,params_mp_nangtheta_,params_mp_nangphi_;
extern int params_mp_ipes_,params_mp_indicesmp_[maxnang*maxnang];
extern double params_mp_angthetah_,params_mp_angthetal_;
extern double params_mp_angphil_,params_mp_angphih_,params_mp_delomega_,params_mp_xango_,params_mp_yango_,params_mp_zango_;
extern bool *params_mp_tgridabso_;       //  array tagging absorbing points
extern double **params_mp_rhoabsoorb_;




//     common fields for the spatial moments

//      common /moment/ ql(kmom)_,

//extern int params_mp_kmom_;
extern int params_mp_nrmom_;
extern double params_mp_qe_[kmom],params_mp_se_[5],params_mp_ajx_,params_mp_ajy_,params_mp_ajz_;
//COMMON /moment/ qe,se,ajx,ajy,ajz,nrmom

// storage for the case of 1ph rotation (see 'phangle')
extern cufftDoubleComplex *params_mp_oldhole_,*params_mp_newhole_;

// storage for base wavefunctions in case of dynamics with exact exchange 
// or in case of computing overlaps with initial state
extern cufftDoubleComplex **params_mp_psisavex_;

//     the energy transmitted from calc-lda to info etc

extern double params_mp_enrear_,params_mp_ecback_,params_mp_ecrho_,params_mp_ecorr_;
extern double params_mp_dt12_,params_mp_sgaus_,params_mp_ekion_,params_mp_energy_;
extern double params_mp_energ2_,params_mp_enerpw_,params_mp_encoulsp_,params_mp_entrop_,params_mp_epot_,params_mp_espnb_,params_mp_esh1_;
extern double params_mp_etot_,params_mp_ekionold_,params_mp_qold2_,params_mp_qold3_,params_mp_qold4_;
extern double params_mp_ekmat_,params_mp_engg_,params_mp_enii_,params_mp_enig_,params_mp_ecrhoimage_;
extern double *params_mp_ekinsp_,*params_mp_evarsp_,*params_mp_evarsp2_,*params_mp_epotsp_;
extern int params_mp_jekion_,params_mp_iquery4_;
#if(fullsic||twostsic)  
extern double **params_mp_hmatrix_;
extern double params_mp_symcon_;
#endif


// dynamic variables of ionic motion

extern int params_mp_nxsg_,params_mp_nysg_,params_mp_nzsg_;      // size of subgrids
extern int params_mp_inewforce_;
extern int params_mp_mzforce_,params_mp_myforce_,params_mp_mxforce_;       // symmetrized forces
extern int params_mp_nrare_,params_mp_nfix_;                    //  Nr. of raregas atoms
extern int *params_mp_nfixed_;                    //  Nr. of fixed ions
extern int params_mp_idielec_;                     //  switch to dielectricum
extern bool *params_mp_tblock_;
extern double *params_mp_cx_,*params_mp_cy_,*params_mp_cz_,*params_mp_cpx_,*params_mp_cpy_,*params_mp_cpz_; 
extern double *params_mp_dcx_,*params_mp_dcy_,*params_mp_dcz_,*params_mp_dcpx_,*params_mp_dcpy_,*params_mp_dcpz_; 
extern double *params_mp_fx_,*params_mp_fy_,*params_mp_fz_,*params_mp_flx_,*params_mp_fly_,*params_mp_flz_;
extern double *params_mp_fprojx_,*params_mp_fprojy_,*params_mp_fprojz_;

//                                      book keeping for LCGO initialization
extern double *params_mp_radini_,*params_mp_ipol_;
extern int *params_mp_initord_;
extern int *params_mp_nmaxst_;



//     parameters for simulated annealinig
extern double *params_mp_eloc_,*params_mp_enoloc_,*params_mp_eion_,*params_mp_eiinew_;
extern double params_mp_cptemp_,params_mp_delpos_,params_mp_ERR_,params_mp_binerg_;
extern double params_mp_errtot_,params_mp_erfac1_,params_mp_trfac1_,params_mp_prfac1_;
extern double params_mp_trfac2_,params_mp_prfac2_,params_mp_errsim_,params_mp_eiontot_,params_mp_enoloctot_,params_mp_facann_;
extern double params_mp_eloctot_,params_mp_errks0_,params_mp_errks_,params_mp_sumvar_,params_mp_sumvar2_;
extern int params_mp_ionsin_,params_mp_nrun_,params_mp_nloop1_,params_mp_nloop2_,params_mp_loop1_;
extern int params_mp_ifall_,params_mp_nyes_,params_mp_ncon_, params_mp_ncsim_,params_mp_iknow_;



//     parameters for external excitations by laser or projectile
extern int params_mp_itft_;
extern double params_mp_tnode_,params_mp_deltat_,params_mp_tpeak_,params_mp_omega_,params_mp_e0_,params_mp_time_,params_mp_tfs_;  
extern double params_mp_e1x_,params_mp_e1y_,params_mp_e1z_,params_mp_e2x_,params_mp_e2y_,params_mp_e2z_,params_mp_phi_;
extern double *params_mp_fl_,params_mp_power_;
extern double params_mp_elaser_;
extern int params_mp_ijel_;
extern double params_mp_acc1old_,params_mp_acc2old_,params_mp_foft1old_,params_mp_foft2old_,params_mp_timeold_;
extern int params_mp_ilas_;
extern double params_mp_fpulseinteg1_;        // integrated pulse for gauge trasnf
extern double params_mp_fpulseinteg2_;        // integrated pulse for gauge trasnf.
extern double params_mp_projcharge_;                   // projectile charge
extern double params_mp_projvelx_,params_mp_projvely_,params_mp_projvelz_;   // projectile velocity
extern double params_mp_projinix_,params_mp_projiniy_,params_mp_projiniz_;   // initial projectile position
                   // impact parameter = min(projinix,projiniy,projiniz)



// workspace for communication
extern double params_mp_rvectmp2_[3],params_mp_rvectmp_[3];
extern int params_mp_iindtmp_[3];

extern int params_mp_num_gpus_;  //total number of gpus on the node
extern int params_mp_mygpu_; //number of the actual gpu used by the node

cudaError_t error;
cudaStream_t *stream1,*stream2,*stream3;

void Check_CUDA_Error(cudaError_t error)
{
	if(error!=cudaSuccess) 
	{
		printf("Cuda error: %s\n",cudaGetErrorString(error));
		exit(-1);
	}
}
