#ifndef Physics_Analysis_Testpi_H
#define Physics_Analysis_Testpi_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"

class Testpi : public Algorithm {

public:
    Testpi(const std::string& name, ISvcLocator* pSvcLocator);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();  

private:

    // Declare r0, z0 cut for charged tracks
    double m_vr0cut;
    double m_vz0cut;
    double m_cthcut;
    
    //Declare energy, dphi, dthe cuts for fake gamma's
    double m_BenergyThreshold;
    double m_EenergyThreshold;
    double m_EenergyTdc;
    double m_gammaAngCut;
                  
    // 
    int m_test4C;
    int m_test5C;

    // 
    int m_checkDedx;
    int m_checkTof;

    NTuple::Tuple*  m_tuple1;     // rhopi 4C
    NTuple::Item<long>    m_run;
    NTuple::Item<long>    m_rec;
    NTuple::Item<long>    m_idxmc;  
    NTuple::Array<long>   m_pdgid; 
    NTuple::Array<long>   m_motheridx;        
    NTuple::Item<double>  m_dpxpip;
    NTuple::Item<double>  m_dpypip;
    NTuple::Item<double>  m_dpzpip;
    NTuple::Item<double>  m_dpxpim;
    NTuple::Item<double>  m_dpypim;
    NTuple::Item<double>  m_dpzpim;
    NTuple::Item<double>  m_dpxpi0;
    NTuple::Item<double>  m_dpypi0;
    NTuple::Item<double>  m_dpzpi0;
    NTuple::Item<double>  m_pi0_recoil;
    NTuple::Item<double>  m_dpxpip_no;
    NTuple::Item<double>  m_dpypip_no;
    NTuple::Item<double>  m_dpzpip_no;
    NTuple::Item<double>  m_dpxpim_no;
    NTuple::Item<double>  m_dpypim_no;
    NTuple::Item<double>  m_dpzpim_no;  
    NTuple::Item<double>  m_piptot_mt;
    NTuple::Item<double>  m_pimtot_mt;
    NTuple::Item<double>  m_pi0tot_mt;
    NTuple::Item<double>  m_mpip_mt;
    NTuple::Item<double>  m_mpim_mt;
    NTuple::Item<double>  m_mpi0_mt;
    NTuple::Item<double>  m_pi0rec;
    NTuple::Item<double>  m_piprec;
    NTuple::Item<double>  m_pimrec;
    NTuple::Item<double>  m_umiss;
    NTuple::Item<double>  m_umiss2;  
    NTuple::Item<long>    m_nch;
    NTuple::Item<long>    m_nneu;  
    NTuple::Item<double>  m_chi1;
    NTuple::Item<double>  m_mpi0;
    NTuple::Item<double>  m_emcpi0;
    NTuple::Item<double>  m_prho0;
    NTuple::Item<double>  m_prhop;
    NTuple::Item<double>  m_prhom;
    NTuple::Item<long>    m_good;
    NTuple::Item<long>    m_gam;
    NTuple::Item<long>    m_pip;
    NTuple::Item<long>    m_pim;
    NTuple::Item<double>  m_2gam;
    NTuple::Item<double>  m_outpi0;
    NTuple::Item<double>  m_outpip;
    NTuple::Item<double>  m_outpim;
    NTuple::Item<double>  m_enpip;
    NTuple::Item<double>  m_dcpip;
    NTuple::Item<double>  m_enpim;
    NTuple::Item<double>  m_dcpim;
    NTuple::Item<double>  m_pipf;
    NTuple::Item<double>  m_pimf;
    NTuple::Item<double>  m_pi0f;
    NTuple::Item<long>    m_nangecc;
    NTuple::Array<double> m_dthec;
    NTuple::Array<double> m_dphic;
    NTuple::Array<double> m_dangc;
    NTuple::Array<double> m_mspippim;
     
    NTuple::Item<double>  m_pmax;
    NTuple::Item<double>  m_ppx;
    NTuple::Item<double>  m_ppy;
    NTuple::Item<double>  m_ppz;
    NTuple::Item<double>  m_costhep;
    NTuple::Item<double>  m_ppxkal;
    NTuple::Item<double>  m_ppykal;
    NTuple::Item<double>  m_ppzkal;
    NTuple::Item<double>  m_mpx;
    NTuple::Item<double>  m_mpy;
    NTuple::Item<double>  m_mpz;
    NTuple::Item<double>  m_costhem;
    NTuple::Item<double>  m_phipin;
    NTuple::Item<double>  m_phip;
    NTuple::Item<double>  m_phimin;
    NTuple::Item<double>  m_phim;
    
    NTuple::Item<double>  m_mpxkal;
    NTuple::Item<double>  m_mpykal;
    NTuple::Item<double>  m_mpzkal;
    NTuple::Item<double>  m_vxpin;
    NTuple::Item<double>  m_vypin;
    NTuple::Item<double>  m_vzpin;
    NTuple::Item<double>  m_vrpin;
    NTuple::Item<double>  m_costhepin;  
    NTuple::Item<double>  m_vxmin;
    NTuple::Item<double>  m_vymin;
    NTuple::Item<double>  m_vzmin;      
    NTuple::Item<double>  m_vrmin;
    NTuple::Item<double>  m_costhemin;  
    
    NTuple::Item<double>  m_vxp;
    NTuple::Item<double>  m_vyp;
    NTuple::Item<double>  m_vzp;
    NTuple::Item<double>  m_vrp;
    NTuple::Item<double>  m_vxm;
    NTuple::Item<double>  m_vym;
    NTuple::Item<double>  m_vzm;      
    NTuple::Item<double>  m_vrm;
    
    NTuple::Item<double>  dangsg;
    NTuple::Item<double>  dthesg;
    NTuple::Item<double>  dphisg;
    NTuple::Item<double>  cosgth1;
    NTuple::Item<double>  cosgth2;
    NTuple::Item<double>  m_txt1;
    NTuple::Item<double>  m_txt2;
    NTuple::Item<double>  m_txt3;
    NTuple::Item<double>  m_txt4;  
    NTuple::Item<double>  m_chi5;
    NTuple::Item<double>  m_kpi0;
    NTuple::Item<double>  m_kpkm;
    NTuple::Item<double>  m_kpp0;
    NTuple::Item<double>  m_kmp0;
    NTuple::Item<double>  cosgve1;
    NTuple::Item<double>  m_pgam2pi1;
    NTuple::Item<double>  m_pgam2pi2;
    NTuple::Item<double>  cosva1;
    NTuple::Item<double>  cosva2;
    NTuple::Item<double>  m_laypi1;
    NTuple::Item<double>  m_hit1; 
    NTuple::Item<double>  m_laypi2;
    NTuple::Item<double>  m_hit2; 
    NTuple::Item<double>  m_anglepm;

    NTuple::Item<long>    m_ngch;
    NTuple::Array<double> m_ptrk;
    NTuple::Array<double> m_chie;
    NTuple::Array<double> m_chimu;
    NTuple::Array<double> m_chipi;
    NTuple::Array<double> m_chik;
    NTuple::Array<double> m_chip;
    NTuple::Array<double> m_probPH;
    NTuple::Array<double> m_normPH;
    NTuple::Array<double> m_PHe;
    NTuple::Array<double> m_PHmu;
    NTuple::Array<double> m_PHpion;
    NTuple::Array<double> m_PHkaon;
    NTuple::Array<double> m_PHproton;  
    NTuple::Array<double> m_PHreso_e;
    NTuple::Array<double> m_PHreso_mu;
    NTuple::Array<double> m_PHreso_pion;
    NTuple::Array<double> m_PHreso_kaon;
    NTuple::Array<double> m_PHreso_proton;
             
    NTuple::Array<double> m_ghit;
    NTuple::Array<double> m_thit;

    NTuple::Item<long>  m_tof1vo_trk1;
    NTuple::Item<long>  m_tof2vo_trk1;
    NTuple::Item<long>  m_tof1vo_trk2;
    NTuple::Item<long>  m_tof2vo_trk2;
                          
    NTuple::Array<long> m_tofvo1;
    NTuple::Array<long> m_tofvo2;

    NTuple::Array<double> m_ph_btof11;
    NTuple::Array<double> m_zhit_btof11;
    NTuple::Array<double> m_cntr_btof11;
    NTuple::Array<double> m_qual_btof11;
    NTuple::Array<double> m_tpi_btof11;
    NTuple::Array<double> m_tk_btof11;
    NTuple::Array<double> m_toff_pi_btof11;
    NTuple::Array<double> m_tsig_pi_btof11;
    
    NTuple::Array<double> m_ph_btof12;
    NTuple::Array<double> m_zhit_btof12;
    NTuple::Array<double> m_cntr_btof12;
    NTuple::Array<double> m_qual_btof12;
    NTuple::Array<double> m_tpi_btof12;
    NTuple::Array<double> m_tk_btof12;
    NTuple::Array<double> m_toff_pi_btof12;
    NTuple::Array<double> m_tsig_pi_btof12;
    
    NTuple::Array<double> m_ph_btof21;
    NTuple::Array<double> m_zhit_btof21;
    NTuple::Array<double> m_cntr_btof21;
    NTuple::Array<double> m_qual_btof21;
    NTuple::Array<double> m_tpi_btof21;
    NTuple::Array<double> m_tk_btof21;
    NTuple::Array<double> m_toff_pi_btof21;
    NTuple::Array<double> m_tsig_pi_btof21;

    NTuple::Array<double> m_ph_btof22;
    NTuple::Array<double> m_zhit_btof22;
    NTuple::Array<double> m_cntr_btof22;
    NTuple::Array<double> m_qual_btof22;
    NTuple::Array<double> m_tpi_btof22;
    NTuple::Array<double> m_tk_btof22;
    NTuple::Array<double> m_toff_pi_btof22;
    NTuple::Array<double> m_tsig_pi_btof22;
    
    NTuple::Array<double> m_ph_btof1;

    NTuple::Array<double> m_zhit_btof1;
    NTuple::Array<double> m_cntr_btof1;
    NTuple::Array<double> m_qual_btof1;
    NTuple::Array<double> m_tpi_btof1;
    NTuple::Array<double> m_tpi_btof1_mea;
    NTuple::Array<double> m_tpi_btof1_exp;
    NTuple::Array<double> m_tpi_btof2_mea;
    NTuple::Array<double> m_tpi_btof2_exp;
    NTuple::Array<double> m_tpi_btof_mea;
    NTuple::Array<double> m_tpi_btof_exp;    
    NTuple::Array<double> m_tk_btof1;
    NTuple::Array<double> m_toff_pi_btof1;
    NTuple::Array<double> m_tsig_pi_btof1;
    
    NTuple::Array<double> m_ph_btof2;
    NTuple::Array<double> m_zhit_btof2;
    NTuple::Array<double> m_cntr_btof2;
    NTuple::Array<double> m_qual_btof2;
    NTuple::Array<double> m_tpi_btof2;
    NTuple::Array<double> m_tk_btof2;
    NTuple::Array<double> m_toff_pi_btof2;
    NTuple::Array<double> m_tsig_pi_btof2;

    NTuple::Array<double> m_ph_btof;
    NTuple::Array<double> m_zhit_btof;
    NTuple::Array<double> m_cntr_btof;
    NTuple::Array<double> m_qual_btof;
    NTuple::Array<double> m_tpi_btof;
    NTuple::Array<double> m_tk_btof;
    NTuple::Array<double> m_toff_pi_btof;
    NTuple::Array<double> m_tsig_pi_btof;

    NTuple::Array<double> m_ph_etof;
    NTuple::Array<double> m_zhit_etof;
    NTuple::Array<double> m_cntr_etof;
    NTuple::Array<double> m_qual_etof;
    NTuple::Array<double> m_tpi_etof;
    NTuple::Array<double> m_tk_etof;
    NTuple::Array<double> m_toff_pi_etof;
    NTuple::Array<double> m_tsig_pi_etof;
    
    NTuple::Array<double> m_dedx_pid;
    NTuple::Array<double> m_tof_pid;
    NTuple::Array<double> m_tof1_pid;
    NTuple::Array<double> m_tof2_pid;
    NTuple::Array<double> m_chi_pid;
    NTuple::Array<double> m_prob_pid;
    NTuple::Array<double> m_ptrk_pid;
    NTuple::Array<double> m_cost_pid;  
    NTuple::Item<long>    m_pnp;
    NTuple::Item<long>    m_pnm;  
    NTuple::Item<long>    m_pnp_dedx;
    NTuple::Item<long>    m_pnm_dedx;    
    NTuple::Item<long>    m_pnp_tof;
    NTuple::Item<long>    m_pnm_tof;

    NTuple::Item<long>     m_pnp_misk;
    NTuple::Item<long>     m_pnm_misk;
    NTuple::Item<long>     m_pnp_misp;
    NTuple::Item<long>     m_pnm_misp;
    NTuple::Item<long>     m_pnp_dedxmisk;
    NTuple::Item<long>     m_pnm_dedxmisk;
    NTuple::Item<long>     m_pnp_dedxmisp;
    NTuple::Item<long>     m_pnm_dedxmisp;
    NTuple::Item<long>     m_pnp_tofmisk;
    NTuple::Item<long>     m_pnm_tofmisk;
    NTuple::Item<long>     m_pnp_tofmisp;
    NTuple::Item<long>     m_pnm_tofmisp;  
    NTuple::Item<double>   m_pi_prob_pip;
    NTuple::Item<double>   m_ki_prob_pip;
    NTuple::Item<double>   m_pr_prob_pip;
    NTuple::Item<double>   m_pi_prob_pim;
    NTuple::Item<double>   m_ki_prob_pim;
    NTuple::Item<double>   m_pr_prob_pim;
    NTuple::Item<double>   m_pi_chi_pip;
    NTuple::Item<double>   m_pi_chi_pim;
    NTuple::Item<double>   m_ki_chi_pip;
    NTuple::Item<double>   m_ki_chi_pim;
    NTuple::Item<double>   m_pr_chi_pip;
    NTuple::Item<double>   m_pr_chi_pim;
    NTuple::Item<double>   m_pi_probdedx_pip;
    NTuple::Item<double>   m_ki_probdedx_pip;
    NTuple::Item<double>   m_pr_probdedx_pip;
    NTuple::Item<double>   m_pi_probdedx_pim;
    NTuple::Item<double>   m_ki_probdedx_pim;
    NTuple::Item<double>   m_pr_probdedx_pim;
    NTuple::Item<double>   m_pi_chidedx_pip;
    NTuple::Item<double>   m_pi_chidedx_pim;
    NTuple::Item<double>   m_ki_chidedx_pip;
    NTuple::Item<double>   m_ki_chidedx_pim;
    NTuple::Item<double>   m_pr_chidedx_pip;
    NTuple::Item<double>   m_pr_chidedx_pim;
    NTuple::Item<double>   m_pi_probtof_pip;
    NTuple::Item<double>   m_ki_probtof_pip;
    NTuple::Item<double>   m_pr_probtof_pip;
    NTuple::Item<double>   m_pi_probtof_pim;
    NTuple::Item<double>   m_ki_probtof_pim;
    NTuple::Item<double>   m_pr_probtof_pim;
    NTuple::Item<double>   m_pi_chitof_pip;
    NTuple::Item<double>   m_pi_chitof_pim;
    NTuple::Item<double>   m_ki_chitof_pip;
    NTuple::Item<double>   m_ki_chitof_pim;
    NTuple::Item<double>   m_pr_chitof_pip;
    NTuple::Item<double>   m_pr_chitof_pim;  
    NTuple::Item<double>   m_chitof_pip;
    NTuple::Item<double>   m_chitof_pim;
    
    NTuple::Item<long>     m_nggneu;
    NTuple::Array<double>  m_numHits;    // Total number of hits
    NTuple::Array<double>  m_secondmoment;
    NTuple::Array<double>  m_x;       //  Shower coordinates and errors
    NTuple::Array<double>  m_y;
    NTuple::Array<double>  m_z;
    NTuple::Array<double>  m_cosemc;   // Shower Counter angles and errors
    NTuple::Array<double>  m_phiemc;
    NTuple::Array<double>  m_energy;  // Total energy observed in Emc
    NTuple::Array<double>  m_eSeed;
    NTuple::Array<double>  m_e3x3; 
    NTuple::Array<double>  m_e5x5; 
    NTuple::Array<double>  m_lat;
    NTuple::Array<double>  m_a20;
    NTuple::Array<double>  m_a42;  

    NTuple::Array<double>  m_dr0;
    NTuple::Array<double>  m_phi00;
    NTuple::Array<double>  m_kappa0;
    NTuple::Array<double>  m_dz0;
    NTuple::Array<double>  m_tanl0;
    NTuple::Array<double>  m_err00;
    NTuple::Array<double>  m_err01;
    NTuple::Array<double>  m_err02;
    NTuple::Array<double>  m_err03;
    NTuple::Array<double>  m_err04;
    NTuple::Array<double>  m_err05;
    NTuple::Array<double>  m_err06;
    NTuple::Array<double>  m_err07;
    NTuple::Array<double>  m_err08;    
    NTuple::Array<double>  m_err09;
    NTuple::Array<double>  m_err10;
    NTuple::Array<double>  m_err11;
    NTuple::Array<double>  m_err12;
    NTuple::Array<double>  m_err13;
    NTuple::Array<double>  m_err14;
    
    NTuple::Array<double>  m_drk0;
    NTuple::Array<double>  m_phik00;
    NTuple::Array<double>  m_kappak0;
    NTuple::Array<double>  m_dzk0;
    NTuple::Array<double>  m_tanlk0;
    NTuple::Array<double>  m_errk00;
    NTuple::Array<double>  m_errk01;
    NTuple::Array<double>  m_errk02;
    NTuple::Array<double>  m_errk03;
    NTuple::Array<double>  m_errk04;
    NTuple::Array<double>  m_errk05;
    NTuple::Array<double>  m_errk06;
    NTuple::Array<double>  m_errk07;
    NTuple::Array<double>  m_errk08;    
    NTuple::Array<double>  m_errk09;
    NTuple::Array<double>  m_errk10;
    NTuple::Array<double>  m_errk11;
    NTuple::Array<double>  m_errk12;
    NTuple::Array<double>  m_errk13;
    NTuple::Array<double>  m_errk14;     
    NTuple::Item<double>   m_kal_stat_pip;
    NTuple::Item<double>   m_kal_stat_pim;          
    NTuple::Item<double>   m_kal_stat0_pip;
    NTuple::Item<double>   m_kal_stat0_pim;
    NTuple::Item<double>   m_kal_stat1_pip; 
    NTuple::Item<double>   m_kal_stat1_pim; 
    NTuple::Item<double>   m_lp_tof;
    NTuple::Item<double>   m_lm_tof;
    NTuple::Item<double>   m_lp_tof_t0;
    NTuple::Item<double>   m_lm_tof_t0;
    NTuple::Item<int>      m_lp_tof_flag; //1 barral; 2 endcap
    NTuple::Item<int>      m_lm_tof_flag; //1 barral; 2 endcap
    NTuple::Item<double>   m_lp_tof_ph;
    NTuple::Item<double>   m_lm_tof_ph;  

};

#endif 
