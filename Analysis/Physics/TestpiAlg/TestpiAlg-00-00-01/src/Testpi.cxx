#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"
#include "McTruth/McParticle.h"
#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "TestpiAlg/Testpi.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"
#include <vector>
using namespace Event;
const double mpi = 0.13957;
const double mk  = 0.493677;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
const double velc = 299.792458;   // tof path unit in mm
const double par_mc[5]={550*0.068461,9.07863,0.00159196,2.18615,1.56239};
//Data
const double par_data[5]= {550.0*0.0154625,34.9635,4.04586e-14,2.39614,7.31565};
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
const HepLorentzVector ecms(0.034,0,0,3.097);
const HepLorentzVector ecms2(0.034,0,0,3.097);
const Hep3Vector u_cms = -ecms2.boostVector();
int Ncounter0;
int Ncut0,Ncut1,Ncut2,Ncut3,Ncut4,Ncut5,Ncut6,Ncut7,Ncut8,Ncut9,Ncut10,Ncut11;

/////////////////////////////////////////////////////////////////////////////

Testpi::Testpi(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator) {
    //Declare the properties  
    declareProperty("Vr0cut", m_vr0cut=1.0);
    declareProperty("Vz0cut", m_vz0cut=10.0);
    declareProperty("Vctcut", m_cthcut=0.93);
    declareProperty("BEnergyThreshold", m_BenergyThreshold=0.025);
    declareProperty("EEnergyThreshold", m_EenergyThreshold=0.05);
    declareProperty("EnergyTdc",   m_EenergyTdc=14);
    declareProperty("GammaAngCut", m_gammaAngCut=10.0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode Testpi::initialize(){
  MsgStream Log(msgSvc(), name());
    Log << MSG::INFO << "in initialize()" << endmsg;
    StatusCode status;
    
    NTuplePtr nt1(ntupleSvc(), "FILE1/fit4c");
    if (nt1) m_tuple1 = nt1;
    else {
        m_tuple1 = ntupleSvc()->book ("FILE1/fit4c", CLID_ColumnWiseTuple, "ks N-Tuple example");
        if (m_tuple1) {
            status = m_tuple1->addItem ("run",             m_run);
            status = m_tuple1->addItem ("rec",             m_rec);     
            status = m_tuple1->addItem("indexmc",          m_idxmc, 0, 100);
            status = m_tuple1->addIndexedItem("pdgid",     m_idxmc, m_pdgid);
            status = m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);
            status = m_tuple1->addItem("mdpxpip",          m_dpxpip);
            status = m_tuple1->addItem("mdpypip",          m_dpypip);
            status = m_tuple1->addItem("mdpzpip",          m_dpzpip);
            status = m_tuple1->addItem("mdpxpim",          m_dpxpim);
            status = m_tuple1->addItem("mdpypim",          m_dpypim);
            status = m_tuple1->addItem("mdpzpim",          m_dpzpim);
            status = m_tuple1->addItem("mdpxpi0",          m_dpxpi0);
            status = m_tuple1->addItem("mdpypi0",          m_dpypi0);
            status = m_tuple1->addItem("mdpzpi0",          m_dpzpi0);
            status = m_tuple1->addItem("pi0_recoil",       m_pi0_recoil);
            status = m_tuple1->addItem("mdpxpip_no",       m_dpxpip_no);
            status = m_tuple1->addItem("mdpypip_no",       m_dpypip_no);
            status = m_tuple1->addItem("mdpzpip_no",       m_dpzpip_no);
            status = m_tuple1->addItem("mdpxpim_no",       m_dpxpim_no);
            status = m_tuple1->addItem("mdpypim_no",       m_dpypim_no);
            status = m_tuple1->addItem("mdpzpim_no",       m_dpzpim_no);
            
            status = m_tuple1->addItem("mpiptot_mt", m_piptot_mt);
            status = m_tuple1->addItem("mpimtot_mt", m_pimtot_mt);
            status = m_tuple1->addItem("mpi0tot_mt", m_pi0tot_mt);
            status = m_tuple1->addItem("mmpip_mt",   m_mpip_mt);  
            status = m_tuple1->addItem("mmpim_mt",   m_mpim_mt);  
            status = m_tuple1->addItem("mmpi0_mt",   m_mpi0_mt);  
            status = m_tuple1->addItem("mpi0rec",    m_pi0rec);   
            status = m_tuple1->addItem("mpiprec",    m_piprec); 
            status = m_tuple1->addItem("mpimrec",    m_pimrec); 
            status = m_tuple1->addItem("mumiss",     m_umiss);
            status = m_tuple1->addItem("mumiss2",    m_umiss2);
            status = m_tuple1->addItem("nch",        m_nch);
            status = m_tuple1->addItem("nneu",       m_nneu);
            status = m_tuple1->addItem("chi1",       m_chi1);
            status = m_tuple1->addItem("mpi0",       m_mpi0);
            status = m_tuple1->addItem("memcpi0",    m_emcpi0);
            status = m_tuple1->addItem("mprho0",     m_prho0);
            status = m_tuple1->addItem("mprhop",     m_prhop);
            status = m_tuple1->addItem("mprhom",     m_prhom);
            status = m_tuple1->addItem("mgood",      m_good);
            status = m_tuple1->addItem("mgam",       m_gam);
            status = m_tuple1->addItem("mpip",       m_pip);
            status = m_tuple1->addItem("mpim",       m_pim);
            status = m_tuple1->addItem("m2gam",      m_2gam);
            status = m_tuple1->addItem("ngch",       m_ngch, 0, 10);
            status = m_tuple1->addItem("nggneu",     m_nggneu,0, 10);
            status = m_tuple1->addItem("moutpi0",    m_outpi0);
            status = m_tuple1->addItem("moutpip",    m_outpip);
            status = m_tuple1->addItem("moutpim",    m_outpim);
            status = m_tuple1->addItem("menpip",     m_enpip);
            status = m_tuple1->addItem("mdcpip",     m_dcpip );
            status = m_tuple1->addItem("menpim",     m_enpim );
            status = m_tuple1->addItem("mdcpim",     m_dcpim );   
            status = m_tuple1->addItem("mpipf",      m_pipf);        
            status = m_tuple1->addItem("mpimf",      m_pimf);        
            status = m_tuple1->addItem("mpi0f",      m_pi0f);
            
            status = m_tuple1->addItem("mpmax",    m_pmax);
            status = m_tuple1->addItem("mppx",     m_ppx);
            status = m_tuple1->addItem("mppy",     m_ppy);
            status = m_tuple1->addItem("mppz",     m_ppz);
            status = m_tuple1->addItem("mcosthep", m_costhep);
            status = m_tuple1->addItem("mppxkal",  m_ppxkal);
            status = m_tuple1->addItem("mppykal",  m_ppykal);
            status = m_tuple1->addItem("mppzkal",  m_ppzkal);
            status = m_tuple1->addItem("mmpx",     m_mpx);
            status = m_tuple1->addItem("mmpy",     m_mpy);
            status = m_tuple1->addItem("mmpz",     m_mpz);
            status = m_tuple1->addItem("mcosthem", m_costhem);
            status = m_tuple1->addItem("m_phipin", m_phipin);
            status = m_tuple1->addItem("m_phip",   m_phip);
            status = m_tuple1->addItem("m_phimin", m_phimin);
            status = m_tuple1->addItem("m_phim",   m_phim);        
            
            status = m_tuple1->addItem("mmpxkal", m_mpxkal);
            status = m_tuple1->addItem("mmpykal", m_mpykal);
            status = m_tuple1->addItem("mmpzkal", m_mpzkal);
            
            status = m_tuple1->addItem("mvxpin",          m_vxpin);
            status = m_tuple1->addItem("mvypin",          m_vypin);
            status = m_tuple1->addItem("mvzpin",          m_vzpin);
            status = m_tuple1->addItem("mvrpin",          m_vrpin);
            status = m_tuple1->addItem("mcosthepin",      m_costhepin);
            status = m_tuple1->addItem("mvxmin",          m_vxmin);
            status = m_tuple1->addItem("mvymin",          m_vymin);
            status = m_tuple1->addItem("mvzmin",          m_vzmin);
            status = m_tuple1->addItem("mvrmin",          m_vrmin);
            status = m_tuple1->addItem("mcosthemin",      m_costhemin);
            status = m_tuple1->addItem("mvxp",            m_vxp);
            status = m_tuple1->addItem("mvyp",            m_vyp);
            status = m_tuple1->addItem("mvzp",            m_vzp);
            status = m_tuple1->addItem("mvrp",            m_vrp);
            status = m_tuple1->addItem("mvxm",            m_vxm);
            status = m_tuple1->addItem("mvym",            m_vym);
            status = m_tuple1->addItem("mvzm",            m_vzm);
            status = m_tuple1->addItem("mvrm",            m_vrm);
            status = m_tuple1->addItem("nangecc",         m_nangecc,0,10);
            status = m_tuple1->addIndexedItem("mdthec",   m_nangecc, m_dthec);
            status = m_tuple1->addIndexedItem("mdphic",   m_nangecc, m_dphic);
            status = m_tuple1->addIndexedItem("mdangc",   m_nangecc, m_dangc);
            status = m_tuple1->addIndexedItem("mspippim", m_nangecc, m_mspippim);
            
            status = m_tuple1->addItem("angsg",     dangsg);
            status = m_tuple1->addItem("thesg",     dthesg);
            status = m_tuple1->addItem("phisg",     dphisg);
            status = m_tuple1->addItem("cosgth1",   cosgth1);
            status = m_tuple1->addItem("cosgth2",   cosgth2);
            status = m_tuple1->addItem("mtxt1",     m_txt1);
            status = m_tuple1->addItem("mtxt2",     m_txt2);
            status = m_tuple1->addItem("mtxt3",     m_txt3);
            status = m_tuple1->addItem("mtxt4",     m_txt4);
            status = m_tuple1->addItem("mchi5",     m_chi5);
            status = m_tuple1->addItem("mkpi0",     m_kpi0);
            status = m_tuple1->addItem("mkpkm",     m_kpkm);
            status = m_tuple1->addItem("mkpp0",     m_kpp0);
            status = m_tuple1->addItem("mkmp0",     m_kmp0);            
            status = m_tuple1->addItem("cosgve1",   cosgve1);            
            status = m_tuple1->addItem("mpgam2pi1", m_pgam2pi1);
            status = m_tuple1->addItem("mpgam2pi2", m_pgam2pi2);
            status = m_tuple1->addItem("cosva1",    cosva1);
            status = m_tuple1->addItem("cosva2",    cosva2);
            status = m_tuple1->addItem("laypi1",    m_laypi1);
            status = m_tuple1->addItem("hit1",      m_hit1);
            status = m_tuple1->addItem("laypi2",    m_laypi2);
            status = m_tuple1->addItem("hit2",      m_hit2);
            status = m_tuple1->addItem("anglepm",   m_anglepm);
            
            status = m_tuple1->addIndexedItem("mptrk",          m_ngch, m_ptrk);
            status = m_tuple1->addIndexedItem("chie",           m_ngch, m_chie);
            status = m_tuple1->addIndexedItem("chimu",          m_ngch,m_chimu);
            status = m_tuple1->addIndexedItem("chipi",          m_ngch,m_chipi);
            status = m_tuple1->addIndexedItem("chik",           m_ngch,m_chik);
            status = m_tuple1->addIndexedItem("chip",           m_ngch,m_chip);
            status = m_tuple1->addIndexedItem("probPH",         m_ngch,m_probPH);
            status = m_tuple1->addIndexedItem("normPH",         m_ngch,m_normPH);
            status = m_tuple1->addIndexedItem("mPHe",           m_ngch,m_PHe);
            status = m_tuple1->addIndexedItem("mPHmu",          m_ngch,m_PHmu);
            status = m_tuple1->addIndexedItem("mPHpion",        m_ngch,m_PHpion);
            status = m_tuple1->addIndexedItem("mPHkaon",        m_ngch,m_PHkaon);
            status = m_tuple1->addIndexedItem("mPHproton",      m_ngch,m_PHproton);
            status = m_tuple1->addIndexedItem("mPHreso_e",      m_ngch,m_PHreso_e);
            status = m_tuple1->addIndexedItem("mPHreso_mu",     m_ngch,m_PHreso_mu);
            status = m_tuple1->addIndexedItem("mPHreso_pion",   m_ngch,m_PHreso_pion);
            status = m_tuple1->addIndexedItem("mPHreso_kaon",   m_ngch,m_PHreso_kaon);
            status = m_tuple1->addIndexedItem("mPHreso_proton", m_ngch,m_PHreso_proton);
            status = m_tuple1->addIndexedItem("ghit",           m_ngch,m_ghit);
            status = m_tuple1->addIndexedItem("thit",           m_ngch,m_thit);
            
            status = m_tuple1->addIndexedItem("ph_btof11",      m_ngch,     m_ph_btof11);
            status = m_tuple1->addIndexedItem("zhit_btof11",    m_ngch,     m_zhit_btof11);
            status = m_tuple1->addIndexedItem("cntr_btof11",    m_ngch,     m_cntr_btof11);
            status = m_tuple1->addIndexedItem("qual_btof11",    m_ngch,     m_qual_btof11);
            status = m_tuple1->addIndexedItem("tpi_btof11",     m_ngch,     m_tpi_btof11);
            status = m_tuple1->addIndexedItem("tk_btof11",      m_ngch,     m_tk_btof11);
            status = m_tuple1->addIndexedItem("toff_pi_btof11", m_ngch,     m_toff_pi_btof11);
            status = m_tuple1->addIndexedItem("tsig_pi_btof11", m_ngch,     m_tsig_pi_btof11);
            
            status = m_tuple1->addIndexedItem("ph_btof12",        m_ngch,      m_ph_btof12);
            status = m_tuple1->addIndexedItem("zhit_btof12",      m_ngch,      m_zhit_btof12);
            status = m_tuple1->addIndexedItem("cntr_btof12",      m_ngch,      m_cntr_btof12);
            status = m_tuple1->addIndexedItem("qual_btof12",      m_ngch,      m_qual_btof12);
            status = m_tuple1->addIndexedItem("tpi_btof12",       m_ngch,      m_tpi_btof12);
            status = m_tuple1->addIndexedItem("tk_btof12",        m_ngch,      m_tk_btof12);
            status = m_tuple1->addIndexedItem("toff_pi_btof12",   m_ngch,      m_toff_pi_btof12);
            status = m_tuple1->addIndexedItem("tsig_pi_btof12",   m_ngch,      m_tsig_pi_btof12);
            
            status = m_tuple1->addIndexedItem("ph_btof21",        m_ngch,     m_ph_btof21);
            status = m_tuple1->addIndexedItem("zhit_btof21",      m_ngch,     m_zhit_btof21);
            status = m_tuple1->addIndexedItem("cntr_btof21",      m_ngch,     m_cntr_btof21);
            status = m_tuple1->addIndexedItem("qual_btof21",      m_ngch,     m_qual_btof21);
            status = m_tuple1->addIndexedItem("tpi_btof21",       m_ngch,     m_tpi_btof21);
            status = m_tuple1->addIndexedItem("tk_btof21",        m_ngch,     m_tk_btof21);
            status = m_tuple1->addIndexedItem("toff_pi_btof21",   m_ngch,     m_toff_pi_btof21);
            status = m_tuple1->addIndexedItem("tsig_pi_btof21",   m_ngch,     m_tsig_pi_btof21);
            
            status = m_tuple1->addIndexedItem("ph_btof22",        m_ngch,     m_ph_btof22);
            status = m_tuple1->addIndexedItem("zhit_btof22",      m_ngch,     m_zhit_btof22);
            status = m_tuple1->addIndexedItem("cntr_btof22",      m_ngch,     m_cntr_btof22);
            status = m_tuple1->addIndexedItem("qual_btof22",      m_ngch,     m_qual_btof22);
            status = m_tuple1->addIndexedItem("tpi_btof22",       m_ngch,     m_tpi_btof22);
            status = m_tuple1->addIndexedItem("tk_btof22",        m_ngch,     m_tk_btof22);
            status = m_tuple1->addIndexedItem("toff_pi_btof22",   m_ngch,     m_toff_pi_btof22);
            status = m_tuple1->addIndexedItem("tsig_pi_btof22",   m_ngch,     m_tsig_pi_btof22);
            
            status = m_tuple1->addIndexedItem("ph_btof1",        m_ngch,      m_ph_btof1);
            status = m_tuple1->addIndexedItem("zhit_btof1",      m_ngch,      m_zhit_btof1);
            status = m_tuple1->addIndexedItem("cntr_btof1",      m_ngch,      m_cntr_btof1);
            status = m_tuple1->addIndexedItem("qual_btof1",      m_ngch,      m_qual_btof1);
            status = m_tuple1->addIndexedItem("tpi_btof1",       m_ngch,      m_tpi_btof1);
            status = m_tuple1->addIndexedItem("tk_btof1",        m_ngch,      m_tk_btof1);
            status = m_tuple1->addIndexedItem("m_tpi_btof1_mea", m_ngch,      m_tpi_btof1_mea);
            status = m_tuple1->addIndexedItem("m_tpi_btof1_exp", m_ngch,      m_tpi_btof1_exp);
            status = m_tuple1->addIndexedItem("m_tpi_btof2_mea", m_ngch,      m_tpi_btof2_mea);
            status = m_tuple1->addIndexedItem("m_tpi_btof2_exp", m_ngch,      m_tpi_btof2_exp);
            status = m_tuple1->addIndexedItem("m_tpi_btof_mea",  m_ngch,      m_tpi_btof_mea);
            status = m_tuple1->addIndexedItem("m_tpi_btof_exp",  m_ngch,      m_tpi_btof_exp);
            status = m_tuple1->addIndexedItem("toff_pi_btof1",   m_ngch,      m_toff_pi_btof1);
            status = m_tuple1->addIndexedItem("tsig_pi_btof1",   m_ngch,      m_tsig_pi_btof1);
            
            status = m_tuple1->addIndexedItem("ph_btof2",       m_ngch,      m_ph_btof2);
            status = m_tuple1->addIndexedItem("zhit_btof2",     m_ngch,      m_zhit_btof2);
            status = m_tuple1->addIndexedItem("cntr_btof2",     m_ngch,      m_cntr_btof2);
            status = m_tuple1->addIndexedItem("qual_btof2",     m_ngch,      m_qual_btof2);
            status = m_tuple1->addIndexedItem("tpi_btof2",      m_ngch,      m_tpi_btof2);
            status = m_tuple1->addIndexedItem("tk_btof2",       m_ngch,      m_tk_btof2);
            status = m_tuple1->addIndexedItem("toff_pi_btof2",  m_ngch,      m_toff_pi_btof2);
            status = m_tuple1->addIndexedItem("tsig_pi_btof2",  m_ngch,      m_tsig_pi_btof2);
            
            status = m_tuple1->addIndexedItem("ph_btof",        m_ngch,      m_ph_btof);
            status = m_tuple1->addIndexedItem("zhit_btof",      m_ngch,      m_zhit_btof);
            status = m_tuple1->addIndexedItem("cntr_btof",      m_ngch,      m_cntr_btof);
            status = m_tuple1->addIndexedItem("qual_btof",      m_ngch,      m_qual_btof);
            status = m_tuple1->addIndexedItem("tpi_btof",       m_ngch,      m_tpi_btof);
            status = m_tuple1->addIndexedItem("tk_btof",        m_ngch,      m_tk_btof);
            status = m_tuple1->addIndexedItem("toff_pi_btof",   m_ngch,      m_toff_pi_btof);
            status = m_tuple1->addIndexedItem("tsig_pi_btof",   m_ngch,      m_tsig_pi_btof);
            
            status = m_tuple1->addIndexedItem("ph_etof",        m_ngch,      m_ph_etof); 
            status = m_tuple1->addIndexedItem("zhit_etof",      m_ngch,      m_zhit_etof);
            status = m_tuple1->addIndexedItem("cntr_etof",      m_ngch,      m_cntr_etof);
            status = m_tuple1->addIndexedItem("qual_etof",      m_ngch,      m_qual_etof);
            status = m_tuple1->addIndexedItem("tpi_etof",       m_ngch,      m_tpi_etof);
            status = m_tuple1->addIndexedItem("tk_etof",        m_ngch,      m_tk_etof);
            status = m_tuple1->addIndexedItem("toff_pi_etof",   m_ngch,      m_toff_pi_etof);
            status = m_tuple1->addIndexedItem("tsig_pi_etof",   m_ngch,      m_tsig_pi_etof);
            
            status = m_tuple1->addIndexedItem("mdedx_pid",      m_ngch,      m_dedx_pid);
            status = m_tuple1->addIndexedItem("mtof_pid",       m_ngch,      m_tof_pid);
            status = m_tuple1->addIndexedItem("mchi_pid",       m_ngch,      m_chi_pid);
            status = m_tuple1->addIndexedItem("mtof1_pid",      m_ngch,      m_tof1_pid);
            status = m_tuple1->addIndexedItem("mtof2_pid",      m_ngch,      m_tof2_pid);
            status = m_tuple1->addIndexedItem("mprob_pid",      m_ngch,      m_prob_pid);
            status = m_tuple1->addIndexedItem("mptrk_pid",      m_ngch,      m_ptrk_pid);
            status = m_tuple1->addIndexedItem("mcost_pid",      m_ngch,      m_cost_pid);
            status = m_tuple1->addItem("mpnp",                  m_pnp);
            status = m_tuple1->addItem("mpnm",                  m_pnm);
            status = m_tuple1->addItem("mpnp_dedx",             m_pnp_dedx);
            status = m_tuple1->addItem("mpnm_dedx",             m_pnm_dedx);
            status = m_tuple1->addItem("mpnp_tof",              m_pnp_tof);
            status = m_tuple1->addItem("mpnm_tof",              m_pnm_tof);
            
            status = m_tuple1->addItem("m_pnp_misk",        m_pnp_misk);
            status = m_tuple1->addItem("m_pnm_misk",        m_pnm_misk);
            status = m_tuple1->addItem("m_pnp_misp",        m_pnp_misp);
            status = m_tuple1->addItem("m_pnm_misp",        m_pnm_misp);
            status = m_tuple1->addItem("m_pnp_dedxmisk",    m_pnp_dedxmisk);
            status = m_tuple1->addItem("m_pnm_dedxmisk",    m_pnm_dedxmisk);
            status = m_tuple1->addItem("m_pnp_dedxmisp",    m_pnp_dedxmisp);
            status = m_tuple1->addItem("m_pnm_dedxmisp",    m_pnm_dedxmisp);
            status = m_tuple1->addItem("m_pnp_tofmisk",     m_pnp_tofmisk);
            status = m_tuple1->addItem("m_pnm_tofmisk",     m_pnm_tofmisk);
            status = m_tuple1->addItem("m_pnp_tofmisp",     m_pnp_tofmisp);
            status = m_tuple1->addItem("m_pnm_tofmisp",     m_pnm_tofmisp);
            status = m_tuple1->addItem("m_pi_prob_pip",     m_pi_prob_pip);      
            status = m_tuple1->addItem("m_ki_prob_pip",     m_ki_prob_pip);
            status = m_tuple1->addItem("m_pr_prob_pip",     m_pr_prob_pip);
            status = m_tuple1->addItem("m_pi_prob_pim",     m_pi_prob_pim);
            status = m_tuple1->addItem("m_ki_prob_pim",     m_ki_prob_pim);
            status = m_tuple1->addItem("m_pr_prob_pim",     m_pr_prob_pim);
            status = m_tuple1->addItem("m_pi_chi_pip",      m_pi_chi_pip);
            status = m_tuple1->addItem("m_pi_chi_pim",      m_pi_chi_pim);
            status = m_tuple1->addItem("m_ki_chi_pip",      m_ki_chi_pip);
            status = m_tuple1->addItem("m_ki_chi_pim",      m_ki_chi_pim);
            status = m_tuple1->addItem("m_pr_chi_pip",      m_pr_chi_pip);
            status = m_tuple1->addItem("m_pr_chi_pim",      m_pr_chi_pim);
            status = m_tuple1->addItem("m_pi_probdedx_pip", m_pi_probdedx_pip);
            status = m_tuple1->addItem("m_ki_probdedx_pip", m_ki_probdedx_pip);
            status = m_tuple1->addItem("m_pr_probdedx_pip", m_pr_probdedx_pip);  
            status = m_tuple1->addItem("m_pi_probdedx_pim", m_pi_probdedx_pim);
            status = m_tuple1->addItem("m_ki_probdedx_pim", m_ki_probdedx_pim);
            status = m_tuple1->addItem("m_pr_probdedx_pim", m_pr_probdedx_pim);
            status = m_tuple1->addItem("m_pi_chidedx_pip",  m_pi_chidedx_pip);
            status = m_tuple1->addItem("m_pi_chidedx_pim",  m_pi_chidedx_pim);
            status = m_tuple1->addItem("m_ki_chidedx_pip",  m_ki_chidedx_pip);
            status = m_tuple1->addItem("m_ki_chidedx_pim",  m_ki_chidedx_pim);
            status = m_tuple1->addItem("m_pr_chidedx_pip",  m_pr_chidedx_pip);
            status = m_tuple1->addItem("m_pr_chidedx_pim",  m_pr_chidedx_pim);
            status = m_tuple1->addItem("m_pi_probtof_pip",  m_pi_probtof_pip);
            status = m_tuple1->addItem("m_ki_probtof_pip",  m_ki_probtof_pip);
            status = m_tuple1->addItem("m_pr_probtof_pip",  m_pr_probtof_pip);
            status = m_tuple1->addItem("m_pi_probtof_pim",  m_pi_probtof_pim);
            status = m_tuple1->addItem("m_ki_probtof_pim",  m_ki_probtof_pim);
            status = m_tuple1->addItem("m_pr_probtof_pim",  m_pr_probtof_pim);
            status = m_tuple1->addItem("m_pi_chitof_pip",   m_pi_chitof_pip);
            status = m_tuple1->addItem("m_pi_chitof_pim",   m_pi_chitof_pim);
            status = m_tuple1->addItem("m_ki_chitof_pip",   m_ki_chitof_pip);
            status = m_tuple1->addItem("m_ki_chitof_pim",   m_ki_chitof_pim);
            status = m_tuple1->addItem("m_pr_chitof_pip",   m_pr_chitof_pip);
            status = m_tuple1->addItem("m_pr_chitof_pim",   m_pr_chitof_pim);
            status = m_tuple1->addItem("m_chitof_pip",      m_chitof_pip);
            status = m_tuple1->addItem("m_chitof_pim",      m_chitof_pim);
            
            status = m_tuple1->addItem("mtof1vo_trk1",      m_tof1vo_trk1);
            status = m_tuple1->addItem("mtof2vo_trk1",      m_tof2vo_trk1);
            status = m_tuple1->addItem("mtof1vo_trk2",      m_tof1vo_trk2);
            status = m_tuple1->addItem("mtof2vo_trk2",      m_tof2vo_trk2);
            
            status = m_tuple1->addIndexedItem("mtofvo1",      m_nggneu,      m_tofvo1);
            status = m_tuple1->addIndexedItem("mtofvo2",      m_nggneu,      m_tofvo2);
            status = m_tuple1->addIndexedItem("numHits",      m_nggneu,      m_numHits);    // Total number of hits
            status = m_tuple1->addIndexedItem("secondmoment", m_nggneu,      m_secondmoment);
            status = m_tuple1->addIndexedItem("mx",           m_nggneu,      m_x);       //  Shower coordinates and errors
            status = m_tuple1->addIndexedItem("my",           m_nggneu,      m_y);
            status = m_tuple1->addIndexedItem("mz",           m_nggneu,      m_z);
            status = m_tuple1->addIndexedItem("cosemc",       m_nggneu,      m_cosemc);   // Shower Counter angles and errors
            status = m_tuple1->addIndexedItem("phiemc",       m_nggneu,      m_phiemc);
            status = m_tuple1->addIndexedItem("energy",       m_nggneu,      m_energy);  // Total energy observed in Emc
            status = m_tuple1->addIndexedItem("eseed",        m_nggneu,      m_eSeed);
            status = m_tuple1->addIndexedItem("me9",          m_nggneu,      m_e3x3); 
            status = m_tuple1->addIndexedItem("me25",         m_nggneu,      m_e5x5); 
            status = m_tuple1->addIndexedItem("mlat",         m_nggneu,      m_lat);
            status = m_tuple1->addIndexedItem("ma20",         m_nggneu,      m_a20);
            status = m_tuple1->addIndexedItem("ma42",         m_nggneu,      m_a42);  
            
            status = m_tuple1->addIndexedItem("m_dr0",        m_ngch,      m_dr0);   
            status = m_tuple1->addIndexedItem("m_phi00",      m_ngch,      m_phi00);
            status = m_tuple1->addIndexedItem("m_kappa0",     m_ngch,      m_kappa0);
            status = m_tuple1->addIndexedItem("m_dz0",        m_ngch,      m_dz0);
            status = m_tuple1->addIndexedItem("m_tanl0",      m_ngch,      m_tanl0);
            status = m_tuple1->addIndexedItem("m_err00",      m_ngch,      m_err00);
            status = m_tuple1->addIndexedItem("m_err01",      m_ngch,      m_err01);
            status = m_tuple1->addIndexedItem("m_err02",      m_ngch,      m_err02);
            status = m_tuple1->addIndexedItem("m_err03",      m_ngch,      m_err03);
            status = m_tuple1->addIndexedItem("m_err04",      m_ngch,      m_err04);
            status = m_tuple1->addIndexedItem("m_err05",      m_ngch,      m_err05);
            status = m_tuple1->addIndexedItem("m_err06",      m_ngch,      m_err06);
            status = m_tuple1->addIndexedItem("m_err07",      m_ngch,      m_err07);
            status = m_tuple1->addIndexedItem("m_err08",      m_ngch,      m_err08);    
            status = m_tuple1->addIndexedItem("m_err09",      m_ngch,      m_err09);
            status = m_tuple1->addIndexedItem("m_err10",      m_ngch,      m_err10);
            status = m_tuple1->addIndexedItem("m_err11",      m_ngch,      m_err11);
            status = m_tuple1->addIndexedItem("m_err12",      m_ngch,      m_err12);
            status = m_tuple1->addIndexedItem("m_err13",      m_ngch,      m_err13);
            status = m_tuple1->addIndexedItem("m_err14",      m_ngch,      m_err14);
            
            status = m_tuple1->addIndexedItem("m_drk0",       m_ngch,      m_drk0);
            status = m_tuple1->addIndexedItem("m_phik00",     m_ngch,      m_phik00);
            status = m_tuple1->addIndexedItem("m_kappak0",    m_ngch,      m_kappak0);
            status = m_tuple1->addIndexedItem("m_dzk0",       m_ngch,      m_dzk0);
            status = m_tuple1->addIndexedItem("m_tanlk0",     m_ngch,      m_tanlk0);
            status = m_tuple1->addIndexedItem("m_errk00",     m_ngch,      m_errk00);
            status = m_tuple1->addIndexedItem("m_errk01",     m_ngch,      m_errk01);
            status = m_tuple1->addIndexedItem("m_errk02",     m_ngch,      m_errk02);
            status = m_tuple1->addIndexedItem("m_errk03",     m_ngch,      m_errk03);
            status = m_tuple1->addIndexedItem("m_errk04",     m_ngch,      m_errk04);
            status = m_tuple1->addIndexedItem("m_errk05",     m_ngch,      m_errk05);
            status = m_tuple1->addIndexedItem("m_errk06",     m_ngch,      m_errk06);
            status = m_tuple1->addIndexedItem("m_errk07",     m_ngch,      m_errk07);
            status = m_tuple1->addIndexedItem("m_errk08",     m_ngch,      m_errk08);    
            status = m_tuple1->addIndexedItem("m_errk09",     m_ngch,      m_errk09);
            status = m_tuple1->addIndexedItem("m_errk10",     m_ngch,      m_errk10);
            status = m_tuple1->addIndexedItem("m_errk11",     m_ngch,      m_errk11);
            status = m_tuple1->addIndexedItem("m_errk12",     m_ngch,      m_errk12);
            status = m_tuple1->addIndexedItem("m_errk13",     m_ngch,      m_errk13);
            status = m_tuple1->addIndexedItem("m_errk14",     m_ngch,      m_errk14);
            
            status = m_tuple1->addItem("kal_stat_pip",    m_kal_stat_pip);
            status = m_tuple1->addItem("kal_stat_pim",    m_kal_stat_pim);
            status = m_tuple1->addItem("kal_stat0_pip",   m_kal_stat0_pip);
            status = m_tuple1->addItem("kal_stat0_pim",   m_kal_stat0_pim);
            status = m_tuple1->addItem("kal_stat1_pip",   m_kal_stat1_pip);
            status = m_tuple1->addItem("kal_stat1_pim",   m_kal_stat1_pim);       
            
            status = m_tuple1->addItem("lptof",           m_lp_tof);
            status = m_tuple1->addItem("lmtof",           m_lm_tof);
            status = m_tuple1->addItem("lptoft0",         m_lp_tof_t0);
            status = m_tuple1->addItem("lmtoft0",         m_lm_tof_t0);
            status = m_tuple1->addItem("lptofflag",       m_lp_tof_flag); //1 barral; 2 endcap
            status = m_tuple1->addItem("lmtofflag",       m_lm_tof_flag); //1 barral; 2 endcap
            status = m_tuple1->addItem("lptof_ph",        m_lp_tof_ph);
            status = m_tuple1->addItem("lmtof_ph",        m_lm_tof_ph);
        }
        else    { 
            Log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1) << endmsg;
            return StatusCode::FAILURE;
        }
    }
    
    //
    //--------end of book--------
    //
    
    Log << MSG::INFO << "successfully return from initialize()" <<endmsg;
    return StatusCode::SUCCESS;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode Testpi::execute() {
  MsgStream Log(msgSvc(), name());
    Log << MSG::INFO << "in execute()" << endreq;
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
    int run=eventHeader->runNumber();
    int event=eventHeader->eventNumber();
    Ncounter0=eventHeader->runNumber();
    m_run = eventHeader->runNumber();  
    m_rec = eventHeader->eventNumber();
    
    //MC information
    HepLorentzVector pip_mctruth,pim_mctruth,pi0_mctruth;
    double par[5]={0.0, 0.0, 0.0, 0.0, 0.0};
    if (eventHeader->runNumber() > 0) {
        par[0] = par_data[0];
        par[1] = par_data[1];
        par[2] = par_data[2];
        par[3] = par_data[3];
        par[4] = par_data[4]; 
    } else {
        par[0] = par_mc[0];
        par[1] = par_mc[1];
        par[2] = par_mc[2];
        par[3] = par_mc[3];
        par[4] = par_mc[4];  
    }
    if (eventHeader->runNumber() < 0) {
        SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");
        int m_numParticle = 0;
        if (!mcParticleCol) {
            std::cout << "Could not retrieve McParticelCol" << std::endl;
            return StatusCode::FAILURE;
        }
        else {
            bool psipDecay = false;
            int rootIndex = -1;
            Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
            for (; iter_mc != mcParticleCol->end(); iter_mc++) {
                if ((*iter_mc)->primaryParticle()) continue;
                if (!(*iter_mc)->decayFromGenerator()) continue;
                if ((*iter_mc)->particleProperty() == 443) {
                    psipDecay = true;
                    rootIndex = (*iter_mc)->trackIndex();
                }
                int pid = (*iter_mc)->particleProperty();
                int pidmo = ((*iter_mc)->mother()).particleProperty();
                HepLorentzVector ip = (*iter_mc)->initialFourMomentum();
                if (pid == 211 && pidmo == 113) {
                    pip_mctruth = ip;
                }
                if (pid == -211 && pidmo == 113) {
                    pim_mctruth = ip;
                }
                if (pid == 111 && pidmo == 443) {
                    pi0_mctruth  = ip;
                }
                if (!psipDecay) continue;
                int mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
                int pdgid = (*iter_mc)->particleProperty();
                m_pdgid[m_numParticle] = pdgid;
                m_motheridx[m_numParticle] = mcidx;
                m_numParticle += 1;    
            }
        }  
        m_idxmc = m_numParticle;
    }

    Ncut0++;
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
    Log << MSG::DEBUG << "run, evtnum = " << run << " , " << event << endreq;
    m_nch  = evtRecEvent->totalCharged();
    m_nneu = evtRecEvent->totalNeutral();
    if (evtRecEvent->totalCharged() > 10) return StatusCode::SUCCESS;
    if (evtRecEvent->totalNeutral() > 50) return StatusCode::SUCCESS;
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(), EventModel::EvtRec::EvtRecTrackCol);
    
    //
    // check x0, y0, z0, r0
    // suggest cut: |z0|<5 && r0<1
    //
    Vint iGood, ipip, ipim, ipnp,ipnm;
    iGood.clear();
    ipip.clear();
    ipim.clear();
    ipnp.clear();
    ipnm.clear();
    Vp4 ppip, ppim , ppnp, ppnm;
    ppip.clear();
    ppim.clear();
    ppnp.clear();
    ppnm.clear();
    Hep3Vector xorigin(0,0,0);
    
    IVertexDbSvc*  vtxsvc;
    Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
    if (vtxsvc->isVertexValid()) { 
        double* dbv = vtxsvc->PrimaryVertex();
        double*  vv = vtxsvc->SigmaPrimaryVertex();
        xorigin.setX(dbv[0]);
        xorigin.setY(dbv[1]);
        xorigin.setZ(dbv[2]);
    }
    
    int nCharge = 0;
    for (int i = 0; i < evtRecEvent->totalCharged(); i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isMdcTrackValid()) continue;
        if (!(*itTrk)->isMdcKalTrackValid()) continue;
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();    
        
        double pch   = mdcTrk->p();
        double x0    = mdcTrk->x();
        double y0    = mdcTrk->y();
        double z0    = mdcTrk->z();
        double phi0  = mdcTrk->helix(1);
        double xv    = xorigin.x();
        double yv    = xorigin.y();
        double Rxy   = fabs((x0 - xv) * cos(phi0) + (y0 - yv) * sin(phi0));
        double m_vx0 = x0;
        double m_vy0 = y0;
        double m_vz0 = z0 - xorigin.z();
        double m_vr0 = Rxy;
        double m_Vct = cos(mdcTrk->theta());
        
        HepVector     a = mdcTrk->helix();
        HepSymMatrix Ea = mdcTrk->err();
        HepPoint3D point0(0., 0., 0.);   // the initial point for MDC recosntruction
        HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
        VFHelix helixip(point0, a, Ea);
        helixip.pivot(IP);
        HepVector vecipa = helixip.a();
        double Rvxy0  = fabs(vecipa[0]);  //the nearest distance to IP in xy plane
        double Rvz0   = vecipa[3];         //the nearest distance to IP in z direction
        double Rvphi0 = vecipa[1];
        
        if (fabs(Rvz0)  >= m_vz0cut) continue;
        if (fabs(Rvxy0) >= m_vr0cut) continue;
        if (fabs(m_Vct) >= m_cthcut) continue;
        
        iGood.push_back(i);
        nCharge += mdcTrk->charge();
    }
  
    //
    // Finish Good Charged Track Selection
    //
    int nGood = iGood.size();
    Log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;
    if ((nGood != 2) || (nCharge!=0)) {
        return StatusCode::SUCCESS;
    }
    Ncut1++;
    
    Vint iGam;
    iGam.clear();
    for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isEmcShowerValid()) continue;
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
        // find the nearest charged track
        double dthe = 200.;
        double dphi = 200.;
        double dang = 200.; 
        for (int j = 0; j < evtRecEvent->totalCharged(); j++) {
            EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + j;
            if (!(*jtTrk)->isExtTrackValid()) continue;
            RecExtTrack *extTrk = (*jtTrk)->extTrack();
            if (extTrk->emcVolumeNumber() == -1) continue;
            Hep3Vector extpos = extTrk->emcPosition();
            double angd = extpos.angle(emcpos);
            double thed = extpos.theta() - emcpos.theta();
            double phid = extpos.deltaPhi(emcpos);
            thed = fmod(thed + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;
            phid = fmod(phid + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;
            if (fabs(thed) < fabs(dthe)) dthe = thed;
            if (fabs(phid) < fabs(dphi)) dphi = phid;
            if (angd < dang) dang = angd;
        }
        if (dang >= 200) continue;
        double eraw = emcTrk->energy();
        dthe = dthe * 180 / (CLHEP::pi);
        dphi = dphi * 180 / (CLHEP::pi);
        dang = dang * 180 / (CLHEP::pi);
        double gamcos  = fabs(cos(emcTrk->theta()));
        double mtdc_ch = emcTrk->time();
        if (!((gamcos < 0.8 && eraw > m_BenergyThreshold) || (gamcos > 0.86 && gamcos < 0.92 && eraw > m_EenergyThreshold))) continue;
        if (mtdc_ch < 0 || mtdc_ch > m_EenergyTdc) continue;
        if (dang < m_gammaAngCut) continue;
        iGam.push_back(i);
    }
  
    //
    // Finish Good Photon Selection
    //
    int nGam = iGam.size();
    
    Log << MSG::DEBUG << "num Good Photon " << nGam  << " , " << evtRecEvent->totalNeutral() << endreq;
    if (nGam < 2) {
        return StatusCode::SUCCESS;
    }
    Ncut2++;
    
    //
    // Assign 4-momentum to each photon
    // 
    Vp4 pGam;
    pGam.clear();
    for (int i = 0; i < nGam; i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i]; 
        RecEmcShower* emcTrk = (*itTrk)->emcShower();
        double eraw = emcTrk->energy();
        double phi  = emcTrk->phi();
        double the  = emcTrk->theta();
        HepLorentzVector ptrk;
        ptrk.setPx(eraw * sin(the) * cos(phi));
        ptrk.setPy(eraw * sin(the) * sin(phi));
        ptrk.setPz(eraw * cos(the));
        ptrk.setE(eraw);
        pGam.push_back(ptrk);
    }
    
    for (int i = 0; i < nGood; i++) {  //for rhopi without PID
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
        if (!(*itTrk)->isMdcTrackValid()) continue;
        if (!(*itTrk)->isMdcKalTrackValid()) continue;
        RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack(); 
        RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
        RecMdcKalTrack::setPidType (RecMdcKalTrack::pion);
        if(mdcTrk->charge() > 0) {
            ipip.push_back(iGood[i]);
            HepLorentzVector ptrk;
            ptrk.setPx(mdcKalTrk->px());
            ptrk.setPy(mdcKalTrk->py());
            ptrk.setPz(mdcKalTrk->pz());
            double p3 = ptrk.mag();
            ptrk.setE(sqrt(p3 * p3 + mpi * mpi));
            ppip.push_back(ptrk);
            m_kal_stat_pip  = mdcKalTrk->stat();
            m_kal_stat0_pip = mdcKalTrk->getStat(0,2);
            m_kal_stat1_pip = mdcKalTrk->getStat(1,2);
        } else {
            ipim.push_back(iGood[i]);
            HepLorentzVector ptrk;
            ptrk.setPx(mdcKalTrk->px());
            ptrk.setPy(mdcKalTrk->py());
            ptrk.setPz(mdcKalTrk->pz());
            double p3 = ptrk.mag();
            ptrk.setE(sqrt(p3 * p3 + mpi * mpi));
            ppim.push_back(ptrk);
            m_kal_stat_pim  = mdcKalTrk->stat();
            m_kal_stat0_pim = mdcKalTrk->getStat(0,2);
            m_kal_stat1_pim = mdcKalTrk->getStat(1,2);
        }
    }// without PID

    int npip = ipip.size();
    int npim = ipim.size();
    Log << MSG::DEBUG << "num of pion " << ipip.size() << "," << ipim.size() << endreq;
    if (npip != 1 || npim != 1) {
        return StatusCode::SUCCESS;
    } 
    Ncut3++;

    //***********************************************//
    // the angle between the two charged tracks      //
    //***********************************************//

    int langcc   = 0;
    double dthec = 200.; 
    double dphic = 200.; 
    double dangc = 200.; 
    for (int i = 0; i < npip; i++) {
        EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + ipip[i] ;
        RecMdcTrack* mdcTrk1 = (*itTrk)->mdcTrack();
        RecMdcKalTrack* mdcKalTrk1 = (*itTrk)->mdcKalTrack();
        Hep3Vector emcpos(ppip[i].px(), ppip[i].py(), ppip[i].pz());
        for (int j = 0; j < npim; j++) {
            EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + ipim[j];
            RecMdcTrack* mdcTrk2 = (*jtTrk)->mdcTrack();
            RecMdcKalTrack* mdcKalTrk2 = (*jtTrk)->mdcKalTrack();
            HepLorentzVector pippim = ppip[i] + ppim[j];
            Hep3Vector extpos(ppim[j].px(), ppim[j].py(), ppim[j].pz());
            double angd = extpos.angle(emcpos);
            double thed = extpos.theta() - emcpos.theta();
            double phid = extpos.deltaPhi(emcpos);
            thed = fmod(thed + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;
            phid = fmod(phid + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;
            m_dthec[langcc] = thed * 180 / (CLHEP::pi);
            m_dphic[langcc] = phid * 180 / (CLHEP::pi);
            m_dangc[langcc] = angd * 180 / (CLHEP::pi);
            m_mspippim[langcc] = pippim.m();
            langcc++;
        }
    }
    m_nangecc = langcc;
    Log << MSG::DEBUG << "liuf_JYZHANG Good Photon "  << endreq;

    //
    // Loop each gamma pair, check ppi0 and pTot 
    //
    double m_m2gg, m_momentpi0;
    HepLorentzVector pTot, p2g;
    //******************************************************//
    //   asign the momentum of MDC and KALFIT               //
    //   the deposite energy of EMC  for charged tracks     //
    //******************************************************// 
    double m_momentpip, m_momentpim, extmot1, extmot2;
    double emcTg1  = 0.0;
    double emcTg2  = 0.0;
    double nlaypi1 = 0.0;
    double nhit1   = 0.0;
    double nlaypi2 = 0.0;
    double nhit2   = 0.0;
    EvtRecTrackIterator itTrk11 = evtRecTrkCol->begin() + ipip[0];
    RecMdcTrack* mdcTrk11 = (*itTrk11)->mdcTrack();
    RecMdcKalTrack* mdcKalTrk11 = (*itTrk11)->mdcKalTrack();
    RecEmcShower* emcTrk11 = (*itTrk11)->emcShower();
    RecMucTrack *mucTrk11 = (*itTrk11)->mucTrack();
    double phi01 = mdcTrk11->helix(1);
    EvtRecTrackIterator itTrk12 = evtRecTrkCol->begin() + ipim[0];
    RecMdcTrack* mdcTrk12 = (*itTrk12)->mdcTrack();
    RecMdcKalTrack* mdcKalTrk12 = (*itTrk12)->mdcKalTrack();
    RecEmcShower* emcTrk12 = (*itTrk12)->emcShower();
    RecMucTrack *mucTrk12 = (*itTrk12)->mucTrack();
    double phi02 = mdcTrk12->helix(1);
    m_vxpin = mdcTrk11->x();
    m_vypin = mdcTrk11->y();
    m_vzpin = mdcTrk11->z()-xorigin.z();
    m_vrpin = fabs((mdcTrk11->x() - xorigin.x()) * cos(phi01) + (mdcTrk11->y() - xorigin.y()) * sin(phi01));
    m_costhepin = cos(mdcTrk11->theta());
    m_momentpip = mdcTrk11->p();
    m_ppx = mdcTrk11->px();
    m_ppy = mdcTrk11->py();
    m_ppz = mdcTrk11->pz();
    m_vxp = mdcKalTrk11->x();
    m_vyp = mdcKalTrk11->y();
    m_vzp = mdcKalTrk11->z() - xorigin.z();
    m_vrp = fabs((mdcKalTrk11->x() - xorigin.x()) * cos(phi01) + (mdcKalTrk11->y() - xorigin.y()) * sin(phi01));
    m_costhep = cos(mdcKalTrk11->theta());
    m_phipin  = mdcTrk11->phi();
    m_phip    = mdcKalTrk11->phi();
    extmot1   = mdcKalTrk11->p();
    m_ppxkal  = mdcKalTrk11->px();
    m_ppykal  = mdcKalTrk11->py();
    m_ppzkal  = mdcKalTrk11->pz();
    m_vxmin = mdcTrk12->x();  
    m_vymin = mdcTrk12->y();  
    m_vzmin = mdcTrk12->z() - xorigin.z();  
    m_vrmin = fabs((mdcTrk12->x() - xorigin.x()) * cos(phi02) + (mdcTrk12->y() - xorigin.y()) * sin(phi02));
    m_costhemin = cos(mdcTrk12->theta());
    m_momentpim = mdcTrk12->p();
    m_mpx = mdcTrk12->px();
    m_mpy = mdcTrk12->py();
    m_mpz = mdcTrk12->pz();
    m_vxm = mdcKalTrk12->x();
    m_vym = mdcKalTrk12->y();
    m_vzm = mdcKalTrk12->z() - xorigin.z();  
    m_vrm = fabs((mdcKalTrk12->x() - xorigin.x()) * cos(phi02) + (mdcKalTrk12->y() - xorigin.y()) * sin(phi02));
    m_costhem = cos(mdcKalTrk12->theta());
    m_phimin  = mdcTrk12->phi();
    m_phim    = mdcKalTrk12->phi();
    extmot2   = mdcKalTrk12->p();
    m_mpxkal  = mdcKalTrk12->px();
    m_mpykal  = mdcKalTrk12->py();
    m_mpzkal  = mdcKalTrk12->pz();

    if ((*itTrk11)->isEmcShowerValid()) {
        emcTg1 = emcTrk11->energy();
    }
    if ((*itTrk12)->isEmcShowerValid()) {
        emcTg2 = emcTrk12->energy();
    }
    if ((*itTrk11)->isMucTrackValid()) {
        nlaypi1 = mucTrk11->numLayers();
        nhit1   = mucTrk11->numHits();
    }
    if ((*itTrk12)->isMucTrackValid()) {
        nlaypi2 = mucTrk12->numLayers();
        nhit2   = mucTrk12->numHits();  
    }
    double tof1vo_t1 = -1;
    double tof2vo_t1 = -1;
    double tof1vo_t2 = -1;
    double tof2vo_t2 = -1;
    if ((*itTrk11)->isExtTrackValid()) {
        RecExtTrack *extTrk11 = (*itTrk11)->extTrack();
        tof1vo_t1 = extTrk11->tof1VolumeNumber();
        tof2vo_t1 = extTrk11->tof2VolumeNumber();
    }
    if ((*itTrk12)->isExtTrackValid()) {
        RecExtTrack *extTrk12 = (*itTrk12)->extTrack();
        tof1vo_t2 = extTrk12->tof1VolumeNumber();
        tof2vo_t2 = extTrk12->tof2VolumeNumber();
    }
    m_tof1vo_trk1 = tof1vo_t1;
    m_tof2vo_trk1 = tof2vo_t1;
    m_tof1vo_trk2 = tof1vo_t2;
    m_tof2vo_trk2 = tof2vo_t2;
    m_laypi1 = nlaypi1;
    m_hit1   = nhit1;  
    m_laypi2 = nlaypi2;
    m_hit2   = nhit2;  

    RecMdcKalTrack *pipTrk = (*(evtRecTrkCol->begin() + ipip[0]))->mdcKalTrack();
    RecMdcKalTrack *pimTrk = (*(evtRecTrkCol->begin() + ipim[0]))->mdcKalTrack();
    WTrackParameter wvpipTrk, wvpimTrk,wvkpTrk, wvkmTrk;
    wvpipTrk = WTrackParameter(mpi, pipTrk->getZHelix(), pipTrk->getZError());
    wvpimTrk = WTrackParameter(mpi, pimTrk->getZHelix(), pimTrk->getZError());
    wvkpTrk  = WTrackParameter(mk, pipTrk->getZHelixK(), pipTrk->getZErrorK());//kaon
    wvkmTrk  = WTrackParameter(mk, pimTrk->getZHelixK(), pimTrk->getZErrorK());//kaon

    //****************************************************//
    //    Test vertex fit                                 //
    //****************************************************//
    HepPoint3D vx(0., 0., 0.);
    HepSymMatrix Evx(3, 0);
    double bx = 1E+6;
    double by = 1E+6;
    double bz = 1E+6;
    Evx[0][0] = bx * bx;
    Evx[1][1] = by * by;
    Evx[2][2] = bz * bz;
    VertexParameter vxpar;
    vxpar.setVx(vx);
    vxpar.setEvx(Evx);
    VertexFit* vtxfit = VertexFit::instance();
    //****************************************************//
    //   if the charged particle is kaon                  //
    //****************************************************//
    double chi5  = 9999.0;
    double jkpi0 = -0.5;
    double jkpkm = 0.0;
    double jkpp0 = 0.0;
    double jkmp0 = 0.0;  
    vtxfit->init();
    vtxfit->AddTrack(0,  wvkpTrk);
    vtxfit->AddTrack(1,  wvkmTrk);
    vtxfit->AddVertex(0, vxpar,0, 1);
    if (vtxfit->Fit(0)) {
        vtxfit->Swim(0);  
        WTrackParameter wkp = vtxfit->wtrk(0);
        WTrackParameter wkm = vtxfit->wtrk(1);
        HepPoint3D xorigin2 = vtxfit->vx(0);
        HepSymMatrix xem2 = vtxfit->Evx(0);
        KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
        //
        //  Apply Kinematic 4C fit
        // 
        for (int i = 0; i < nGam - 1; i++) {
            RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i]))->emcShower();
            for(int j = i + 1; j < nGam; j++) {
                RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[j]))->emcShower();
                kmfit->init();
                kmfit->setEspread(0.00085);
                kmfit->setBeamPosition(xorigin2);
                kmfit->setVBeamPosition(xem2);  
                kmfit->AddTrack(0, wkp);
                kmfit->AddTrack(1, wkm);
                kmfit->AddTrack(2, 0.0, g1Trk);
                kmfit->AddTrack(3, 0.0, g2Trk);
                kmfit->AddFourMomentum(0, ecms);
                bool oksq = kmfit->Fit();
                if (oksq) {
                    double chi2 = kmfit->chisq();
                    if (chi2 < chi5) {
                        HepLorentzVector kpi0 = kmfit->pfit(2) + kmfit->pfit(3);
                        HepLorentzVector kpkm = kmfit->pfit(0) + kmfit->pfit(1);
                        HepLorentzVector kpp0 = kmfit->pfit(0) + kpi0;
                        HepLorentzVector kmp0 = kmfit->pfit(1) + kpi0;
                        chi5  = chi2;
                        jkpi0 = kpi0.m();
                        jkpkm = kpkm.m();
                        jkpp0 = kpp0.m();
                        jkmp0 = kmp0.m();
                    }
                }
            }
        }
    } // vetex is true

    //****************************************************//
    //   if the charged particles are pions for real      //
    //****************************************************//
    vtxfit->init();
    vtxfit->AddTrack(0,  wvpipTrk);
    vtxfit->AddTrack(1,  wvpimTrk);
    vtxfit->AddVertex(0, vxpar,0, 1);
    if (!vtxfit->Fit(0)) return SUCCESS;
    vtxfit->Swim(0);
    WTrackParameter wpip = vtxfit->wtrk(0);
    WTrackParameter wpim = vtxfit->wtrk(1);
    HepPoint3D xorigin1 = vtxfit->vx(0);
    HepSymMatrix xem1 = vtxfit->Evx(0);
    KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
    
    //
    //  Apply Kinematic 4C fit
    // 
    double chisq = 9999.;
    int ig1 = -1;
    int ig2 = -1;
    HepLorentzVector gg1, gg2;    
    for (int i = 0; i < nGam - 1; i++) {
            RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i]))->emcShower();
            for (int j = i + 1; j < nGam; j++) {
            RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[j]))->emcShower();
            kmfit->init();
            kmfit->setEspread(0.00085);
            kmfit->setBeamPosition(xorigin1);
            kmfit->setVBeamPosition(xem1);
            kmfit->AddTrack(0, wpip);
            kmfit->AddTrack(1, wpim);
            kmfit->AddTrack(2, 0.0, g1Trk);
            kmfit->AddTrack(3, 0.0, g2Trk);
            kmfit->AddFourMomentum(0, ecms);
            bool oksq = kmfit->Fit();
            if (oksq) {
                double chi2 = kmfit->chisq();
                if (chi2 < chisq) {
                    chisq = chi2;
                    ig1 = iGam[i];
                    ig2 = iGam[j];
                    gg1 = pGam[i];
                    gg2 = pGam[j];
                }
            }
        }
    }
    
    p2g = gg1 + gg2;
    m_pmax = gg1.e()>gg2.e()?gg1.e():gg2.e();
    m_m2gg = p2g.m();
    m_momentpi0 = sqrt(p2g.px() * p2g.px() + p2g.py() * p2g.py() + p2g.pz() * p2g.pz());  
    Log << MSG::DEBUG << " chisq for 4c " << chisq << endreq;
    Ncut4++;
    if (chisq > 200) {
        return StatusCode::SUCCESS; 
    }
    Ncut5++;

    // select charge track and nneu track
    Vint jGood;
    jGood.clear();
    jGood.push_back(ipip[0]);  
    jGood.push_back(ipim[0]);
    m_ngch = jGood.size();
    Vint jGgam;
    jGgam.clear();
    jGgam.push_back(ig1);
    jGgam.push_back(ig2);
    m_nggneu = jGgam.size();
    Vp4 jpGood;
    jpGood.clear();
    jpGood.push_back(ppip[0]);
    jpGood.push_back(ppim[0]);
    Vp4 jpGam;
    jpGam.clear();
    jpGam.push_back(gg1);
    jpGam.push_back(gg2);
    //*******************************************************//
    //   the angle between pmiss and the plane of two angle  //
    //*******************************************************//
    HepLorentzVector pEmiss(0.0, 0.0, 0.0, 0.0);
    double pg2x, pg2y, pg2z;
    pEmiss = ecms - (ppip[0] + ppim[0]);
    double umiss  = pEmiss.e() - pEmiss.rho();
    double umiss2 = pEmiss.e() * pEmiss.e() - pEmiss.rho() * pEmiss.rho();
    m_umiss  = umiss;
    m_umiss2 = umiss2;
    Hep3Vector pTmiss(pEmiss.px(), pEmiss.py(), pEmiss.pz());
    Hep3Vector gT1(gg1.px(), gg1.py(), gg1.pz());
    Hep3Vector gT2(gg2.px(), gg2.py(), gg2.pz());
    Hep3Vector pPmgg = gT1.cross(gT2);
    double angsg = pTmiss.angle(pPmgg);
    double thesg = pTmiss.theta() - pPmgg.theta();
    double phisg = pTmiss.deltaPhi(pPmgg);
    thesg = fmod(thesg + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;
    phisg = fmod(phisg + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) - CLHEP::pi;
    dthesg = thesg * 180 / (CLHEP::pi);
    dphisg = phisg * 180 / (CLHEP::pi);
    dangsg = angsg * 180 / (CLHEP::pi);
    RecEmcShower *g1Trk = (*(evtRecTrkCol->begin()+ig1))->emcShower();
    RecEmcShower *g2Trk = (*(evtRecTrkCol->begin()+ig2))->emcShower();
    kmfit->init();
    kmfit->setEspread(0.00085);
    kmfit->setBeamPosition(xorigin1);
    kmfit->setVBeamPosition(xem1);
    kmfit->AddTrack(0, wpip);
    kmfit->AddTrack(1, wpim);
    kmfit->AddTrack(2, 0.0, g1Trk);
    kmfit->AddTrack(3, 0.0, g2Trk);
    kmfit->AddFourMomentum(0, ecms);
    bool oksq = kmfit->Fit();
    if (!oksq) return SUCCESS;
    HepLorentzVector ppi0  = kmfit->pfit(2) + kmfit->pfit(3);
    HepLorentzVector prho0 = kmfit->pfit(0) + kmfit->pfit(1);
    HepLorentzVector prhop = kmfit->pfit(0) + ppi0;
    HepLorentzVector prhom = kmfit->pfit(1) + ppi0;
    HepLorentzVector pgam2pi1 = prho0 + kmfit->pfit(2);
    HepLorentzVector pgam2pi2 = prho0 + kmfit->pfit(3);
    HepLorentzVector pgam11 = kmfit->pfit(2);
    HepLorentzVector pgam12 = kmfit->pfit(3);
    HepLorentzVector emcpi0 = jpGam[0] + jpGam[1];
    m_chi1 = kmfit->chisq();
    m_mpi0 = ppi0.m();
    m_emcpi0 = emcpi0.m();
    m_prho0 = prho0.m();
    m_prhop = prhop.m();
    m_prhom = prhom.m();
    m_good  = nGood;
    m_gam   = nGam;
    m_pip   = npip;
    m_pim   = npim;
    m_2gam  = m_m2gg;
    m_outpi0 = m_momentpi0;
    m_outpip = m_momentpip;
    m_outpim = m_momentpim;
    m_enpip  = emcTg1;
    m_dcpip  = extmot1;
    m_enpim  = emcTg2;
    m_dcpim  = extmot2;
    m_pipf = kmfit->pfit(0).rho();
    m_pimf = kmfit->pfit(1).rho();
    m_pi0f = ppi0.rho();
    m_chi5 = chi5;
    m_kpi0 = jkpi0;
    m_kpkm = jkpkm;
    m_kpp0 = jkpp0;
    m_kmp0 = jkmp0;
    m_pgam2pi1 = pgam2pi1.m();
    m_pgam2pi2 = pgam2pi2.m();
    cosva1 = kmfit->pfit(2).rho();
    cosva2 = kmfit->pfit(3).rho();

    //
    // check dedx infomation
    //
    for (int ii = 0; ii < m_ngch; ii++) {
        // dedx
        m_ptrk[ii]  = 9999.0;
        m_chie[ii]  = 9999.0;
        m_chimu[ii] = 9999.0;
        m_chipi[ii] = 9999.0;
        m_chik[ii] = 9999.0;
        m_chip[ii] = 9999.0;
        m_ghit[ii] = 9999.0;
        m_thit[ii] = 9999.0;
        m_probPH[ii] = 9999.0;
        m_normPH[ii] = 9999.0;
        m_PHe[ii] = 9999.0;
        m_PHmu[ii] = 9999.0;
        m_PHpion[ii] = 9999.0;
        m_PHkaon[ii] = 9999.0;
        m_PHproton[ii] = 9999.0;
        m_PHreso_e[ii] = 9999.0;
        m_PHreso_mu[ii] = 9999.0;
        m_PHreso_pion[ii] = 9999.0;
        m_PHreso_kaon[ii] = 9999.0;
        m_PHreso_proton[ii] = 9999.0;
        //endtof
        m_tofvo1[ii] = -1;
        m_tofvo2[ii] = -1;
        m_ph_btof11[ii]    = 9999.0;
        m_zhit_btof11[ii]  = 9999.0;
        m_cntr_btof11[ii]  = 9999.0;
        m_qual_btof11[ii]  = 9999.0;
        m_tpi_btof11[ii]   = 9999.0;
        m_tk_btof11[ii]    = 9999.0;
        m_toff_pi_btof11[ii]   =  9999.0;
        m_tsig_pi_btof11[ii]   = 9999.0;
        m_ph_btof12[ii]   = 9999.0;
        m_zhit_btof12[ii]  = 9999.0;
        m_cntr_btof12[ii]  = 9999.0;
        m_qual_btof12[ii]  = 9999.0;
        m_tpi_btof12[ii]   = 9999.0;
        m_tk_btof12[ii]    = 9999.0;
        m_toff_pi_btof12[ii]   =  9999.0;
        m_tsig_pi_btof12[ii]   = 9999.0;
        m_ph_btof21[ii]   = 9999.0;
        m_zhit_btof21[ii]  = 9999.0;
        m_cntr_btof21[ii]  = 9999.0;
        m_qual_btof21[ii]  = 9999.0;
        m_tpi_btof21[ii]   = 9999.0;
        m_tk_btof21[ii]    = 9999.0;
        m_toff_pi_btof21[ii]   =  9999.0;
        m_tsig_pi_btof21[ii]   = 9999.0;
        m_ph_btof22[ii]   = 9999.0;
        m_zhit_btof22[ii]  = 9999.0;
        m_cntr_btof22[ii]  = 9999.0;
        m_qual_btof22[ii]  = 9999.0;
        m_tpi_btof22[ii]   = 9999.0;
        m_tk_btof22[ii]    = 9999.0;
        m_toff_pi_btof22[ii]   =  9999.0;
        m_tsig_pi_btof22[ii]   = 9999.0;
        m_ph_btof1[ii]   = 9999.0;
        m_zhit_btof1[ii]  = 9999.0;
        m_cntr_btof1[ii]  = 9999.0;
        m_qual_btof1[ii]  = 9999.0;
        m_tpi_btof1[ii]   = 9999.0;
        m_tk_btof1[ii]    = 9999.0;
        m_toff_pi_btof1[ii]   =  9999.0;
        m_tsig_pi_btof1[ii]   = 9999.0;
        m_ph_btof2[ii]   = 9999.0;
        m_zhit_btof2[ii]  = 9999.0;
        m_cntr_btof2[ii]  = 9999.0;
        m_qual_btof2[ii]  = 9999.0;
        m_tpi_btof2[ii]   = 9999.0;
        m_tk_btof2[ii]    = 9999.0;
        m_toff_pi_btof2[ii]   =  9999.0;
        m_tsig_pi_btof2[ii]   = 9999.0;
        m_ph_btof[ii]   = 9999.0;
        m_zhit_btof[ii]  = 9999.0;
        m_cntr_btof[ii]  = 9999.0;
        m_qual_btof[ii]  = 9999.0;
        m_tpi_btof[ii]   = 9999.0;
        m_tk_btof[ii]    = 9999.0;
        m_toff_pi_btof[ii]   =  9999.0;
        m_tsig_pi_btof[ii]   = 9999.0;
        m_ph_etof[ii]   = 9999.0;
        m_zhit_etof[ii]  = 9999.0;
        m_cntr_etof[ii]  = 9999.0;
        m_qual_etof[ii]  = 9999.0;
        m_tpi_etof[ii]   = 9999.0;
        m_tk_etof[ii]    = 9999.0;
        m_toff_pi_etof[ii]   =  9999.0;
        m_tsig_pi_etof[ii]   = 9999.0;
        //pid
        m_dedx_pid[ii] = 9999.0;
        m_tof_pid[ii] = 9999.0;
        m_chi_pid[ii] = 9999.0;
        m_tof1_pid[ii] = 9999.0;
        m_tof2_pid[ii] = 9999.0;
        m_prob_pid[ii] = 9999.0; 
        m_ptrk_pid[ii] = 9999.0;  
        m_cost_pid[ii] = 9999.0;
    }
    double bethe_B[5] = {0., 0., 0., 0., 0.};
    double ex_sigma[5] = {0., 0., 0., 0., 0.};
    int indx0 = 0;
    for (int i = 0; i < m_ngch; i++) {
        EvtRecTrackIterator  itTrk = evtRecTrkCol->begin() + jGood[i];
        if (!(*itTrk)->isMdcTrackValid()) continue;
        if (!(*itTrk)->isMdcDedxValid()) continue;
        RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
        RecMdcDedx* dedxTrk = (*itTrk)->mdcDedx();
        RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
        RecMdcKalTrack::setPidType (RecMdcKalTrack::pion);
        m_ptrk[indx0] = mdcKalTrk->p();
        m_chie[indx0] = dedxTrk->chiE();
        m_chimu[indx0] = dedxTrk->chiMu();
        m_chipi[indx0] = dedxTrk->chiPi();
        m_chik[indx0] = dedxTrk->chiK();
        m_chip[indx0] = dedxTrk->chiP();
        m_ghit[indx0] = dedxTrk->numGoodHits();
        m_thit[indx0] = dedxTrk->numTotalHits();
        m_probPH[indx0] = dedxTrk->probPH();
        m_normPH[indx0] = dedxTrk->normPH();
        double ppde = mdcKalTrk->p();
        double sig_the = sin(mdcTrk->theta());
        double Nohit = dedxTrk->numGoodHits();
        for (Int_t it = 0; it < 5; it++) {
            double beta_G = ppde/xmass[it];
            //  -----  beta square  ------
            double beta = beta_G/sqrt(1+(beta_G)*(beta_G));
            double betterm = par[1]-log(par[2]+pow(1/beta_G,par[4]));
            bethe_B[it] = par[0]/pow(beta,par[3])*betterm-par[0];
            double f_betagamma=0.0;
            double g_sinth=0.0;
            double h_nhit=0.0;
            int i_t0=1;
            if (eventHeader->runNumber() > 0) {
                //data
                f_betagamma = 1.41792e+01*pow(beta_G,-2.82457 ) + 2.91780e+01;
                g_sinth = ( 1.59138e-02*sig_the*sig_the+ 4.25527e-02) / (1.59138e-02*0.8034*0.8034+ 4.25527e-02) ;
                h_nhit = ( 1.87903e-04*Nohit*Nohit-1.09517e-02*Nohit+2.15604e-01) / ( 1.87903e-04*26.84*26.84-1.09517e-02*26.84+ 2.15604e-01);
                i_t0 = 1;
            } else {
                //mc
                f_betagamma = 2.82701e+01*pow(beta_G,-2.58827 ) + 2.32476e+01;
                g_sinth = ( -1.84295e-04*sig_the*sig_the+ 4.80499e-02) / (-1.84295e-04*0.7459*0.7459+ 4.80499e-02);
                h_nhit = ( 3.97405e-04*Nohit*Nohit-2.24273e-02*Nohit+3.63487e-01) / ( 3.97405e-04*29.04*29.04-2.24273e-02*29.04+ 3.63487e-01);
                i_t0 = 1;
            }
            ex_sigma[it] = f_betagamma* g_sinth * h_nhit * i_t0;
        }
        m_PHe[indx0]=bethe_B[0];
        m_PHmu[indx0]=bethe_B[1];
        m_PHpion[indx0]=bethe_B[2];
        m_PHkaon[indx0]=bethe_B[3];
        m_PHproton[indx0]=bethe_B[4];
        m_PHreso_e[indx0]=ex_sigma[0];
        m_PHreso_mu[indx0]=ex_sigma[1];
        m_PHreso_pion[indx0]=ex_sigma[2];
        m_PHreso_kaon[indx0]=ex_sigma[3];
        m_PHreso_proton[indx0]=ex_sigma[4];
        indx0++;
    }

    //
    // check TOF infomation
    //
    int indx1=0;
    for (int i = 0; i < m_ngch; i++) {
        EvtRecTrackIterator  itTrk = evtRecTrkCol->begin() + jGood[i];
        if(!(*itTrk)->isMdcTrackValid()) continue;
        if(!(*itTrk)->isTofTrackValid()) continue;
        RecMdcTrack * mdcTrk = (*itTrk)->mdcTrack();
        SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
        double ptrk = mdcTrk->p();
        SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
        for(;iter_tof != tofTrkCol.end(); iter_tof++ ) { 
            TofHitStatus* hitStatus = new TofHitStatus;
            hitStatus->setStatus( (*iter_tof)->status() );    
            double path=(*iter_tof)->path(); // ? 
            double tof  = (*iter_tof)->tof();
            double ph   = (*iter_tof)->ph();
            double rhit = (*iter_tof)->zrhit();
            double qual = 0.0 + (*iter_tof)->quality();
            double cntr = 0.0 + (*iter_tof)->tofID();
            double texp[5];
            double toff[5];
            double tsig[5];
            for(int j = 0; j < 5; j++) {
                texp[j] = (*iter_tof)->texp(j);
                toff[j] = (*iter_tof)->toffset(j);
            }
            if( !(hitStatus->is_raw()) ) {  // have tof information
                if( hitStatus->is_barrel() ) {  // barrel TOF
                    if( hitStatus->is_readout() ) { // single end
                        if( hitStatus->layer()==1 ) { // inner layer
                            if( hitStatus->is_east() ) { // east end
                                m_ph_btof11[indx1]   = ph;
                                m_zhit_btof11[indx1]  = rhit;
                                m_cntr_btof11[indx1]  = cntr;
                                m_qual_btof11[indx1]  = qual;
                                m_tpi_btof11[indx1]   = tof - toff[2]-texp[2];
                                m_tk_btof11[indx1]    = tof - toff[3]-texp[3];
                                m_toff_pi_btof11[indx1]   =  (*iter_tof)->toffset(2);
                                m_tsig_pi_btof11[indx1]   = (*iter_tof)->sigma(2);
                            }
                            else {                       // west end
                                m_ph_btof12[indx1]   = ph;
                                m_zhit_btof12[indx1]  = rhit;
                                m_cntr_btof12[indx1]  = cntr;
                                m_qual_btof12[indx1]  = qual;
                                m_tpi_btof12[indx1]   = tof - toff[2]-texp[2];
                                m_tk_btof12[indx1]    = tof - toff[3]-texp[3];
                                m_toff_pi_btof12[indx1]   =  (*iter_tof)->toffset(2);
                                m_tsig_pi_btof12[indx1]   = (*iter_tof)->sigma(2);
                            }
                        }
                        else if( hitStatus->layer()==2  ) { // outer layer
                            if( hitStatus->is_east() ) { // east end
                                m_ph_btof21[indx1]   = ph;
                                m_zhit_btof21[indx1]  = rhit;
                                m_cntr_btof21[indx1]  = cntr;
                                m_qual_btof21[indx1]  = qual;
                                m_tpi_btof21[indx1]   = tof - toff[2]-texp[2];
                                m_tk_btof21[indx1]    = tof - toff[3]-texp[3];
                                m_toff_pi_btof21[indx1]   =  (*iter_tof)->toffset(2);
                                m_tsig_pi_btof21[indx1]   = (*iter_tof)->sigma(2);
                            }
                            else {                       // west end
                                m_ph_btof22[indx1]   = ph;
                                m_zhit_btof22[indx1]  = rhit;
                                m_cntr_btof22[indx1]  = cntr;
                                m_qual_btof22[indx1]  = qual;
                                m_tpi_btof22[indx1]   = tof - toff[2]-texp[2];
                                m_tk_btof22[indx1]    = tof - toff[3]-texp[3];
                                m_toff_pi_btof22[indx1]   =  (*iter_tof)->toffset(2);
                                m_tsig_pi_btof22[indx1]   = (*iter_tof)->sigma(2);
                            }
                        }
                    }
                    else {
                        if( hitStatus->is_counter() ) { // single layer
                            if( hitStatus->layer()==1 ) { // inner layer
                                m_ph_btof1[indx1]   = ph;
                                m_zhit_btof1[indx1]  = rhit;
                                m_cntr_btof1[indx1]  = cntr;
                                m_qual_btof1[indx1]  = qual;
                                m_tpi_btof1_mea[indx1]   = tof - toff[2];
                                m_tpi_btof1_exp[indx1]   = texp[2];
                                m_tpi_btof1[indx1]   = tof - toff[2]-texp[2];
                                m_tk_btof1[indx1]    = tof - toff[3]-texp[3];
                                m_toff_pi_btof1[indx1]   =  (*iter_tof)->toffset(2);
                                m_tsig_pi_btof1[indx1]   = (*iter_tof)->sigma(2);
                            }
                            else if( hitStatus->layer()==2  ) { // outer layer
                                m_ph_btof2[indx1]   = ph;
                                m_zhit_btof2[indx1]  = rhit;
                                m_cntr_btof2[indx1]  = cntr;
                                m_qual_btof2[indx1]  = qual;
                                m_tpi_btof2_mea[indx1]   = tof - toff[2];
                                m_tpi_btof2_exp[indx1]   = texp[2];
                                m_tpi_btof2[indx1]   = tof - toff[2]-texp[2];
                                m_tk_btof2[indx1]    = tof - toff[3]-texp[3];
                                m_toff_pi_btof2[indx1]   =  (*iter_tof)->toffset(2);
                                m_tsig_pi_btof2[indx1]   = (*iter_tof)->sigma(2);
                            }
                        }
                        else {
                            if( hitStatus->is_cluster() ) { // double layer
                                m_ph_btof[indx1]   = ph;
                                m_zhit_btof[indx1]  = rhit;
                                m_cntr_btof[indx1]  = cntr;
                                m_qual_btof[indx1]  = qual;
                                m_tpi_btof_mea[indx1]   = tof - toff[2];
                                m_tpi_btof_exp[indx1]   = texp[2];
                                m_tpi_btof[indx1]   = tof - toff[2]-texp[2];
                                m_tk_btof[indx1]    = tof - toff[3]-texp[3];
                                m_toff_pi_btof[indx1]   =  (*iter_tof)->toffset(2);
                                m_tsig_pi_btof[indx1]   = (*iter_tof)->sigma(2);
                            }
                        }
                    }
                }
                else { // endcap TOF
                    if( hitStatus->is_readout() ) { // single end
                        m_ph_etof[indx1]   = ph;
                        m_zhit_etof[indx1]  = rhit;
                        m_cntr_etof[indx1]  = cntr;
                        m_qual_etof[indx1]  = qual;
                        m_tpi_etof[indx1]   = tof - toff[2]-texp[2];
                        m_tk_etof[indx1]    = tof - toff[3]-texp[3];
                        m_toff_pi_etof[indx1]   =  (*iter_tof)->toffset(2);
                        m_tsig_pi_etof[indx1]   = (*iter_tof)->sigma(2);
                    }
                }
            }
            delete hitStatus;
        } 
        indx1++;
    } // loop all charged track

    //
    // Assign 4-momentum to each charged track
    //
    Vint ipnp_misk,ipnm_misk,ipnp_misp,ipnm_misp;
    ipnp_misk.clear();
    ipnm_misk.clear();
    ipnp_misp.clear();
    ipnm_misp.clear();
    double pi_chi_pip=9999.0;
    double ki_chi_pip=9999.0;
    double pr_chi_pip=9999.0;
    double pi_chi_pim=9999.0;
    double ki_chi_pim=9999.0;
    double pr_chi_pim=9999.0;
    double pi_prob_pip=-9999.0;
    double ki_prob_pip=-9999.0;
    double pr_prob_pip=-9999.0;
    double pi_prob_pim=-9999.0;
    double ki_prob_pim=-9999.0;
    double pr_prob_pim=-9999.0;
    int indx2=0;
    ParticleID *pid = ParticleID::instance();
    for(int i = 0; i < m_ngch; i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + jGood[i]; 
        pid->init();
        pid->setMethod(pid->methodProbability());
        pid->setChiMinCut();
        pid->setRecTrack(*itTrk);
        pid->usePidSys(pid->useDedx() | pid->useTofCorr()); // use PID sub-system
        pid->identify(pid->onlyPionKaonProton());
        pid->calculate();
        if(!(pid->IsPidInfoValid())) continue;
        RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
        m_dedx_pid[indx2] = pid->chiDedx(2);
        m_tof_pid[indx2]  = pid->chiTof(2);
        m_tof1_pid[indx2] = pid->chiTof1(2);
        m_tof2_pid[indx2] = pid->chiTof2(2);
        m_chi_pid[indx2]= pid->chi(2);
        m_prob_pid[indx2] = pid->probPion();
        RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
        RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);//PID can set to electron, muon, pion, kaon and proton;The default setting is pion
        m_ptrk_pid[indx2] = mdcKalTrk->p();
        m_cost_pid[indx2] = cos(mdcTrk->theta());
        if(mdcTrk->charge() >0 && pid->probPion() > pid->probKaon()&& pid->probPion() >pid->probProton()) {
            ipnp.push_back(jGood[i]);
            HepLorentzVector ptrk;
            ptrk.setPx(mdcKalTrk->px());
            ptrk.setPy(mdcKalTrk->py());
            ptrk.setPz(mdcKalTrk->pz());
            double p3 = ptrk.mag();
            ptrk.setE(sqrt(p3*p3+mpi*mpi));
            ppnp.push_back(ptrk);
        }     //plus charge with with PID
        if(mdcTrk->charge() <0 && pid->probPion() > pid->probKaon()&& pid->probPion() >pid->probProton()) {
            ipnm.push_back(jGood[i]);
            HepLorentzVector ptrk;
            ptrk.setPx(mdcKalTrk->px());
            ptrk.setPy(mdcKalTrk->py());
            ptrk.setPz(mdcKalTrk->pz());
            double p3 = ptrk.mag();
            ptrk.setE(sqrt(p3*p3+mpi*mpi));
            ppnm.push_back(ptrk);
        }     //minus charge with with PID
        if(mdcTrk->charge() >0 && pid->probKaon()> pid->probPion() && pid->probKaon() >pid->probProton()) {
            ipnp_misk.push_back(jGood[i]);
        }     //plus charge with with PID
        if(mdcTrk->charge() <0 && pid->probKaon() >pid->probPion() && pid->probKaon() >pid->probProton() ) {
            ipnm_misk.push_back(jGood[i]);
        }     //minus charge with with PID
        if(mdcTrk->charge() >0 && pid->probProton() >pid->probPion() && pid->probProton() >pid->probKaon()) {
            ipnp_misp.push_back(jGood[i]);
        }     //plus charge with with PID
        if(mdcTrk->charge() <0 && pid->probProton() >pid->probPion() && pid->probProton() >pid->probKaon()) {
            ipnm_misp.push_back(jGood[i]);
        }     //minus charge with with PID
        if(mdcTrk->charge() >0) {
            pi_chi_pip=pid->chi(2);
            ki_chi_pip=pid->chi(3);
            pr_chi_pip=pid->chi(4);
            pi_prob_pip=pid->probPion();
            ki_prob_pip=pid->probKaon();
            pr_prob_pip=pid->probProton();
        }
        if(mdcTrk->charge() <0) {
            pi_chi_pim=pid->chi(2);
            ki_chi_pim=pid->chi(3);
            pr_chi_pim=pid->chi(4);
            pi_prob_pim=pid->probPion();
            ki_prob_pim=pid->probKaon();
            pr_prob_pim=pid->probProton();
        }
        indx2++;
    }
    m_pnp = ipnp.size();
    m_pnm = ipnm.size();
    m_pnp_misk = ipnp_misk.size();
    m_pnm_misk = ipnm_misk.size();
    m_pnp_misp = ipnp_misp.size();
    m_pnm_misp = ipnm_misp.size();
    m_pi_prob_pip=pi_prob_pip;
    m_ki_prob_pip=ki_prob_pip;
    m_pr_prob_pip=pr_prob_pip;
    m_pi_prob_pim=pi_prob_pim;
    m_ki_prob_pim=ki_prob_pim;
    m_pr_prob_pim=pr_prob_pim;
    m_pi_chi_pip=pi_chi_pip;
    m_pi_chi_pim=pi_chi_pim;
    m_ki_chi_pip=ki_chi_pip;
    m_ki_chi_pim=ki_chi_pim;
    m_pr_chi_pip=pr_chi_pip;
    m_pr_chi_pim=pr_chi_pim;   //new

    Vint ipnp_dedx,ipnm_dedx,ipnp_dedxmisk,ipnm_dedxmisk,ipnp_dedxmisp,ipnm_dedxmisp;
    Vint ipnp_tof,ipnm_tof,ipnp_tofmisk,ipnm_tofmisk,ipnp_tofmisp,ipnm_tofmisp;   
    ipnp_dedx.clear();
    ipnm_dedx.clear();
    ipnp_dedxmisk.clear();
    ipnm_dedxmisk.clear();
    ipnp_dedxmisp.clear();
    ipnm_dedxmisp.clear();
    ipnp_tof.clear();
    ipnm_tof.clear();
    ipnp_tofmisk.clear();
    ipnm_tofmisk.clear();
    ipnp_tofmisp.clear();
    ipnm_tofmisp.clear();
    
    //
    // pid use dedx
    //
    double pi_chidedx_pip=9999.0;
    double pi_chidedx_pim=9999.0;
    double ki_chidedx_pip=9999.0;
    double ki_chidedx_pim=9999.0;
    double pr_chidedx_pip=9999.0;
    double pr_chidedx_pim=9999.0;
    double pi_probdedx_pip=-9999.0;
    double ki_probdedx_pip=-9999.0;
    double pr_probdedx_pip=-9999.0;
    double pi_probdedx_pim=-9999.0;
    double ki_probdedx_pim=-9999.0;
    double pr_probdedx_pim=-9999.0;
    
    int indx12=0;
    for(int i = 0; i < m_ngch; i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + jGood[i]; 
        pid->init();
        pid->setMethod(pid->methodProbability());
        pid->setChiMinCut(10);
        pid->setChiMinCut();
        pid->setRecTrack(*itTrk);
        pid->usePidSys(pid->useDedx()); // use PID sub-system
        pid->identify(pid->onlyPionKaonProton());
        pid->calculate();
        if(!(pid->IsPidInfoValid())) continue;
        RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
        RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
        RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);//PID can set to electron, muon, pion, kaon and proton;The default setting is pion
        if(mdcTrk->charge() >0 && pid->probPion() > pid->probKaon()&& pid->probPion() >pid->probProton()) {
            ipnp_dedx.push_back(jGood[i]);
        }     //plus charge with with PID
        if(mdcTrk->charge() <0 && pid->probPion() > pid->probKaon()&& pid->probPion() >pid->probProton()) {
            ipnm_dedx.push_back(jGood[i]);
        }     //minus charge with with PID
        if(mdcTrk->charge() >0 && pid->probKaon()> pid->probPion() && pid->probKaon() >pid->probProton()) {
            ipnp_dedxmisk.push_back(jGood[i]);
        }     //plus charge with with PID
        if(mdcTrk->charge() <0 && pid->probKaon() >pid->probPion() && pid->probKaon() >pid->probProton() ) {
            ipnm_dedxmisk.push_back(jGood[i]);
        }     //minus charge with with PID
        if(mdcTrk->charge() >0 && pid->probProton() >pid->probPion() && pid->probProton() >pid->probKaon()) {
            ipnp_dedxmisp.push_back(jGood[i]);
        }     //plus charge with with PID
        if(mdcTrk->charge() <0 && pid->probProton() >pid->probPion() && pid->probProton() >pid->probKaon()) {
            ipnm_dedxmisp.push_back(jGood[i]);
        }     //minus charge with with PID
        if(mdcTrk->charge() >0) {
            pi_chidedx_pip=pid->chi(2);
            ki_chidedx_pip=pid->chi(3);
            pr_chidedx_pip=pid->chi(4);
            pi_probdedx_pip=pid->probPion();
            ki_probdedx_pip=pid->probKaon();
            pr_probdedx_pip=pid->probProton();
        }
        if(mdcTrk->charge() <0) {
            pi_chidedx_pim=pid->chi(2);
            ki_chidedx_pim=pid->chi(3);
            pr_chidedx_pim=pid->chi(4);
            pi_probdedx_pim=pid->probPion();
            ki_probdedx_pim=pid->probKaon();
            pr_probdedx_pim=pid->probProton();
        }
    indx12++;
    }

    m_pnp_dedx  = ipnp_dedx.size();
    m_pnm_dedx  = ipnm_dedx.size();
    m_pnp_dedxmisk  = ipnp_dedxmisk.size();
    m_pnm_dedxmisk  = ipnm_dedxmisk.size();
    m_pnp_dedxmisp  = ipnp_dedxmisp.size();
    m_pnm_dedxmisp  = ipnm_dedxmisp.size();
    m_pi_probdedx_pip=pi_probdedx_pip;
    m_ki_probdedx_pip=ki_probdedx_pip;
    m_pr_probdedx_pip=pr_probdedx_pip;
    m_pi_probdedx_pim=pi_probdedx_pim;
    m_ki_probdedx_pim=ki_probdedx_pim;
    m_pr_probdedx_pim=pr_probdedx_pim;
    m_pi_chidedx_pip=pi_chidedx_pip;
    m_pi_chidedx_pim=pi_chidedx_pim;
    m_ki_chidedx_pip=ki_chidedx_pip;
    m_ki_chidedx_pim=ki_chidedx_pim;
    m_pr_chidedx_pip=pr_chidedx_pip;
    m_pr_chidedx_pim=pr_chidedx_pim;   //new

    //
    // pid use tof
    //
    int indx13=0;
    double chitof_pip=9999.0;
    double chitof_pim=9999.0;
    double pi_chitof_pip=9999.0;
    double pi_chitof_pim=9999.0;
    double ki_chitof_pip=9999.0;
    double ki_chitof_pim=9999.0;
    double pr_chitof_pip=9999.0;
    double pr_chitof_pim=9999.0;
    double pi_probtof_pip=-9999.0;
    double ki_probtof_pip=-9999.0;
    double pr_probtof_pip=-9999.0;
    double pi_probtof_pim=-9999.0;
    double ki_probtof_pim=-9999.0;
    double pr_probtof_pim=-9999.0;

    for(int i = 0; i < m_ngch; i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + jGood[i]; 
        if(!(*itTrk)->isExtTrackValid()) {
            m_tofvo1[indx13]=-1;
            m_tofvo2[indx13]=-1;
        } else {
            RecExtTrack *extTrk = (*itTrk)->extTrack();
            m_tofvo1[indx13]=extTrk->tof1VolumeNumber();
            m_tofvo2[indx13]=extTrk->tof2VolumeNumber();
        }
        indx13++;
        pid->init();
        pid->setMethod(pid->methodProbability());
        pid->setChiMinCut();
        pid->setRecTrack(*itTrk);
        pid->usePidSys(pid->useTofCorr()); // use PID sub-system
        pid->identify(pid->onlyPion() | pid->onlyKaon());    // seperater Pion/Kaon
        pid->identify(pid->onlyPionKaonProton());
        pid->calculate();
        if(!(pid->IsPidInfoValid())) continue;
        RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
        RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
        RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);//PID can set to electron, muon, pion, kaon and proton;The default setting is pion
        double pidpion1=pid->chiTofCorr(2);
        double pidkaon1=pid->chiTofCorr(3);
        double pidpron1=pid->chiTofCorr(4); 
        if(mdcTrk->charge() >0 && pid->probPion() > pid->probKaon()&& pid->probPion() >pid->probProton()) {
            ipnp_tof.push_back(jGood[i]);
            chitof_pip=pidpion1;
        }     //plus charge with with PID
        if(mdcTrk->charge() <0 && pid->probPion() > pid->probKaon()&& pid->probPion() >pid->probProton() ) {
            ipnm_tof.push_back(jGood[i]);
            chitof_pim=pidpion1;
        }     //minus charge with with PID
        if(mdcTrk->charge() >0 && pid->probKaon()> pid->probPion() && pid->probKaon() >pid->probProton()) {
            ipnp_tofmisk.push_back(jGood[i]);
        }     //plus charge with with PID
        if(mdcTrk->charge() <0 && pid->probKaon() >pid->probPion() && pid->probKaon() >pid->probProton() ) {
            ipnm_tofmisk.push_back(jGood[i]);
        }     //minus charge with with PID
        if(mdcTrk->charge() >0 && pid->probProton() >pid->probPion() && pid->probProton() >pid->probKaon()) {
            ipnp_tofmisp.push_back(jGood[i]);
        }     //plus charge with with PID
        if(mdcTrk->charge() <0 && pid->probProton() >pid->probPion() && pid->probProton() >pid->probKaon()) {
            ipnm_tofmisp.push_back(jGood[i]);
        }     //minus charge with with PID
        if(mdcTrk->charge() >0) {
            pi_chitof_pip=pidpion1;
            ki_chitof_pip=pidkaon1;
            pr_chitof_pip=pidpron1;
            pi_probtof_pip=pid->probPion();
            ki_probtof_pip=pid->probKaon();
            pr_probtof_pip=pid->probProton();
        }
        if(mdcTrk->charge() <0) {
            pi_chitof_pim=pidpion1;
            ki_chitof_pim=pidkaon1;
            pr_chitof_pim=pidpron1;
            pi_probtof_pim=pid->probPion();
            ki_probtof_pim=pid->probKaon();
            pr_probtof_pim=pid->probProton();
        }
    }

    m_pi_probtof_pip=pi_probtof_pip;
    m_ki_probtof_pip=ki_probtof_pip;
    m_pr_probtof_pip=pr_probtof_pip;
    m_pi_probtof_pim=pi_probtof_pim;
    m_ki_probtof_pim=ki_probtof_pim;
    m_pr_probtof_pim=pr_probtof_pim;
    m_chitof_pip=chitof_pip;
    m_chitof_pim=chitof_pim;
    m_pi_chitof_pip=pi_chitof_pip;
    m_pi_chitof_pim=pi_chitof_pim;
    m_ki_chitof_pip=ki_chitof_pip;
    m_ki_chitof_pim=ki_chitof_pim;
    m_pr_chitof_pip=pr_chitof_pip;
    m_pr_chitof_pim=pr_chitof_pim;   //new
    m_pnp_tof  = ipnp_tof.size();
    m_pnm_tof  = ipnm_tof.size();
    m_pnp_tofmisk  = ipnp_tofmisk.size();
    m_pnm_tofmisk  = ipnm_tofmisk.size();
    m_pnp_tofmisp  = ipnp_tofmisp.size();
    m_pnm_tofmisp  = ipnm_tofmisp.size();
    double tof_l=-999,tof_t0=-999,tof_ph=-999;
    int tofflag=0;
    double tof_lm=-999,tof_t0m=-999,tof_phm=-999;
    int tofflagm=0;
    for(int i = 0; i < m_ngch; i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + jGood[i]; 
        if(!(*itTrk)->isMdcTrackValid()) continue;
        if (!(*itTrk)->isMdcKalTrackValid()) continue;
        RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack(); 
        RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
        RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);
        if(mdcTrk->charge() >0) {
            if(!(*itTrk)->isTofTrackValid()) {tof_l=-999;
            } else {
                SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
                SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
                for( ; iter_tof != tofTrkCol.end(); iter_tof++ ) {
                TofHitStatus *status = new TofHitStatus; 
                status->setStatus((*iter_tof)->status());
                    if( !( status->is_raw())){
                        if( status->is_cluster()){  //1-2-3-4
                            tof_l = (*iter_tof)->tof();
                            tof_t0 = (*iter_tof)->t0();
                            tof_ph   = (*iter_tof)->ph();
                            if( status->is_barrel()){
                                tofflag=1;
                            } else {
                                tofflag = 2;
                            }
                        }
                    }    
                }// tof loop finish!
            }//isTofTrackValid judgement finish!
        }
        if(mdcTrk->charge() <0) {
            if(!(*itTrk)->isTofTrackValid()) {tof_lm=-999;
            } else {
                SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
                SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
                for( ; iter_tof != tofTrkCol.end(); iter_tof++ ) {
                    TofHitStatus *status = new TofHitStatus; 
                    status->setStatus((*iter_tof)->status());
                    if( !( status->is_raw())){
                        if( status->is_cluster()){  //1-2-3-4
                            tof_lm = (*iter_tof)->tof();
                            tof_t0m = (*iter_tof)->t0();
                            tof_phm   = (*iter_tof)->ph();
                            if( status->is_barrel()){
                                tofflagm=1;
                            } else {
                                tofflagm = 2;
                            }
                        }
                    }    
                }// tof loop finish!
            }//isTofTrackValid judgement finish!
        }
    }
    m_lp_tof=tof_l;
    m_lp_tof_t0=tof_t0;
    m_lp_tof_flag=tofflag;
    m_lp_tof_ph=tof_ph;
    m_lm_tof=tof_lm;
    m_lm_tof_t0=tof_t0m;
    m_lm_tof_flag=tofflagm;
    m_lm_tof_ph=tof_phm;

    int iphoton = 0;
    for (int i=0; i<m_nggneu; i++) {
        EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + jGgam[i];
        if (!(*itTrk)->isEmcShowerValid()) continue; 
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        m_numHits[iphoton] = emcTrk->numHits();
        m_secondmoment[iphoton] = emcTrk->secondMoment();
        m_x[iphoton] = emcTrk->x();
        m_y[iphoton]= emcTrk->y();
        m_z[iphoton]= emcTrk->z();
        m_cosemc[iphoton]   = cos(emcTrk->theta());
        m_phiemc[iphoton]   = emcTrk->phi();      
        m_energy[iphoton]   = emcTrk->energy();
        m_eSeed[iphoton]    = emcTrk->eSeed(); 
        m_e3x3[iphoton]     = emcTrk->e3x3();  
        m_e5x5[iphoton]     = emcTrk->e5x5();  
        m_lat[iphoton]      = emcTrk->latMoment();
        m_a20[iphoton]      = emcTrk->a20Moment();
        m_a42[iphoton]      = emcTrk->a42Moment();
        iphoton++;
    }

    int idex5 = 0;
    for(int i = 0; i < m_ngch; i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + jGood[i]; 
        if(!(*itTrk)->isMdcTrackValid()) continue;
        if (!(*itTrk)->isMdcKalTrackValid()) continue;
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
        RecMdcKalTrack::setPidType  (RecMdcKalTrack::pion);
        m_dr0[idex5]=mdcTrk->helix(0);
        m_phi00[idex5]=mdcTrk->helix(1);
        m_kappa0[idex5]=mdcTrk->helix(2);
        m_dz0[idex5]=mdcTrk->helix(3);
        m_tanl0[idex5]=mdcTrk->helix(4);
        m_err00[idex5]=mdcTrk->err()[0][0];
        m_err01[idex5]=mdcTrk->err()[0][1];
        m_err02[idex5]=mdcTrk->err()[0][2];
        m_err03[idex5]=mdcTrk->err()[0][3];
        m_err04[idex5]=mdcTrk->err()[0][4];
        m_err05[idex5]=mdcTrk->err()[1][1];
        m_err06[idex5]=mdcTrk->err()[1][2];
        m_err07[idex5]=mdcTrk->err()[1][3];
        m_err08[idex5]=mdcTrk->err()[1][4];
        m_err09[idex5]=mdcTrk->err()[2][2];
        m_err10[idex5]=mdcTrk->err()[2][3];
        m_err11[idex5]=mdcTrk->err()[2][4];
        m_err12[idex5]=mdcTrk->err()[3][3];
        m_err13[idex5]=mdcTrk->err()[3][4];
        m_err14[idex5]=mdcTrk->err()[4][4];
        m_drk0[idex5]=mdcKalTrk->helix()[0];
        m_phik00[idex5]=mdcKalTrk->helix()[1];
        m_kappak0[idex5]=mdcKalTrk->helix()[2];
        m_dzk0[idex5]=mdcKalTrk->helix()[3];
        m_tanlk0[idex5]=mdcKalTrk->helix()[4];
        m_errk00[idex5]=mdcKalTrk->err()[0][0];
        m_errk01[idex5]=mdcKalTrk->err()[0][1];
        m_errk02[idex5]=mdcKalTrk->err()[0][2];
        m_errk03[idex5]=mdcKalTrk->err()[0][3];
        m_errk04[idex5]=mdcKalTrk->err()[0][4];
        m_errk05[idex5]=mdcKalTrk->err()[1][1];
        m_errk06[idex5]=mdcKalTrk->err()[1][2];
        m_errk07[idex5]=mdcKalTrk->err()[1][3];
        m_errk08[idex5]=mdcKalTrk->err()[1][4];
        m_errk09[idex5]=mdcKalTrk->err()[2][2];
        m_errk10[idex5]=mdcKalTrk->err()[2][3];
        m_errk11[idex5]=mdcKalTrk->err()[2][4];
        m_errk12[idex5]=mdcKalTrk->err()[3][3];
        m_errk13[idex5]=mdcKalTrk->err()[3][4];
        m_errk14[idex5]=mdcKalTrk->err()[4][4];
        idex5++;
    }
	m_tuple1->write();

    HepLorentzVector pip_mc,pim_mc,pi0_mc;
    pip_mc=ppip[0];
    pim_mc=ppim[0];
    pi0_mc=emcpi0;
    HepLorentzVector delpip,delpim,delpi0;
    delpip=pip_mctruth-pip_mc;
    delpim=pim_mctruth-pim_mc;
    delpi0=pi0_mctruth-pi0_mc;
    m_dpxpip=delpip.px();
    m_dpypip=delpip.py();
    m_dpzpip=delpip.pz();
    m_dpxpim=delpim.px();
    m_dpypim=delpim.py();
    m_dpzpim=delpim.pz();
    m_dpxpi0=delpi0.px();
    m_dpypi0=delpi0.py();
    m_dpzpi0=delpi0.pz();
    m_pi0_recoil=pi0_mc.m();
    m_dpxpip_no=pip_mctruth.px()-mdcTrk11->px();
    m_dpypip_no=pip_mctruth.py()-mdcTrk11->py();
    m_dpzpip_no=pip_mctruth.pz()-mdcTrk11->pz();
    m_dpxpim_no=pim_mctruth.px()-mdcTrk12->px();
    m_dpypim_no=pim_mctruth.py()-mdcTrk12->py();
    m_dpzpim_no=pim_mctruth.pz()-mdcTrk12->pz();
    m_piptot_mt=pip_mctruth.rho();
    m_pimtot_mt=pim_mctruth.rho();
    m_pi0tot_mt=pi0_mctruth.rho();
    m_mpip_mt=pip_mctruth.m();
    m_mpim_mt=pim_mctruth.m(); 
    m_mpi0_mt=pi0_mctruth.m();
    m_pi0rec=(ecms-pip_mc-pim_mc).m();
    m_piprec=(ecms-pi0_mc-pim_mc).m();
    m_pimrec=(ecms-pi0_mc-pip_mc).m();
    Ncut6++;

    if(fabs(ppi0.m()-0.135)>=0.015) return SUCCESS;
    Ncut7++;
    if(kmfit->chisq()>=chi5) return SUCCESS;
    Ncut8++;
    if(kmfit->chisq()>=100) return SUCCESS;
    Ncut9++;
    if(nGam!=2)  return SUCCESS;
    if(emcTg1/extmot1>=0.8)  return SUCCESS;
    Ncut10++;
    if(emcTg2/extmot2>=0.8)  return SUCCESS;
    Ncut11++;

    return StatusCode::SUCCESS;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode Testpi::finalize() {
    cout<<"total number:         "<<Ncut0<<endl;
    cout<<"nGood==2, nCharge==0: "<<Ncut1<<endl;
    cout<<"nGam>=2:              "<<Ncut2<<endl;
    cout<<"no Pass Pid:          "<<Ncut3<<endl;
    cout<<"pre_vertex fit        "<<Ncut4<<endl;
    cout<<"Pass 4C:              "<<Ncut5<<endl;
    cout<<"initial  event Cut:   "<<Ncut6<<endl;
    cout<<"event Cut_pi0:        "<<Ncut7<<endl;
    cout<<"event Cut_chikk:      "<<Ncut8<<endl;
    cout<<"event Cut_chisq<100:  "<<Ncut9<<endl;
    cout<<"event Cute/p<0.8:     "<<Ncut10<<endl;
    cout<<"event Cute/p<0.8_2:   "<<Ncut11<<endl;
    MsgStream Log(msgSvc(), name());
    Log << MSG::INFO << "in finalize()" << endmsg;
    return StatusCode::SUCCESS;
}

