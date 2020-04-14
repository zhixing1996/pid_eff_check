// Package EvtPreSlt
// Author Bo Zheng 
// it is used to do the initial selection of Bhabha, Dimu, Digam events
// how to use: 
// 1: add  " use EvtPreSlt  EvtPreSlt-*  Analysis " in "requirement"
// 2: add  " #include "EvtPreSlt/EvtPreSlt.h" in the head file
// 4: then you can use the function and the variable defined in this package 

#ifndef EvtPreSlt_selection_H 
#define EvtPreSlt_selection_H

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Service.h"
#include "EventModel/Event.h"
#include "EventModel/EventModel.h"
#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "EvtRecEvent/EvtRecVeeVertex.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "VertexFit/WTrackParameter.h"
#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "ParticleID/ParticleID.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "EvtPreSlt/tools.h"
#include "DstEvent/TofHitStatus.h"

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector< std::vector<int> > VVint;
typedef std::vector< std::vector<double> > VVdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<WTrackParameter> VWTrkPara;

const double PI = 3.1415927;
const double M_Pi0 = 0.1349766;
const double mass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272};

class  EvtPreSlt {
    public:
        ~EvtPreSlt(){;}
        void PiPiPiPiSel(HepLorentzVector , Vint );
        void PiPiPi0Pi0Sel(HepLorentzVector , Vint , Vint );
        void PiPiPiPiPiPiSel(HepLorentzVector , Vint );
        void PiPiPiPiPi0Pi0Sel(HepLorentzVector , Vint , Vint);
        void KKKKSel(HepLorentzVector , Vint );
        void KKPi0Pi0Sel(HepLorentzVector , Vint , Vint );
        void KKKKKKSel(HepLorentzVector , Vint );
        void KKKKPi0Pi0Sel(HepLorentzVector , Vint , Vint );
        static EvtPreSlt * instance();
        bool isPiPiPiPiSelValid() {return PiPiPiPiSelValid;}
        bool isPiPiPi0Pi0SelValid() {return PiPiPi0Pi0SelValid;}
        bool isPiPiPiPiPiPiSelValid() {return PiPiPiPiPiPiSelValid;}
        bool isPiPiPiPiPi0Pi0SelValid() {return PiPiPiPiPi0Pi0SelValid;}
        bool isKKKKSelValid() {return KKKKSelValid;}
        bool isKKPi0Pi0SelValid() {return KKPi0Pi0SelValid;}
        bool isKKKKKKSelValid() {return KKKKKKSelValid;}
        bool isKKKKPi0Pi0SelValid() {return KKKKPi0Pi0SelValid;}
        Vint I_pi() {return i_pi;}
        Vint Charge_pi() {return charge_pi;}
        Vint GoodHits_pi() {return goodHits_pi;}
        Vint N_hits_pi() {return n_hits_pi;}
        Vdouble P_pi() {return p_pi;}
        Vdouble Cos_pi() {return cos_pi;}
        Vdouble Phi_pi() {return phi_pi;}
        Vdouble ProbPH_pi() {return probPH_pi;}
        Vdouble Beta_pi() {return beta_pi;}
        VVdouble Chi_pi() {return chi_pi;}
        VVdouble Prob_pi() {return prob_pi;}
        VVdouble DEdx_hits_pi() {return dEdx_hits_pi;}
        VVdouble Dt_pi() {return dt_pi;}
        Vint I_K() {return i_K;}
        Vint Charge_K() {return charge_K;}
        Vint GoodHits_K() {return goodHits_K;}
        Vint N_hits_K() {return n_hits_K;}
        Vdouble P_K() {return p_K;}
        Vdouble Cos_K() {return cos_K;}
        Vdouble Phi_K() {return phi_K;}
        Vdouble ProbPH_K() {return probPH_K;}
        Vdouble Beta_K() {return beta_K;}
        VVdouble Chi_K() {return chi_K;}
        VVdouble Prob_K() {return prob_K;}
        VVdouble DEdx_hits_K() {return dEdx_hits_K;}
        VVdouble Dt_K() {return dt_K;}
        double Chi2_vf() {return chi2_vf;}
        double Chi2_kf() {return chi2_kf;}
  
    private:
        EvtPreSlt() {;}
        static EvtPreSlt *m_pointer;
        bool PiPiPiPiSelValid;
        bool PiPiPi0Pi0SelValid;
        bool PiPiPiPiPiPiSelValid;
        bool PiPiPiPiPi0Pi0SelValid;
        bool KKKKSelValid;
        bool KKPi0Pi0SelValid;
        bool KKKKKKSelValid;
        bool KKKKPi0Pi0SelValid;
        Vint i_pi, charge_pi, goodHits_pi, n_hits_pi;
        Vdouble p_pi, cos_pi, phi_pi, probPH_pi, chi_pi_sub, prob_pi_sub, dEdx_hits_pi_sub, dt_pi_sub, beta_pi;
        VVdouble chi_pi, prob_pi, dEdx_hits_pi, dt_pi;
        VWTrkPara vwtrkpara_pi;
        Vint i_K, charge_K, goodHits_K, n_hits_K;
        Vdouble p_K, cos_K, phi_K, probPH_K, chi_K_sub, prob_K_sub, dEdx_hits_K_sub, dt_K_sub, beta_K;
        VVdouble chi_K, prob_K, dEdx_hits_K, dt_K;
        VWTrkPara vwtrkpara_K;
        VertexParameter birth;
        double chi2_vf, chi2_kf;
};

#endif
