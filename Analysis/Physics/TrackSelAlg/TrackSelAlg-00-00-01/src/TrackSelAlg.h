#ifndef Physics_Analysis_TrackSelAlg_H
#define Physics_Analysis_TrackSelAlg_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GaudiKernel/LoadFactoryEntries.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/PropertyMgr.h"
#include "DatabaseSvc/IDatabaseSvc.h"

// don't know their functions
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "EventModel/Event.h"
/////////////////////////////

#include "EventModel/EventModel.h"
#include "EventModel/EventHeader.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "MdcRecEvent/RecMdcKalTrack.h"

#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/Helix.h"
#include "VertexFit/WTrackParameter.h"
#include "ParticleID/ParticleID.h"

#include "TrkPreSlt/TrkPreSlt.h"
#include "EvtPreSlt/EvtPreSlt.h"

#include "TMath.h"
#include <vector>
#include <TTree.h>

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
using CLHEP::Hep3Vector;
using CLHEP::HepLorentzVector;
typedef std::vector<int> Vint;
typedef std::vector<double> Vdouble;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<WTrackParameter> VWTrkPara;
typedef std::vector< std::vector<int> > VVint;
typedef std::vector< std::vector<double> > VVdouble;

class TrackSelAlg : public Algorithm {

    public:
        TrackSelAlg(const std::string& name, ISvcLocator* pSvcLocator);
        StatusCode initialize();
        StatusCode execute();
        StatusCode finalize();  

    private:
        std::vector<int> m_Modes;
        bool m_isMonteCarlo;
        bool m_pid;
        bool m_debug;
        double m_BarrelEnergyCut;
        double m_EndcapEnergyCut;
        double m_AngleCut;
        int m_EmcTdcMinCut;
        int m_EmcTdcMaxCut;
        double m_KsVxyCut;
        double m_KsVzCut;
        double m_VxyCut;
        double m_VzCut;
        double m_CosThetaCut;  

        // Ntuple1 info
        NTuple::Tuple* m_tuple1;
        NTuple::Item<int> m_runNo;
        NTuple::Item<int> m_evtNo;
        NTuple::Item<int> m_mode;
        NTuple::Item<double> m_beamE;
        NTuple::Item<double> m_chi2_vf;
        NTuple::Item<double> m_chi2_kf;
        NTuple::Item<int> m_n_pi;
        NTuple::Matrix<double> m_trk_pi;
        NTuple::Item<int> m_n_K;
        NTuple::Matrix<double> m_trk_K;
        NTuple::Item<bool> m_valid_PiPiPiPi;
        NTuple::Item<bool> m_valid_PiPiPi0Pi0;
        NTuple::Item<bool> m_valid_PiPiPiPiPiPi;
        NTuple::Item<bool> m_valid_PiPiPiPiPi0Pi0;
        NTuple::Item<bool> m_valid_KKKK;
        NTuple::Item<bool> m_valid_KKPi0Pi0;
        NTuple::Item<bool> m_valid_KKKKKK;
        NTuple::Item<bool> m_valid_KKKKPi0Pi0;

        // common info
        int runNo, evtNo;
        double beamE;
        Hep3Vector xorigin;
        Vint i_pi, change_pi, goodHits_pi;
        Vdouble p_pi, cos_pi, phi_pi, probPH_pi, beta_pi;
        VVdouble chi_pi, prob_pi, dt_pi;
        Vint i_K, change_K, goodHits_K;
        Vdouble p_K, cos_K, phi_K, probPH_K, beta_K;
        VVdouble chi_K, prob_K, dt_K;

        // function
        double ECMS(int runNo);
};
#endif
