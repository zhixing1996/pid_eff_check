// Package TrkPreSlt
// Author Bo Zheng 
// it is used to do the initial selection of good track and good gamma, also cut the vertex of the charged track
// how to use: 
// 1: add  " use TrkPreSlt  TrkPreSlt-*  Analysis " in "requirement"
// 2: add  " #include "TrkPreSlt/TrkPreSlt.h" in the head file
// 3: add  " public TrkPreSlt " in the head file
// 4: then you can direct use the function and the variable defined in this package 

#ifndef TrkPreSlt_selection_H 
#define TrkPreSlt_selection_H
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventTag/EventTagSvc.h"
#include "McTruth/McParticle.h"
#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "TrkPreSlt/TrkPreSlt.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#include <vector>
#include "VertexFit/KinematicFit.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/IVertexDbSvc.h"
#include "ParticleID/ParticleID.h"

typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<Hep3Vector> Vp3;
typedef std::vector<double> Vdouble;

const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};

class TrkPreSlt {
    public:
        ~TrkPreSlt(){;}
        static TrkPreSlt* instance();
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // the function used to select the good charged track
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void TrkPreSelection(Hep3Vector xorigin); // if no input vr/vz cut, the initial value is vr0cut=1.0,vz0cut=5.0;
        void TrkPreSelection(Hep3Vector xorigin, double m_VxycCut, double m_VzCut, double m_CosThetaCut); // cut vlaue of vr, vz,costheta
        void TrkPreSelection(Hep3Vector xorigin, double m_vr0cut, double m_vz0cut, double m_costhetacut, Vint iTrkUsed); // cut vlaue of vr, vz,costheta, and thd id not used
        void TrkSlt_RejectCosmic(Hep3Vector xorigin, double m_vr0cut, double m_vz0cut, double m_costhetacut); // cut vlaue of vr, vz,costheta
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // the function used to select the good Gam
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void GamSelection(double m_energyThreshold, double m_gammaAngleCut); // cut for the energy of gamma and the angle between the gamma and the nearest charged track;
        void GamSelection(double m_energyThreshold, double m_gammaThetaCut, double m_gammaPhiCut); // cut for the energy of gamma and the phi, the theta between the gamma and the nearest charged track
        void GamSelection(double m_BarrelEnergyCut, double m_EndcapEnergyCut, double m_gammaAngleCut, int m_TdcMinCut, int m_TdcMaxCut);
        void GamSelection(double m_BarrelEnergyCut, double m_EndcapEnergyCut, double m_gammaAngleCut, int m_TdcMinCut, int m_TdcMaxCut, HepLorentzVector p4, double m_AngleOutVolume);
        void GamSelection(double m_BarrelEnergyCut, double m_EndcapEnergyCut, double m_gammaAngleCut, int m_TdcMinCut, int m_TdcMaxCut, Vint iGamUsed); // cut for the energy of barrel EMC, endcap EMC,  angle, TDC 
        void GamSelection(Vint iGood, double m_BarrelEnergyCut, double m_EndcapEnergyCut, double m_gammaAngleCut, int m_TdcMinCut, int m_TdcMaxCut); // to selecte the gammas with angle cut with good tracks
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // the function used to cut the vertex of the charged track
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        bool Vertex(EvtRecTrackIterator itTrk, Hep3Vector xorigin, double m_vr0cut, double m_vz0cut, double m_CosThetaCut);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // the function to get return value of good charged track selection
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        double GetRxy(EvtRecTrackIterator itTrk, Hep3Vector xorigin); // the return value is the Rxy of the track
        double GetVz(EvtRecTrackIterator itTrk, Hep3Vector xorigin); // the return value is the Vz of the track
        double GetAngle(RecEmcShower* emcTrk); 
        HepLorentzVector GetP4OfRun(int runNo);
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // this is the function to get the four momentum of the run with input the run number
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        Hep3Vector GetXorigin(int runNo);
        Vint igood() {return iGood;}
        Vint icosmic() {return iCosmic;}
        Vint igam() {return iGam;} // the return value of good Gam selection
 
  private:
        TrkPreSlt() {;}
        static TrkPreSlt *m_pointer;
        HepLorentzVector m_p4Lab;
        Hep3Vector m_xorigin;
        Vint iGood, iGam, iSltTot, iCosmic;
        bool m_VertexValid;
        double m_Rxy, m_Vz, m_dang;
};
#endif
