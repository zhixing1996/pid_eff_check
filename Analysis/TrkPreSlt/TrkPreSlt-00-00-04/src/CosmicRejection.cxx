#include "TrkPreSlt/TrkPreSlt.h"

void TrkPreSlt::TrkSlt_RejectCosmic(Hep3Vector xorigin, double m_vr0cut,double m_vz0cut,double m_costhetacut) {
    IDataProviderSvc* eventSvc = NULL;
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    if (evtRecEvent->totalCharged() > 40)  return;
    Vdouble Vtof,Vphi;
    Vint    Vchrg;
    Vp3     P3Trk;
    Vtof.clear();
    Vphi.clear();
    Vchrg.clear();
    P3Trk.clear();
    iSltTot.clear();
    iGood.clear();
    iCosmic.clear();
    double Tof = 999.;
    double TexpPro = 999.;
    for (int i = 0; i < evtRecEvent->totalCharged(); i++) {
        IDataProviderSvc* eventSvc = NULL;
        Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
        SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isMdcTrackValid()) continue;
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        double pch = mdcTrk->p();
        if (!(*itTrk)->isTofTrackValid()) continue;
        SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
        SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
        for ( ; iter_tof != tofTrkCol.end(); iter_tof++) {
            TofHitStatus *status = new TofHitStatus;
            status->setStatus((*iter_tof)->status());
            if ((status->is_cluster())) {
                double tofpath_tem = (*iter_tof)->path();
                TexpPro = tofpath_tem*xmass[4]/pch;
                Tof = (*iter_tof)->tof();
            }
            delete status;
        }
        if (Tof > TexpPro + 1 || Tof < 2) continue;
        iSltTot.push_back(i);
        P3Trk.push_back(mdcTrk->p3());
        Vtof.push_back(Tof);
        Vphi.push_back(mdcTrk->helix(1));
        Vchrg.push_back(mdcTrk->charge());
    }
    int nSltTot = iSltTot.size();
    for (int i1 = 0; i1 < nSltTot-1; i1++) { 
        for (int i2 = i1+1; i2<nSltTot; i2++) {
            if (fabs(Vtof[i1] - Vtof[i2]) < 4.) continue;
            if (Vchrg[i1] + Vchrg[i2] != 0) continue;
            double delphi = P3Trk[i1].deltaPhi(P3Trk[i2])*180./3.14159265;
            double delphimdc = (Vphi[i1]-Vphi[i2])*180./3.14159265-180.;
            if (fabs(delphimdc) > 5.) continue;
            iCosmic.push_back(iSltTot[i1]);
            iCosmic.push_back(iSltTot[i2]);
        }
    }
    int nCosmic = iCosmic.size();
    Tof = 999.;
    TexpPro = 999.;    
    bool cosmic;
    for (int i = 0; i < evtRecEvent->totalCharged(); i++) {
        cosmic = false;
        for (int ic = 0; ic < nCosmic; ic++) {
            if (i==iCosmic[ic]) cosmic = true;
        }
        if (cosmic) continue;
        IDataProviderSvc* eventSvc = NULL;
        Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
        SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isMdcTrackValid()) continue;
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        Hep3Vector p3trk = mdcTrk->p3();
        double pch = mdcTrk->p();
        double x0 = mdcTrk->x();
        double y0 = mdcTrk->y();
        double z0 = mdcTrk->z();
        double phi0 = mdcTrk->helix(1);
        double xv = xorigin.x();
        double yv = xorigin.y();
        double Rxy = fabs((x0 - xv)*cos(phi0) + (y0 - yv)*sin(phi0));
        double m_costheta = cos(mdcTrk->theta());
        if (pch>3.0) continue;
        if (fabs(z0) >= m_vz0cut) continue;
        if (Rxy >= m_vr0cut) continue;
        if (fabs(m_costheta) > m_costhetacut) continue;
        if ((*itTrk)->isTofTrackValid()) {
            SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
            SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();    
            for ( ; iter_tof != tofTrkCol.end(); iter_tof++) {
            TofHitStatus *status = new TofHitStatus;
            status->setStatus((*iter_tof)->status());
            if ((status->is_cluster())) {
                double tofpath_tem = (*iter_tof)->path();
                TexpPro = tofpath_tem*xmass[4]/pch;
                Tof = (*iter_tof)->tof();
            }
            delete status;
            }
        }
        if (Tof > TexpPro + 1 || Tof < 2) continue;
        iGood.push_back(i);
    }
}
