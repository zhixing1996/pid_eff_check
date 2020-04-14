#include "TrkPreSlt/TrkPreSlt.h"
                                                                                                          
double TrkPreSlt::GetRxy(EvtRecTrackIterator itTrk, Hep3Vector xorigin) {
    m_Rxy = -99.;
    if ((*itTrk)->isMdcTrackValid()) {
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        double x0 = mdcTrk->x();
        double y0 = mdcTrk->y();
        double z0 = mdcTrk->z();
        double phi0 = mdcTrk->helix(1);
        double xv = xorigin.x();
        double yv = xorigin.y();
        double m_vz = z0 - xorigin.z();
        double m_vr = fabs((x0 - xv)*cos(phi0) + (y0 - yv)*sin(phi0));
        m_Rxy = m_vr;
    }
    return m_Rxy; 
} 

double TrkPreSlt::GetVz(EvtRecTrackIterator itTrk, Hep3Vector xorigin) {
    m_Vz = -99.;
    if ((*itTrk)->isMdcTrackValid()) {
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        double z0 = mdcTrk->z();
        double m_vz = z0-xorigin.z();
        m_Vz = m_vz;
    }
    return m_Vz; 
}

double TrkPreSlt::GetAngle(RecEmcShower *emcTrk) {
    double dthe = 200.;
    double dphi = 200.;
    double dang = 200.;
    m_dang = dang;
    IDataProviderSvc* eventSvc = NULL;
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
    Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
    // find the nearest charged track
    if (evtRecEvent->totalCharged() > 0) {
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
        dthe = dthe*180/(CLHEP::pi);
        dphi = dphi*180/(CLHEP::pi);
        dang = dang*180/(CLHEP::pi);
    }
    m_dang = dang;
    return m_dang;
}
