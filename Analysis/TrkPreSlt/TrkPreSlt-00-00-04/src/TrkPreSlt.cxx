#include "TrkPreSlt/TrkPreSlt.h"
TrkPreSlt* TrkPreSlt::m_pointer = 0;
TrkPreSlt* TrkPreSlt::instance() {
    if (!m_pointer) m_pointer = new TrkPreSlt();
    return m_pointer;
}

void TrkPreSlt::TrkPreSelection(Hep3Vector xorigin, double m_vr0cut, double m_vz0cut, double m_costhetacut) {
    iGood.clear();
    int nCharge = 0;
    IDataProviderSvc* eventSvc = NULL;
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
    if (evtRecEvent->totalCharged() > 50)  return;
    for (int i = 0; i < evtRecEvent->totalCharged(); i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isMdcTrackValid()) continue;
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        double pch = mdcTrk->p();
        double x0 = mdcTrk->x();
        double y0 = mdcTrk->y();
        double z0 = mdcTrk->z();
        double phi0 = mdcTrk->helix(1);
        double xv = xorigin.x();
        double yv = xorigin.y();
        double Rxy = fabs((x0 - xv)*cos(phi0) + (y0 - yv)*sin(phi0));
        double m_costheta = cos(mdcTrk->theta());
        if (pch > 5.0) continue;
        if (fabs(z0) >= m_vz0cut) continue;
        if (Rxy >= m_vr0cut) continue;
        if (fabs(m_costheta) > m_costhetacut) continue;
        double tof_tem = 99.;
        if ((*itTrk)->isTofTrackValid()) {
            SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
            SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
            for ( ; iter_tof != tofTrkCol.end(); iter_tof++) {
                TofHitStatus *status = new TofHitStatus;
                status->setStatus((*iter_tof)->status());
                if((status->is_cluster())) {
                    tof_tem = (*iter_tof)->tof();
                }
            delete status;
            }
        }
        if (tof_tem < 0.) continue;
        iGood.push_back(i);
        nCharge += mdcTrk->charge();
    }
}

void TrkPreSlt::TrkPreSelection(Hep3Vector xorigin, double m_vr0cut, double m_vz0cut, double m_costhetacut, Vint iTrkUsed) {
    iGood.clear();
    int nCharge = 0;
    IDataProviderSvc* eventSvc = NULL;
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
    if (evtRecEvent->totalCharged() > 50) return;
    bool iused;
    for (int i = 0; i < evtRecEvent->totalCharged(); i++) {
        iused = false;
        for (int iu = 0; iu < iTrkUsed.size(); iu++) {
            if( i == iTrkUsed[iu]) iused = true;
        }
        if (iused) continue;
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isMdcTrackValid()) continue;
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        double pch = mdcTrk->p();
        double x0 = mdcTrk->x();
        double y0 = mdcTrk->y();
        double z0 = mdcTrk->z();
        double phi0 = mdcTrk->helix(1);
        double xv = xorigin.x();
        double yv = xorigin.y();
        double Rxy = fabs((x0 - xv)*cos(phi0) + (y0 - yv)*sin(phi0));
        double m_costheta = cos(mdcTrk->theta());
        if (pch > 5.0) continue;    
        if (fabs(z0) >= m_vz0cut) continue;
        if (Rxy >= m_vr0cut) continue;
        if (fabs(m_costheta) > m_costhetacut) continue;
        double tof_tem = 99.;
        if ((*itTrk)->isTofTrackValid()) {
            SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
            SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
            for ( ; iter_tof != tofTrkCol.end(); iter_tof++) {
                TofHitStatus *status = new TofHitStatus;
                status->setStatus((*iter_tof)->status());
                if ((status->is_cluster())) {
                    tof_tem = (*iter_tof)->tof();
                }
                delete status;
            }
        }
        if (tof_tem < 0.) continue;
        iGood.push_back(i);
    }
}

void TrkPreSlt::GamSelection(double m_energyThreshold, double m_gammaAngleCut) {
    iGam.clear();
    IDataProviderSvc* eventSvc = NULL;
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc, "/Event/EventHeader");
    int irun = eventHeader->runNumber();
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
    if (evtRecEvent->totalNeutral() > 20) return;
    for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isEmcShowerValid()) continue;
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
        double dthe = 200.;
        double dphi = 200.;
        double dang = 200.;
        double eraw = 0.0;
        double costheta = cos(emcTrk->theta());
        double phi = emcTrk->phi();
        bool hotchannel = false;
        if (abs(irun) >= 11615 && abs(irun) <= 11655) {
            if (fabs(phi - 0.4480) < 0.0002 && fabs(costheta - 0.2855) < 0.0002) hotchannel = true;
        }
        if (hotchannel) continue;
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
            if (dang >= 200) continue;
            dthe = dthe*180/(CLHEP::pi);
            dphi = dphi*180/(CLHEP::pi);
            dang = dang*180/(CLHEP::pi);
            if (dang<m_gammaAngleCut) continue;
        }
        eraw = emcTrk->energy();
        if (eraw < m_energyThreshold) continue;
        iGam.push_back(i);
    }
}

void TrkPreSlt::GamSelection(double m_energyThreshold, double m_gammaThetaCut, double m_gammaPhiCut) {
    iGam.clear();
    IDataProviderSvc* eventSvc = NULL;
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc, "/Event/EventHeader");
    int irun = eventHeader->runNumber();
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc,  EventModel::EvtRec::EvtRecTrackCol);
    if (evtRecEvent->totalNeutral() > 50) return;
    for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isEmcShowerValid()) continue;
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
        double dthe = 200.;
        double dphi = 200.;
        double dang = 200.;
        double eraw = 0.;
        double costheta = cos(emcTrk->theta());
        double phi = emcTrk->phi();
        bool hotchannel = false;
        if (abs(irun) >= 11615 && abs(irun) <= 11655) {
            if (fabs(phi - 0.4480) < 0.0002 && fabs(costheta - 0.2855) < 0.0002) hotchannel = true;
        }
        if (hotchannel) continue;
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
            if (dang >= 200) continue;
            dthe = dthe*180/(CLHEP::pi);
            dphi = dphi*180/(CLHEP::pi);
            dang = dang*180/(CLHEP::pi);
            if ((fabs(dthe) < m_gammaThetaCut) && (fabs(dphi)<m_gammaPhiCut) ) continue;
        }
        eraw = emcTrk->energy();
        if (eraw < m_energyThreshold) continue;
        iGam.push_back(i);
    }
}

void TrkPreSlt::GamSelection(double m_BarrelEnergyCut, double m_EndcapEnergyCut, double m_gammaAngleCut, int m_TdcMinCut, int m_TdcMaxCut) {
    iGam.clear();
    IDataProviderSvc* eventSvc = NULL;
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc, "/Event/EventHeader");
    int irun = eventHeader->runNumber();
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc,  EventModel::EvtRec::EvtRecTrackCol);
    if (evtRecEvent->totalNeutral()>50)  return;
    for (int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isEmcShowerValid()) continue;
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        double eraw = 0.;
        eraw = emcTrk->energy();
        double costheta = cos(emcTrk->theta());
        double phi = emcTrk->phi();
        bool hotchannel = false;
        bool EnergyCut = false;
        if (abs(irun) >= 11615 && abs(irun) <= 11655) {
            if (fabs(phi - 0.4480) < 0.0002 && fabs(costheta - 0.2855) < 0.0002) hotchannel = true;
        }
        if (hotchannel) continue;
        if (fabs(costheta) < 0.81) {
            if (eraw > m_BarrelEnergyCut) EnergyCut = true;
        }
        if (fabs(costheta) < 0.92 && fabs(costheta) > 0.86) {
            if (eraw > m_EndcapEnergyCut) EnergyCut = true;
        }
        if (!EnergyCut) continue;
        double eraw0, gtime0;
        if (i == 0) {
            eraw0 = emcTrk->energy();
            gtime0 = emcTrk->time();
        }
        double gammatime = emcTrk->time();
        double dtime = gammatime - gtime0;
        if (evtRecEvent->totalCharged() == 0) {
            if (dtime < -10 || dtime > 10) continue;
        }
        Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
        double dthe = 200.;
        double dphi = 200.;
        double dang = 200.;
        if (evtRecEvent->totalCharged() > 0) {
            int emctdc = emcTrk->time();
            if (emctdc > m_TdcMaxCut || emctdc < m_TdcMinCut) continue;
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
            dthe = dthe*180/(CLHEP::pi);
            dphi = dphi*180/(CLHEP::pi);
            dang = dang*180/(CLHEP::pi);
            if (dang < m_gammaAngleCut) continue;
        }
        iGam.push_back(i);
    }
}

void TrkPreSlt::GamSelection(double m_BarrelEnergyCut, double m_EndcapEnergyCut, double m_gammaAngleCut, int m_TdcMinCut, int m_TdcMaxCut, HepLorentzVector p4, double m_AngleOutVolume) {
    iGam.clear();
    IDataProviderSvc* eventSvc = NULL;
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc, "/Event/EventHeader");
    int irun = eventHeader->runNumber();
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc,  EventModel::EvtRec::EvtRecTrackCol);
    if (evtRecEvent->totalNeutral() > 50) return;
    for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isEmcShowerValid()) continue;
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        double eraw = 0.;
        eraw = emcTrk->energy();
        double costheta = cos(emcTrk->theta());
        double the = emcTrk->theta();
        double phi = emcTrk->phi();
        bool hotchannel = false;
        bool EnergyCut = false;
        if (abs(irun) >= 11615 && abs(irun) <= 11655) {
            if (fabs(phi - 0.4480) < 0.0002 && fabs(costheta - 0.2855) < 0.0002) hotchannel = true;
        }
        if (hotchannel) continue;
        if (fabs(costheta) < 0.81) {
            if (eraw > m_BarrelEnergyCut) EnergyCut = true;
        }
        if (fabs(costheta) < 0.92 && fabs(costheta) > 0.86) {
            if (eraw > m_EndcapEnergyCut) EnergyCut = true;
        }
        if (!EnergyCut) continue;
        double eraw0, gtime0;
        if (i == 0) {
            eraw0 = emcTrk->energy();
            gtime0 = emcTrk->time();
        }
        double gammatime = emcTrk->time();
        double dtime = gammatime - gtime0;
        if (evtRecEvent->totalCharged() == 0) {
            if (dtime < -10 || dtime > 10) continue;
        }
        HepLorentzVector p4gam;
        p4gam.setPx(eraw*sin(the)*cos(phi));
        p4gam.setPy(eraw*sin(the)*sin(phi));
        p4gam.setPz(eraw*cos(the));
        p4gam.setE(eraw);
        if (p4gam.vect().angle(p4.vect())*180./3.141593 < m_AngleOutVolume) continue;
        Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
        double dthe = 200.;
        double dphi = 200.;
        double dang = 200.;
        if (evtRecEvent->totalCharged()>0) {
            int emctdc = emcTrk->time();
            if (emctdc > m_TdcMaxCut || emctdc < m_TdcMinCut) continue;
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
            if (dang>=200) continue;
            dthe = dthe*180/(CLHEP::pi);
            dphi = dphi*180/(CLHEP::pi);
            dang = dang*180/(CLHEP::pi);
            if (dang < m_gammaAngleCut) continue;
        }
        iGam.push_back(i);
    }
}

void TrkPreSlt::GamSelection(double m_BarrelEnergyCut, double m_EndcapEnergyCut, double m_gammaAngleCut, int m_TdcMinCut, int m_TdcMaxCut, Vint iGamUsed) {
    iGam.clear();
    IDataProviderSvc* eventSvc = NULL;
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc, "/Event/EventHeader");
    int irun = eventHeader->runNumber();
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
    if (evtRecEvent->totalNeutral() > 50) return;
    bool iused = false;
    for (int i = evtRecEvent->totalCharged(); i< evtRecEvent->totalTracks(); i++) {
        iused = false;
        for (int iu = 0; iu < iGamUsed.size(); iu++){
            if (i == iGamUsed[iu]) iused=true;
        }
        if (iused) continue;    
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isEmcShowerValid()) continue;
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        double eraw = 0.;
        eraw = emcTrk->energy();
        double costheta = cos(emcTrk->theta());
        double phi = emcTrk->phi();
        bool hotchannel = false;
        bool EnergyCut = false;
        if (abs(irun) >= 11615 && abs(irun) <= 11655) {
            if (fabs(phi - 0.4480) < 0.0002 && fabs(costheta - 0.2855) < 0.0002) hotchannel = true;
        }
        if (hotchannel) continue;
        if (fabs(costheta) < 0.81) {
            if (eraw > m_BarrelEnergyCut) EnergyCut = true;
        }
        if (fabs(costheta) < 0.92 && fabs(costheta) > 0.86) {
            if (eraw > m_EndcapEnergyCut) EnergyCut = true;
        }
        if (!EnergyCut) continue;
        double eraw0, gtime0;
        if (i == 0) {
            eraw0 = emcTrk->energy();
            gtime0 = emcTrk->time();
        }
        double gammatime = emcTrk->time();
        double dtime = gammatime - gtime0;
        if (evtRecEvent->totalCharged() == 0) {
            if (dtime < -10 || dtime > 10) continue;
        }
        Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
        double dthe = 200.;
        double dphi = 200.;
        double dang = 200.;
        if (evtRecEvent->totalCharged() > 0) {
            int emctdc = emcTrk->time();
            if (emctdc > m_TdcMaxCut || emctdc < m_TdcMinCut) continue;
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
            dthe = dthe*180/(CLHEP::pi);
            dphi = dphi*180/(CLHEP::pi);
            dang = dang*180/(CLHEP::pi);
            if (dang < m_gammaAngleCut) continue;
        }
        iGam.push_back(i);
    }
}

void TrkPreSlt::GamSelection(Vint iGood, double m_BarrelEnergyCut, double m_EndcapEnergyCut, double m_gammaAngleCut, int m_TdcMinCut, int m_TdcMaxCut) {
    iGam.clear();
    IDataProviderSvc* eventSvc = NULL;
    Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc, "/Event/EventHeader");
    int irun = eventHeader->runNumber();
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
    if (evtRecEvent->totalNeutral() > 50) return;
    for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
        if (!(*itTrk)->isEmcShowerValid()) continue;
        RecEmcShower *emcTrk = (*itTrk)->emcShower();
        double eraw = 0.;
        eraw = emcTrk->energy();
        double costheta = cos(emcTrk->theta());
        double phi = emcTrk->phi();
        bool hotchannel = false;
        if (abs(irun) >= 11615 && abs(irun) <= 11655) {
            if (fabs(phi - 0.4480) < 0.0002 && fabs(costheta - 0.2855) < 0.0002) hotchannel = true;
        }
        if (hotchannel) continue;
        if (fabs(costheta) < 0.81) {
            if (eraw < m_BarrelEnergyCut) continue;
        }
        if (fabs(costheta) < 0.92 && fabs(costheta) > 0.86) {
            if (eraw < m_EndcapEnergyCut) continue;
        }
        Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
        double dthe = 200.;
        double dphi = 200.;
        double dang = 200.;
        if (iGood.size() > 0) {
            int emctdc = emcTrk->time();
            if (emctdc > m_TdcMaxCut || emctdc < m_TdcMinCut) continue;
        }
        double eraw0, gtime0;
        if (iGood.size() == 0) {
            eraw0 = emcTrk->energy();
            gtime0 = emcTrk->time();
        }
        double gammatime = emcTrk->time();
        double dtime = gammatime - gtime0;
        if (evtRecEvent->totalCharged() == 0) {
            if (dtime < -10 || dtime > 10) continue;
        }
        if (iGood.size() > 0) {
            for (int j = 0; j < iGood.size(); j++) {
                EvtRecTrackIterator jtTrk = evtRecTrkCol->begin() + iGood[j];
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
            dthe = dthe*180/(CLHEP::pi);
            dphi = dphi*180/(CLHEP::pi);
            dang = dang*180/(CLHEP::pi);
            if (dang < m_gammaAngleCut) continue;
        }
        iGam.push_back(i);
    }
}

bool TrkPreSlt::Vertex(EvtRecTrackIterator itTrk, Hep3Vector xorigin, double m_vr0cut, double m_vz0cut, double m_CosThetaCut) {
    m_VertexValid = false;
    if ((*itTrk)->isMdcTrackValid()) {
        RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
        double pch = mdcTrk->p();
        double x0 = mdcTrk->x();
        double y0 = mdcTrk->y();
        double z0 = mdcTrk->z();
        double phi0 = mdcTrk->helix(1);
        double xv = xorigin.x();
        double yv = xorigin.y();
        double m_vz = z0 - xorigin.z();
        double m_vr = fabs((x0 - xv)*cos(phi0) + (y0 - yv)*sin(phi0));
        double m_costheta = fabs(cos(mdcTrk->theta()));
        if ((m_vr < m_vr0cut) && (fabs(m_vz) < m_vz0cut) && (m_costheta < m_CosThetaCut)) {
            m_VertexValid = true;
        }
    }
    return m_VertexValid;
}
