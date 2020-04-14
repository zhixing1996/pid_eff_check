#include "EvtPreSlt/tools.h"

// Vertex fit
double vertex_fit(VWTrkPara &vwtrkpara, VertexParameter &birth) {
    VertexParameter vxpar;
    Hep3Vector xorigin(0,0,0);
    HepSymMatrix Exorigin(3,0); // error matrix
    double bx = 1E+6, by = 1E+6, bz = 1E+6;
    Exorigin[0][0] = bx*bx; Exorigin[1][1] = by*by; Exorigin[2][2] = bz*bz;
    vxpar.setVx(xorigin); vxpar.setEvx(Exorigin); // vx: vertex
    VertexFit * vtxfit = VertexFit::instance();
    vtxfit->init();
    for (int i = 0; i < vwtrkpara.size(); i++) vtxfit->AddTrack(i, vwtrkpara[i]);
    std::vector<int> trkId(vwtrkpara.size(), 0);
    for (int i = 0; i < vwtrkpara.size(); i++) trkId[i] = i;
    vxpar.setVx(xorigin); vxpar.setEvx(Exorigin); // vx: vertex
    vtxfit->AddVertex(0, vxpar, trkId);
    if (!vtxfit->Fit(0)) return -999.;
    double chi2_vf = vtxfit->chisq(0);
    vtxfit->Swim(0);
    for (int i = 0; i < vwtrkpara.size(); i++) vwtrkpara[i] = vtxfit->wtrk(i);
    if (chi2_vf == 1) std::cout << "==========  VF chi2: " << chi2_vf << std::endl;
    birth = vtxfit->vpar(0);
    return chi2_vf;
}

// Kinematic fit
double kinematic_fit(const VWTrkPara &vwtrkpara, const HepLorentzVector &p4) {
    KinematicFit * kmfit = KinematicFit::instance();
    kmfit->init();
    for (int i = 0; i < vwtrkpara.size(); i++) kmfit->AddTrack(i, vwtrkpara[i]);
    kmfit->AddFourMomentum(0, p4);
    if (!kmfit->Fit()) return -999.;
    return kmfit->chisq();
}

// Kinematic fit with two pi0
double kinematic_fit_two_pi0(const VWTrkPara &vwtrkpara, const HepLorentzVector &p4, Vint iGam) {
    IDataProviderSvc* eventSvc = NULL;
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc, "/Event/EventHeader");
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
    if (!evtRecTrkCol) return -999.;
    KinematicFit * kmfit = KinematicFit::instance();
    double chi2_min = 999.;
    double chi2_kf = 999.;
    for (int i = 0; i < iGam.size(); i++) {
        RecEmcShower *g1Trk = (*(evtRecTrkCol->begin() + iGam[i]))->emcShower();
        for (int j = 0 ; j < iGam.size(); j++) {
            if (i == j) continue;
            RecEmcShower *g2Trk = (*(evtRecTrkCol->begin() + iGam[j]))->emcShower();
            for (int k = 0; k < iGam.size(); k++) {
                if (k == j || k == i) continue;
                RecEmcShower *g3Trk = (*(evtRecTrkCol->begin() + iGam[k]))->emcShower();
                for (int l = 0; l < iGam.size(); l++) {
                    if (l == i || l == j || l == k) continue;
                    RecEmcShower *g4Trk = (*(evtRecTrkCol->begin() + iGam[l]))->emcShower();
                    kmfit->init();
                    for (int m = 0; m < vwtrkpara.size(); m++) kmfit->AddTrack(m, vwtrkpara[m]);
                    kmfit->AddTrack(vwtrkpara.size(), 0.0, g1Trk);
                    kmfit->AddTrack(vwtrkpara.size() + 1, 0.0, g2Trk);
                    kmfit->AddTrack(vwtrkpara.size() + 2, 0.0, g3Trk);
                    kmfit->AddTrack(vwtrkpara.size() + 3, 0.0, g4Trk);
                    kmfit->AddResonance(0, M_Pi0, vwtrkpara.size(), vwtrkpara.size() + 1);
                    kmfit->AddResonance(1, M_Pi0, vwtrkpara.size() + 2, vwtrkpara.size() + 3);
                    kmfit->AddFourMomentum(3, p4);
                    bool oksq = kmfit->Fit();
                    if (oksq) {
                        chi2_kf = kmfit->chisq();
                        if (chi2_kf < chi2_min) {
                            chi2_min  = chi2_kf;
                        }
                    }
                }
            }
        }
    }
    return chi2_kf;
}
