#include "EvtPreSlt/EvtPreSlt.h"

void EvtPreSlt::KKKKKKSel(HepLorentzVector p4Lab, Vint iGood) {
    KKKKKKSelValid = false;
    if (iGood.size() != 6) return ;
    IDataProviderSvc* eventSvc = NULL;
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc, "/Event/EventHeader");
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    if (!evtRecEvent) return ;
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
    if (!evtRecTrkCol) return ;
    int n_Kp = 0, n_Km = 0;
    i_K.clear();
    charge_K.clear();
    p_K.clear();
    goodHits_K.clear();
    cos_K.clear();
    phi_K.clear();
    probPH_K.clear();
    chi_K_sub.clear();
    prob_K_sub.clear();
    chi_K.clear();
    prob_K.clear();
    n_hits_K.clear();
    dEdx_hits_K_sub.clear();
    dt_K_sub.clear();
    beta_K.clear();
    dEdx_hits_K.clear();
    dt_K.clear();
    vwtrkpara_K.clear();
    ParticleID *pid = ParticleID::instance();
    pid->init();
    pid->setMethod(pid->methodProbability());
    pid->setChiMinCut(4);
    for (int i = 0; i < iGood.size(); i++) {
        EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
        if (!(*itTrk)->isMdcTrackValid()) return ;
        if (!(*itTrk)->isMdcKalTrackValid()) return ;
        RecMdcTrack* mdcTrk = (*itTrk)->mdcTrack();
        RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();
        RecMdcKalTrack::setPidType (RecMdcKalTrack::kaon);
        if (mdcTrk->charge() > 0) {
            n_Kp++;
            i_K.push_back(iGood[i]);
            charge_K.push_back(mdcKalTrk->charge());
            p_K.push_back(mdcKalTrk->p());
            cos_K.push_back(mdcKalTrk->pz()/mdcKalTrk->p());
            phi_K.push_back(mdcKalTrk->phi());
            if ((*itTrk)->isMdcDedxValid()) {
                // save MDC information
                RecMdcDedx* dedxtrk = (*itTrk)->mdcDedx();
                probPH_K.push_back(dedxtrk->probPH());
                goodHits_K.push_back(dedxtrk->numGoodHits());
                DedxHitRefVec v_dedx_hit = dedxtrk->getVecDedxHits();
                SmartRefVector<RecMdcDedxHit>::iterator iter_dedxhit = v_dedx_hit.begin();
                int k(0);
                for (int j = 0; j < v_dedx_hit.size(); j++) {
                    Identifier mdcid = (*(iter_dedxhit + j))->getMdcId();
                    int layid = MdcID::layer(mdcid);
                    if (layid < 4) continue;
                    dEdx_hits_K_sub.push_back((*(iter_dedxhit + j))->getDedx());
                    k++;
                }
                n_hits_K.push_back(k);
                dEdx_hits_K.push_back(dEdx_hits_K_sub);
                // save chi and probability from dE/dx
                pid->setRecTrack(*itTrk);
                pid->usePidSys(pid->useDedx()); // use PID sub-system
                pid->identify(pid->all() | pid->onlyElectron() | pid->onlyMuon() | pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());
                pid->calculate();
                for (int j = 0; j < 5; j++) {
                    chi_K_sub.push_back(dedxtrk->chi(j));
                    if (!(pid->IsPidInfoValid())) return ;
                    prob_K_sub.push_back(pid->prob(j));
                }
                chi_K.push_back(chi_K_sub);
                prob_K.push_back(prob_K_sub);

                // save TOF information
                bool isCluster = false;
                if ((*itTrk)->isTofTrackValid()) {
                    SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
                    SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
                    for ( ; iter_tof != tofTrkCol.end(); iter_tof++) {
                        TofHitStatus *status = new TofHitStatus;
                        status->setStatus((*iter_tof)->status());
                        if (status->is_cluster()) {
                            isCluster = true;
                            dt_K_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(0));
                            dt_K_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(1));
                            dt_K_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(2));
                            dt_K_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(3));
                            dt_K_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(4));
                            beta_K.push_back((*iter_tof)->beta());
                            dt_K.push_back(dt_K_sub);
                        }
                    }
                }
            }
            vwtrkpara_K.push_back(WTrackParameter(mass[3], mdcKalTrk->getZHelixK(), mdcKalTrk->getZErrorK()));
        } else {
            n_Km++;
            i_K.push_back(iGood[i]);
            charge_K.push_back(mdcKalTrk->charge());
            p_K.push_back(mdcKalTrk->p());
            cos_K.push_back(mdcKalTrk->pz()/mdcKalTrk->p());
            phi_K.push_back(mdcKalTrk->phi());
            if ((*itTrk)->isMdcDedxValid()) {
                // save MDC information
                RecMdcDedx* dedxtrk = (*itTrk)->mdcDedx();
                probPH_K.push_back(dedxtrk->probPH());
                goodHits_K.push_back(dedxtrk->numGoodHits());
                DedxHitRefVec v_dedx_hit = dedxtrk->getVecDedxHits();
                SmartRefVector<RecMdcDedxHit>::iterator iter_dedxhit = v_dedx_hit.begin();
                int k(0);
                for (int j = 0; j < v_dedx_hit.size(); j++) {
                    Identifier mdcid = (*(iter_dedxhit + j))->getMdcId();
                    int layid = MdcID::layer(mdcid);
                    if (layid < 4) continue;
                    dEdx_hits_K_sub.push_back((*(iter_dedxhit + j))->getDedx());
                    k++;
                }
                n_hits_K.push_back(k);
                dEdx_hits_K.push_back(dEdx_hits_K_sub);
                // save chi and probability from dE/dx
                pid->setRecTrack(*itTrk);
                pid->usePidSys(pid->useDedx()); // use PID sub-system
                pid->identify(pid->all() | pid->onlyElectron() | pid->onlyMuon() | pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());
                pid->calculate();
                for (int j = 0; j < 5; j++) {
                    chi_K_sub.push_back(dedxtrk->chi(j));
                    if (!(pid->IsPidInfoValid())) return ;
                    prob_K_sub.push_back(pid->prob(j));
                }
                chi_K.push_back(chi_K_sub);
                prob_K.push_back(prob_K_sub);

                // save TOF information
                bool isCluster = false;
                if ((*itTrk)->isTofTrackValid()) {
                    SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
                    SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
                    for ( ; iter_tof != tofTrkCol.end(); iter_tof++) {
                        TofHitStatus *status = new TofHitStatus;
                        status->setStatus((*iter_tof)->status());
                        if (status->is_cluster()) {
                            isCluster = true;
                            dt_K_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(0));
                            dt_K_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(1));
                            dt_K_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(2));
                            dt_K_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(3));
                            dt_K_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(4));
                            beta_K.push_back((*iter_tof)->beta());
                            dt_K.push_back(dt_K_sub);
                        }
                    }
                }
            }
            vwtrkpara_K.push_back(WTrackParameter(mass[3], mdcKalTrk->getZHelixK(), mdcKalTrk->getZErrorK()));
        }
    }
    if (n_Kp != 3 || n_Km != 3) return ;
    chi2_vf = vertex_fit(vwtrkpara_K, birth);
    if (chi2_vf == -999.) return ;
    chi2_kf = kinematic_fit(vwtrkpara_K, p4Lab);
    if (chi2_kf == -999.) return ;
    KKKKKKSelValid = true;
    return ;
}
