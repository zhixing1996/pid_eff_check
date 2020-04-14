#include "EvtPreSlt/EvtPreSlt.h"

void EvtPreSlt::PiPiPi0Pi0Sel(HepLorentzVector p4Lab, Vint iGood, Vint iGam) {
    PiPiPi0Pi0SelValid = false;
    if (iGood.size() != 2) return ;
    if (iGam.size() < 4) return ;
    IDataProviderSvc* eventSvc = NULL;
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc, "/Event/EventHeader");
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc, EventModel::EvtRec::EvtRecEvent);
    if (!evtRecEvent) return ;
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc, EventModel::EvtRec::EvtRecTrackCol);
    if (!evtRecTrkCol) return ;
    int n_pip = 0, n_pim = 0;
    i_pi.clear();
    charge_pi.clear();
    p_pi.clear();
    goodHits_pi.clear();
    cos_pi.clear();
    phi_pi.clear();
    probPH_pi.clear();
    chi_pi_sub.clear();
    prob_pi_sub.clear();
    chi_pi.clear();
    prob_pi.clear();
    n_hits_pi.clear();
    dEdx_hits_pi_sub.clear();
    dt_pi_sub.clear();
    beta_pi.clear();
    dEdx_hits_pi.clear();
    dt_pi.clear();
    vwtrkpara_pi.clear();
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
        RecMdcKalTrack::setPidType (RecMdcKalTrack::pion);
        if (mdcTrk->charge() > 0) {
            n_pip++;
            i_pi.push_back(iGood[i]);
            charge_pi.push_back(mdcKalTrk->charge());
            p_pi.push_back(mdcKalTrk->p());
            cos_pi.push_back(mdcKalTrk->pz()/mdcKalTrk->p());
            phi_pi.push_back(mdcKalTrk->phi());
            if ((*itTrk)->isMdcDedxValid()) {
                // save MDC information
                RecMdcDedx* dedxtrk = (*itTrk)->mdcDedx();
                probPH_pi.push_back(dedxtrk->probPH());
                goodHits_pi.push_back(dedxtrk->numGoodHits());
                DedxHitRefVec v_dedx_hit = dedxtrk->getVecDedxHits();
                SmartRefVector<RecMdcDedxHit>::iterator iter_dedxhit = v_dedx_hit.begin();
                int k(0);
                for (int j = 0; j < v_dedx_hit.size(); j++) {
                    Identifier mdcid = (*(iter_dedxhit + j))->getMdcId();
                    int layid = MdcID::layer(mdcid);
                    if (layid < 4) continue;
                    dEdx_hits_pi_sub.push_back((*(iter_dedxhit + j))->getDedx());
                    k++;
                }
                n_hits_pi.push_back(k);
                dEdx_hits_pi.push_back(dEdx_hits_pi_sub);
                // save chi and probability from dE/dx
                pid->setRecTrack(*itTrk);
                pid->usePidSys(pid->useDedx()); // use PID sub-system
                pid->identify(pid->all() | pid->onlyElectron() | pid->onlyMuon() | pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());
                pid->calculate();
                for (int j = 0; j < 5; j++) {
                    chi_pi_sub.push_back(dedxtrk->chi(j));
                    if (!(pid->IsPidInfoValid())) return ;
                    prob_pi_sub.push_back(pid->prob(j));
                }
                chi_pi.push_back(chi_pi_sub);
                prob_pi.push_back(prob_pi_sub);

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
                            dt_pi_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(0));
                            dt_pi_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(1));
                            dt_pi_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(2));
                            dt_pi_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(3));
                            dt_pi_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(4));
                            beta_pi.push_back((*iter_tof)->beta());
                            dt_pi.push_back(dt_pi_sub);
                        }
                    }
                }
            }
            vwtrkpara_pi.push_back(WTrackParameter(mass[3], mdcKalTrk->getZHelixK(), mdcKalTrk->getZErrorK()));
        } else {
            n_pim++;
            i_pi.push_back(iGood[i]);
            charge_pi.push_back(mdcKalTrk->charge());
            p_pi.push_back(mdcKalTrk->p());
            cos_pi.push_back(mdcKalTrk->pz()/mdcKalTrk->p());
            phi_pi.push_back(mdcKalTrk->phi());
            if ((*itTrk)->isMdcDedxValid()) {
                // save MDC information
                RecMdcDedx* dedxtrk = (*itTrk)->mdcDedx();
                probPH_pi.push_back(dedxtrk->probPH());
                goodHits_pi.push_back(dedxtrk->numGoodHits());
                DedxHitRefVec v_dedx_hit = dedxtrk->getVecDedxHits();
                SmartRefVector<RecMdcDedxHit>::iterator iter_dedxhit = v_dedx_hit.begin();
                int k(0);
                for (int j = 0; j < v_dedx_hit.size(); j++) {
                    Identifier mdcid = (*(iter_dedxhit + j))->getMdcId();
                    int layid = MdcID::layer(mdcid);
                    if (layid < 4) continue;
                    dEdx_hits_pi_sub.push_back((*(iter_dedxhit + j))->getDedx());
                    k++;
                }
                n_hits_pi.push_back(k);
                dEdx_hits_pi.push_back(dEdx_hits_pi_sub);
                // save chi and probability from dE/dx
                pid->setRecTrack(*itTrk);
                pid->usePidSys(pid->useDedx()); // use PID sub-system
                pid->identify(pid->all() | pid->onlyElectron() | pid->onlyMuon() | pid->onlyPion() | pid->onlyKaon() | pid->onlyProton());
                pid->calculate();
                for (int j = 0; j < 5; j++) {
                    chi_pi_sub.push_back(dedxtrk->chi(j));
                    if (!(pid->IsPidInfoValid())) return ;
                    prob_pi_sub.push_back(pid->prob(j));
                }
                chi_pi.push_back(chi_pi_sub);
                prob_pi.push_back(prob_pi_sub);

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
                            dt_pi_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(0));
                            dt_pi_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(1));
                            dt_pi_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(2));
                            dt_pi_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(3));
                            dt_pi_sub.push_back((*iter_tof)->tof() - (*iter_tof)->texp(4));
                            beta_pi.push_back((*iter_tof)->beta());
                            dt_pi.push_back(dt_pi_sub);
                        }
                    }
                }
            }
            vwtrkpara_pi.push_back(WTrackParameter(mass[3], mdcKalTrk->getZHelixK(), mdcKalTrk->getZErrorK()));
        }
    }
    if (n_pip != 1 || n_pim != 1) return ;
    chi2_vf = vertex_fit(vwtrkpara_pi, birth);
    if (chi2_vf == -999.) return ;
    chi2_kf = kinematic_fit_two_pi0(vwtrkpara_pi, p4Lab, iGam);
    if (chi2_kf == -999.) return ;
    PiPiPi0Pi0SelValid = true;
    return ;
}
