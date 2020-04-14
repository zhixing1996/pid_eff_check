// -*- C++ -*- //
//
// Description: e+e- --> 2(pi+pi-), pi+pi-2pi0, 3(pi+pi-), 2(pi+pi-)2pi0, 2(K+K-), K+K-2pi0, 3(K+K-), 2(K+K-)2pi0
//              mode ==  1,         2,          3,         4,             5,       6,        7,       8
//
// Original Author:  Maoqiang JING <jingmq@ihep.ac.cn>
//          Created: [2020-04-10 Fri 17:03] 
//          Inspired by Suyu XIAO's code 
// 
//

#include "TrackSelAlg.h"

// 
// module declare
//


VertexFit * vtxfit = VertexFit::instance();
KinematicFit * kmfit = KinematicFit::instance();

DECLARE_ALGORITHM_FACTORY( TrackSelAlg )
DECLARE_FACTORY_ENTRIES( TrackSelAlg ) {
     
         DECLARE_ALGORITHM( TrackSelAlg );
}

LOAD_FACTORY_ENTRIES( TrackSelAlg )

TrackSelAlg::TrackSelAlg(const std::string& name, ISvcLocator* pSvcLocator) :
	Algorithm(name, pSvcLocator) {
    m_Modes.push_back(0);
    m_Modes.push_back(1);
    m_Modes.push_back(2);
    m_Modes.push_back(3);
    m_Modes.push_back(4);
    m_Modes.push_back(5);
    m_Modes.push_back(6);
    m_Modes.push_back(7);
    declareProperty("DMode", m_Modes);
    declareProperty("IsMonteCarlo", m_isMonteCarlo = false);
    declareProperty("UsePID", m_pid = false);
    declareProperty("Debug", m_debug = false);
    declareProperty("BarrelEnergyCut", m_BarrelEnergyCut = 0.025);
    declareProperty("EndcapEnergyCut", m_EndcapEnergyCut = 0.05);
    declareProperty("AngleCut", m_AngleCut = 10.0);
    declareProperty("EmcTdcMinCut", m_EmcTdcMinCut = 0);
    declareProperty("EmcTdcMaxCut", m_EmcTdcMaxCut = 14);
    declareProperty("VxyCut", m_VxyCut = 1.);
    declareProperty("VzCut", m_VzCut = 10.0);
    declareProperty("CosThetaCut", m_CosThetaCut = 0.93);

}

StatusCode TrackSelAlg::initialize() {
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << ">>>>>>> in initialize()" << endmsg;

    StatusCode status;

    NTuplePtr nt1(ntupleSvc(), "FILE1/track");
    if (nt1) m_tuple1 = nt1;
    else {
        m_tuple1 = ntupleSvc()->book("FILE1/track", CLID_ColumnWiseTuple, "Tracks used for pid efficiency calculation");
        if (m_tuple1) {
            status = m_tuple1->addItem("runNo", m_runNo);
            status = m_tuple1->addItem("evtNo", m_evtNo);
            status = m_tuple1->addItem("mode", m_mode);
            status = m_tuple1->addItem("beamE", m_beamE);
            status = m_tuple1->addItem("n_pi", m_n_pi, 0, 6);
            status = m_tuple1->addIndexedItem("trk_pi", m_n_pi, 22, m_trk_pi);
            status = m_tuple1->addItem("n_K", m_n_K, 0, 6);
            status = m_tuple1->addIndexedItem("trk_K", m_n_K, 22, m_trk_K);
            status = m_tuple1->addItem("chi2_vf", m_chi2_vf);
            status = m_tuple1->addItem("chi2_kf", m_chi2_kf);
            status = m_tuple1->addItem("valid_PiPiPiPi", m_valid_PiPiPiPi);
            status = m_tuple1->addItem("valid_PiPiPi0Pi0", m_valid_PiPiPi0Pi0);
            status = m_tuple1->addItem("valid_PiPiPiPiPiPi", m_valid_PiPiPiPiPiPi);
            status = m_tuple1->addItem("valid_PiPiPiPiPi0Pi0", m_valid_PiPiPiPiPi0Pi0);
            status = m_tuple1->addItem("valid_KKKK", m_valid_KKKK);
            status = m_tuple1->addItem("valid_KKPi0Pi0", m_valid_KKPi0Pi0);
            status = m_tuple1->addItem("valid_KKKKKK", m_valid_KKKKKK);
            status = m_tuple1->addItem("valid_KKKKPi0Pi0", m_valid_KKKKPi0Pi0);
        }
        else {
            log << MSG::ERROR << "Cannot book N-tuple:" << long(m_tuple1) << endmsg;
            return StatusCode::FAILURE;
        }
    }

    log << MSG::INFO << "successfully return from initialize()" << endmsg;
    return StatusCode::SUCCESS;
}

StatusCode TrackSelAlg::execute() {
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "in execute()" << endreq;

    // grt common info
    SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
    if (!eventHeader) return StatusCode::FAILURE;
    m_runNo = eventHeader->runNumber();
    m_evtNo = eventHeader->eventNumber();

    if (m_debug) {
        std::cout << "**************************************" << std::endl;
        std::cout << "run = "  <<  m_runNo  <<  ", event = " << m_evtNo << std::endl;
        std::cout << "**************************************" << std::endl;
    }

    // read beam energy
    int runNo = m_runNo;
    beamE = -999.;
    if (fabs(runNo) > 0) {
        char stmt1[400];
        snprintf(stmt1, 1024,
        				"select BER_PRB, BPR_PRB "
        				"from RunParams where run_number = %d", runNo);
        DatabaseRecordVector res;
        IDatabaseSvc* dbsvc;
        //read db use service
        Gaudi::svcLocator()->service("DatabaseSvc", dbsvc, true);
        int row_no = dbsvc->query("run", stmt1, res);
        if (row_no != 0) {
            DatabaseRecord* records = res[0];
            double E_E = 0, E_P = 0;
            E_E = records->GetDouble("BER_PRB");
            E_P = records->GetDouble("BPR_PRB");
            beamE = (E_E + E_P);
        }
    }
    HepLorentzVector p4Lab(0.011*beamE, 0, 0, beamE);
    if (m_debug) std::cout << "Beam energy: " << beamE << std::endl;

    IVertexDbSvc* vtxsvc;
    Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
    if (vtxsvc->isVertexValid()) {
        double* dbv = vtxsvc->PrimaryVertex();
        double* vv = vtxsvc->SigmaPrimaryVertex();
        xorigin.set(dbv[0],dbv[1],dbv[2]);
    }

    if (m_debug) std::cout << "Vertex origin (x, y, z): (" << xorigin.x() << ", " << xorigin.y() << ", " << xorigin.z() << ")" << std::endl;

    m_valid_PiPiPiPi = false;
    m_valid_PiPiPi0Pi0 = false;
    m_valid_PiPiPiPiPiPi = false;
    m_valid_PiPiPiPiPi0Pi0 = false;
    m_valid_KKKK = false;
    m_valid_KKPi0Pi0 = false;
    m_valid_KKKKKK = false;
    m_valid_KKKKPi0Pi0 = false;
    TrkPreSlt* trkslt = TrkPreSlt::instance();
    trkslt->GamSelection(m_BarrelEnergyCut, m_EndcapEnergyCut, m_AngleCut, m_EmcTdcMinCut, m_EmcTdcMaxCut);
    trkslt->TrkPreSelection(xorigin, m_VxyCut, m_VzCut, m_CosThetaCut);
    if (m_debug) std::cout << "After selection N(charged): " << trkslt->igood().size() << ", N(photon): " << trkslt->igam().size() << std::endl;
    EvtPreSlt* evtslt = EvtPreSlt::instance();

    // mode 1: e+e- --> 2(pi+pi-)
    m_runNo = eventHeader->runNumber();
    m_evtNo = eventHeader->eventNumber();
    m_beamE = beamE;
    m_mode = 1;
    evtslt->PiPiPiPiSel(p4Lab, trkslt->igood());
    i_pi = evtslt->I_pi();
    change_pi = evtslt->Charge_pi();
    p_pi = evtslt->P_pi();
    cos_pi = evtslt->Cos_pi();
    phi_pi = evtslt->Phi_pi();
    chi_pi = evtslt->Chi_pi();
    prob_pi = evtslt->Prob_pi();
    probPH_pi = evtslt->ProbPH_pi();
    goodHits_pi = evtslt->GoodHits_pi();
    dt_pi = evtslt->Dt_pi();
    beta_pi = evtslt->Beta_pi();
    m_n_pi = i_pi.size();
    for (int i = 0; i < m_n_pi; i++) {
        for (int j = 0; j < 5; j++) m_trk_pi[i][j] = prob_pi[i][j];
        for (int j = 0; j < 5; j++) m_trk_pi[i][5 + j] = chi_pi[i][j];
        m_trk_pi[i][10] = change_pi[i];
        m_trk_pi[i][11] = p_pi[i];
        m_trk_pi[i][12] = cos_pi[i];
        m_trk_pi[i][13] = phi_pi[i];
        m_trk_pi[i][14] = probPH_pi[i];
        m_trk_pi[i][15] = goodHits_pi[i];
        for (int j = 0; j < 5; j++) m_trk_pi[i][16 + j] = dt_pi[i][j];
        m_trk_pi[i][21] = beta_pi[i];
    }
    m_chi2_vf = evtslt->Chi2_vf();
    m_chi2_kf = evtslt->Chi2_kf();
    m_valid_PiPiPiPi = evtslt->isPiPiPiPiSelValid();
    if (trkslt->igam().size() != 0) m_valid_PiPiPiPi = false;
    if (m_debug) std::cout << "Mode : e+e- --> 2(pi+pi-), N(pi): " << m_n_pi << ", N(photon): " << trkslt->igam().size() << ", chi2(vf): " << m_chi2_vf << ", chi2(kf): " << m_chi2_kf << std::endl;
    if (m_debug) {
        if (m_valid_PiPiPiPi) std::cout << "----> Status: successful!" << std::endl;
        else std::cout << "----> Status: failed!" << std::endl;
    }
    m_tuple1->write();

    // mode 2: e+e- --> pi+pi-2pi0
    m_runNo = eventHeader->runNumber();
    m_evtNo = eventHeader->eventNumber();
    m_beamE = beamE;
    m_mode = 2;
    evtslt->PiPiPi0Pi0Sel(p4Lab, trkslt->igood(), trkslt->igam());
    i_pi = evtslt->I_pi();
    change_pi = evtslt->Charge_pi();
    p_pi = evtslt->P_pi();
    cos_pi = evtslt->Cos_pi();
    phi_pi = evtslt->Phi_pi();
    chi_pi = evtslt->Chi_pi();
    prob_pi = evtslt->Prob_pi();
    probPH_pi = evtslt->ProbPH_pi();
    goodHits_pi = evtslt->GoodHits_pi();
    dt_pi = evtslt->Dt_pi();
    beta_pi = evtslt->Beta_pi();
    m_n_pi = i_pi.size();
    for (int i = 0; i < m_n_pi; i++) {
        for (int j = 0; j < 5; j++) m_trk_pi[i][j] = prob_pi[i][j];
        for (int j = 0; j < 5; j++) m_trk_pi[i][5 + j] = chi_pi[i][j];
        m_trk_pi[i][10] = change_pi[i];
        m_trk_pi[i][11] = p_pi[i];
        m_trk_pi[i][12] = cos_pi[i];
        m_trk_pi[i][13] = phi_pi[i];
        m_trk_pi[i][14] = probPH_pi[i];
        m_trk_pi[i][15] = goodHits_pi[i];
        for (int j = 0; j < 5; j++) m_trk_pi[i][16 + j] = dt_pi[i][j];
        m_trk_pi[i][21] = beta_pi[i];
    }
    m_chi2_vf = evtslt->Chi2_vf();
    m_chi2_kf = evtslt->Chi2_kf();
    m_valid_PiPiPi0Pi0 = evtslt->isPiPiPi0Pi0SelValid();
    if (m_debug) std::cout << "Mode : e+e- --> pi+pi-2pi0, N(pi): " << m_n_pi << ", N(photon): " << trkslt->igam().size() << ", chi2(vf): " << m_chi2_vf << ", chi2(kf): " << m_chi2_kf << std::endl;
    if (m_debug) {
        if (m_valid_PiPiPi0Pi0) std::cout << "----> Status: successful!" << std::endl;
        else std::cout << "----> Status: failed!" << std::endl;
    }
    m_tuple1->write();

    // mode 3: e+e- --> 3(pi+pi-)
    m_runNo = eventHeader->runNumber();
    m_evtNo = eventHeader->eventNumber();
    m_beamE = beamE;
    m_mode = 3;
    evtslt->PiPiPiPiPiPiSel(p4Lab, trkslt->igood());
    i_pi = evtslt->I_pi();
    change_pi = evtslt->Charge_pi();
    p_pi = evtslt->P_pi();
    cos_pi = evtslt->Cos_pi();
    phi_pi = evtslt->Phi_pi();
    chi_pi = evtslt->Chi_pi();
    prob_pi = evtslt->Prob_pi();
    probPH_pi = evtslt->ProbPH_pi();
    goodHits_pi = evtslt->GoodHits_pi();
    dt_pi = evtslt->Dt_pi();
    beta_pi = evtslt->Beta_pi();
    m_n_pi = i_pi.size();
    for (int i = 0; i < m_n_pi; i++) {
        for (int j = 0; j < 5; j++) m_trk_pi[i][j] = prob_pi[i][j];
        for (int j = 0; j < 5; j++) m_trk_pi[i][5 + j] = chi_pi[i][j];
        m_trk_pi[i][10] = change_pi[i];
        m_trk_pi[i][11] = p_pi[i];
        m_trk_pi[i][12] = cos_pi[i];
        m_trk_pi[i][13] = phi_pi[i];
        m_trk_pi[i][14] = probPH_pi[i];
        m_trk_pi[i][15] = goodHits_pi[i];
        for (int j = 0; j < 5; j++) m_trk_pi[i][16 + j] = dt_pi[i][j];
        m_trk_pi[i][21] = beta_pi[i];
    }
    m_chi2_vf = evtslt->Chi2_vf();
    m_chi2_kf = evtslt->Chi2_kf();
    m_valid_PiPiPiPiPiPi = evtslt->isPiPiPiPiPiPiSelValid();
    if (trkslt->igam().size() != 0) m_valid_PiPiPiPiPiPi = false;
    if (m_debug) std::cout << "Mode : e+e- --> 3(pi+pi-), N(pi): " << m_n_pi << ", N(photon): " << trkslt->igam().size() << ", chi2(vf): " << m_chi2_vf << ", chi2(kf): " << m_chi2_kf << std::endl;
    if (m_debug) {
        if (m_valid_PiPiPiPiPiPi) std::cout << "----> Status: successful!" << std::endl;
        else std::cout << "----> Status: failed!" << std::endl;
    }
    m_tuple1->write();

    // mode 4: e+e- --> 2(pi+pi-)2pi0
    m_runNo = eventHeader->runNumber();
    m_evtNo = eventHeader->eventNumber();
    m_beamE = beamE;
    m_mode = 4;
    evtslt->PiPiPiPiPi0Pi0Sel(p4Lab, trkslt->igood(), trkslt->igam());
    i_pi = evtslt->I_pi();
    change_pi = evtslt->Charge_pi();
    p_pi = evtslt->P_pi();
    cos_pi = evtslt->Cos_pi();
    phi_pi = evtslt->Phi_pi();
    chi_pi = evtslt->Chi_pi();
    prob_pi = evtslt->Prob_pi();
    probPH_pi = evtslt->ProbPH_pi();
    goodHits_pi = evtslt->GoodHits_pi();
    dt_pi = evtslt->Dt_pi();
    beta_pi = evtslt->Beta_pi();
    m_n_pi = i_pi.size();
    for (int i = 0; i < m_n_pi; i++) {
        for (int j = 0; j < 5; j++) m_trk_pi[i][j] = prob_pi[i][j];
        for (int j = 0; j < 5; j++) m_trk_pi[i][5 + j] = chi_pi[i][j];
        m_trk_pi[i][10] = change_pi[i];
        m_trk_pi[i][11] = p_pi[i];
        m_trk_pi[i][12] = cos_pi[i];
        m_trk_pi[i][13] = phi_pi[i];
        m_trk_pi[i][14] = probPH_pi[i];
        m_trk_pi[i][15] = goodHits_pi[i];
        for (int j = 0; j < 5; j++) m_trk_pi[i][16 + j] = dt_pi[i][j];
        m_trk_pi[i][21] = beta_pi[i];
    }
    m_chi2_vf = evtslt->Chi2_vf();
    m_chi2_kf = evtslt->Chi2_kf();
    m_valid_PiPiPiPiPi0Pi0 = evtslt->isPiPiPiPiPi0Pi0SelValid();
    if (m_debug) std::cout << "Mode : e+e- --> 2(pi+pi-)2pi0, N(pi): " << m_n_pi << ", N(photon): " << trkslt->igam().size() << ", chi2(vf): " << m_chi2_vf << ", chi2(kf): " << m_chi2_kf << std::endl;
    if (m_debug) {
        if (m_valid_PiPiPiPiPi0Pi0) std::cout << "----> Status: successful!" << std::endl;
        else std::cout << "----> Status: failed!" << std::endl;
    }
    m_tuple1->write();

    // mode 5: e+e- --> 2(K+K-)
    m_runNo = eventHeader->runNumber();
    m_evtNo = eventHeader->eventNumber();
    m_beamE = beamE;
    m_mode = 5;
    evtslt->KKKKSel(p4Lab, trkslt->igood());
    i_K = evtslt->I_K();
    change_K = evtslt->Charge_K();
    p_K = evtslt->P_K();
    cos_K = evtslt->Cos_K();
    phi_K = evtslt->Phi_K();
    chi_K = evtslt->Chi_K();
    prob_K = evtslt->Prob_K();
    probPH_K = evtslt->ProbPH_K();
    goodHits_K = evtslt->GoodHits_K();
    dt_K = evtslt->Dt_K();
    beta_K = evtslt->Beta_K();
    m_n_K = i_K.size();
    for (int i = 0; i < m_n_K; i++) {
        for (int j = 0; j < 5; j++) m_trk_K[i][j] = prob_K[i][j];
        for (int j = 0; j < 5; j++) m_trk_K[i][5 + j] = chi_K[i][j];
        m_trk_K[i][10] = change_K[i];
        m_trk_K[i][11] = p_K[i];
        m_trk_K[i][12] = cos_K[i];
        m_trk_K[i][13] = phi_K[i];
        m_trk_K[i][14] = probPH_K[i];
        m_trk_K[i][15] = goodHits_K[i];
        for (int j = 0; j < 5; j++) m_trk_K[i][16 + j] = dt_K[i][j];
        m_trk_K[i][21] = beta_K[i];
    }
    m_chi2_vf = evtslt->Chi2_vf();
    m_chi2_kf = evtslt->Chi2_kf();
    m_valid_KKKK = evtslt->isKKKKSelValid();
    if (trkslt->igam().size() != 0) m_valid_KKKK = false;
    if (m_debug) std::cout << "Mode : e+e- --> 2(K+K-), N(K): " << m_n_K << ", N(photon): " << trkslt->igam().size() << ", chi2(vf): " << m_chi2_vf << ", chi2(kf): " << m_chi2_kf << std::endl;
    if (m_debug) {
        if (m_valid_KKKK) std::cout << "----> Status: successful!" << std::endl;
        else std::cout << "----> Status: failed!" << std::endl;
    }
    m_tuple1->write();

    // mode 6: e+e- --> K+K-2pi0
    m_runNo = eventHeader->runNumber();
    m_evtNo = eventHeader->eventNumber();
    m_beamE = beamE;
    m_mode = 6;
    evtslt->KKPi0Pi0Sel(p4Lab, trkslt->igood(), trkslt->igam());
    i_K = evtslt->I_K();
    change_K = evtslt->Charge_K();
    p_K = evtslt->P_K();
    cos_K = evtslt->Cos_K();
    phi_K = evtslt->Phi_K();
    chi_K = evtslt->Chi_K();
    prob_K = evtslt->Prob_K();
    probPH_K = evtslt->ProbPH_K();
    goodHits_K = evtslt->GoodHits_K();
    dt_K = evtslt->Dt_K();
    beta_K = evtslt->Beta_K();
    m_n_K = i_K.size();
    for (int i = 0; i < m_n_K; i++) {
        for (int j = 0; j < 5; j++) m_trk_K[i][j] = prob_K[i][j];
        for (int j = 0; j < 5; j++) m_trk_K[i][5 + j] = chi_K[i][j];
        m_trk_K[i][10] = change_K[i];
        m_trk_K[i][11] = p_K[i];
        m_trk_K[i][12] = cos_K[i];
        m_trk_K[i][13] = phi_K[i];
        m_trk_K[i][14] = probPH_K[i];
        m_trk_K[i][15] = goodHits_K[i];
        for (int j = 0; j < 5; j++) m_trk_K[i][16 + j] = dt_K[i][j];
        m_trk_K[i][21] = beta_K[i];
    }
    m_chi2_vf = evtslt->Chi2_vf();
    m_chi2_kf = evtslt->Chi2_kf();
    m_valid_KKPi0Pi0 = evtslt->isKKPi0Pi0SelValid();
    if (m_debug) std::cout << "Mode : e+e- --> K+K-2pi0, N(K): " << m_n_K << ", N(photon): " << trkslt->igam().size() << ", chi2(vf): " << m_chi2_vf << ", chi2(kf): " << m_chi2_kf << std::endl;
    if (m_debug) {
        if (m_valid_KKPi0Pi0) std::cout << "----> Status: successful!" << std::endl;
        else std::cout << "----> Status: failed!" << std::endl;
    }
    m_tuple1->write();

    // mode 7: e+e- --> 3(K+K-)
    m_runNo = eventHeader->runNumber();
    m_evtNo = eventHeader->eventNumber();
    m_beamE = beamE;
    m_mode = 7;
    evtslt->KKKKKKSel(p4Lab, trkslt->igood());
    i_K = evtslt->I_K();
    change_K = evtslt->Charge_K();
    p_K = evtslt->P_K();
    cos_K = evtslt->Cos_K();
    phi_K = evtslt->Phi_K();
    chi_K = evtslt->Chi_K();
    prob_K = evtslt->Prob_K();
    probPH_K = evtslt->ProbPH_K();
    goodHits_K = evtslt->GoodHits_K();
    dt_K = evtslt->Dt_K();
    beta_K = evtslt->Beta_K();
    m_n_K = i_K.size();
    for (int i = 0; i < m_n_K; i++) {
        for (int j = 0; j < 5; j++) m_trk_K[i][j] = prob_K[i][j];
        for (int j = 0; j < 5; j++) m_trk_K[i][5 + j] = chi_K[i][j];
        m_trk_K[i][10] = change_K[i];
        m_trk_K[i][11] = p_K[i];
        m_trk_K[i][12] = cos_K[i];
        m_trk_K[i][13] = phi_K[i];
        m_trk_K[i][14] = probPH_K[i];
        m_trk_K[i][15] = goodHits_K[i];
        for (int j = 0; j < 5; j++) m_trk_K[i][16 + j] = dt_K[i][j];
        m_trk_K[i][21] = beta_K[i];
    }
    m_chi2_vf = evtslt->Chi2_vf();
    m_chi2_kf = evtslt->Chi2_kf();
    m_valid_KKKKKK = evtslt->isKKKKKKSelValid();
    if (trkslt->igam().size() != 0) m_valid_KKKKKK = false;
    if (m_debug) std::cout << "Mode : e+e- --> 3(K+K-), N(K): " << m_n_K << ", N(photon): " << trkslt->igam().size() << ", chi2(vf): " << m_chi2_vf << ", chi2(kf): " << m_chi2_kf << std::endl;
    if (m_debug) {
        if (m_valid_KKKKKK) std::cout << "----> Status: successful!" << std::endl;
        else std::cout << "----> Status: failed!" << std::endl;
    }
    m_tuple1->write();

    // mode 8: e+e- --> 2(K+K-)2pi0
    m_runNo = eventHeader->runNumber();
    m_evtNo = eventHeader->eventNumber();
    m_beamE = beamE;
    m_mode = 8;
    evtslt->KKKKPi0Pi0Sel(p4Lab, trkslt->igood(), trkslt->igam());
    i_K = evtslt->I_K();
    change_K = evtslt->Charge_K();
    p_K = evtslt->P_K();
    cos_K = evtslt->Cos_K();
    phi_K = evtslt->Phi_K();
    chi_K = evtslt->Chi_K();
    prob_K = evtslt->Prob_K();
    probPH_K = evtslt->ProbPH_K();
    goodHits_K = evtslt->GoodHits_K();
    dt_K = evtslt->Dt_K();
    beta_K = evtslt->Beta_K();
    m_n_K = i_K.size();
    for (int i = 0; i < m_n_K; i++) {
        for (int j = 0; j < 5; j++) m_trk_K[i][j] = prob_K[i][j];
        for (int j = 0; j < 5; j++) m_trk_K[i][5 + j] = chi_K[i][j];
        m_trk_K[i][10] = change_K[i];
        m_trk_K[i][11] = p_K[i];
        m_trk_K[i][12] = cos_K[i];
        m_trk_K[i][13] = phi_K[i];
        m_trk_K[i][14] = probPH_K[i];
        m_trk_K[i][15] = goodHits_K[i];
        for (int j = 0; j < 5; j++) m_trk_K[i][16 + j] = dt_K[i][j];
        m_trk_K[i][21] = beta_K[i];
    }
    m_chi2_vf = evtslt->Chi2_vf();
    m_chi2_kf = evtslt->Chi2_kf();
    m_valid_KKKKPi0Pi0 = evtslt->isKKKKPi0Pi0SelValid();
    if (m_debug) std::cout << "Mode : e+e- --> 2(K+K-)2pi0, N(K): " << m_n_K << ", N(photon): " << trkslt->igam().size() << ", chi2(vf): " << m_chi2_vf << ", chi2(kf): " << m_chi2_kf << std::endl;
    if (m_debug) {
        if (m_valid_KKKKPi0Pi0) std::cout << "----> Status: successful!" << std::endl;
        else std::cout << "----> Status: failed!" << std::endl;
    }
    m_tuple1->write();

}

StatusCode TrackSelAlg::finalize() {
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << ">>>>>>> in finalize()" << endmsg;
    return StatusCode::SUCCESS;
}

double TrackSelAlg::ECMS(int runNo) {
    if (runNo >= 59163 && runNo <= 59573) return 4.12848;
    else if (runNo >= 59574 && runNo <= 59896) return 4.15744;
    else if (runNo >= 59902 && runNo <= 60363) return 4.28788;
    else if (runNo >= 60364 && runNo <= 60805) return 4.31205;
    else if (runNo >= 60808 && runNo <= 61242) return 4.33739;
    else if (runNo >= 61249 && runNo <= 61762) return 4.37737;
    else if (runNo >= 61763 && runNo <= 62285) return 4.39645;
    else if (runNo >= 62286 && runNo <= 62823) return 4.43624;
    else if (runNo >= 63075 && runNo <= 63515) return 4.6260;
    else if (runNo >= 63516 && runNo <= 63715) return 4.6400;
    else if (runNo >= 63718 && runNo <= 63852) return 4.6600;
    else if (runNo >= 63867 && runNo <= 63888) return 4.6800;
    else return 999.;
}
