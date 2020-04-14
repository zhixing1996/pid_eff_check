#include "TrkPreSlt/TrkPreSlt.h"

HepLorentzVector TrkPreSlt::GetP4OfRun(int runNo) {
    // first round @ 3.686
    if (fabs(runNo) >= 8093 && fabs(runNo) <= 9025)   { m_p4Lab[0] = 0.0393; m_p4Lab[1] = -0.001; m_p4Lab[2] = 0.0043; m_p4Lab[3] = 3.68631; }
    // first round @ 3.650 
    if (fabs(runNo) >= 9613 && fabs(runNo) <=9779)    { m_p4Lab[0] = 0.040; m_p4Lab[1] = -0.001; m_p4Lab[2] = 0.004; m_p4Lab[3] = 3.6502; }    
    //first round @ 3.097
    if (fabs(runNo) >= 9819 && fabs(runNo) < 10878)   { m_p4Lab[0] = 0.034067; m_p4Lab[1] = -0.000; m_p4Lab[2] = 0.000; m_p4Lab[3] = 3.097; }
    //first round @ 3.773 
    if (fabs(runNo) >= 11414 && fabs(runNo) <= 14604) { m_p4Lab[0] = 0.04249; m_p4Lab[1] = 0.0007; m_p4Lab[2] = 0.001; m_p4Lab[3] = 3.7730; }
    //second round @ 3.097
    if (fabs(runNo) >= 27255 && fabs(runNo) < 28293)  { m_p4Lab[0] = 0.034067; m_p4Lab[1] = -0.000; m_p4Lab[2] = 0.000; m_p4Lab[3] = 3.097; }
    // second round @ 3.650 
    if (fabs(runNo) >= 33725 && fabs(runNo) <= 33772) { m_p4Lab[0] = 0.040; m_p4Lab[1] = -0.001; m_p4Lab[2] = 0.004; m_p4Lab[3] = 3.6502; }    
    if (fabs(runNo) >= 35227 && fabs(runNo) <= 36213) { m_p4Lab[0] = 0.044; m_p4Lab[1] = -0.001; m_p4Lab[2] = 0.004; m_p4Lab[3] = 4.600; }    
    return m_p4Lab;
}
