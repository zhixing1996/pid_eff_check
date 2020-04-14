#include "TrkPreSlt/TrkPreSlt.h"

Hep3Vector TrkPreSlt::GetXorigin(int runNo) {
    if(runNo>=8093&&runNo<=10878)  m_xorigin.set(0.164816,-0.257954,-0.330269);
    if(runNo==11414)               m_xorigin.set(0.866032,-0.231788,0.318504);
    if(runNo==11418)               m_xorigin.set(0.853652,-0.230285,0.154119);
    if(runNo==11420)               m_xorigin.set(0.853506,-0.229947,-0.0657398);
    if(runNo>=11428&&runNo<=12091) m_xorigin.set(0.853443,-0.228519,-0.0854736);
    if(runNo>=12092&&runNo<=12103) m_xorigin.set(0.891239,-0.228108,-0.17582);
    if(runNo>=12104&&runNo<=12145) m_xorigin.set(0.918063,-0.226754,-0.184516);
    if(runNo>=12146&&runNo<=12326) m_xorigin.set(0.937982,-0.227077,0.215204);
    if(runNo>=12327              ) m_xorigin.set(0.942276,-0.227446,0.232697);
    return m_xorigin;
}
