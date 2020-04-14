#ifndef TOOLS_HEAD
#define TOOLS_HEAD

#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/SecondVertexFit.h"
#include "ParticleID/ParticleID.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "EvtPreSlt/EvtPreSlt.h"

using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<WTrackParameter> VWTrkPara;

extern VertexFit * vtxfit;
extern KinematicFit * kmfit;

// Vertex fit
double vertex_fit(VWTrkPara &vwtrkpara, VertexParameter &birth);
// Kinematic fit
double kinematic_fit(const VWTrkPara &vwtrkpara, const HepLorentzVector &p4);
// Kinematic fit with two pi0
double kinematic_fit_two_pi0(const VWTrkPara &vwtrkpara, const HepLorentzVector &p4, Vint iGam);

#endif
