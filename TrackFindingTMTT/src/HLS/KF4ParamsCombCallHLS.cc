#include "L1Trigger/TrackFindingTMTT/interface/HLS/KF4ParamsCombCallHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanUpdateHLS.h"

#include "L1Trigger/TrackFindingTMTT/interface/TP.h"
#include "L1Trigger/TrackFindingTMTT/interface/StubCluster.h"
#include "L1Trigger/TrackFindingTMTT/interface/kalmanState.h"

#include "FWCore/Utilities/interface/Exception.h"

#include <vector>
#include <iostream>

namespace TMTT {

//=== Update KF helix params with this stub.
//=== (Override KF state updator in L1KalmanComb with version suitable for HLS).

const kalmanState* KF4ParamsCombCallHLS::kalmanUpdate( unsigned skipped, unsigned layer, const StubCluster *stubCluster, const kalmanState &stateIn, const TP *tpa ) {

  //  cout.setf(ios::scientific, ios::floatfield); // Get useful debug printout ...
  cout.unsetf(ios::floatfield); // Get useful debug printout ...
  cout.precision(8);
  
  // Get digitisation multipliers.
  rMult_   = pow(2, getSettings()->rtBits() )   / (getSettings()->rtRange());
  zMult_   = pow(2, getSettings()->zBits() )    / (getSettings()->zRange()); 
  phiMult_ = pow(2, getSettings()->phiSBits() ) / (getSettings()->phiSRange()); 
  // Multiplier of (phiMult/rMult) for helix param "inv2R" simplifies the KF maths, as explained in 
  //https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/demonstrator/specifications/demonstrator2_formats_working_doc.docx ,
  inv2R_Mult_ = (phiMult_/rMult_);
  // Reference radius in r-phi plane.
  chosenRofPhi_ = getSettings()->chosenRofPhi();
  // Number of eta sectors.
  numEtaRegions_ = getSettings()->numEtaRegions();
  // Get digitised stub info
  KalmanHLS::StubHLS stubDigi = this->getDigiStub(stubCluster, &stateIn);

  // Get digitised KF state info
  KalmanHLS::KFstateHLS stateInDigi = this->getDigiStateIn(skipped, layer, &stateIn);

  // Call HLS code to add stub to helix state.
  KalmanHLS::KFstateHLS stateOutDigi;
  KalmanHLS::ExtraOutHLS extraOut;
  KalmanHLS::KalmanUpdateHLS(stubDigi, stateInDigi, stateOutDigi, extraOut);

  // Store the extra info provided by the HLS updator about whether the state passes cuts.
  extraOut_ = extraOut;

  // Convert digitized ourput KF state to floating point.
  const kalmanState* newState = this->getStateOut(&stateIn, stubCluster, stateOutDigi);

  return newState;
}

//=== Get digital stub that the KF VHDL injects into the KF state updater (Maxeller/HLS)

KalmanHLS::StubHLS KF4ParamsCombCallHLS::getDigiStub(const StubCluster* stubCluster, const kalmanState* state) {
  // Get digitised stub(s) making up stub cluster.
  const vector<const Stub*> stubs = stubCluster->stubs();
  if (stubs.size() != 1) throw cms::Exception("KF4ParamsCombCallHLS: Can't cope with StubCluster that doesn't contain a single stub")<<stubs.size()<<endl;

  const DigitalStub& digiStub = stubs[0]->digitalStub();
 
  KalmanHLS::StubHLS stubDigi;
  // KF uses stub r, not rT. 
  stubDigi.r    = digiStub.iDigi_Rt() + std::round(rMult_*chosenRofPhi_); 
  // In octant format, KF VHDL z digitisation multiplier is same as r multiplier, with BSZ1 = BSZ bits
  // assigned to integer part of z in both VHDL & HLS. The following statement assigns stubDigi.z to
  // digiStub.iDigi_Z_KF().
  // In nonant format, KF VHDL z multiplier is factor 2 smaller than r multiplier, with BSZ1 = BSZ - 1
  // bits assigned to integer part of z in VHDL, and BSZ bits assigned to it in HLS. Following statement
  // therefore doubles the effective z multiplier.
  enum {BSZ = KalmanHLS::BSZ, BSZ1 = KalmanHLS::BSZ1};
  stubDigi.z.range(BSZ1-1, 0) = AP_FIXED(BSZ1,BSZ1)(digiStub.iDigi_Z_KF()).range(BSZ1-1, 0); 
  stubDigi.phiS = digiStub.iDigi_PhiS();

#ifdef IRT_DEBUG
  if (state->candidate().getMatchedTP() != nullptr) {
    unsigned int iPhiSec = state->candidate().iPhiSec();
    float phiSec = 2.*M_PI * (0.5 + float(iPhiSec)) / float(getSettings()->numPhiSectors()) - M_PI;
    float phiStubOff = reco::deltaPhi(digiStub.phi(), digiStub.phiS());
    cout<<"KF4 sector phi check "<<phiSec<<" "<<phiStubOff<<endl;
    cout<<"KF4 input stub: float (r,phi) = "<<stubs[0]->r()<<" "<<reco::deltaPhi(stubs[0]->phi(),phiSec)<<endl;
    cout<<"KF4 input stub: digi (r,phi) = "<<double(stubDigi.r)/rMult_<<" "<<double(stubDigi.phiS)/phiMult_<<endl;
  }
#endif

#ifdef IRT_DEBUG
  cout<<"PS MODULE "<<stubs[0]->psModule()<<endl;;
#endif

  stubDigi.valid = true;

  return stubDigi;
}

//=== Get digitised KF state info that the KF VHDL injects into the KF state updater (Maxeller/HLS)

KalmanHLS::KFstateHLS KF4ParamsCombCallHLS::getDigiStateIn(unsigned int skipped, unsigned int layer, const kalmanState* state) {
  // Calculate factors to convert floating point helix params to digitized ones.
  // Based on constants & functions named *HWU* in
//https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/cgn/src/formats/Constants.maxj .

  KalmanHLS::KFstateHLS stateDigi;

  // Cell HT found track in, with (mbin,cbin) centred on zero.
  pair<unsigned int, unsigned int> htCell = state->candidate().getCellLocationHT();
  // Digitized HT cell info must be centred on zero. (See L1fittedTrack::digitalTrack()).
  stateDigi.mBin = htCell.first   - floor(getSettings()->houghNbinsPt()/2); 
  stateDigi.cBin = htCell.second  - floor(getSettings()->houghNbinsPhi()/2);

  // Fitted helix parameters.
  vector<double> helixParams = state->xa();
  double inv2R = helixParams[0]; // Half inverse radius of curvature.
  double phi0  = helixParams[1]; // Measured with respect to centre of phi sector.
  double tanL  = helixParams[2];
  double z0    = helixParams[3];
  TMatrixD cov = state->pxxa();
  double cov_inv2R_inv2R = cov[0][0];
  double cov_phi0_phi0   = cov[1][1];
  double cov_tanL_tanL   = cov[2][2];
  double cov_z0_z0       = cov[3][3];
  double cov_inv2R_phi0  = cov[0][1];
  double cov_tanL_z0     = cov[2][3];  

  // Digitize helix parameters, using multiplication factors in 
  // Demonstrator Units" section of https://twiki.cern.ch/twiki/pub/CMS/CMSUKTrkTrig/KFequations.pdf .

  // Using the specified multipliers, which are related to those used for stubs (see getDigiStub()) simplifies the maths.
  // The helix params need better granularity than the stubs, which is accomodated by using additional bits after the 
  // decimal point in the digitized numbers (profitting from Maxeller/HLS ability to handle floats).

  stateDigi.inv2R = inv2R*inv2R_Mult_; // See inv2RToHWU() in Maxeller code in above web link.
  stateDigi.phi0  = phi0*phiMult_;

#ifdef IRT_DEBUG
  if (state->candidate().getMatchedTP() != nullptr) {
    unsigned int iPhiSec = state->candidate().iPhiSec();
    float phiSec = 2.*M_PI * (0.5 + float(iPhiSec)) / float(getSettings()->numPhiSectors()) - M_PI;
    cout<<"KF4 Input track (float): q/pt = "<<inv2R/(0.5*getSettings()->invPtToInvR())<<" phi0 = "<<phi0<<endl;
    cout<<"KF4 Input track (digi): q/pt = "<<double(stateDigi.inv2R)/inv2R_Mult_/(0.5*getSettings()->invPtToInvR())<<" phi0 = "<<double(stateDigi.phi0)/phiMult_<<endl;
    cout<<"          truth: q/pt = "<<state->candidate().getMatchedTP()->qOverPt()<<" phi0 = "<<reco::deltaPhi(state->candidate().getMatchedTP()->phi0(),phiSec)<<endl;
  }
#endif

  stateDigi.tanL  = tanL;              // Multiplier is 1 for tanL.
  // Multiplier of rMult instead of zMult here simplifies maths, like for stub z in getDigiStub().
  stateDigi.z0    = z0*rMult_;         // See rMmToHWU() in Maxeller code.

  // Check digitisation range of helix parameters is sufficient.
  KalmanHLS::CHECK_AP::checkCalc("helix0", stateDigi.inv2R, inv2R*inv2R_Mult_, 9.9e9, 0.0001);
  KalmanHLS::CHECK_AP::checkCalc("helix1", stateDigi.phi0 , phi0*phiMult_    , 9.9e9, 0.001);
  KalmanHLS::CHECK_AP::checkCalc("helix2", stateDigi.tanL , tanL             , 9.9e9, 0.001);
  KalmanHLS::CHECK_AP::checkCalc("helix3", stateDigi.z0   , z0*rMult_        , 9.9e9, 0.1);

  // Multipliers for covariance matrix follow from those of helix params.
  stateDigi.cov_00 = cov_inv2R_inv2R * inv2R_Mult_ * inv2R_Mult_;
  stateDigi.cov_11 = cov_phi0_phi0 * phiMult_ * phiMult_;
  stateDigi.cov_22 = cov_tanL_tanL;
  stateDigi.cov_33 = cov_z0_z0 * rMult_ * rMult_;
  stateDigi.cov_01 = cov_inv2R_phi0 * rMult_ * phiMult_;
  stateDigi.cov_23 = cov_tanL_z0 * rMult_;

  // Check digitisation range of covariance matrix is sufficient.
  KalmanHLS::CHECK_AP::checkCalc("C00_old", stateDigi.cov_00, cov_inv2R_inv2R * inv2R_Mult_ * inv2R_Mult_);
  KalmanHLS::CHECK_AP::checkCalc("C11_old", stateDigi.cov_11, cov_phi0_phi0   * phiMult_    * phiMult_);
  KalmanHLS::CHECK_AP::checkCalc("C22_old", stateDigi.cov_22, cov_tanL_tanL);
  KalmanHLS::CHECK_AP::checkCalc("C33_old", stateDigi.cov_33, cov_z0_z0       * rMult_      * rMult_);
  KalmanHLS::CHECK_AP::checkCalc("C01_old", stateDigi.cov_01, cov_inv2R_phi0  * rMult_      * phiMult_);
  KalmanHLS::CHECK_AP::checkCalc("C23_old", stateDigi.cov_23, cov_tanL_z0     * rMult_);

  stateDigi.chiSquared = state->chi2();

  // This is the KF layer that we are currently looking for stubs in, incremented by L1KalmanComb::doKF(), which in any eta region increases from 0-7 as a particle goes through each layer in turn.
  stateDigi.layerID        = layer;
  // This is the number of skipped layers assuming we find a stub in the layer currently being searched.
  stateDigi.nSkippedLayers = skipped;

  stateDigi.candidateID = 0; // Not used by KF updator.
  stateDigi.eventID     = 0; // Not used by KF updator.

  unsigned int iEtaReg = state->candidate().iEtaReg(); // Although this comes from the state, it is actually the eta region of the stub.
  // This is encoded in tortuous way copied from Maxeller code (lines 127-133).
//https://svnweb.cern.ch/cern/wsvn/UK-TrackTrig/firmware/trunk/cactusupgrades/projects/tracktrigger/kalmanfit/firmware/hdl/packages/stubs.vhd
  if (iEtaReg < numEtaRegions_/2) {
    stateDigi.etaSectorID = (numEtaRegions_/2 - 1) - iEtaReg; // Count eta regions increasing away from theta = 90 degrees.
    stateDigi.etaSectorZsign = true; // Eta sector in -ve z half of Tracker.
  } else {
    stateDigi.etaSectorID = iEtaReg - numEtaRegions_/2;
    stateDigi.etaSectorZsign = false;
  }

  stateDigi.valid = true;

  return stateDigi;
}

//=== Convert digitized ourput KF state to floating point.

const kalmanState* KF4ParamsCombCallHLS::getStateOut(const kalmanState* stateIn, const StubCluster* stubCluster, const KalmanHLS::KFstateHLS& stateOutDigi) {
  // Convert digitized helix state to floating point one.
  // Also copy some info directly from input floating point to output floating point state, if unchanged.

  // Fill arguments of L1KalmanComb::mkState(), which is called to make a kalmanState object.
  const L1track3D& candidate = stateIn->candidate();
  unsigned int n_skipped     = stateOutDigi.nSkippedLayers;
  unsigned int kLayer_next   = stateOutDigi.layerID; // Unchanged by KF updator.
  unsigned int layerId       = stateIn->layerId();
  const kalmanState* last_state = stateIn; 
  // Factors to convert digitized helix params to floating ones are inverse of those in getDigiStateIn().
  vector<double> x(4); // helix params
  x[0] = double(stateOutDigi.inv2R) / inv2R_Mult_;
  x[1] = double(stateOutDigi.phi0) / phiMult_;
  x[2] = double(stateOutDigi.tanL);
  x[3] = double(stateOutDigi.z0) / rMult_;
  TMatrixD pxx(4,4); // helix covariance matrix
  pxx[0][0] = double(stateOutDigi.cov_00) / (inv2R_Mult_ * inv2R_Mult_);
  pxx[1][1] = double(stateOutDigi.cov_11) / (phiMult_ * phiMult_);
  pxx[2][2] = double(stateOutDigi.cov_22);
  pxx[3][3] = double(stateOutDigi.cov_33) / (rMult_ * rMult_);
  pxx[0][1] = double(stateOutDigi.cov_01) / (rMult_ * phiMult_);
  pxx[1][0] = pxx[0][1];
  pxx[2][3] = double(stateOutDigi.cov_23) / (rMult_);
  pxx[3][2] = pxx[2][3];
  TMatrixD K(4,2); // KF gain matrix - don't provide, as can't be used?
  TMatrixD dcov(2,2); // Stub (phi,z) position covariance matrix - don't provide as can't be used?
  const StubCluster* stubcl = stubCluster;
  double chi2 = stateOutDigi.chiSquared;

  const kalmanState* ks = this->mkState(candidate, n_skipped, kLayer_next, layerId, last_state,
				        x, pxx, K, dcov, stubcl, chi2);

#ifdef IRT_DEBUG
  if (ks->candidate().getMatchedTP() != nullptr) {
    unsigned int iPhiSec = ks->candidate().iPhiSec();
    float phiSec = 2.*M_PI * (0.5 + float(iPhiSec)) / float(getSettings()->numPhiSectors()) - M_PI;
    double inv2R = x[0]; // Half inverse radius of curvature.
    double phi0  = x[1]; // Measured with respect to centre of phi sector.
    double tanL  = x[2];
    double z0    = x[3];
    cout<<"KF4 Output track (digi): q/pt = "<<inv2R/inv2R_Mult_/(0.5*getSettings()->invPtToInvR())<<" phi0 = "<<phi0/phiMult_<<endl;
    cout<<"          truth: q/pt = "<<ks->candidate().getMatchedTP()->qOverPt()<<" phi0 = "<<reco::deltaPhi(ks->candidate().getMatchedTP()->phi0(),phiSec)<<endl;
  }
#endif

  return ks;
}

//=== This is identical to version in KF4ParamsComb, deciding if a state passes cuts,
//=== except that it also checks the cut decisions produced by the HLS KalmanUpdate.

bool KF4ParamsCombCallHLS::isGoodState( const kalmanState &state ) const
{

#ifdef IRT_DEBUG
  bool trueTP = (state.candidate().getMatchedTP() != nullptr);
  bool vTrueTP = (trueTP && state.candidate().getPurity() < 0.99 && state.candidate().pt() > 3);
  cout<<"GOOD CHECK "<<trueTP<<" "<<vTrueTP<<endl;
#endif

  // Check if this state passes cuts using standard C++ KF code.
  bool goodState = KF4ParamsComb::isGoodState(state);

  // Check if this state passes cuts using HLS KF xode.
  bool goodState_HLS = extraOut_.z0Cut && extraOut_.ptCut && extraOut_.chiSquaredCut && extraOut_.sufficientPScut;

  unsigned nStubLayers = state.nStubLayers();
  double pt=fabs( getSettings()->invPtToInvR() / (2*state.xa()[INV2R]) ); 
  double z0=fabs( state.xa()[Z0] ); 

  // Check if this state passes Algo50 duplicate removal cuts.
  // (Approximate since it negelects "rescue" step of Algo50).
  //bool algo50_HLS = true;
  if (nStubLayers > 3) {
    // algo50_HLS = (extraOut_.consistent && extraOut_.sectorCut);
    // Debug printout to debug Algo50. To be used with "Cheat" cfg option. 
    // And print fittedTrack.consistent() from TMTrackProducer.cc.
    //
    // cout<<"algo50: HT consistent="<<extraOut_.consistent<<" sector consistent="<<extraOut_.sectorCut
    //   <<" mbin = "<<(state.candidate().getCellLocationHT().first-floor(getSettings()->houghNbinsPt()/2))
    //   <<" vs "<<extraOut_.mBinHelix
    //   <<" cbin = "<<state.candidate().getCellLocationHT().second-floor(getSettings()->houghNbinsPhi()/2)
    //   <<" vs "<<extraOut_.cBinHelix<<endl;
  }

  // Check if HLS & C++ KF agree ...

  if (goodState && not goodState_HLS) {
    cout<<"ERROR: KF HLS incorrectly rejected state "<<nStubLayers<<" "<<extraOut_.z0Cut<<" "<<extraOut_.ptCut<<" "<<extraOut_.chiSquaredCut<<" "<<extraOut_.sufficientPScut<<endl;
    cout<<"chi2="<<state.chi2()<<" pt="<<pt<<" 1/2R="<<state.xa()[INV2R]<<endl;
  } else if (not goodState && goodState_HLS) {
    cout<<"ERROR: KF HLS incorrectly kept state "<<nStubLayers<<" "<<extraOut_.z0Cut<<" "<<extraOut_.ptCut<<" "<<extraOut_.chiSquaredCut<<" "<<extraOut_.sufficientPScut<<endl;
#ifdef IRT_DEBUG
  } else {
    cout<<"KF HLS CUTS GOOD "<<nStubLayers<<" "<<goodState_HLS<<" algo50 "<<algo50_HLS<<endl;
#endif
  }

#ifdef IRT_DEBUG
  if (trueTP) cout<<"L1K GOOD STATE CHECK: vtrue="<<vTrueTP<<" good="<<goodState_HLS<<" nLay="<<nStubLayers<<" chi2="<<state.chi2()<<" pt="<<pt<<" z0="<<state.xa()[Z0]<<" skip="<<state.nSkippedLayers()<<endl;
#endif

  //return goodState;
  return goodState_HLS;
  //return goodState_HLS && algo50_HLS; // includes approximate algo50 duplicate removal, without rescue step.
}

}
