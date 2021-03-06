#ifndef __HLSconstants__
#define __HLSconstants__

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSutilities.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/StubHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KFstateHLS.h"
#else
#include "HLSutilities.h"
#include "StubHLS.h"
#include "KFstateHLS.h"
#endif

#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

//--- Number of helix parameters for track fit.

static const unsigned int N_HELIX_PAR = 5;

//-- Copied from Maxeller code Constants.maxj

// Digitisation multipliers (from data format doc).
// KF uses same multiplier for r as for stubs in DTC, but one extra bit to accomodate larger range,
// since KF measures r w.r.t. beamline. And it uses r multiplier for z too.
static const float rMult = pow(2.,BSR-1)/103.1103;
static const float phiMult = pow(2.,BSP)/0.698131700;
static const float rphiMult = rMult*phiMult;
static const float inv2R_Mult = (phiMult/rMult);

// Beam spot length & reference radii w.r.t. beamline.
static const float beamSpotLength= 15.0;
static const float chosenRofPhi_flt = 61.273;
static const StubHLS::TR chosenRofPhi = chosenRofPhi_flt*rMult;
static const float chosenRofZ_flt = 50.0;
static const StubHLS::TR chosenRofZ = chosenRofZ_flt*rMult;

static const float bField = 3.8112;
static const float invPtToInvR = bField*(2.9979e8/1.0e11);
static const float minPt_HT = 3.; // Range of Hough transform
static const float invRmin_HT = invPtToInvR*(1./minPt_HT);

static const float kalmanMultScatTerm = 0.00075; // Same as cfg param of same name in CMSSW TMTT code.

// Phi sectors
static const float TWO_PI = 2*3.14159265;
static const int numPhiSectors = 18;
static const float phiSectorWidth = TWO_PI / numPhiSectors;

// Bit shift *_bitShift to calculate HT cell from digitized (phi, invR) of helix params.
// Chosen such that pow(2,+-shift) = (dcBin_digi, dmBin_digi) calculated below.
// (where sign diff is because in KalmanUpdate.cc, one is used to shift right & one to shift left).
// Adjust if you change bits assigned to stubs.
enum {phiToCbin_bitShift = 7, invRToMbin_bitShift = 4};
enum {BCH=BH1-phiToCbin_bitShift, BMH=BH0+invRToMbin_bitShift};

// Size of HT cells
static const int numPhiBins = 64; // Bins in HT
static const int numPtBins = 36;  // Nonants
//static const int numPtBins = 32;  // Octants
static const AP_INT(BCH) minPhiBin = -numPhiBins/2; // BCH & BMH should be larger than BHT_C & BHT_M to monitor overflow.
static const AP_INT(BCH) maxPhiBin =  numPhiBins/2 - 1;
static const AP_INT(BMH) minPtBin  = -numPtBins/2;
static const AP_INT(BMH) maxPtBin  =  numPtBins/2 - 1;

/*
static const float dcBin = numPhiBins / phiSectorWidth; 
static const float dmBin = numPtBins / (invRmin_HT); 
static const float dcBin_digi = dcBin/phiMult; // = 1 / pow(2,7)
static const float dmBin_digi = dmBin/inv2R_Mult; // = pow(2,2+BEX)
*/

// Eta sector boundaries in z at reference radius (assumed symmetric).
// (As this is complex, ROM initialization fails unless stored in a class ...)

class EtaBoundaries {
public:
  enum {nSec=9};

  EtaBoundaries() {
    static const float eta[nSec+1] = {0.0, 0.31, 0.61, 0.89, 1.16, 1.43, 1.7, 1.95, 2.16, 2.4};
    for (unsigned int i = 0; i <= nSec; i++) {
      float zAtRefR = chosenRofZ_flt/tan(2 * atan(exp(-eta[i])));
      z_[i] = rMult*zAtRefR;
    }
  }

public:
  StubHLS::TZ  z_[nSec+1];
};

// Also failed in VHDL
//static const EtaBoundaries etaBoundaries;

//--- Cuts to select acceptable fitted track states.
//--- (Array vs #stubs on track, where element 0 is never used).

// Pt or 1/2R cut.
static const float ptCut_flt_tight = 2.95; // Smaller than HT cut to allow for resolution during KF fit.
static const float ptCut_flt_loose = 2.90;
static const float inv2Rcut_flt_tight = 0.5*invPtToInvR*(1./ptCut_flt_tight);
static const float inv2Rcut_flt_loose = 0.5*invPtToInvR*(1./ptCut_flt_loose);
static const KFstateHLS<5>::TR inv2Rcut_tight = inv2R_Mult*inv2Rcut_flt_tight;
static const KFstateHLS<5>::TR inv2Rcut_loose = inv2R_Mult*inv2Rcut_flt_loose;
static const KFstateHLS<5>::TR inv2Rcut_none = (1 << (BH0 - 1)) - 1;
static const KFstateHLS<5>::TR inv2Rcut[] = {inv2Rcut_none, inv2Rcut_none, inv2Rcut_loose, inv2Rcut_loose, inv2Rcut_tight, inv2Rcut_tight, inv2Rcut_tight};

// z0 cut
static const KFstateHLS<5>::TZ z0Cut_tight = rMult*beamSpotLength; // r multiplier used for z in KF. 
static const KFstateHLS<5>::TZ z0Cut_none = (1 << (BH3-1)) - 1; 
static const KFstateHLS<5>::TZ z0Cut[] = {z0Cut_none, z0Cut_none,  z0Cut_tight, z0Cut_tight, z0Cut_tight, z0Cut_tight, z0Cut_tight}; 

// d0 cut
static const float d0Cut_flt_tight = 5.;
static const float d0Cut_flt_loose = 10.;
static const KFstateHLS<5>::TD d0Cut_tight = rphiMult*d0Cut_flt_tight;
static const KFstateHLS<5>::TD d0Cut_loose = rphiMult*d0Cut_flt_loose;
static const KFstateHLS<5>::TD d0Cut_none = (1 << (BH4-1)) - 1;
static const KFstateHLS<5>::TD d0Cut[] = {d0Cut_none, d0Cut_none, d0Cut_none, d0Cut_loose, d0Cut_tight, d0Cut_tight, d0Cut_tight};

// Chi2 cut (assuming digitisation multiplier of unity).
static const KFstateHLS<5>::TCHI chi2Cut_none = (1 << BCHI) - 1;
static const KFstateHLS<5>::TCHI chi2Cut[] = {chi2Cut_none, chi2Cut_none, 10, 30, 80, 120, 160}; 

#ifdef CMSSW_GIT_HASH
}
}
#endif

#endif
