///=== This is the base class for all the track fit algorithms

///=== Written by: Alexander D. Morton and Sioni Summers

#include "L1Trigger/TrackFindingTMTT/interface/TrackFitGeneric.h"
#include "L1Trigger/TrackFindingTMTT/interface/Settings.h"
#include "L1Trigger/TrackFindingTMTT/interface/Utility.h"
#include "L1Trigger/TrackFindingTMTT/interface/ChiSquared4ParamsApprox.h"
#include "L1Trigger/TrackFindingTMTT/interface/KF4ParamsComb.h"
#include "L1Trigger/TrackFindingTMTT/interface/KF5ParamsComb.h"
#include "L1Trigger/TrackFindingTMTT/interface/SimpleLR.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <map> 
#include <new>
 
namespace TMTT {

//=== Set configuration parameters.
 
TrackFitGeneric::TrackFitGeneric( const Settings* settings, const string &fitterName ) : settings_(settings), fitterName_(fitterName), nDupStubs_(0) {
}
 
 
//=== Fit a track candidate obtained from the Hough Transform.
//=== Specify which phi sector and eta region it is in.
 
L1fittedTrack TrackFitGeneric::fit(const L1track3D& l1track3D,  unsigned int iPhiSec, unsigned int iEtaReg) {
  return L1fittedTrack (settings_, l1track3D, l1track3D.getStubs(), 0, 0, 0, 0, 0, 999999., 0, iPhiSec, iEtaReg);
}
 
TrackFitGeneric* TrackFitGeneric::create(std::string fitter, const Settings* settings) {
    if (fitter.compare("ChiSquared4ParamsApprox")==0) {
	return new ChiSquared4ParamsApprox(settings, 4);
    } else if (fitter.compare("KF4ParamsComb")==0) {
	return new KF4ParamsComb(settings, 4, fitter );
    } else if (fitter.compare("KF5ParamsComb")==0) {
	return new KF5ParamsComb(settings, 5, fitter );
    } else if (fitter.compare("SimpleLR")==0) {
      return new SimpleLR(settings);
    } else {
      throw cms::Exception("TrackFitGeneric: ERROR you requested unknown track fitter")<<fitter<<endl;
    }
} 

}
