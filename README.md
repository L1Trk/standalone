# standalone
not a fork of CMSSW

Recipes coming soon



To run the hybrid code:

cmsRun L1TrackNtupleMaker_cfg.py

To change between hybrid and tracklet recompile FPGAConstants.hh with the following lines changed:

#define USEHYBRID
static bool doKF=true; //true if using KF (and USEHYBRID uncommented)

To run TMTT:

cmsRun tmtt_tf_analysis_cfg.py
