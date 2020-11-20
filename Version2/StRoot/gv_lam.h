#ifndef gvlam_def
#define gvlam_def

/// C++ headers
#include <iostream>

/// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
//#include "TVector3.hh"
#include <assert.h>
#include <TTree.h>
#include <vector>
#include <TLorentzVector.h>
#include <deque>

/// PicoDst headers
#include "./StPicoDstReader.h"
#include "./StPicoDst.h"
#include "./StPicoEvent.h"
#include "./StPicoTrack.h"
#include "./StPicoBTofHit.h"
#include "./StPicoBTowHit.h"
//#include "./StPicoHelix.h"
#include "./StPicoPhysicalHelix.h"
#include "./StPicoBTofPidTraits.h"
#include "./StRefMultCorr/StRefMultCorr.h"
#include "./StRefMultCorr/CentralityMaker.h"
#include "./particle.h"
#include "./event.h"

/// Analysis headers
#include "./StV0Type.h"
#include "./StDcaService.h"
#include "../structs.h"

namespace gv_lam
{
    //Global Variables
    //Structure 1 - conditions for analysis
    extern struct condition curr;
    
    //Structure 3 - Event Cuts
    extern struct event_cuts cut_eve;
    
    //Structure 4 - Constants related to V0
    extern struct v0_const constants_v0;
    
    //Structure 5 - Track Cuts
    extern struct track_cuts cut_trk;
    
    //Structure 7 - daughter vectors to find V0
    extern struct v0_dau_vectors *mDau1;
    extern struct v0_dau_vectors *mDau2;
    
    extern struct TRK track1, track2;
    
    extern struct lambda_info linfo;
    
    extern StPicoTrack *ptrack;
    extern StPicoTrack *gtrack;
    
    extern StPicoTrack *track_tmp;
    
    extern struct crucial_v0calc curr_trk;
    
    extern StRefMultCorr *refmultCorrUtil;
    
    extern int trk1_ids[5000];
    extern int trk2_ids[5000];
    
    extern std::deque<v0_dau_vectors> buffer_protons[40];
    extern struct v0_dau_vectors buffer_protons_tmp;
}

#endif
