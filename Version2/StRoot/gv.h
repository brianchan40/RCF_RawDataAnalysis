#ifndef gv_def
#define gv_def

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
//#include "./StRefMultCorr.h"
//#include "./CentralityMaker.h"
#include "./StRefMultCorr/StRefMultCorr.h"
#include "./StRefMultCorr/CentralityMaker.h"
#include "./particle.h"
#include "./event.h"

/// Analysis headers
#include "./StV0Type.h"
#include "./StDcaService.h"
#include "../structs.h"

namespace gv
{
    extern StPicoEvent *picoEvent;
    extern StPicoDst *picoDst;

    //Structure 2 - particles
    extern particle *p_lambda;
    extern particle *p_proton;
    extern particle *p_all;
    extern particle *p_dau;

    //Structure 3 - Event details
    extern event *details;

    //Global Variables
    extern TTree *corr_tree;
    extern struct histograms_created all_hists;
    extern struct histograms_xi all_hists_xi;
    
    extern TTree *mXiTree;

    extern int check;
    
    extern StPicoTrack *track_forp;
}

#endif
