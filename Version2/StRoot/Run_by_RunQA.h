#ifndef Run_byRunQA_def
#define Run_byRunQA_def

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
#include <TObject.h>
#include <deque>
#include "TProfile.h"

/// PicoDst headers
#include "./StPicoDstReader.h"
#include "./StPicoDst.h"
#include "./StPicoEvent.h"
#include "./StPicoTrack.h"
#include "./StPicoBTofHit.h"
#include "./StPicoBTowHit.h"
#include "./StPicoHelix.h"
#include "./StPicoPhysicalHelix.h"
#include "./StPicoBTofPidTraits.h"
#include "./StRefMultCorr/StRefMultCorr.h"
#include "./StRefMultCorr/CentralityMaker.h"
//#include "./StRefMultCorr.h"
//#include "./CentralityMaker.h"
#include "./particle.h"
#include "./event.h"

/// Analysis headers
#include "./StV0Type.h"
#include "./StDcaService.h"
#include "../structs.h"
#include "./gv.h"
#include "./gv_lam.h"

// Global Variables
Int_t Run_num = 0;
double pt_event = 0., eta_event = 0., phi_event = 0., dca_event = 0., nHitsFit_event = 0., nHitsDedx_event = 0.;
int number_tracks_counted = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Initialize Function
//void init_Lambda();
//Make Function
bool Run_by_Run();
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Supporting Functions for Initialize
//void init_macro();//setting condition for analysis
//void init_eventcuts();
//void init_v0const();
//void init_trkcuts();
void init_hist_run_by_run();
//void init_tree();
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Supporting Functions for Make
bool event_cut_passed_run_by_run();//function to determine whether or not they passed the event cuts, and also filling QA plots
//void set_centralitydependent_cuts(int cen);
//void store_eventDetails();//function to store event details, and QA plots related to that
//void init_dau_vectors();

////////////////// Beginning of V0's //////////////////
void v0_selection_run_by_run();
///// Primary Tracks //////
//bool v0_selection_primary();
//bool primary_cuts_decision();
//void event_plane_storage(); //selection = 0 (primary) 1 (daughter)
///// Global Tracks //////
bool v0_dau_selection_run_by_run();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Initialize Histograms Here
TProfile *mean_pt;
TProfile *mean_eta;
TProfile *mean_phi;
TProfile *mean_dca;
TProfile *mean_nHitsFit;
TProfile *mean_nHitsDedx;
TProfile *vz_dist;
TProfile *vxy_dist;
TProfile *refmult_dist;

#endif
