#ifndef LambdaMaker_def
#define LambdaMaker_def

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
bool mix_events;
int vz_index;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Initialize Function
void init_Lambda();
//Make Function
bool make_Lambda(int cen);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Supporting Functions for Initialize
void init_macro();//setting condition for analysis
void init_eventcuts();
void init_v0const();
void init_trkcuts();
void init_hist();
//void init_tree();
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////\
//Supporting Functions for Make
bool event_cut_passed(int cen);//function to determine whether or not they passed the event cuts, and also filling QA plots
void set_centralitydependent_cuts(int cen);
void store_eventDetails();//function to store event details, and QA plots related to that
void init_dau_vectors();

////////////////// Beginning of V0's //////////////////
void v0_selection();
///// Global Tracks //////
bool v0_dau_selection();
void v0_dau_selection_QAPlots();
void v0_dau_selection_sigma();
bool v0_dau_selection_trkcuts();
void event_plane_storage_1();
bool v0_dau_selection_finalized();
bool filling_mDau1();
bool filling_mDau2();
///// Lambda Reconstruction //////
void reconstruct_V0();
void init_track1();
bool getinfofortracks1(int i);
void init_track2();
bool getinfofortracks2(int i);
void calc_Dcatrks();
void not_DcaAlgoLong();
bool cuts_relationsoftwotracks();
bool cut_v0mass();
bool cut_v0rdotp();
bool cut_v0decaylen();
void fill_finalplots_v0_dau();
//void calc_p_p1(TVector3 p_p, TVector3 xv0);
//void calc_p_p2(TVector3 p_p, TVector3 xv0);
bool final_data_organization(int i, int j);
void fill_tree(int i, int j);
void fill_hist(int temp_index_cent);
int determine_pt_bins(float lambda_pt_val);

//void clearing_everything();

///// Mix Events //////
void mix_events_check();
void mix_events_fillbuffer();

#endif
