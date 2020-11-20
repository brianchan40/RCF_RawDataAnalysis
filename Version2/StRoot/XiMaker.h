#ifndef XiMaker_def
#define XiMaker_def

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
#include "TVector3.h"
#include <assert.h>
#include <TTree.h>
#include <vector>
#include <TLorentzVector.h>
#include "TObject.h"

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
#include "./particle.h"
#include "./event.h"

/// Analysis headers
#include "./StV0Type.h"
#include "./StXiType.h"
#include "./StDcaService.h"
#include "../structs.h"
#include "./gv.h"

//// Initializing ////
void init_xi();
void init_macro_xi();
void init_const_xi();
void init_param_xi();
void init_hist_xi();
void init_tree_xi();
//// End of Initializing ////


//// Making /////
void make_Xi();
void init_dau_vectors_xi();
void storeEventDetails_xi();
////// Beginning of xi_selection() ////////
void xi_selection();
void xi_selection_QAPlots(StPicoTrack *track, TVector3 p, TVector3 origin);
void Bachelor_selection(StPicoTrack *track);
////// End of xi_selection() ////////
////// Beginning of xi_reconstruction() ////////
void xi_reconstruction();
void helix_bachelor();
bool calc_dca_xi(int j);
bool calc_ximass();
bool tof_record();
void fill_mXiDst(int nXi, int i);

#endif
