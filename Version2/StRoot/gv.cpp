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
#include "./gv.h"

namespace gv{
	
	//Structure 2 - particles
    particle *p_lambda = new particle[5000];
    particle *p_all = new particle[5000];
    particle *p_dau = new particle[5000];
    particle *p_proton = new particle[5000];

    //Structure 3 - Event details
    event *details = new event();

	TTree *corr_tree = new TTree("corr_tree", "Tree for Correlations");
    TTree *mXiTree = new TTree("XiPicoDst", "XiPicoDst from StXiMaker");

    int check = 1;
    
    struct histograms_created all_hists;
    struct histograms_xi all_hists_xi;
    
    StPicoEvent *picoEvent;
    StPicoDst *picoDst;
    
    StPicoTrack *track_forp = new StPicoTrack();
}
