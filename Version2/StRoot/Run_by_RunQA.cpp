/*
 * Author: Brian Chan
 * Date: February 27, 2019
 *
 * Lambdamaker tries to reconstruct lambdas for further analysis
 *
 **/


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
#include "TProfile.h"

/// gv::picoDst headers
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
#include "./Run_by_RunQA.h"
#include "./gv.h"
#include "./gv_lam.h"

using namespace std;

void init_hist_run_by_run(){
    mean_pt = new TProfile("mean_pt", "Mean Pt by Run", 37980, 70, 38050, 0, 5, "");
    mean_eta = new TProfile("mean_eta", "Mean Eta by Run", 37980, 70, 38050, -1, 1, "");
    mean_phi = new TProfile("mean_phi", "Mean Phi by Run", 37980, 70, 38050, -1, 1, "");
    mean_dca = new TProfile("mean_dca", "Mean DCA by Run", 37980, 70, 38050, 0, 10, "");
    mean_nHitsFit = new TProfile("mean_nHitsFit", "Mean nHitsFit by Run", 37980, 70, 38050, 0, 50, "");
    mean_nHitsDedx = new TProfile("mean_nHitsDedx", "Mean nHitsDedx by Run", 37980, 70, 38050,0, 50, "");
    vz_dist = new TProfile("vz_dist", "Vz by Run", 37980, 70, 38050, -130, 130, "");
    vxy_dist = new TProfile("vxy_dist", "Vxy by Run", 37980, 70, 38050, 0, 3, "");
    refmult_dist = new TProfile("refmult_dist", "Refmult by Run", 37980, 70, 38050, 0, 800, "");

}

/// Make FUNCTION
bool Run_by_Run()
{
    /// Event Cuts & gv::details Stored
    if(!event_cut_passed_run_by_run()) return false;

    //*******************************************/
    /// Track analysis

    // Global & Primary Tracks Initialized
    v0_selection_run_by_run();

    if(gv_lam::curr.debug_reconstruct == true) std::cout << "50" << endl;

    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Supporting Functions for Make

//function to determine whether or not they passed the event cuts, and also filling QA plots
bool event_cut_passed_run_by_run()
{
    //Trigger Information
    if (!gv::picoEvent->isTrigger(610001) && !gv::picoEvent->isTrigger(610011) && !gv::picoEvent->isTrigger(610021) && !gv::picoEvent->isTrigger(610031) && !gv::picoEvent->isTrigger(610041) && !gv::picoEvent->isTrigger(610051)) return false;

    //Centrality definition BEGINS
    Run_num = gv::picoEvent->runId();
    Int_t Event_num = gv::picoEvent->eventId();
    Int_t refmult = gv::picoEvent->grefMult();
    int Ntofmatch = gv::picoEvent->nBTOFMatch();
    double zdcx = gv::picoEvent->ZDCx();
    
    // For refmult
    gv_lam::refmultCorrUtil  = CentralityMaker::instance()->getRefMultCorr() ;
    
    // You need to specify the run number you are going to process
    gv_lam::refmultCorrUtil->init(Run_num);

    // Call initEvent(const UShort_t RefMult, const Double_t z) function
    // event-by-event at the beginning before using any other functions
    gv_lam::refmultCorrUtil->initEvent(refmult,gv::picoEvent->primaryVertex().z(),zdcx);
    //if(gv_lam::refmultCorrUtil->isBadRun(Run_num)) return false;
    //if(!gv_lam::refmultCorrUtil->passnTofMatchRefmultCut(1.*refmult, 1.*Ntofmatch)) return false;

    //   see StRefMultCorr.h for the definition of centrality bins
    const Int_t cent9  = gv_lam::refmultCorrUtil->getCentralityBin9() ;
    
    //Duplicate Events? - Begins
    if ( gv::picoDst->numberOfTracks() == gv_lam::cut_eve.mTestNTrack && gv_lam::cut_eve.mEventsProcessed != 0 && gv_lam::cut_eve.mTestVZ != 0 &&  gv::picoEvent->primaryVertex().z() == gv_lam::cut_eve.mTestVZ )
    {
        std::cout << gv_lam::cut_eve.mEventsProcessed << " " << "seems a duplicated event!" << std::endl;
        return false;
    }
    gv_lam::cut_eve.mTestVZ = gv::picoEvent->primaryVertex().z();
    gv_lam::cut_eve.mTestNTrack = gv::picoDst->numberOfTracks();
    //Duplicate Events? - Ends

    double vpdVz = -990.0;
    vpdVz = gv::picoEvent->vzVpd();
    
    gv::details->Magn = gv::picoEvent->bField();
    
    Run_num = Run_num - 19130000;
    
    vz_dist->Fill(Run_num, gv::picoEvent->primaryVertex().z());
    vxy_dist->Fill(Run_num, sqrt(fabs((gv::picoEvent->primaryVertex().x()) * (gv::picoEvent->primaryVertex().x()) + (gv::picoEvent->primaryVertex().y()) * (gv::picoEvent->primaryVertex().y()))));
    refmult_dist->Fill(Run_num, refmult);
    
    return true;
}

////////////////// Beginning of V0's //////////////////
void v0_selection_run_by_run()
{
    Int_t nTracks = gv::picoDst->numberOfTracks();

    for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
    {
        StPicoTrack *track = gv::picoDst->track(iTrk);
        if(!track) continue;
        
        number_tracks_counted++;

        /// Global Tracks Analysis
        gv_lam::gtrack = track;
        if(!v0_dau_selection_run_by_run()) continue;

        track = NULL;
    }
    mean_pt->Fill(Run_num, pt_event/(double)number_tracks_counted);
    mean_eta->Fill(Run_num, eta_event/(double)number_tracks_counted);
    mean_phi->Fill(Run_num, phi_event/(double)number_tracks_counted);
    mean_dca->Fill(Run_num, dca_event/(double)number_tracks_counted);
    mean_nHitsFit->Fill(Run_num, nHitsFit_event/(double)number_tracks_counted);
    mean_nHitsDedx->Fill(Run_num, nHitsDedx_event/(double)number_tracks_counted);
    
    pt_event = 0., eta_event = 0., phi_event = 0., dca_event = 0., nHitsFit_event = 0., nHitsDedx_event = 0.;
    number_tracks_counted = 0;
}

///////////////////// Beginning of Global Track Analysis ///////////////////////
bool v0_dau_selection_run_by_run()
{
    nHitsFit_event += gv_lam::gtrack->nHitsFit();
    nHitsDedx_event += gv_lam::gtrack->nHitsDedx();
    
    StPicoPhysicalHelix helix = gv_lam::gtrack->helix(gv::details->Magn); //inner helix. good for dca to PV.
    TVector3 p = helix.momentum(gv::details->Magn * kilogauss);  //momentum at origin
    double pt = p.Perp();
    
    pt_event += pt;
    eta_event += p.Eta();
    phi_event += p.Phi();
    dca_event += gv_lam::gtrack->gDCA(gv::picoEvent->primaryVertex().x(), gv::picoEvent->primaryVertex().y(), gv::picoEvent->primaryVertex().z());
    
    return true;
}
