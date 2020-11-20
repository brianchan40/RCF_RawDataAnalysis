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
#include "./LambdaMaker.h"
#include "./gv.h"
#include "./gv_lam.h"

using namespace std;

/// INITIALIZATION
void init_Lambda()
{
    if(gv_lam::curr.debug_Event == true) std::cout << "got into init_Lambda()" << std::endl;
    init_macro(); //Initializing Conditions for Code
    if(gv_lam::curr.debug_Event == true) std::cout << "init_macro();" << std::endl;
    init_eventcuts(); //initializing event cuts
    if(gv_lam::curr.debug_Event == true) std::cout << "init_macro()1;" << std::endl;
    init_v0const(); //initializing constants related to V0 based on V0 types.
    if(gv_lam::curr.debug_Event == true) std::cout << "init_macro()2;" << std::endl;
    init_trkcuts(); //initializing Track cuts based on the V0 type
    if(gv_lam::curr.debug_Event == true) std::cout << "init_macro()3;" << std::endl;
    init_hist(); //Create Histograms
    if(gv_lam::curr.debug_Event == true) std::cout << "init_macro()4;" << std::endl;
    //init_tree();
}

/// Make FUNCTION
bool make_Lambda(int cen)
{
    gv::check = 3;

    if(gv_lam::curr.debug_Event == true) std::cout << "Make_Lambda" << endl;
    
    /// Event Cuts & gv::details Stored
    if(!event_cut_passed(cen)) return false;

    if(gv_lam::curr.debug_Event == true) std::cout << "3" << endl;

    store_eventDetails();

    //*******************************************/
    /// Track analysis

    // Global & Primary Tracks Initialized
    init_dau_vectors();

    if(gv_lam::curr.debug_priTrack == true) std::cout << "10" << endl;

    // Finding V0's
    v0_selection();

    if(gv_lam::curr.debug_reconstruct == true) std::cout << "50" << endl;

    return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Supporting Functions for Initialize

//setting condition for analysis
void init_macro()
{
    gv_lam::curr.mV0Type = kAntiLambda;
    gv_lam::curr.mRotate = true;
    mix_events = false;

    gv_lam::curr.mSameSignPlus = false;
    gv_lam::curr.mSameSignMinus = false;
    gv_lam::curr.mDcaAlgoLong = true;
    gv_lam::curr.mDumpNull = false;

    //choosing what to debug
    gv_lam::curr.debug_Event = false;
    gv_lam::curr.debug_priTrack = false;
    gv_lam::curr.debug_gTrack = false;
    gv_lam::curr.debug_gTrack2 = false;
    gv_lam::curr.debug_reconstruct = false;
}

void init_eventcuts()
{
    //statistic information
    gv_lam::cut_eve.mEventsProcessed = 0;                 //  Number of Events read and processed

    //some diagnosing variables
    gv_lam::cut_eve.mTestVZ = 0;
    gv_lam::cut_eve.mTestNTrack = 0;

    gv_lam::cut_eve.cutAbsVertexZLeEq  = 70.;

    vz_index = 0;
}

void init_v0const()
{
    if(gv_lam::curr.mV0Type == kLambda || gv_lam::curr.mV0Type == kAntiLambda)
    {
        // for Lambda and AntiLambda
        gv_lam::constants_v0.mMass1      = 0.93827; // mass of proton
        gv_lam::constants_v0.mMass2      = 0.13957; // mass of pion
        gv_lam::constants_v0.mMassV0     = 1.115684;// mass of Lambda

        if(gv_lam::curr.mV0Type == kLambda)
        {
            gv_lam::constants_v0.mCharge1  = 1;
            gv_lam::constants_v0.mCharge2  = -1;
        }
        else
        {
            gv_lam::constants_v0.mCharge1  = -1;
            gv_lam::constants_v0.mCharge2  = 1;
        }
        //do not setup the cut values here. those are parameters

        //parameters for StDcaService.cxx
        kShiftConnect = 0.3;
        kShiftContain = 0.3;
    }
    else if(gv_lam::curr.mV0Type == kKs)
    {
        // for Ks
        gv_lam::constants_v0.mMass1     = 0.13957;     // mass of pion+
        gv_lam::constants_v0.mMass2     = 0.13957;     // mass of pion-
        gv_lam::constants_v0.mMassV0    = 0.49768;     // mass of K0s

        gv_lam::constants_v0.mCharge1  = 1;
        gv_lam::constants_v0.mCharge2  = -1;
        //do not setup the cut values here. those are parameters

        //parameters for StDcaService.cxx
        kShiftConnect = 0.3;
        kShiftContain = 0.3;
    }
    else
    {
        // for photon
        assert(gv_lam::curr.mV0Type == kPhoton);
        gv_lam::constants_v0.mMass1  = 0.51099907e-3;
        gv_lam::constants_v0.mMass2  = 0.51099907e-3;
        gv_lam::constants_v0.mMassV0 = 0;

        gv_lam::constants_v0.mCharge1  = 1;
        gv_lam::constants_v0.mCharge2  = -1;

        if(gv_lam::curr.mSameSignPlus)
        {
            gv_lam::constants_v0.mCharge1    = 1;
            gv_lam::constants_v0.mCharge2    = 1;
        }

        if(gv_lam::curr.mSameSignMinus)
        {
            gv_lam::constants_v0.mCharge1    = -1;
            gv_lam::constants_v0.mCharge2    = -1;
        }
        kShiftConnect = 2;
        kShiftContain = 2;
    }
}

void init_trkcuts()
{

    if(gv_lam::curr.mV0Type == kPhoton)
    {
        gv_lam::cut_trk.cutNHitsGr = 10;
        gv_lam::cut_trk.cutPtGrEq1 = 0.006;
        gv_lam::cut_trk.cutPtGrEq2 = 0.006;

        gv_lam::cut_trk.cutAbsNSigma1Le = 4.;
        gv_lam::cut_trk.cutAbsNSigma2Le = 4.;
        gv_lam::cut_trk.cutDca1GrEq  = 0.7;
        gv_lam::cut_trk.cutDca2GrEq  = 1.0;

        gv_lam::cut_trk.cutDca1to2LeEq = 1.2;
        gv_lam::cut_trk.cutV0MassWidthLeEq = 0.2;
        gv_lam::cut_trk.cutDcaV0Le    = 0;
        gv_lam::cut_trk.cutV0DecLenGrEq = 3.0;
        gv_lam::cut_trk.cutDau1Dau2Ang3DLe = 0.2;
        gv_lam::cut_trk.cutDau1Dau2DipAngDiffLe = 0.05;
    }
    else if(gv_lam::curr.mV0Type == kLambda || gv_lam::curr.mV0Type == kAntiLambda)
    {
        gv_lam::cut_trk.cutNHitsGr = 15;
        //if(gv_lam::curr.mV0Type == kLambda)
        //{
        gv_lam::cut_trk.cutPtGrEq1  = 0.2; //0.4
        gv_lam::cut_trk.cutPtGrEq2  = 0.15;
        //}
        //else
        //{
        //gv_lam::cut_trk.cutPtGrEq1  = 0.2;
        //gv_lam::cut_trk.cutPtGrEq2  = 0.2;
        //}

        gv_lam::cut_trk.cutAbsNSigma1Le = 4.;
        gv_lam::cut_trk.cutAbsNSigma2Le = 4.;
        //gv_lam::cut_trk.cutDca1GrEq  = 0.5; //0.5
        //gv_lam::cut_trk.cutDca2GrEq  = 1.0; //1.5
        //gv_lam::cut_trk.cutDca1GrEq  = 2.0; //0.5
        //gv_lam::cut_trk.cutDca2GrEq  = 2.0; //1.5

        //gv_lam::cut_trk.cutDca1to2LeEq = 1.2; //1.0, 0.8
        gv_lam::cut_trk.cutV0MassWidthLeEq = 0.07;
        gv_lam::cut_trk.cutV0rdotpGr  = 0.;
        gv_lam::cut_trk.cutDcaV0Le    = 3.5; // 1.0
        //gv_lam::cut_trk.cutV0DecLenGrEq = 3.0;
    }
    else
    {
        gv_lam::cut_trk.cutNHitsGr = 15;
        gv_lam::cut_trk.cutPtGrEq1 = 0.15;
        gv_lam::cut_trk.cutPtGrEq2 = 0.15;

        gv_lam::cut_trk.cutAbsNSigma1Le = 4.;
        gv_lam::cut_trk.cutAbsNSigma2Le = 4.;
        gv_lam::cut_trk.cutDca1GrEq  = 1.0;
        gv_lam::cut_trk.cutDca2GrEq  = 1.0;

        gv_lam::cut_trk.cutDca1to2LeEq = 1.0;
        gv_lam::cut_trk.cutV0MassWidthLeEq = 0.07;
        gv_lam::cut_trk.cutV0rdotpGr  = 0.;
        gv_lam::cut_trk.cutDcaV0Le    = 1.0;
        gv_lam::cut_trk.cutV0DecLenGrEq = 3.0;
    }

    gv_lam::cut_trk.nsigma = 0;
}

void init_hist()
{
    if(gv_lam::curr.debug_Event == true) std::cout << "init_hist" << endl;
    const Int_t nbins = 100   ;

    //QA for events
    gv::all_hists.hNPrimVertex  = new TH1F( "PrimVertex", "Number of Primary Vertex", 10, 0.0, 10.0 ) ;

    //std::cout << "init_macro()5;" << std::endl;
    gv::all_hists.hVertexZ  = new TH1F( "VertexZ", "Event Vertex Z Position", nbins * 4, -100.0, 100.0 ) ;
    gv::all_hists.hNRefMult  = new TH1F( "RefMult", "Reference Multiplicity", 1000, 0.0, 1000.0 ) ;
    gv::all_hists.hSelectNRefMultM0  = new TH1F( "SelectRefMultM0", "Reference Multiplicity of selected events", 1000, 0.0, 1000.0 ) ;
    gv::all_hists.hSelectNRefMultM1  = new TH1F( "SelectRefMultM1", "Number of Lambdas in Each Centrality", 20, 0.0, 10.0 ) ;
    gv::all_hists.hSelectNRefMultM2  = new TH1F( "SelectRefMultM2", "Number of Tracks in Each Centrality", 20, 0.0, 10.0 ) ;

    gv::all_hists.hSelectNRefMultFtpcE  = new TH1F( "SelectRefMultFtpcE", "Reference Multiplicity of selected events (FTPC East)", 1000, 0.0, 1000.0 ) ;
    gv::all_hists.hSelectBbcEAdcSum = new TH1F("SelectBbcEAdcSum", "BBC East Adc Sum", 5000, 0.0, 5000.0);

    if(gv_lam::curr.debug_Event == true) std::cout << "done with QA for events" << endl;
    
    //QA for global tracks
    gv::all_hists.hPtRaw  = new TH1F( "PtRaw", "Transverse Momentum for all particles", nbins * 6, 0.0, 30.0 ) ;
    gv::all_hists.hEtaRaw  = new TH1F( "EtaRaw", "Eta for all particles", nbins * 5, -2, 2 ) ;
    gv::all_hists.hPhiRaw  = new TH1F( "PhiRaw", "Phi for all particles", nbins * 10, -TMath::Pi(), TMath::Pi() ) ;
    gv::all_hists.hEta  = new TH1F( "Eta", "Eta for selected particles", nbins * 5, -2, 2 ) ;
    gv::all_hists.hPhi  = new TH1F( "Phi", "Phi for selected particles", nbins * 10, -TMath::Pi(), TMath::Pi() ) ;
    gv::all_hists.hPhiLowPt  = new TH1F( "PhiLowPt", "Phi for selected particles", nbins * 10, -TMath::Pi(), TMath::Pi() ) ;
    gv::all_hists.hPhiHighPt  = new TH1F( "PhiHighPt", "Phi for selected particles", nbins * 10, -TMath::Pi(), TMath::Pi() ) ;
    gv::all_hists.hDedxP  = new TH2F( "DedxP", "dEdx for selected particles", nbins * 3, -3, 3, nbins * 3, 0, 3e-5) ;
    gv::all_hists.hNSigmaPion  = new TH1F( "nSigmaPion", "nSigmaPion for selected particles", nbins * 2, -10, 10 ) ;
    gv::all_hists.hNSigmaProton  = new TH1F( "nSigmaProton", "nSigmaProton for selected particles", nbins * 2, -10, 10 ) ;
    gv::all_hists.hNSigmaKaon  = new TH1F( "nSigmaKaon", "nSigmaKaon for selected particles", nbins * 2, -10, 10 ) ;

    gv::all_hists.hPtDiff  = new TH1F( "ptdiff", "Reference Multiplicity", 100, 0.0, 0.02 ) ;
    gv::all_hists.hOrDiff  = new TH1F( "ordiff", "Reference Multiplicity", 100, 0.0, 0.2 ) ;
    gv::all_hists.hPDiff  = new TH1F( "pdiff", "Reference Multiplicity", 100, 0.0, 0.4 ) ;
    gv::all_hists.hDcaDiff  = new TH1F( "dcadiff", "Reference Multiplicity", 100, -0.4, 0.4 ) ;
    gv::all_hists.hVertexZDiff = new TH1F("VertexZDiff", "Reference Multiplicity", 200, -10, 10);

    gv::all_hists.charge_dist = new TH1D("charge_dist", "Distribution of charges", 5, -2.5, 2.5);

    gv::all_hists.lambda_pt = new TH1F("lambda_pt", "Pt of Lambda", 300, 0, 15);

    gv::all_hists.primary_vertex = new TH3F("primary_vertex", "Primary Vertex", 200, -10, 10, 200, -10, 10, 200, -10, 10);

    if(gv_lam::curr.debug_Event == true) std::cout << "done with QA for global tracks" << endl;
    
    //QA for V0's
    gv::all_hists.hInvMass = new TH1F("V0Mass", "V0 Inv. Mass", 200, gv_lam::constants_v0.mMassV0 - gv_lam::cut_trk.cutV0MassWidthLeEq, gv_lam::constants_v0.mMassV0 + gv_lam::cut_trk.cutV0MassWidthLeEq) ;

    for(int i = 0; i < 9; i++)
    {
        for(int j = 0; j < 17; j++)
        {
            gv::all_hists.V0Mass_cent[i][j] = new TH1F(TString::Format("V0Mass_%d_%d", i, j), TString::Format("V0 Inv. Mass for Centrality-%d Pt Bin - %d", i, j), 200, gv_lam::constants_v0.mMassV0 - gv_lam::cut_trk.cutV0MassWidthLeEq, gv_lam::constants_v0.mMassV0 + gv_lam::cut_trk.cutV0MassWidthLeEq);
            //V0Mass_pt_cent[i] = new TH1F(TString::Format("V0Mass_%d", i), TString::Format("V0 Inv. Mass for Pt-%d", i), 200, mMassV0-cutV0MassWidthLeEq, mMassV0+cutV0MassWidthLeEq);
        }
    }

    gv::all_hists.cent_check_9 = new TH1I("cent_check_9", "Centrality Check for 9", 20, 0, 20);
    gv::all_hists.cent_check_16 = new TH1I("cent_check_16", "Centrality Check for 16", 20, 0, 20);

    if(gv_lam::curr.debug_Event == true) std::cout << "done with QA for V0's" << endl;
    
    //Deciding Cuts
    gv::all_hists.mom_track_gpt = new TH1F( "Pt", "Transverse Momentum for particles", 600, 0.0, 5.0 ) ;
    gv::all_hists.nHitsFitDist = new TH1F("nHitsFitDist", "nHitsFit Distribution", 500, 0, 100);
    gv::all_hists.nsigmadist = new TH1D("nsigmadist", "nsigma Distribution", 400, -20, 20);
    gv::all_hists.DCA_Dist = new TH1F("DCA_Dist", "DCA Distribution", 200, 0, 10);
    gv::all_hists.dca1to2_Dist = new TH1F("dca1to2_Dist", "dca1to2 Distribution", 200, 0, 10);
    gv::all_hists.v0decaylengthQA = new TH1D("v0decaylengthQA", "v0 decaylength QA", 200, 0, 5);
    gv::all_hists.dcav0toPVQA = new TH1D("dcav0toPVQA", "dca v0 to PV", 180, 0, 3);
    gv::all_hists.entireV0Mass = new TH1D("entireV0Mass", "V0Mass before 0.07 Cut", 200, gv_lam::constants_v0.mMassV0 - gv_lam::cut_trk.cutV0MassWidthLeEq, gv_lam::constants_v0.mMassV0 + gv_lam::cut_trk.cutV0MassWidthLeEq);
    gv::all_hists.transversevertexmag = new TH1D("transversevertexmag", "Magnitude of Transverse Vertex", 400, 0, 10);
    gv::all_hists.pseudoRapidityQA1 = new TH1D("pseudoRapidityQA1", "pseudoRapidity for op1", 400, -2.5, 2.5);
    gv::all_hists.pseudoRapidityQA2 = new TH1D("pseudoRapidityQA2", "pseudoRapidity for op2", 400, -2.5, 2.5);
    gv::all_hists.RapidityQA = new TH1D("RapidityQA", "Rapidity for W0", 400, -2.5, 2.5);
    gv::all_hists.lambdarapidity = new TH1D("lambdarapidity", "Rapidity of Lambda", 400, -2.5, 2.5);
    gv::all_hists.trigger_IDs = new TH1I("trigger_IDs", "Trigger IDs", 700000, 0, 700000);
    cout << "trigger_IDs = " << gv::all_hists.trigger_IDs << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Supporting Functions for Make

//function to determine whether or not they passed the event cuts, and also filling QA plots
bool event_cut_passed(int cen)
{
    //Trigger Information
    std::vector<unsigned int> triggers = gv::picoEvent->triggerIds();
    
    int size_triggers = triggers.size();
    
    if(gv_lam::curr.debug_Event == true) std::cout << "Triggers?" << endl;
    
    for(int i = 0; i < size_triggers; i++)
    {
        if(gv_lam::curr.debug_Event == true) std::cout << "i = " << triggers[i] << endl;
        gv::all_hists.trigger_IDs->Fill(triggers[i]);
    }
    triggers.clear();
    
    if(gv_lam::curr.debug_Event == true) std::cout << "Not Triggers" << endl;
    //if(gv_lam::curr.debug_Event == true) std::cout << gv::picoEvent->triggerIds().at(0) << ", " << gv::picoEvent->triggerIds().at(1) << endl;
    //if (!gv::picoEvent->isTrigger(350003) && !gv::picoEvent->isTrigger(350013) && !gv::picoEvent->isTrigger(350023) && !gv::picoEvent->isTrigger(350033) && !gv::picoEvent->isTrigger(350043)) return false;
    //if (!gv::picoEvent->isTrigger(440001) && !gv::picoEvent->isTrigger(440004) && !gv::picoEvent->isTrigger(440005) && !gv::picoEvent->isTrigger(440015)) return false;
    //if (!gv::picoEvent->isTrigger(450050) && !gv::picoEvent->isTrigger(450060) && !gv::picoEvent->isTrigger(450015) && !gv::picoEvent->isTrigger(450025)) return false;
    //if (gv::picoEvent->isTrigger(450010) || gv::picoEvent->isTrigger(450020) || (!gv::picoEvent->isTrigger(450013) && !gv::picoEvent->isTrigger(450023))) return false;
    if (!gv::picoEvent->isTrigger(610001) && !gv::picoEvent->isTrigger(610011) && !gv::picoEvent->isTrigger(610021) && !gv::picoEvent->isTrigger(610031) && !gv::picoEvent->isTrigger(610041) && !gv::picoEvent->isTrigger(610051)) return false;

    if(gv_lam::curr.debug_Event == true) std::cout << "11" << endl;
    //Centrality definition BEGINS
    Int_t Run_num = gv::picoEvent->runId();
    Int_t Event_num = gv::picoEvent->eventId();
    //Int_t refmult = gv::picoEvent->refMult();
    Int_t refmult = gv::picoEvent->grefMult();
    int Ntofmatch = gv::picoEvent->nBTOFMatch();
    double zdcx = gv::picoEvent->ZDCx();
    
	if(gv_lam::curr.debug_Event == true) std::cout << "12" << endl;

    // For refmult
    gv_lam::refmultCorrUtil  = CentralityMaker::instance()->getRefMultCorr() ;
    //gv_lam::refmultCorrUtil  = CentralityMaker::instance()->getgRefMultCorr_VpdMB30();
    /*if(!gv_lam::refmultCorrUtil){
        gv_lam::refmultCorrUtil  = CentralityMaker::instance()->getgRefMultCorr_P16id();
        gv_lam::refmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
        gv_lam::refmultCorrUtil->readScaleForWeight("/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14_P16id.txt");
    }*/
    
    //gv_lam::refmultCorrUtil  = CentralityMaker::instance()->getgRefMultCorr_VpdMBnoVtx();
    
    // You need to specify the run number you are going to process
    
    //cout << "Run_num = " << Run_num << endl;
    //cout << "gv_lam::refmultCorrUtil->getScaleForWeight(); = " << gv_lam::refmultCorrUtil->getScaleForWeight() << endl;
    if(gv_lam::curr.debug_Event == true) std::cout << Run_num << endl;
    gv_lam::refmultCorrUtil->init(Run_num);

    // Call initEvent(const UShort_t RefMult, const Double_t z) function
    // event-by-event at the beginning before using any other functions
    if(gv_lam::curr.debug_Event == true) std::cout << "Run_num" << endl;
    gv_lam::refmultCorrUtil->initEvent(refmult,gv::picoEvent->primaryVertex().z(),zdcx);
    if(gv_lam::curr.debug_Event == true) std::cout << Run_num << endl;
    if(gv_lam::refmultCorrUtil->isBadRun(Run_num)) return false;
    if(!gv_lam::refmultCorrUtil->passnTofMatchRefmultCut(1.*refmult, 1.*Ntofmatch)) return false;

    if(gv_lam::curr.debug_Event == true) std::cout << "13" << endl;
    //   see StRefMultCorr.h for the definition of centrality bins
    //const Int_t cent16 = gv_lam::refmultCorrUtil->getCentralityBin16() ;
    const Int_t cent9  = gv_lam::refmultCorrUtil->getCentralityBin9() ;
    
    if(gv_lam::curr.mV0Type == kLambda || gv_lam::curr.mV0Type == kAntiLambda) set_centralitydependent_cuts(cent9);

    gv::all_hists.cent_check_9->Fill(cent9);
    //gv::all_hists.cent_check_16->Fill(cent16);

    if(gv_lam::curr.debug_Event == true) std::cout << "6" << endl;

    gv::all_hists.primary_vertex->Fill(gv::picoEvent->primaryVertex().x(), gv::picoEvent->primaryVertex().y(), gv::picoEvent->primaryVertex().z());
    if ( fabs(gv::picoEvent->primaryVertex().x()) < 1e-5 && fabs(gv::picoEvent->primaryVertex().y()) < 1e-5 && fabs(gv::picoEvent->primaryVertex().z()) < 1e-5 )  return false;//Cut on vertext: no-vertex got reported as (0,0,0)

    //Duplicate Events? - Begins
    if ( gv::picoDst->numberOfTracks() == gv_lam::cut_eve.mTestNTrack && gv_lam::cut_eve.mEventsProcessed != 0 && gv_lam::cut_eve.mTestVZ != 0 &&  gv::picoEvent->primaryVertex().z() == gv_lam::cut_eve.mTestVZ )
    {
        std::cout << gv_lam::cut_eve.mEventsProcessed << " " << "seems a duplicated event!" << std::endl;
        return false;
    }
    gv_lam::cut_eve.mTestVZ = gv::picoEvent->primaryVertex().z();
    gv_lam::cut_eve.mTestNTrack = gv::picoDst->numberOfTracks();
    //Duplicate Events? - Ends

    if(gv_lam::curr.debug_Event == true) std::cout << "1" << endl;

    //QA Plots
    gv::all_hists.hVertexZ -> Fill(gv::picoEvent->primaryVertex().z()) ; // Make histogram of the vertex Z distribution
    gv::all_hists.hNRefMult -> Fill( gv::picoEvent->refMult() );
    gv::all_hists.transversevertexmag->Fill(fabs((gv::picoEvent->primaryVertex().x()) * (gv::picoEvent->primaryVertex().x()) + (gv::picoEvent->primaryVertex().y()) * (gv::picoEvent->primaryVertex().y())));

    //cut on vertex z
    if (fabs(gv::picoEvent->primaryVertex().z()) > gv_lam::cut_eve.cutAbsVertexZLeEq ) return false ;

    //Cut on transverse vertex
    if(fabs((gv::picoEvent->primaryVertex().x()) * (gv::picoEvent->primaryVertex().x()) + (gv::picoEvent->primaryVertex().y()) * (gv::picoEvent->primaryVertex().y())) > 2.0) return false;

    if(gv_lam::curr.debug_Event == true) std::cout << "2" << endl;

    //More QA Plots
    gv::all_hists.hSelectNRefMultM0 -> Fill( gv::picoEvent->refMult() );//this is an ESSENTIAL histogram to record the total number of events for certain centrality. always make sure it is filled AFTER event selection!

    gv::all_hists.hSelectNRefMultFtpcE -> Fill( gv::picoEvent->refMultFtpcEast() );
    
    double vpdVz = -990.0;
    vpdVz = gv::picoEvent->vzVpd();
    
    if(fabs(vpdVz - gv::picoEvent->primaryVertex().z()) >= 4) return false;

    return true;
}

void set_centralitydependent_cuts(int cen){
    if(cen <= 1){
        gv_lam::cut_trk.cutDca1to2LeEq = 1.2; //1.0, 0.8
        gv_lam::cut_trk.cutV0DecLenGrEq = 2.5;
        gv_lam::cut_trk.cutDca1GrEq  = 0.1; //0.5
        gv_lam::cut_trk.cutDca2GrEq  = 0.8; //1.5
    }
    else if(cen <= 3){
        gv_lam::cut_trk.cutDca1to2LeEq = 1.1; //1.0, 0.8
        gv_lam::cut_trk.cutV0DecLenGrEq = 3.0;
        gv_lam::cut_trk.cutDca1GrEq  = 0.2; //0.5
        gv_lam::cut_trk.cutDca2GrEq  = 1.0; //1.5
    }
    else if(cen <= 5){
        gv_lam::cut_trk.cutDca1to2LeEq = 1.0; //1.0, 0.8
        gv_lam::cut_trk.cutV0DecLenGrEq = 3.5;
        gv_lam::cut_trk.cutDca1GrEq  = 0.3; //0.5
        gv_lam::cut_trk.cutDca2GrEq  = 1.3; //1.5
    }
    else {
        gv_lam::cut_trk.cutDca1to2LeEq = 0.9; //1.0, 0.8
        gv_lam::cut_trk.cutV0DecLenGrEq = 4.0;
        gv_lam::cut_trk.cutDca1GrEq  = 0.4; //0.5
        gv_lam::cut_trk.cutDca2GrEq  = 1.5; //1.5
    }
}

//function to store event gv::details, and QA plots related to that
void store_eventDetails()
{
    // find primary vertex from VPD vertex and TOF matched global tracks beam projection.
    // if there is no VpdVz, skip events; if there are no TOF primary tracks, skip events.
    // condition for TOF primary tracks, |DCA_xy| < 1cm; |Zproj - vpzvz|<6 cm;
    // the xy coordinates of pv will be beam line, z is from dca piont of closest primary tracks.
    double vpdVz = -990.0;
    vpdVz = gv::picoEvent->vzVpd();

    double bbce = gv::picoEvent->bbcEastRate();
    double bbcw = gv::picoEvent->bbcWestRate();
    double bbcx = gv::picoEvent->BBCx();
    double zdcx = gv::picoEvent->ZDCx();

    //Magnetic Field
    gv::details->Magn = gv::picoEvent->bField();

    if(gv_lam::curr.debug_Event == true) std::cout << "4" << endl;

    //Centrality definition BEGINS
    Int_t Run_num = gv::picoEvent->runId();
    Int_t Event_num = gv::picoEvent->eventId();
    //Int_t refmult = gv::picoEvent->refMult();
    Int_t refmult = gv::picoEvent->grefMult();

    if(gv_lam::curr.debug_Event == true) std::cout << "5" << endl;

    //Centrality definition ENDS

    if(gv_lam::curr.debug_Event == true) std::cout << "7" << endl;

    gv::details->PVtxz = gv::picoEvent->primaryVertex().z();
    gv::details->PVtxx = gv::picoEvent->primaryVertex().x();
    gv::details->PVtxy = gv::picoEvent->primaryVertex().y();
    gv::details->VPDvz = vpdVz;
    gv::details->Run = Run_num;
    gv::details->EventID = Event_num;
    gv::details->Eweight = gv_lam::refmultCorrUtil->getWeight();
    //gv::details->TOFMult = gv::picoEvent->btofTrayMultiplicity();
    gv::details->TOFMult = gv::picoEvent->nBTOFMatch();
    gv::details->RefMult = refmult;
    const Int_t cent9  = gv_lam::refmultCorrUtil->getCentralityBin9() ;
    gv::details->cent = cent9;

    if(gv_lam::curr.debug_Event == true) std::cout << "8" << endl;

    //Prepare for Track analysis
    gv::details->nLambda = 0;
    gv::details->num_trk = 0;
    gv::details->n_Proton = 0;

    //cout << "Did I even get here?" << endl;
    //gv_lam::refmultCorrUtil->clear();

    if(gv_lam::curr.debug_Event == true) std::cout << "9" << endl;
}

void init_dau_vectors()
{
    //gv_lam::mDau1 = new v0_dau_vectors();
    //gv_lam::mDau2 = new v0_dau_vectors();

    //std::cout << gv_lam::mDau1->mDauDcaVec.size() << endl;
    //std::cout << gv_lam::mDau1->nPrimary << endl;

    gv_lam::mDau1->mDauDcaVec.clear();
    gv_lam::mDau2->mDauDcaVec.clear();
    gv_lam::mDau1->mDauVec.clear();
    gv_lam::mDau2->mDauVec.clear();

    memset(gv_lam::mDau1->PrimaryTrackID, 0, sizeof(gv_lam::mDau1->PrimaryTrackID));
    memset(gv_lam::mDau1->PrimaryTrackPx, 0, sizeof(gv_lam::mDau1->PrimaryTrackPx));
    memset(gv_lam::mDau1->PrimaryTrackPy, 0, sizeof(gv_lam::mDau1->PrimaryTrackPy));
    memset(gv_lam::mDau1->PrimaryTrackPz, 0, sizeof(gv_lam::mDau1->PrimaryTrackPz));
    gv_lam::mDau1->nPrimary = 0;
}

////////////////// Beginning of V0's //////////////////
void v0_selection()
{
    if(gv_lam::curr.debug_priTrack == true) std::cout << "v0_selection" << endl;
    Int_t nTracks = gv::picoDst->numberOfTracks();

    for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
    {
        if(gv_lam::curr.debug_priTrack == true) std::cout << "trk loop" << endl;

        StPicoTrack *track = gv::picoDst->track(iTrk);
        if(!track) continue;

        if(gv_lam::curr.debug_priTrack == true) std::cout << "11" << endl;

        /// Primary Tracks Analysis
        /*if(track->isPrimary())
        {
            gv_lam::ptrack = track;
            if(!v0_selection_primary()) continue;
        }*/

        /// Global Tracks Analysis
        //else if(!track->isPrimary())
        //{
            gv_lam::ptrack = track;
            gv_lam::gtrack = track;
            gv::all_hists.hSelectNRefMultM2->Fill(gv::details->cent);
            if(!v0_dau_selection()) continue;
        //}

        //else std::cout << "Something really strange has happened... not primary and not not primary" << std::endl;

        //delete track;
        track = NULL;
    }
    
    gv_lam::mDau1->mag_field = gv::details->Magn;
    gv_lam::mDau2->mag_field = gv::details->Magn;

    //std::cout << gv_lam::mDau1->mDauVec.size() << ", " << gv_lam::mDau2->mDauVec.size() << endl;

    assert(gv_lam::mDau1->mDauVec.size() == gv_lam::mDau1->mDauDcaVec.size() && gv_lam::mDau2->mDauVec.size() == gv_lam::mDau2->mDauDcaVec.size());

    if(gv_lam::curr.debug_gTrack2) std::cout << "35" << endl;

    // for mix_events background analysis: replace mDau1 with buffer protons
    if(mix_events) mix_events_check();

    //reconstruct V0
    reconstruct_V0();

    // for mix_events background analysis: fill buffer protons with original mDau1 protons
    if(mix_events) mix_events_fillbuffer();
}

///////////////////// Beginning of Global Track Analysis ///////////////////////
bool v0_dau_selection()
{

//    if(primary_cuts_decision() && gv_lam::ptrack->isPrimary()){
//        gv_lam::mDau1->PrimaryTrackID[gv_lam::mDau1->nPrimary] = gv_lam::ptrack->id();
//        gv_lam::mDau1->nPrimary++;
//    }

    v0_dau_selection_QAPlots();

    unsigned short nHitsFit = gv_lam::gtrack->nHitsFit();
    short charge = 0;
    charge = gv_lam::gtrack->charge();

    if(gv_lam::curr.debug_gTrack == true) std::cout << "18" << endl;

    double nsigmapion = gv_lam::gtrack->nSigmaPion();
    double nsigmaproton = gv_lam::gtrack->nSigmaProton();
    double nsigmakaon = gv_lam::gtrack->nSigmaKaon();
    double nsigmaelectron = gv_lam::gtrack->nSigmaElectron();
    double dedx = gv_lam::gtrack->dEdx();

    if(gv_lam::curr.debug_gTrack == true) std::cout << "19" << endl;

    StPicoPhysicalHelix helix = gv_lam::gtrack->helix(gv::details->Magn); //inner helix. good for dca to PV.
    TVector3 p = helix.momentum(gv::details->Magn * kilogauss);  //momentum at origin
    TVector3 origin = helix.origin();  //origin of helix
    double pt = p.Perp();

    if(gv_lam::curr.debug_gTrack == true) std::cout << "20" << endl;

    v0_dau_selection_sigma();

    //OK. let's cut tracks
    if(!v0_dau_selection_trkcuts()) return false;

    if(gv_lam::curr.debug_gTrack == true) std::cout << "25" << endl;

    /// storing particle gv::details for reconstructing Event Plane
    event_plane_storage_1();

    if(!v0_dau_selection_finalized()) return false;

    if(gv_lam::curr.debug_gTrack2) std::cout << "32" << endl;

    return true;
}
void v0_dau_selection_QAPlots()
{
    gv::all_hists.hPtRaw -> Fill( gv_lam::gtrack->gPt() ) ; //at dca to PV, for a global track, this value is useless. anyway, the pt value is supposed to be the same anywhere.

    short charge = 0;
    charge = gv_lam::gtrack->charge();

    gv::all_hists.charge_dist->Fill(charge);

    StPicoPhysicalHelix helix = gv_lam::gtrack->helix(gv::details->Magn); //inner helix. good for dca to PV.
    TVector3 p = helix.momentum(gv::details->Magn * kilogauss);  //momentum at origin
    TVector3 origin = helix.origin();  //origin of helix

    //some checks.
    gv::all_hists.hPtDiff -> Fill( gv_lam::gtrack->gPt() - p.Perp() ) ;
    gv::all_hists.hOrDiff -> Fill( (gv_lam::gtrack->origin() - origin).Mag() ) ;
    gv::all_hists.hPDiff -> Fill( fabs(gv_lam::gtrack->gMom().Mag() - p.Mag()) ) ;
    //comments: there are difference between the values above. But they seem to be acceptably small!

    unsigned short nHitsFit = gv_lam::gtrack->nHitsFit();
    gv::all_hists.nHitsFitDist->Fill(nHitsFit);
}

void v0_dau_selection_sigma()
{
    short charge = 0;
    charge = gv_lam::gtrack->charge();

    if(charge == gv_lam::constants_v0.mCharge1)
    {
        if(gv_lam::curr.mV0Type == kLambda || gv_lam::curr.mV0Type == kAntiLambda) gv_lam::cut_trk.nsigma = gv_lam::gtrack->nSigmaProton();
        else if(gv_lam::curr.mV0Type == kKs) gv_lam::cut_trk.nsigma = gv_lam::gtrack->nSigmaPion();
        else gv_lam::cut_trk.nsigma = gv_lam::gtrack->nSigmaElectron();

        if(gv_lam::curr.debug_gTrack == true) std::cout << "21" << endl;
    }
    else
    {
        if(gv_lam::curr.mV0Type != kPhoton)  gv_lam::cut_trk.nsigma = gv_lam::gtrack->nSigmaPion();
        else gv_lam::cut_trk.nsigma = gv_lam::gtrack->nSigmaElectron();

        if(gv_lam::curr.debug_gTrack == true) std::cout << "22" << endl;
    }
}

bool v0_dau_selection_trkcuts()
{
    unsigned short nHitsFit = gv_lam::gtrack->nHitsFit();
    short charge = 0;
    charge = gv_lam::gtrack->charge();

    if(((float)nHitsFit) / ((float)gv_lam::gtrack->nHitsMax()) < 0.52) return false;//Added on after the meeting

    if((nHitsFit <= gv_lam::cut_trk.cutNHitsGr) || (abs(charge) != 1)) return false;

    if(gv_lam::curr.debug_gTrack == true) std::cout << "23" << endl;

    gv::all_hists.hNSigmaPion->Fill(gv_lam::gtrack->nSigmaPion());
    gv::all_hists.hNSigmaProton->Fill(gv_lam::gtrack->nSigmaProton());
    gv::all_hists.hNSigmaKaon->Fill(gv_lam::gtrack->nSigmaKaon());

    StPicoPhysicalHelix helix = gv_lam::gtrack->helix(gv::details->Magn); //inner helix. good for dca to PV.
    TVector3 p = helix.momentum(gv::details->Magn * kilogauss);  //momentum at origin
    gv::all_hists.hDedxP->Fill(p.Mag()*charge, gv_lam::gtrack->dEdx());

    gv::all_hists.mom_track_gpt->Fill(p.Perp());

    if(gv_lam::curr.debug_gTrack == true) std::cout << "24" << endl;

    gv::all_hists.nsigmadist->Fill(gv_lam::cut_trk.nsigma);

    return true;
}

void event_plane_storage_1()
{
    TVector3 primary_p;

    StPicoPhysicalHelix helix = gv_lam::gtrack->helix(gv::details->Magn); //inner helix. good for dca to PV.
    primary_p = helix.momentum(gv::details->Magn * kilogauss);  //momentum at origin
    double pathlength = helix.pathLength(gv::picoEvent->primaryVertex(), false); // do scan periods. NOTE: the default is false. this is not necessary for tracks with pt>0.15GeV/c
    TVector3 dca = helix.at(pathlength) - gv::picoEvent->primaryVertex();

    int index = gv_lam::gtrack->bTofPidTraitsIndex();
    int tofflag;
    if(index >= 0) tofflag = (gv::picoDst->btofPidTraits(index))->btofMatchFlag();

    if(gv_lam::curr.debug_gTrack == true) std::cout << "26" << endl;
    
    if(primary_p.Pt() > 2.0 || primary_p.Pt() < 0.15) return;
    if(dca.Mag() > 2) return;
    if(primary_p.Eta() > 1.0 || primary_p.Eta() < -1.0) return;

    //std::cout << "gv::details->num_trk = " << gv::details->num_trk << endl;
    gv::p_all[gv::details->num_trk].px = primary_p.X();
    //std::cout << "gv::p_all[gv::details->num_trk].Px() = " << gv::p_all[gv::details->num_trk].Px() << endl;
    gv::p_all[gv::details->num_trk].py = primary_p.Y();
    gv::p_all[gv::details->num_trk].pz = primary_p.Z() ;
    gv::p_all[gv::details->num_trk].Charge = gv_lam::gtrack->charge();
    gv::p_all[gv::details->num_trk].nSigmaProton = gv_lam::gtrack->nSigmaProton();
    gv::p_all[gv::details->num_trk].trk_id = gv_lam::gtrack->id();

    if(gv_lam::gtrack->isPrimary()) gv::p_all[gv::details->num_trk].prim = 1;
    else gv::p_all[gv::details->num_trk].prim = 0;

    gv::p_all[gv::details->num_trk].dcaglobal = dca.Mag();
    
    //if(primary_p.Perp() >= 2) std::cout << "perp mom = " << primary_p.Perp() << ", while prim = " << gv::p_all[gv::details->num_trk].prim << endl;

    if(gv_lam::gtrack->isTofTrack()) gv::p_all[gv::details->num_trk].isTofTrack = 1;
    else gv::p_all[gv::details->num_trk].isTofTrack = 0;

    if(index > 0)
    {
        gv::p_all[gv::details->num_trk].TOFflag = tofflag;
        
        double tof = (gv::picoDst->btofPidTraits(index))->btof();
        double BtofYLocal = (gv::picoDst->btofPidTraits(index))->btofYLocal();
        gv::p_all[gv::details->num_trk].ToF = tof;
        gv::p_all[gv::details->num_trk].BTOFYLocal = BtofYLocal;
        
        TVector3 pp = helix.momentum(gv::details->Magn * kilogauss);
        double beta = (gv::picoDst->btofPidTraits(index))->btofBeta();
        float mass2 = pp.Mag2() * (1.0 / beta / beta - 1.0);
        gv::p_all[gv::details->num_trk].Mass2 = mass2;

    }
    else
    {
        gv::p_all[gv::details->num_trk].TOFflag = -1;
    }

    if(gv_lam::curr.debug_gTrack == true) std::cout << "27" << endl;

    gv::details->num_trk++;
}

bool v0_dau_selection_finalized()
{
    short charge = 0;
    charge = gv_lam::gtrack->charge();
    StPicoPhysicalHelix helix = gv_lam::gtrack->helix(gv::details->Magn); //inner helix. good for dca to PV.
    TVector3 p = helix.momentum(gv::details->Magn * kilogauss);  //momentum at origin

    if(charge == gv_lam::constants_v0.mCharge1 && fabs(gv_lam::cut_trk.nsigma) < gv_lam::cut_trk.cutAbsNSigma1Le)
    {
        if((gv_lam::curr.mV0Type == kLambda || gv_lam::curr.mV0Type == kAntiLambda) && ((p.Perp() < gv_lam::cut_trk.cutPtGrEq1) || (p.Perp() > 2.0))) return false;
        else if(p.Perp() < gv_lam::cut_trk.cutPtGrEq) return false;

        //if((fabs(gv_lam::gtrack->nSigmaPion()) >= fabs(gv_lam::cut_trk.nsigma)) || (fabs(gv_lam::cut_trk.nsigma) < 0.5)){
        if(!filling_mDau1()) return false;
        //}
    }

    if(charge == gv_lam::constants_v0.mCharge2 && fabs(gv_lam::cut_trk.nsigma) < gv_lam::cut_trk.cutAbsNSigma2Le)
    {
        if((p.Perp() < gv_lam::cut_trk.cutPtGrEq2)) return false;
        if(!filling_mDau2()) return false;
    }

    return true;
}

bool filling_mDau1()
{
    if(gv_lam::curr.debug_gTrack == true) std::cout << "28" << endl;
    //record the first daughter
    //fill the vector
    StPicoPhysicalHelix helix = gv_lam::gtrack->helix(gv::details->Magn); //inner helix. good for dca to PV.

    TVector3 pv = gv::picoEvent->primaryVertex();

    double pathlength = helix.pathLength(gv::picoEvent->primaryVertex(), false); // do scan periods. NOTE: the default is false. this is not necessary for tracks with pt>0.15GeV/c
    TVector3 dca = helix.at(pathlength) - gv::picoEvent->primaryVertex();

    gv::all_hists.DCA_Dist->Fill(dca.Mag());
    if(dca.Mag() < gv_lam::cut_trk.cutDca1GrEq) return false;

    gv_lam::mDau1->mDauDcaVec.push_back(dca.Mag());
    gv_lam::mDau1->mDauVec.push_back(gv_lam::gtrack);

    //std::cout << mDau.mDauCharge.size() << ", " << mDau.mDauCharge[mDau.mDauCharge.size() - 1] << endl;
    if(gv_lam::curr.debug_gTrack2 == true) std::cout << "29" << endl;

    return true;
}

bool filling_mDau2()
{
    if(gv_lam::curr.debug_gTrack == true) std::cout << "30" << endl;
    //record the first daughter
    //fill the vector
    StPicoPhysicalHelix helix = gv_lam::gtrack->helix(gv::details->Magn); //inner helix. good for dca to PV.

    TVector3 pv = gv::picoEvent->primaryVertex();

    //rotate transverse coordinates and momenta for background estimation
    if(gv_lam::curr.mRotate)
    {
        TVector3 p1 = helix.momentum(gv::details->Magn * kilogauss);  //momentum at origin
        TVector3 x1 = helix.origin();    //origin
        p1.SetX(-p1.x());
        p1.SetY(-p1.y());
        x1.SetX(-(x1.x() - pv.x()) + pv.x());
        x1.SetY(-(x1.y() - pv.y()) + pv.y());
        StPicoPhysicalHelix helixtmp(p1, x1, gv::details->Magn * kilogauss, gv_lam::gtrack->charge());
        helix = helixtmp;
    }

    double pathlength = helix.pathLength(gv::picoEvent->primaryVertex(), false); // do scan periods. NOTE: the default is false. this is not necessary for tracks with pt>0.15GeV/c
    TVector3 dca = helix.at(pathlength) - gv::picoEvent->primaryVertex();

    gv::all_hists.DCA_Dist->Fill(dca.Mag());
    if(dca.Mag() < gv_lam::cut_trk.cutDca2GrEq) return false;

    gv_lam::mDau2->mDauDcaVec.push_back(dca.Mag());
    gv_lam::mDau2->mDauVec.push_back(gv_lam::gtrack);

    //std::cout << mDau.mDauCharge.size() << ", " << mDau.mDauCharge[mDau.mDauCharge.size() - 1] << endl;
    if(gv_lam::curr.debug_gTrack2 == true) std::cout << "31" << endl;

    return true;
}
///////////////////// End of Global Track Analysis ///////////////////////

///////////////////// Beginning of Reconstruction of Lambdas ///////////////////////
void reconstruct_V0()
{
    for(unsigned int i = 0; i < gv_lam::mDau1->mDauVec.size(); i++)
    {
        init_track1();

        if(gv_lam::curr.debug_reconstruct) std::cout << "110" << endl;

        bool holder = getinfofortracks1(i); //continue;

        if(gv_lam::curr.debug_reconstruct) std::cout << "33" << endl;

        bool duplicates = false;
        for(int k = 0; k < gv::details->nLambda / 2; k++)
        {
            if(gv_lam::trk1_ids[k] == gv_lam::track1.track->id())
            {
                duplicates = true;
                break;
            }
        }
        if(duplicates) continue;

        for(unsigned int j = 0; j < gv_lam::mDau2->mDauVec.size(); j++)
        {
            //get pion track info here
            init_track2();
            holder = getinfofortracks2(j); //continue;

            if(gv_lam::curr.debug_reconstruct) std::cout << "34" << endl;

            if(gv_lam::track2.track->id() == gv_lam::track1.track->id())continue;  //for same sign
            if((gv_lam::curr.mSameSignPlus || gv_lam::curr.mSameSignMinus) && j <= i)continue; //avoid double counting in s.s.

            for(int k1 = 0; k1 < gv::details->nLambda / 2; k1++)
            {
                if(gv_lam::trk2_ids[k1] == gv_lam::track2.track->id())
                {
                    duplicates = true;
                    break;
                }
            }
            if(duplicates) continue;

            if(gv_lam::curr.debug_reconstruct) std::cout << "35" << endl;

            calc_Dcatrks();

            if(!cuts_relationsoftwotracks()) continue;

            if(gv_lam::curr.debug_reconstruct == true) std::cout << "46" << endl;

            //calculate p1_p and p2_p
            /*TVector3 p1_p, p2_p;
            calc_p_p1(p1_p, curr_trk.xv0);
            calc_p_p2(p2_p, curr_trk.xv0);*/
            //Not sure what this is for

            if(gv_lam::curr.debug_reconstruct == true) std::cout << "47" << endl;

            //Calculate Pt and split into different bins of pt and centrality
            if(!final_data_organization(i, j)) continue;

            if(gv_lam::curr.debug_reconstruct == true) std::cout << "51" << endl;
        }
    }
}

void init_track1()
{
    if(gv_lam::curr.debug_reconstruct) std::cout << "111" << endl;
    //gv_lam::track1.track->Clear();
    //gv_lam::track1.helix.Clear();
    gv_lam::track1.p.SetXYZ(0, 0, 0);
    if(gv_lam::curr.debug_reconstruct) std::cout << "112" << endl;
    gv_lam::track1.pt = 0;
    gv_lam::track1.dca = 0;
    gv_lam::track1.tofflag = 0;
    gv_lam::track1.tof = 0;
    gv_lam::track1.tofpos.SetXYZ(0, 0, 0);
}

bool getinfofortracks1(int i)
{
    if(gv_lam::curr.debug_reconstruct) std::cout << gv_lam::mDau1->mDauDcaVec[i] << endl;

    gv_lam::track1.track = gv_lam::mDau1->mDauVec[i];
    gv_lam::track1.dca = gv_lam::mDau1->mDauDcaVec[i];

    if(!mix_events){
        if(gv_lam::curr.debug_reconstruct) std::cout << "101" << endl;
    
        gv_lam::track1.helix = gv_lam::track1.track->helix(gv::details->Magn); //inner helix. good for dca to PV.

        TVector3 pv = gv::picoEvent->primaryVertex();

        if(gv_lam::curr.debug_reconstruct) std::cout << "102" << endl;

        gv_lam::track1.p = gv_lam::track1.helix.momentum(gv::details->Magn * kilogauss);  //momentum at origin
        //TVector3 origin1 = helix1.origin();  //origin of helix
        gv_lam::track1.pt = gv_lam::track1.p.Perp();

        //record TOF information...
    
        int index = gv_lam::track1.track->bTofPidTraitsIndex();
        if(gv_lam::curr.debug_reconstruct) std::cout << "103" << endl;
        if(index < 0) return false;
        if(gv_lam::curr.debug_reconstruct) std::cout << "index = " << index << endl;

        gv_lam::track1.tofflag = (gv::picoDst->btofPidTraits(index))->btofMatchFlag();
        if(gv_lam::curr.debug_reconstruct) std::cout << "103a" << index << endl;
        gv_lam::track1.tof = (gv::picoDst->btofPidTraits(index))->btof();
        if(gv_lam::curr.debug_reconstruct) std::cout << "103b" << index << endl;
        gv_lam::track1.tofpos = (gv::picoDst->btofPidTraits(index))->btofHitPos();
    }
    else if(mix_events){
        
        if(gv_lam::curr.debug_reconstruct) std::cout << "101" << endl;
        
        gv_lam::track1.helix = gv_lam::track1.track->helix(gv_lam::mDau1->mag_field); //inner helix. good for dca to PV.
        
        TVector3 pv = gv::picoEvent->primaryVertex();
        
        if(gv_lam::curr.debug_reconstruct) std::cout << "102" << endl;
        
        gv_lam::track1.p = gv_lam::track1.helix.momentum(gv_lam::mDau1->mag_field * kilogauss);  //momentum at origin
        gv_lam::track1.pt = gv_lam::track1.p.Perp();
        
    }

    if(gv_lam::curr.debug_reconstruct) std::cout << "104" << endl;

    return true;
}

void init_track2()
{
    gv_lam::track2.p.SetXYZ(0, 0, 0);
    gv_lam::track2.pt = 0;
    gv_lam::track2.dca = 0;
    gv_lam::track2.tofflag = 0;
    gv_lam::track2.tof = 0;
    gv_lam::track2.tofpos.SetXYZ(0, 0, 0);
}

bool getinfofortracks2(int i)
{
    gv_lam::track2.track = gv_lam::mDau2->mDauVec[i];
    gv_lam::track2.dca = gv_lam::mDau2->mDauDcaVec[i];

    gv_lam::track2.helix = gv_lam::track2.track->helix(gv::details->Magn); //inner helix. good for dca to PV.

    TVector3 pv = gv::picoEvent->primaryVertex();
    if(gv_lam::curr.mRotate)
    {
        TVector3 tp1 = gv_lam::track2.helix.momentum(gv::details->Magn * kilogauss);  //momentum at origin
        TVector3 tx1 = gv_lam::track2.helix.origin();    //origin
        tp1.SetX(-tp1.x());
        tp1.SetY(-tp1.y());
        tx1.SetX(-(tx1.x() - pv.x()) + pv.x());
        tx1.SetY(-(tx1.y() - pv.y()) + pv.y());
        StPicoPhysicalHelix helixtmp(tp1, tx1, gv::details->Magn * kilogauss, gv_lam::track2.track->charge());
        gv_lam::track2.helix = helixtmp;
    }

    gv_lam::track2.p = gv_lam::track2.helix.momentum(gv::details->Magn * kilogauss);  //momentum at origin
    //TVector3 origin1 = helix1.origin();  //origin of helix
    gv_lam::track2.pt = gv_lam::track2.p.Perp();

    //record TOF information...
    int index = gv_lam::track2.track->bTofPidTraitsIndex();
    if(index < 0) return false;

    gv_lam::track2.tofflag = (gv::picoDst->btofPidTraits(index))->btofMatchFlag();
    gv_lam::track2.tof = (gv::picoDst->btofPidTraits(index))->btof();
    gv_lam::track2.tofpos = (gv::picoDst->btofPidTraits(index))->btofHitPos();

    return true;
}

void calc_Dcatrks()
{
    if(!gv_lam::curr.mDcaAlgoLong) not_DcaAlgoLong();
    //StHelix method above is VERY SLOW. use it only for checking the consistency of
    //long's code
    else gv_lam::curr_trk.dca1to2 = closestDistance(gv_lam::track1.helix, gv_lam::track2.helix, gv::details->Magn, gv::picoEvent->primaryVertex(), gv_lam::curr_trk.xv0, gv_lam::curr_trk.op1, gv_lam::curr_trk.op2);

    gv_lam::linfo.x = gv_lam::curr_trk.xv0.X();
    gv_lam::linfo.y = gv_lam::curr_trk.xv0.Y();
    gv_lam::linfo.z = gv_lam::curr_trk.xv0.Z();
}

void not_DcaAlgoLong()
{
    pair<double, double> tmps = gv_lam::track1.helix.pathLengths(gv_lam::track2.helix);

    gv_lam::curr_trk.ox1 = gv_lam::track1.helix.at(tmps.first);
    gv_lam::curr_trk.ox2 = gv_lam::track2.helix.at(tmps.second);
    gv_lam::curr_trk.dca1to2 = (gv_lam::curr_trk.ox1 - gv_lam::curr_trk.ox2).Mag();
    gv_lam::curr_trk.xv0 = (gv_lam::curr_trk.ox1 + gv_lam::curr_trk.ox2) * 0.5;
    gv_lam::curr_trk.op1 = gv_lam::track1.helix.momentumAt(tmps.first, gv::details->Magn * kilogauss);
    gv_lam::curr_trk.op2 = gv_lam::track2.helix.momentumAt(tmps.second, gv::details->Magn * kilogauss);

}

bool cuts_relationsoftwotracks()
{
    //cut on dca1to2
    if(gv_lam::curr.debug_reconstruct == true) std::cout << "38" << endl;

    gv::all_hists.dca1to2_Dist->Fill(gv_lam::curr_trk.dca1to2);
    if(gv_lam::curr_trk.dca1to2 > gv_lam::cut_trk.cutDca1to2LeEq) return false;

    //cut on v0 mass
    if(!cut_v0mass()) return false;

    if(gv_lam::curr.debug_reconstruct == true) std::cout << "39" << endl;

    //r dot p for v0. cut on it. should be larger than 0.
    if(!cut_v0rdotp()) return false;

    if(gv_lam::curr.debug_reconstruct == true) std::cout << "43" << endl;

    //cut on decay length
    if(!cut_v0decaylen()) return false;

    if(gv_lam::curr.debug_reconstruct == true) std::cout << "44" << endl;

    //cut on sinth, or theta
    //double sinth = (xv0toPV.Cross(pv0)).Mag() / xv0toPV.Mag() / pv0.Mag();
    //double theta = atan2(sinth, rdotp / xv0toPV.Mag() / pv0.Mag()); //theta range: from 0 to pi

    fill_finalplots_v0_dau();

    return true;
}

bool cut_v0mass()
{
    double oe1 = sqrt(gv_lam::curr_trk.op1.Mag2() + gv_lam::constants_v0.mMass1 * gv_lam::constants_v0.mMass1);
    double oe2 = sqrt(gv_lam::curr_trk.op2.Mag2() + gv_lam::constants_v0.mMass2 * gv_lam::constants_v0.mMass2);

    TLorentzVector dau1_fvec, dau2_fvec;
    dau1_fvec.SetPxPyPzE(gv_lam::curr_trk.op1.x(), gv_lam::curr_trk.op1.y(), gv_lam::curr_trk.op1.z(), oe1);
    dau2_fvec.SetPxPyPzE(gv_lam::curr_trk.op2.x(), gv_lam::curr_trk.op2.y(), gv_lam::curr_trk.op2.z(), oe2);
    gv_lam::curr_trk.Lambda_fvec = dau1_fvec + dau2_fvec;

    double v0mass = gv_lam::curr_trk.Lambda_fvec.M();

    gv_lam::linfo.mass = v0mass;
    if(fabs(gv_lam::linfo.mass - gv_lam::constants_v0.mMassV0) > gv_lam::cut_trk.cutV0MassWidthLeEq) return false;

    /*// Check for Photon Conversion
    double altMass1 = 0.000511;
    double altMass2 = 0.000511;

    double altoe1 = sqrt(gv_lam::curr_trk.op1.Mag2() + altMass1 * altMass1);
    double altoe2 = sqrt(gv_lam::curr_trk.op2.Mag2() + altMass2 * altMass2);

    TLorentzVector altdau1_fvec, altdau2_fvec;
    altdau1_fvec.SetPxPyPzE(gv_lam::curr_trk.op1.x(), gv_lam::curr_trk.op1.y(), gv_lam::curr_trk.op1.z(), altoe1);
    altdau2_fvec.SetPxPyPzE(gv_lam::curr_trk.op2.x(), gv_lam::curr_trk.op2.y(), gv_lam::curr_trk.op2.z(), altoe2);
    TLorentzVector altLambda_fvec = altdau1_fvec + altdau2_fvec;

    double altv0mass = altLambda_fvec.M();

    if(altv0mass < 0.005) return false;
    // End of Check for Photon Conversion*/

    return true;
}

bool cut_v0rdotp()
{
    TVector3 pv0 = gv_lam::curr_trk.op1 + gv_lam::curr_trk.op2;
    TVector3 xv0toPV = gv_lam::curr_trk.xv0 - gv::picoEvent->primaryVertex();

    //helix of v0: straight line
    StPicoPhysicalHelix helixv0(pv0, gv_lam::curr_trk.xv0, 0, 0);

    //pthead, ptarm cut
    double pthead1 = gv_lam::curr_trk.op1.Dot(pv0) / pv0.Mag();
    double pthead2 = gv_lam::curr_trk.op2.Dot(pv0) / pv0.Mag();
    double ptarm   = sqrt(gv_lam::curr_trk.op1.Mag2() - pthead1 * pthead1);

    //forward decay cut
    double ang1 = gv_lam::curr_trk.op1.x() * xv0toPV.x() + gv_lam::curr_trk.op1.y() * xv0toPV.y();
    double ang2 = gv_lam::curr_trk.op2.x() * xv0toPV.x() + gv_lam::curr_trk.op2.y() * xv0toPV.y();

    //r dot p for v0. cut on it. should be larger than 0.
    double rdotp = xv0toPV.Dot(pv0) ;
    if(gv_lam::curr.mV0Type != kPhoton && rdotp <= gv_lam::cut_trk.cutV0rdotpGr) return false;

    //calculate v0 to PV dca. v0 carry no charge. straight line. cut on dca
    double dcav0toPV = rdotp * rdotp / pv0.Mag2();
    dcav0toPV = sqrt( xv0toPV.Mag2() - dcav0toPV);
    gv::all_hists.dcav0toPVQA->Fill(dcav0toPV);
    if(dcav0toPV >= gv_lam::cut_trk.cutDcaV0Le) return false;
    gv_lam::linfo.dcaV0toPv = dcav0toPV;

    return true;
}

bool cut_v0decaylen()
{
    TVector3 xv0toPV = gv_lam::curr_trk.xv0 - gv::picoEvent->primaryVertex();

    double v0decaylength = xv0toPV.Mag();
    gv::all_hists.v0decaylengthQA->Fill(v0decaylength);
    //gv::all_hists.entireV0Mass->Fill(gv_lam::linfo.mass);
    if(v0decaylength < gv_lam::cut_trk.cutV0DecLenGrEq) return false;
    if((gv_lam::curr.mV0Type == kPhoton && acos(gv_lam::curr_trk.op1.Dot(gv_lam::curr_trk.op2) / gv_lam::curr_trk.op1.Mag() / gv_lam::curr_trk.op2.Mag()) > gv_lam::cut_trk.cutDau1Dau2Ang3DLe)) return false;
    if ((gv_lam::curr.mV0Type == kPhoton && fabs(gv_lam::track1.helix.dipAngle() - gv_lam::track2.helix.dipAngle()) > gv_lam::cut_trk.cutDau1Dau2DipAngDiffLe)) return false;

    return true;
}

void fill_finalplots_v0_dau()
{
    TVector3 pv0 = gv_lam::curr_trk.op1 + gv_lam::curr_trk.op2;
    gv::all_hists.RapidityQA->Fill(log( (sqrt(gv_lam::linfo.mass * gv_lam::linfo.mass + pv0.Mag2()) + pv0.Z()) / sqrt(gv_lam::linfo.mass * gv_lam::linfo.mass + pv0.Perp2())));
    gv::all_hists.pseudoRapidityQA1->Fill(gv_lam::curr_trk.op1.PseudoRapidity());
    gv::all_hists.pseudoRapidityQA2->Fill(gv_lam::curr_trk.op2.PseudoRapidity());
    //gv::all_hists.hInvMass->Fill(gv_lam::linfo.mass);

}

/*void calc_p_p1(TVector3 p_p, TVector3 xv0)
{
    p_p.SetXYZ(0, 0, 0);
    int p_id = 0;
    for (int nPri = 0; nPri < gv_lam::mDau1->nPrimary; nPri++)
    {
        if(gv_lam::track1.track->id() == gv_lam::mDau1->PrimaryTrackID[nPri])
        {
            p_id = 1;
            p_p.SetXYZ(gv_lam::mDau1->PrimaryTrackPx[nPri], gv_lam::mDau1->PrimaryTrackPy[nPri], gv_lam::mDau1->PrimaryTrackPz[nPri]);
        }
    }

    double tofpathlen = -999.;
    if(gv_lam::track1.tofflag > 0) tofpathlen = gv_lam::track1.helix.pathLength(gv_lam::track1.tofpos) - gv_lam::track1.helix.pathLength(xv0);
}

void calc_p_p2(TVector3 p_p, TVector3 xv0)
{
    p_p.SetXYZ(0, 0, 0);
    int p_id = 0;
    for (int nPri = 0; nPri < gv_lam::mDau1->nPrimary; nPri++)
    {
        if(gv_lam::track2.track->id() == gv_lam::mDau1->PrimaryTrackID[nPri])
        {
            p_id = 1;
            p_p.SetXYZ(gv_lam::mDau1->PrimaryTrackPx[nPri], gv_lam::mDau1->PrimaryTrackPy[nPri], gv_lam::mDau1->PrimaryTrackPz[nPri]);
        }
    }

    double tofpathlen = -999.;
    if(gv_lam::track2.tofflag > 0) tofpathlen = gv_lam::track2.helix.pathLength(gv_lam::track2.tofpos) - gv_lam::track2.helix.pathLength(xv0);
}*/

bool final_data_organization(int i, int j)
{
    int temp_index_cent = -1;
    TVector3 pv0 = gv_lam::curr_trk.op1 + gv_lam::curr_trk.op2;
    if(gv::details->cent >= 4) temp_index_cent = gv::details->cent - 2;
    else temp_index_cent = floor(gv::details->cent / 2.0);
    if((temp_index_cent < 0) /*|| (fabs(pv0.PseudoRapidity()) > 1)*/) return false;

    //filling daughter partilces p to remove from reaction plane reconstruction, and lambda information into tree
    gv::all_hists.hInvMass->Fill(gv_lam::linfo.mass);
    
    //if(gv_lam::curr_trk.Lambda_fvec.Perp() < 3.0){
    if((gv_lam::curr_trk.Lambda_fvec.Perp() <= 2.2) && (gv_lam::curr_trk.Lambda_fvec.Perp() >= 0.5)){
        if(gv_lam::curr.debug_reconstruct == true) std::cout << "(int)floor((gv_lam::curr_trk.Lambda_fvec.Perp()-0.5)/0.1) = " << (int)floor((gv_lam::curr_trk.Lambda_fvec.Perp()-0.5)/0.1) << endl;
        if(gv_lam::curr.debug_reconstruct == true) std::cout << "gv_lam::curr_trk.Lambda_fvec.Perp() = " << gv_lam::curr_trk.Lambda_fvec.Perp() << endl;
        if(gv_lam::curr.debug_reconstruct == true) std::cout << "gv::details->cent = " << gv::details->cent << endl;
        if((int)floor((gv_lam::curr_trk.Lambda_fvec.Perp()-0.5)/0.1) == 17) {
            if(gv_lam::curr.debug_reconstruct == true) cout << "forced to be 16" << endl;
            gv::all_hists.V0Mass_cent[gv::details->cent][16]->Fill(gv_lam::linfo.mass);
        }
        else {
            if(gv_lam::curr.debug_reconstruct == true) cout << "yay not force " << (int)floor((gv_lam::curr_trk.Lambda_fvec.Perp()-0.5)/0.1) << endl;
            //gv::all_hists.V0Mass_cent[gv::details->cent][6]->Fill(1);
            gv::all_hists.V0Mass_cent[gv::details->cent][(int)floor((gv_lam::curr_trk.Lambda_fvec.Perp()-0.5)/0.1)]->Fill(gv_lam::linfo.mass);
        }
    }
    
    /*if(gv_lam::curr.mRotate){
        if(!((gv_lam::linfo.mass < 1.09) || (gv_lam::linfo.mass > 1.14))){
            if(!((gv_lam::linfo.mass < 1.125) && (gv_lam::linfo.mass > 1.105))){
                fill_tree(i, j);
            }
        }
    }
    else{*/
        if((gv_lam::curr.mV0Type == kLambda) || (gv_lam::curr.mV0Type == kAntiLambda)){
		if(!((gv_lam::linfo.mass < 1.113) || (gv_lam::linfo.mass > 1.119))){
//            	if(((gv_lam::linfo.mass < 1.1) && (gv_lam::linfo.mass > 1.109)) || ((gv_lam::linfo.mass < 1.15) && (gv_lam::linfo.mass > 1.14))){
			fill_tree(i, j);
        	}
	}
	else if(gv_lam::curr.mV0Type == kKs){
		if(!((gv_lam::linfo.mass < 0.49) || (gv_lam::linfo.mass > 0.505))){
                        fill_tree(i, j);
                }
	}
    //}

    //fill histograms with centrality and momentum sorted
    fill_hist(temp_index_cent);

    if(gv_lam::curr.debug_reconstruct == true) std::cout << "49" << endl;

    return true;
}

void fill_tree(int i, int j)
{
    gv::p_lambda[gv::details->nLambda / 2].px = gv_lam::curr_trk.Lambda_fvec.Px();
    gv::p_lambda[gv::details->nLambda / 2].py = gv_lam::curr_trk.Lambda_fvec.Py();
    gv::p_lambda[gv::details->nLambda / 2].pz = gv_lam::curr_trk.Lambda_fvec.Pz();
    gv::p_lambda[gv::details->nLambda / 2].mass = gv_lam::linfo.mass;
    gv::p_lambda[gv::details->nLambda / 2].dcaglobal = gv_lam::linfo.dcaV0toPv;
    gv::p_lambda[gv::details->nLambda / 2].x = gv_lam::linfo.x;
    gv::p_lambda[gv::details->nLambda / 2].y = gv_lam::linfo.y;
    gv::p_lambda[gv::details->nLambda / 2].z = gv_lam::linfo.z;

    if(gv_lam::curr.mV0Type == kLambda){
	gv::p_lambda[gv::details->nLambda / 2].Charge = 1;
    }
    else if(gv_lam::curr.mV0Type == kAntiLambda){
        gv::p_lambda[gv::details->nLambda / 2].Charge = -1;
    }

    gv_lam::trk1_ids[gv::details->nLambda / 2] = gv_lam::track1.track->id();
    gv_lam::trk2_ids[gv::details->nLambda / 2] = gv_lam::track2.track->id();

    gv::p_dau[gv::details->nLambda].px = gv_lam::curr_trk.op1.x();
    gv::p_dau[gv::details->nLambda].py = gv_lam::curr_trk.op1.y();
    gv::p_dau[gv::details->nLambda].pz = gv_lam::curr_trk.op1.z();
    gv::p_dau[gv::details->nLambda].Charge = gv_lam::mDau1->mDauVec[i]->charge();
    gv::p_dau[gv::details->nLambda].trk_id = gv_lam::track1.track->id();

    gv::details->nLambda++;

    gv::p_dau[gv::details->nLambda].px = gv_lam::curr_trk.op2.x();
    gv::p_dau[gv::details->nLambda].py = gv_lam::curr_trk.op2.y();
    gv::p_dau[gv::details->nLambda].pz = gv_lam::curr_trk.op2.z();
    gv::p_dau[gv::details->nLambda].Charge = gv_lam::mDau2->mDauVec[j]->charge();
    gv::p_dau[gv::details->nLambda].trk_id = gv_lam::track2.track->id();

    gv::details->nLambda++;
    
    gv::all_hists.hSelectNRefMultM1->Fill(gv::details->cent);
    gv::all_hists.entireV0Mass->Fill(gv_lam::linfo.mass);
}

void fill_hist(int temp_index_cent)
{
    float lambda_pt_val = gv_lam::curr_trk.Lambda_fvec.Pt();

    gv::all_hists.lambda_pt->Fill(lambda_pt_val);
    gv::all_hists.lambdarapidity->Fill(gv_lam::curr_trk.Lambda_fvec.PseudoRapidity());

    int pt_bin = determine_pt_bins(lambda_pt_val);

    //if(pt_bin >= 0) gv::all_hists.V0Mass_pt_cent[0][temp_index_cent] -> Fill(gv_lam::linfo.mass);
}

int determine_pt_bins(float lambda_pt_val)
{
    int bin_pt = -1;

    if(lambda_pt_val > 0.4 && lambda_pt_val <= 2.0) bin_pt = ceil((lambda_pt_val - 0.4) / 0.2) - 1;
    else if(lambda_pt_val > 2.0 && lambda_pt_val <= 2.6) bin_pt = ceil((lambda_pt_val - 2.0) / 0.3) + 7;
    else if(lambda_pt_val > 2.6 && lambda_pt_val <= 4.2) bin_pt = ceil((lambda_pt_val - 2.6) / 0.4) + 9;

    return bin_pt;
}
///////////////////// End of Reconstruction of Lambdas ///////////////////////

///// Beginning of Mix Events //////
void mix_events_check()
{

    vz_index = int((fabs(gv::picoEvent->primaryVertex().z()) + 50) / 2.5);
    deque<v0_dau_vectors>::iterator pro;

    // Fill tmp with mDau1
    for(unsigned int i = 0; i < gv_lam::mDau1->mDauVec.size(); i++){
        gv_lam::track_tmp = (StPicoTrack*)gv_lam::mDau1->mDauVec[i]->Clone();
        gv_lam::buffer_protons_tmp.mDauDcaVec.push_back(gv_lam::mDau1->mDauDcaVec[i]);
        gv_lam::buffer_protons_tmp.mDauVec.push_back(gv_lam::track_tmp);
    }
    gv_lam::buffer_protons_tmp.mag_field = gv_lam::mDau1->mag_field;

    // Reset mDau1
    gv_lam::mDau1->mDauDcaVec.clear();
    gv_lam::mDau1->mDauVec.clear();

    memset(gv_lam::mDau1->PrimaryTrackID, 0, sizeof(gv_lam::mDau1->PrimaryTrackID));
    memset(gv_lam::mDau1->PrimaryTrackPx, 0, sizeof(gv_lam::mDau1->PrimaryTrackPx));
    memset(gv_lam::mDau1->PrimaryTrackPy, 0, sizeof(gv_lam::mDau1->PrimaryTrackPy));
    memset(gv_lam::mDau1->PrimaryTrackPz, 0, sizeof(gv_lam::mDau1->PrimaryTrackPz));
    gv_lam::mDau1->nPrimary = 0;
    gv_lam::mDau1->mag_field = 0;

    // fill mDau1 with buffer protons
    for(pro = gv_lam::buffer_protons[vz_index].begin(); pro != gv_lam::buffer_protons[vz_index].end(); ++pro)
    {
        for(unsigned int i = 0; i < pro->mDauVec.size(); i++)
        {
            gv_lam::mDau1->mDauDcaVec.push_back(pro->mDauDcaVec[i]);
            gv_lam::mDau1->mDauVec.push_back(pro->mDauVec[i]);
        }
        gv_lam::mDau1->mag_field = pro->mag_field;
    }
}

void mix_events_fillbuffer()
{
    if(gv_lam::buffer_protons_tmp.mDauVec.size())
    {
        //cout << "gv_lam::buffer_protons_tmp.mDauVec.size() = " << gv_lam::buffer_protons_tmp.mDauVec.size() << endl;
        gv_lam::buffer_protons[vz_index].push_back(gv_lam::buffer_protons_tmp);
        while(gv_lam::buffer_protons[vz_index].size() > 10) gv_lam::buffer_protons[vz_index].pop_front();
        //cout << "vz_index = " << vz_index << endl;
    }
    
    // Reset tmp buffer
    gv_lam::buffer_protons_tmp.mDauDcaVec.clear();
    gv_lam::buffer_protons_tmp.mDauVec.clear();
    
    memset(gv_lam::buffer_protons_tmp.PrimaryTrackID, 0, sizeof(gv_lam::buffer_protons_tmp.PrimaryTrackID));
    memset(gv_lam::buffer_protons_tmp.PrimaryTrackPx, 0, sizeof(gv_lam::buffer_protons_tmp.PrimaryTrackPx));
    memset(gv_lam::buffer_protons_tmp.PrimaryTrackPy, 0, sizeof(gv_lam::buffer_protons_tmp.PrimaryTrackPy));
    memset(gv_lam::buffer_protons_tmp.PrimaryTrackPz, 0, sizeof(gv_lam::buffer_protons_tmp.PrimaryTrackPz));
    gv_lam::buffer_protons_tmp.nPrimary = 0;

}
