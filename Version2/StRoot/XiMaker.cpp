/*
 * Author: Brian Chan
 * Date: February 25, 2019
 *
 * XiMaker.C is the new Maker-less version of StXiMaker. Produces Xi particles, and preform physics analysis.
 *
 **/

/// This is needed for calling standalone classes
#define _VANILLA_ROOT_

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

/// picoDst headers
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
#include "./XiMaker.h"
#include "./gv_xi.h"
#include "./gv.h"

using namespace std;


void init_xi()
{
    init_macro_xi();
    init_const_xi();
    init_param_xi();
    init_hist_xi();
    init_tree_xi();
}

void make_Xi()
{
    init_dau_vectors_xi();

    // Do some cleaning here
    gv_xi::mXiDst.nxi = 0;
    if(gv::details->nLambda == 0) return;

    storeEventDetails_xi();

    gv_xi::xi_calc.pv.SetX(gv_xi::mXiDst.primvertexX);
    gv_xi::xi_calc.pv.SetY(gv_xi::mXiDst.primvertexY);
    gv_xi::xi_calc.pv.SetZ(gv_xi::mXiDst.primvertexZ);

    xi_selection(); // selecting tracks for analysis

    //std::cout << "gv_xi::mDau1->mDauVec.size() = " << gv_xi::mDau1->mDauVec.size() << ", gv_xi::mDau1->mDauDcaVec.size() = " << gv_xi::mDau1->mDauDcaVec.size() << endl;
    assert(gv_xi::mDau1->mDauVec.size() == gv_xi::mDau1->mDauDcaVec.size());

    xi_reconstruction();

    gv_xi::cut_eve.mEventsProcessed++;

}

////// Initializing ///////
void init_macro_xi()
{
    gv_xi::curr.mXiType = kXi;  //Lambda as default!

    gv_xi::curr.mRotate = true;
    gv_xi::curr.mDcaAlgoLong = true; //use LONG's dca method as default

    gv_xi::cut_eve.mEventsProcessed = 0     ;                    // Zero the Number of Events processed by the maker
    gv_xi::cut_eve.mTestNTrack = 0;
    gv_xi::cut_eve.mTestVZ = 0;
}
void init_const_xi()
{
    //initialize the constant for different Xi types.

    if(gv_xi::curr.mXiType == kXi || gv_xi::curr.mXiType == kAntiXi)
    {
        // for Xi or AntiXi
        gv_xi::constants_v0.mMass1      = 0.93827; // mass of proton or pbar
        gv_xi::constants_v0.mMass2      = 0.13957; // mass of pion- or pion+
        gv_xi::constants_v0.mMassV0     = 1.115684;// mass of Lambda or AntiLambda
        gv_xi::constants_v0.mMassBachelor = 0.13957; //mass of pion- or pion+
        gv_xi::constants_v0.mMassXi = 1.32131;   //mass of Xi or AntiXi

        if(gv_xi::curr.mXiType == kXi)
        {
            gv_xi::constants_v0.mCharge1    = 1;
            gv_xi::constants_v0.mCharge2    = -1;
            gv_xi::constants_v0.mChargeBachelor = -1;
            gv_xi::constants_v0.mChargeXi = -1;
        }
        else
        {
            gv_xi::constants_v0.mCharge1    = -1;
            gv_xi::constants_v0.mCharge2    = 1;
            gv_xi::constants_v0.mChargeBachelor = 1;
            gv_xi::constants_v0.mChargeXi = 1;
        }
        //do not setup the cut values here. those are parameters

        //parameters for StDcaService.cxx
        kShiftConnect = 0.3;
        kShiftContain = 0.3;
    }
    else
    {
        assert(gv_xi::curr.mXiType == kOmega || gv_xi::curr.mXiType == kAntiOmega);
        // for Omega or AntiOmega
        gv_xi::constants_v0.mMass1      = 0.93827; // mass of proton or pbar
        gv_xi::constants_v0.mMass2      = 0.13957; // mass of pion- or pion+
        gv_xi::constants_v0.mMassV0     = 1.115684;// mass of Lambda or AntiLambda
        gv_xi::constants_v0.mMassBachelor = 0.493677; //mass of K- or K+
        gv_xi::constants_v0.mMassXi = 1.67245;   //mass of Omega or AntiOmega

        if(gv_xi::curr.mXiType == kOmega)
        {
            gv_xi::constants_v0.mCharge1    = 1;
            gv_xi::constants_v0.mCharge2    = -1;
            gv_xi::constants_v0.mChargeBachelor = -1;
            gv_xi::constants_v0.mChargeXi = -1;
        }
        else
        {
            gv_xi::constants_v0.mCharge1    = -1;
            gv_xi::constants_v0.mCharge2    = 1;
            gv_xi::constants_v0.mChargeBachelor = 1;
            gv_xi::constants_v0.mChargeXi = 1;
        }
        //do not setup the cut values here. those are parameters

        //parameters for StDcaService.cxx
        kShiftConnect = 0.3;
        kShiftContain = 0.3;
    }
}

void init_param_xi()
{
    //setup the cut values here. do not hard-code them in ::Make()

    gv_xi::cut_trk.cutNHitsGr = 15;
    gv_xi::cut_trk.cutPtGrEq = 0.15;

    gv_xi::cut_trk.cutAbsNSigma1Le = 4.;
    gv_xi::cut_trk.cutDca1GrEq  = 0.5; //1.0

    gv_xi::cut_trk.cutV0DcaGrEq = 0; //0.3
    gv_xi::cut_trk.cutV0MassWidthLeEq = 0.02;

    gv_xi::cut_trk.cutDca1to2LeEq = 1.0; //0.8
    gv_xi::cut_trk.cutXiMassWidthLeEq = 0.05;
    //cutDauPtArmLeEq = 0.3;
    //cutAbsDausPtShoulderDiffLeEq = 1.2;
    //cutDau1DecAngGr = 0.;
    gv_xi::cut_trk.cutXirdotpGr  = 0;
    gv_xi::cut_trk.cutXircrosspLeEq = 0.5; //0.2
    //cutDcaXiLe    = 3.3;
    gv_xi::cut_trk.cutXiDecLenGrEq = 2.0; //2.0
}

void init_hist_xi()
{
    // Create Histograms
    // there is no better way to set QA histograms. do not use histogram arrays or vectors.
    // it is not useful. there are no need to operate all the histograms at the same time.

    const Int_t    nbins    =  100   ;

    gv::all_hists_xi.hVertexZ  = new TH1F( "Vertex", "Event Vertex Z Position", nbins, -25.0, 25.0 ) ;
    gv::all_hists_xi.hPt  = new TH1F( "Pt_xi", "Transverse Momentum for all particles", nbins, 0.0, 10.0 ) ;
    gv::all_hists_xi.hNPrimVertex  = new TH1F( "PrimVertex_xi", "Number of Primary Vertex", 10, 0.0, 10.0 ) ;
    gv::all_hists_xi.hNRefMult  = new TH1F( "RefMult_xi", "Reference Multiplicity", 200, 0.0, 1000.0 ) ;
    gv::all_hists_xi.hPtDiff  = new TH1F( "ptdiff_xi", "Reference Multiplicity", 100, 0.0, 0.02 ) ;
    gv::all_hists_xi.hOrDiff  = new TH1F( "ordiff_xi", "Reference Multiplicity", 100, 0.0, 0.2 ) ;
    gv::all_hists_xi.hPDiff  = new TH1F( "pdiff_xi", "Reference Multiplicity", 100, 0.0, 0.4 ) ;
    gv::all_hists_xi.hDcaDiff  = new TH1F( "dcadiff_xi", "Reference Multiplicity", 100, -0.4, 0.4 ) ;
    gv::all_hists_xi.hInvMass  = new TH1F( "invmass", "Xi Inv. Mass", 100, gv_xi::constants_v0.mMassXi - gv_xi::cut_trk.cutXiMassWidthLeEq, gv_xi::constants_v0.mMassXi + gv_xi::cut_trk.cutXiMassWidthLeEq ) ;
    gv::all_hists_xi.hInvMassSTD  = new TH1F( "invmassSTD", "STD Xi Inv. Mass", 100, gv_xi::constants_v0.mMassXi - gv_xi::cut_trk.cutXiMassWidthLeEq, gv_xi::constants_v0.mMassXi + gv_xi::cut_trk.cutXiMassWidthLeEq ) ;
}

void init_tree_xi()
{
    gv::mXiTree->Branch("runnumber", &gv_xi::mXiDst.runnumber, "runnumber/I");
    gv::mXiTree->Branch("evtnumber", &gv_xi::mXiDst.evtnumber, "evtnumber/I");
    gv::mXiTree->Branch("nrefmult", &gv_xi::mXiDst.nrefmult, "nrefmult/I");
    gv::mXiTree->Branch("primvertexX", &gv_xi::mXiDst.primvertexX, "primvertexX/F");
    gv::mXiTree->Branch("primvertexY", &gv_xi::mXiDst.primvertexY, "primvertexY/F");
    gv::mXiTree->Branch("primvertexZ", &gv_xi::mXiDst.primvertexZ, "primvertexZ/F");
    gv::mXiTree->Branch("magn", &gv_xi::mXiDst.magn, "magn/F");
    gv::mXiTree->Branch("nxi", &gv_xi::mXiDst.nxi, "nxi/I");

    gv::mXiTree->Branch("v0mass", gv_xi::mXiDst.v0mass, "v0mass[nxi]/F");
    //gv::mXiTree->Branch("v0pt", gv_xi::mXiDst.v0pt, "v0pt[nxi]/F");
    //gv::mXiTree->Branch("v0rapidity", gv_xi::mXiDst.v0rapidity, "v0rapidity[nxi]/F");
    //gv::mXiTree->Branch("v0eta", gv_xi::mXiDst.v0eta, "v0eta[nxi]/F");
    gv::mXiTree->Branch("v0x", gv_xi::mXiDst.v0x, "v0x[nxi]/F");
    gv::mXiTree->Branch("v0y", gv_xi::mXiDst.v0y, "v0y[nxi]/F");
    gv::mXiTree->Branch("v0z", gv_xi::mXiDst.v0z, "v0z[nxi]/F");
    gv::mXiTree->Branch("v0px", gv_xi::mXiDst.v0px, "v0px[nxi]/F");
    gv::mXiTree->Branch("v0py", gv_xi::mXiDst.v0py, "v0py[nxi]/F");
    gv::mXiTree->Branch("v0pz", gv_xi::mXiDst.v0pz, "v0pz[nxi]/F");
    //gv::mXiTree->Branch("v0declen", gv_xi::mXiDst.v0declen, "v0declen[nxi]/F");
    gv::mXiTree->Branch("v0dca", gv_xi::mXiDst.v0dca, "v0dca[nxi]/F");
    //gv::mXiTree->Branch("v0dca2d", gv_xi::mXiDst.v0dca2d, "v0dca2d[nxi]/F");
    gv::mXiTree->Branch("v0pathlen", gv_xi::mXiDst.v0pathlen, "v0pathlen[nxi]/F");

    gv::mXiTree->Branch("dau1id", gv_xi::mXiDst.dau1id, "dau1id[nxi]/I");
    gv::mXiTree->Branch("dau1dca", gv_xi::mXiDst.dau1dca, "dau1dca[nxi]/F");
    //gv::mXiTree->Branch("dau1dca2d", gv_xi::mXiDst.dau1dca2d, "dau1dca2d[nxi]/F");
    //gv::mXiTree->Branch("dau1nhits", gv_xi::mXiDst.dau1nhits, "dau1nhits[nxi]/I");
    //gv::mXiTree->Branch("dau1dedx", gv_xi::mXiDst.dau1dedx, "dau1dedx[nxi]/F");
    //gv::mXiTree->Branch("dau1nsigma", gv_xi::mXiDst.dau1nsigma, "dau1nsigma[nxi]/F");
    //gv::mXiTree->Branch("dau1eta", gv_xi::mXiDst.dau1eta, "dau1eta[nxi]/F");
    //gv::mXiTree->Branch("dau1pt", gv_xi::mXiDst.dau1pt, "dau1pt[nxi]/F");
    gv::mXiTree->Branch("dau1px", gv_xi::mXiDst.dau1px, "dau1px[nxi]/F");
    gv::mXiTree->Branch("dau1py", gv_xi::mXiDst.dau1py, "dau1py[nxi]/F");
    gv::mXiTree->Branch("dau1pz", gv_xi::mXiDst.dau1pz, "dau1pz[nxi]/F");
    //gv::mXiTree->Branch("dau1tpc", gv_xi::mXiDst.dau1tpc, "dau1tpc[nxi]/I");
    //gv::mXiTree->Branch("dau1ssd", gv_xi::mXiDst.dau1ssd, "dau1ssd[nxi]/I");
    //gv::mXiTree->Branch("dau1svt", gv_xi::mXiDst.dau1svt, "dau1svt[nxi]/I");
    //gv::mXiTree->Branch("dau1tofflag", gv_xi::mXiDst.dau1tofflag, "dau1tofflag[nxi]/I");
    //gv::mXiTree->Branch("dau1tof", gv_xi::mXiDst.dau1tof, "dau1tof[nxi]/F");
    //gv::mXiTree->Branch("dau1pathlen", gv_xi::mXiDst.dau1pathlen, "dau1pathlen[nxi]/F");

    gv::mXiTree->Branch("dau2id", gv_xi::mXiDst.dau2id, "dau2id[nxi]/I");
    gv::mXiTree->Branch("dau2dca", gv_xi::mXiDst.dau2dca, "dau2dca[nxi]/F");
    //gv::mXiTree->Branch("dau2dca2d", gv_xi::mXiDst.dau2dca2d, "dau2dca2d[nxi]/F");
    //gv::mXiTree->Branch("dau2nhits", gv_xi::mXiDst.dau2nhits, "dau2nhits[nxi]/I");
    //gv::mXiTree->Branch("dau2dedx", gv_xi::mXiDst.dau2dedx, "dau2dedx[nxi]/F");
    //gv::mXiTree->Branch("dau2nsigma", gv_xi::mXiDst.dau2nsigma, "dau2nsigma[nxi]/F");
    //gv::mXiTree->Branch("dau2eta", gv_xi::mXiDst.dau2eta, "dau2eta[nxi]/F");
    //gv::mXiTree->Branch("dau2pt", gv_xi::mXiDst.dau2pt, "dau2pt[nxi]/F");
    gv::mXiTree->Branch("dau2px", gv_xi::mXiDst.dau2px, "dau2px[nxi]/F");
    gv::mXiTree->Branch("dau2py", gv_xi::mXiDst.dau2py, "dau2py[nxi]/F");
    gv::mXiTree->Branch("dau2pz", gv_xi::mXiDst.dau2pz, "dau2pz[nxi]/F");
    //gv::mXiTree->Branch("dau2tpc", gv_xi::mXiDst.dau2tpc, "dau2tpc[nxi]/I");
    //gv::mXiTree->Branch("dau2ssd", gv_xi::mXiDst.dau2ssd, "dau2ssd[nxi]/I");
    //gv::mXiTree->Branch("dau2svt", gv_xi::mXiDst.dau2svt, "dau2svt[nxi]/I");
    //gv::mXiTree->Branch("dau2tofflag", gv_xi::mXiDst.dau2tofflag, "dau2tofflag[nxi]/I");
    //gv::mXiTree->Branch("dau2tof", gv_xi::mXiDst.dau2tof, "dau2tof[nxi]/F");
    //gv::mXiTree->Branch("dau2pathlen", gv_xi::mXiDst.dau2pathlen, "dau2pathlen[nxi]/F");

    //gv::mXiTree->Branch("dca1to2", gv_xi::mXiDst.dca1to2, "dca1to2[nxi]/F");

    gv::mXiTree->Branch("bachid", gv_xi::mXiDst.bachid, "bachid[nxi]/I");
    gv::mXiTree->Branch("bachdca", gv_xi::mXiDst.bachdca, "bachdca[nxi]/F");
    gv::mXiTree->Branch("bachdca2d", gv_xi::mXiDst.bachdca2d, "bachdca2d[nxi]/F");
    gv::mXiTree->Branch("bachnhits", gv_xi::mXiDst.bachnhits, "bachnhits[nxi]/I");
    gv::mXiTree->Branch("bachdedx", gv_xi::mXiDst.bachdedx, "bachdedx[nxi]/F");
    gv::mXiTree->Branch("bachnsigma", gv_xi::mXiDst.bachnsigma, "bachnsigma[nxi]/F");
    gv::mXiTree->Branch("bacheta", gv_xi::mXiDst.bacheta, "bacheta[nxi]/F");
    gv::mXiTree->Branch("bachpt", gv_xi::mXiDst.bachpt, "bachpt[nxi]/F");
    gv::mXiTree->Branch("bachpx", gv_xi::mXiDst.bachpx, "bachpx[nxi]/F");
    gv::mXiTree->Branch("bachpy", gv_xi::mXiDst.bachpy, "bachpy[nxi]/F");
    gv::mXiTree->Branch("bachpz", gv_xi::mXiDst.bachpz, "bachpz[nxi]/F");
    gv::mXiTree->Branch("bachtpc", gv_xi::mXiDst.bachtpc, "bachtpc[nxi]/I");
    gv::mXiTree->Branch("bachssd", gv_xi::mXiDst.bachssd, "bachssd[nxi]/I");
    gv::mXiTree->Branch("bachsvt", gv_xi::mXiDst.bachsvt, "bachsvt[nxi]/I");
    gv::mXiTree->Branch("bachtofflag", gv_xi::mXiDst.bachtofflag, "bachtofflag[nxi]/I");
    gv::mXiTree->Branch("bachtof", gv_xi::mXiDst.bachtof, "bachtof[nxi]/F");
    gv::mXiTree->Branch("bachpathlen", gv_xi::mXiDst.bachpathlen, "bachpathlen[nxi]/F");

    gv::mXiTree->Branch("dcav0tobach", gv_xi::mXiDst.dcav0tobach, "dcav0tobach[nxi]/F");
    //gv::mXiTree->Branch("stdcav0tobach",gv_xi::mXiDst.stdcav0tobach,"stdcav0tobach[nxi]/F");

    gv::mXiTree->Branch("ximass", gv_xi::mXiDst.ximass, "ximass[nxi]/F");
    gv::mXiTree->Branch("xipt", gv_xi::mXiDst.xipt, "xipt[nxi]/F");
    gv::mXiTree->Branch("xirapidity", gv_xi::mXiDst.xirapidity, "xirapidity[nxi]/F");
    gv::mXiTree->Branch("xieta", gv_xi::mXiDst.xieta, "xieta[nxi]/F");
    gv::mXiTree->Branch("xix", gv_xi::mXiDst.xix, "xix[nxi]/F");
    gv::mXiTree->Branch("xiy", gv_xi::mXiDst.xiy, "xiy[nxi]/F");
    gv::mXiTree->Branch("xiz", gv_xi::mXiDst.xiz, "xiz[nxi]/F");
    gv::mXiTree->Branch("xipx", gv_xi::mXiDst.xipx, "xipx[nxi]/F");
    gv::mXiTree->Branch("xipy", gv_xi::mXiDst.xipy, "xipy[nxi]/F");
    gv::mXiTree->Branch("xipz", gv_xi::mXiDst.xipz, "xipz[nxi]/F");
    gv::mXiTree->Branch("xideclen", gv_xi::mXiDst.xideclen, "xideclen[nxi]/F");
    gv::mXiTree->Branch("xidca", gv_xi::mXiDst.xidca, "xidca[nxi]/F");
    gv::mXiTree->Branch("xidca2d", gv_xi::mXiDst.xidca2d, "xidca2d[nxi]/F");
    gv::mXiTree->Branch("xisinth", gv_xi::mXiDst.xisinth, "xisinth[nxi]/F");
    gv::mXiTree->Branch("xipathlen", gv_xi::mXiDst.xipathlen, "xipathlen[nxi]/F");
}
////// End of Initializing ///////

////// Beginning of Make //////
void init_dau_vectors_xi()
{
    //gv_xi::mDau1 = new v0_dau_vectors();
    //gv_xi::mDau2 = new v0_dau_vectors();
    
    gv_xi::mDau1->mDauDcaVec.clear();
    gv_xi::mDau2->mDauDcaVec.clear();
    gv_xi::mDau1->mDauVec.clear();
    gv_xi::mDau2->mDauVec.clear();
    
    memset(gv_xi::mDau1->PrimaryTrackID, 0, sizeof(gv_xi::mDau1->PrimaryTrackID));
    memset(gv_xi::mDau1->PrimaryTrackPx, 0, sizeof(gv_xi::mDau1->PrimaryTrackPx));
    memset(gv_xi::mDau1->PrimaryTrackPy, 0, sizeof(gv_xi::mDau1->PrimaryTrackPy));
    memset(gv_xi::mDau1->PrimaryTrackPz, 0, sizeof(gv_xi::mDau1->PrimaryTrackPz));
    gv_xi::mDau1->nPrimary = 0;
    gv_xi::mDau2->nPrimary = 0;
}

void storeEventDetails_xi()
{
    gv_xi::mXiDst.runnumber = gv::details->Run;
    gv_xi::mXiDst.evtnumber = gv::details->EventID;
    gv_xi::mXiDst.nrefmult = gv::details->RefMult;
    gv_xi::mXiDst.primvertexX = gv::details->PVtxx;
    gv_xi::mXiDst.primvertexY = gv::details->PVtxy;
    gv_xi::mXiDst.primvertexZ = gv::details->PVtxz;
    // Do 'event' analysis based on event data

    // Record some information...
    Double_t magn = gv::picoEvent->bField();
    //Double_t magn = muEvent->magneticField();   //checked! the same as above
    gv_xi::mXiDst.magn = magn;
}

////// Beginning of xi_selection() ////////
void xi_selection()
{
    Int_t nTracks = gv::picoDst->numberOfTracks();

    for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
    {
        StPicoTrack *track = gv::picoDst->track(iTrk);
        if(!track) continue;

        if(track->isPrimary()) continue; //Selecting only the global tracks

        unsigned short nHits = track->nHits();
        unsigned short nHitsFit = track->nHitsFit();
        short charge = track->charge();
        double nsigmapion = track->nSigmaPion();
        double nsigmaproton = track->nSigmaProton();
        double nsigmakaon = track->nSigmaKaon();
        double dedx = track->dEdx();

        if(gv_xi::curr.mXiType == kXi || gv_xi::curr.mXiType == kAntiXi) gv_xi::cut_trk.nsigma = nsigmapion;
        else gv_xi::cut_trk.nsigma = nsigmakaon;

        StPicoPhysicalHelix helix = track->helix(gv_xi::mXiDst.magn); //inner helix. good for dca to PV.
        TVector3 p = helix.momentum(gv_xi::mXiDst.magn * kilogauss);  //momentum at origin
        TVector3 origin = helix.origin();  //origin of helix
        double pt = p.Perp();

        //QA Plots
        xi_selection_QAPlots(track, p, origin);

        //OK. let's cut tracks
        if(nHits <= gv_xi::cut_trk.cutNHitsGr) continue;
        if(abs(charge) != 1) continue;
        if(pt < gv_xi::cut_trk.cutPtGrEq) continue; //should be larger. like 0.15 or 0.2

        // Recording the Kaon/Pion Track for Omega/Xi Reconstruction
        if(charge == gv_xi::constants_v0.mChargeBachelor && fabs(gv_xi::cut_trk.nsigma) < gv_xi::cut_trk.cutAbsNSigma1Le) Bachelor_selection(track);
    }
}

void xi_selection_QAPlots(StPicoTrack *track, TVector3 p, TVector3 origin)
{
    //some checks.
    gv::all_hists_xi.hPt -> Fill( track->gPt() ) ; //at dca to PV, for a global track, this value is useless. anyway, the pt value is supposed to be the same anywhere.
    gv::all_hists_xi.hPtDiff -> Fill( track->gPt() - p.Perp() ) ;
    gv::all_hists_xi.hOrDiff -> Fill( (track->origin() - origin).Mag() ) ;
    gv::all_hists_xi.hPDiff -> Fill( fabs(track->gMom().Mag() - p.Mag()) ) ;
    //comments: there are difference between the values above. But they seem to be acceptably small!
}

void Bachelor_selection(StPicoTrack *track)
{
    //record the bachelor
    //fill the vector
    StPicoPhysicalHelix helix = track->helix(gv_xi::mXiDst.magn);

    //rotate transverse coordinates and momenta for background estimation
    if(gv_xi::curr.mRotate)
    {
        TVector3 p1 = helix.momentum(gv_xi::mXiDst.magn * kilogauss);  //momentum at origin
        TVector3 x1 = helix.origin();    //origin
        //p1.setX(-p1.x());
        Double_t ang = 1.0 * 3.1415926 / 5;
        p1.SetX(cos(ang) * p1.x() - sin(ang) * p1.y());
        p1.SetY(sin(ang) * p1.x() + cos(ang) * p1.y());

        Double_t x_change = cos(ang) * (x1.x() - gv_xi::xi_calc.pv.x()) - sin(ang) * (x1.y() - gv_xi::xi_calc.pv.y());
        Double_t y_change = sin(ang) * (x1.x() - gv_xi::xi_calc.pv.x()) + cos(ang) * (x1.y() - gv_xi::xi_calc.pv.y());
        x1.SetX(x_change + gv_xi::xi_calc.pv.x());
        x1.SetY(y_change + gv_xi::xi_calc.pv.y());
        StPicoPhysicalHelix helixtmp(p1, x1, gv_xi::mXiDst.magn * kilogauss, track->charge());
        helix = helixtmp;
    }

    double pathlength = helix.pathLength(gv::picoEvent->primaryVertex(), false); // do scan periods. NOTE: the default is false. this is not necessary for tracks with pt>0.15GeV/c
    TVector3 dca = helix.at(pathlength) - gv::picoEvent->primaryVertex();

    if(dca.Mag() < gv_xi::cut_trk.cutDca1GrEq) return;

    gv_xi::mDau1->mDauDcaVec.push_back(dca.Mag());
    gv_xi::mDau1->mDauVec.push_back(track);
    
    return;
}
////// End of xi_selection() ////////

////// Beginning of xi_reconstruction() ////////
void xi_reconstruction()
{
    int nXi = 0;
    //int one = 0;
    //int two = 0;
    //int three = 0;
    
    //std::cout << "total in this event: " << ((gv::details->nLambda) / 2) * gv_xi::mDau1->mDauVec.size() << endl;
    
    for(int i = 0; i < (gv::details->nLambda) / 2; i++)
    {
        //get v0 info here
        //cut them before in LambdaMaker

        if(gv::p_lambda[i].dcaglobal < gv_xi::cut_trk.cutV0DcaGrEq)continue;
        if(fabs(gv::p_lambda[i].mass - gv_xi::constants_v0.mMassV0) > gv_xi::cut_trk.cutV0MassWidthLeEq)continue;

        TVector3 xv0tmp(gv::p_lambda[i].x, gv::p_lambda[i].y, gv::p_lambda[i].z),
                 pv0tmp(gv::p_lambda[i].px, gv::p_lambda[i].py, gv::p_lambda[i].pz);
        gv_xi::xi_calc.xv0 = xv0tmp;
        gv_xi::xi_calc.pv0 = pv0tmp;
        //  cout<<"get "<<i<<"th v0 fine"<<endl;

        for(unsigned int j = 0; j < gv_xi::mDau1->mDauVec.size(); j++)
        {
            //get pion track info here
            gv_xi::xi_trk.track = gv_xi::mDau1->mDauVec[j];
            if(gv_xi::xi_trk.track->id() == gv::p_dau[2 * i + 1].trk_id)continue;

            helix_bachelor();
            if(!calc_dca_xi(j)) continue;
            //one++;
            if(!calc_ximass()) continue;
            //two++;

            //record TOF information...
            if(!tof_record()) continue;
            //three++;

            fill_mXiDst(nXi, i);
            nXi ++;

            gv::all_hists_xi.hInvMass->Fill(gv_xi::xi_calc.ximass);
        }
        //std::cout << "one = " << one << ", two = " << two << ", three = " << three << endl;
    }

    gv_xi::mXiDst.nxi = nXi;
}

void helix_bachelor()
{
    gv_xi::xi_trk.helix = gv_xi::xi_trk.track->helix(gv_xi::mXiDst.magn); //inner helix. good for dca to PV.

    if(gv_xi::curr.mRotate)
    {
        TVector3 tp1 = gv_xi::xi_trk.helix.momentum(gv_xi::mXiDst.magn * kilogauss);  //momentum at origin
        TVector3 tx1 = gv_xi::xi_trk.helix.origin();    //origin
        Double_t ang = 1.0 * 3.1415926 / 5;
        tp1.SetX(cos(ang) * tp1.x() - sin(ang) * tp1.y());
        tp1.SetY(sin(ang) * tp1.x() + cos(ang) * tp1.y());

        Double_t x_change = cos(ang) * (tx1.x() - gv_xi::xi_calc.pv.x()) - sin(ang) * (tx1.y() - gv_xi::xi_calc.pv.y());
        Double_t y_change = sin(ang) * (tx1.x() - gv_xi::xi_calc.pv.x()) + cos(ang) * (tx1.y() - gv_xi::xi_calc.pv.y());
        tx1.SetX(x_change + gv_xi::xi_calc.pv.x());
        tx1.SetY(y_change + gv_xi::xi_calc.pv.y());

        StPicoPhysicalHelix helixtmp(tp1, tx1, gv_xi::mXiDst.magn * kilogauss, gv_xi::xi_trk.track->charge());
        gv_xi::xi_trk.helix = helixtmp;
    }

    StPicoPhysicalHelix helixv0tmp(gv_xi::xi_calc.pv0 * (1.0/gv_xi::xi_calc.pv0.Perp()) * 100, gv_xi::xi_calc.xv0, gv_xi::mXiDst.magn * kilogauss, 1);
    gv_xi::xi_trk.helixv0 = helixv0tmp;
}

bool calc_dca_xi(int j)
{
    gv_xi::xi_calc.dca1 = gv_xi::mDau1->mDauDcaVec[j];
    gv_xi::xi_calc.p1 = gv_xi::xi_trk.helix.momentum(gv_xi::mXiDst.magn * kilogauss);  //momentum at origin
    //cut them before in track selection

    //StHelix dca algorithm is FAST in case of helix to staight line.
    //We will use it as default algorithm,
    //if you want to try out Long's method, switch on
    //mDcaAlgoLong in the constructor.
    if(!gv_xi::curr.mDcaAlgoLong)
    {
        //v0 is a straight line. to simulate it, we boost its transverse momentum to 100GeV/c.
        //about 700 meters transverse radii.
        //ALWAYS call helix1's pathLengths method, treat helixv0 as a parameter.
        //Otherwise, bachelor helix1's periods will be scaned, the background will be HUGE!
        //StPhysicalHelixD helixv0(pv0/pv0.perp()*100, xv0, magn*kilogauss, 1);
        //LOG_QA << 1.0/helixv0.curvature()<<endm;
        pair<double, double> tmps = gv_xi::xi_trk.helix.pathLengths(gv_xi::xi_trk.helixv0);
        TVector3 ox1 = gv_xi::xi_trk.helix.at(tmps.first);
        TVector3 ox2 = gv_xi::xi_trk.helixv0.at(tmps.second);

        gv_xi::xi_calc.dca1tov0 = (ox1 - ox2).Mag();
        gv_xi::xi_calc.xxi = (ox1 + ox2) * 0.5;
        gv_xi::xi_calc.op1 = gv_xi::xi_trk.helix.momentumAt(tmps.first, gv_xi::mXiDst.magn * kilogauss);
    }

    //LONG's two helices dca method, should also work. also need to construct v0 helix first.
    //this method work well for smaller v0 pt (less than 100GeV). but has some difference
    //in some rare cases. should not be used for serious Xi analysis.
    //the following lines are for debug only!
    if(gv_xi::curr.mDcaAlgoLong)
    {
        gv_xi::xi_calc.dca1tov0 = closestDistance(gv_xi::xi_trk.helix, gv_xi::xi_trk.helixv0, gv_xi::mXiDst.magn, gv_xi::xi_calc.pv, gv_xi::xi_calc.xxi, gv_xi::xi_calc.op1, gv_xi::xi_calc.op2);
    }

    //cut on dca1to2
    if(gv_xi::xi_calc.dca1tov0 > gv_xi::cut_trk.cutDca1to2LeEq) return false;
    //LOG_QA<<stdca1tov0<<" "<<longdca1tov0<<" "<<dca1tov0<<endm;

    return true;
}

bool calc_ximass()
{
    gv_xi::xi_calc.oe1 = sqrt(gv_xi::xi_calc.op1.Mag2() + gv_xi::constants_v0.mMassBachelor * gv_xi::constants_v0.mMassBachelor);
    gv_xi::xi_calc.oev0 = sqrt(gv_xi::xi_calc.pv0.Mag2() + gv_xi::constants_v0.mMassV0 * gv_xi::constants_v0.mMassV0);
    gv_xi::xi_calc.ximass = sqrt(gv_xi::constants_v0.mMassBachelor * gv_xi::constants_v0.mMassBachelor + gv_xi::constants_v0.mMassV0 * gv_xi::constants_v0.mMassV0 + 2.*gv_xi::xi_calc.oe1 * gv_xi::xi_calc.oev0 - 2.*gv_xi::xi_calc.op1.Dot(gv_xi::xi_calc.pv0));
    //cut on ximass
    if(fabs(gv_xi::xi_calc.ximass - gv_xi::constants_v0.mMassXi) > gv_xi::cut_trk.cutXiMassWidthLeEq) return false;

    gv_xi::xi_calc.pxi = gv_xi::xi_calc.op1 + gv_xi::xi_calc.pv0;
    gv_xi::xi_calc.xxitoPV = gv_xi::xi_calc.xxi - gv_xi::xi_calc.pv;

    //forward decay cut
    gv_xi::xi_calc.ang1 = (gv_xi::xi_calc.pxi.X() * gv_xi::xi_calc.xxitoPV.X() + gv_xi::xi_calc.pxi.Y() * gv_xi::xi_calc.xxitoPV.Y()) / gv_xi::xi_calc.pxi.Perp() / gv_xi::xi_calc.xxitoPV.Perp();

    //r dot p for xi. cut on it. should be larger than 0.
    gv_xi::xi_calc.rdotp = gv_xi::xi_calc.xxitoPV.Dot(gv_xi::xi_calc.pxi);
    if(gv_xi::xi_calc.rdotp / gv_xi::xi_calc.xxitoPV.Mag() / gv_xi::xi_calc.pxi.Mag() <= gv_xi::cut_trk.cutXirdotpGr)return false;
    if(gv_xi::xi_calc.pv0.Dot(gv_xi::xi_calc.xv0 - gv_xi::xi_calc.xxi) <= gv_xi::cut_trk.cutXirdotpGr)return false;
    gv_xi::xi_calc.pro_cross = gv_xi::xi_calc.xxitoPV.Cross(gv_xi::xi_calc.pxi);
    if(gv_xi::xi_calc.pro_cross.Mag() / gv_xi::xi_calc.pxi.Mag() / gv_xi::xi_calc.xxitoPV.Mag() > 0.15 ) return false;

    //calculate xi to PV dca. cut on dca
    StPicoPhysicalHelix helixxitmp(gv_xi::xi_calc.pxi, gv_xi::xi_calc.xxi, gv_xi::mXiDst.magn * kilogauss, gv_xi::constants_v0.mChargeXi);
    gv_xi::xi_trk.helixxi = helixxitmp;
    gv_xi::xi_calc.pathlength = gv_xi::xi_trk.helixxi.pathLength(gv_xi::xi_calc.pv, false); // do scan periods. NOTE: the default is false. this is not necessary for tracks with pt>0.15GeV/c
    gv_xi::xi_calc.dca = gv_xi::xi_trk.helixxi.at(gv_xi::xi_calc.pathlength) - gv_xi::xi_calc.pv;
    gv_xi::xi_calc.dcaxitoPV = gv_xi::xi_calc.dca.Mag();

    //cut on decay length
    gv_xi::xi_calc.xidecaylength = gv_xi::xi_calc.xxitoPV.Mag();
    if(gv_xi::xi_calc.xidecaylength < gv_xi::cut_trk.cutXiDecLenGrEq) return false;

    //cut on sinth, or theta
    gv_xi::xi_calc.sinth = (gv_xi::xi_calc.xxitoPV.Cross(gv_xi::xi_calc.pxi)).Mag() / gv_xi::xi_calc.xxitoPV.Mag() / gv_xi::xi_calc.pxi.Mag();
    double theta = atan2(gv_xi::xi_calc.sinth, gv_xi::xi_calc.rdotp / gv_xi::xi_calc.xxitoPV.Mag() / gv_xi::xi_calc.pxi.Mag()); //theta range: from 0 to pi
    if(gv_xi::xi_calc.sinth > gv_xi::cut_trk.cutXircrosspLeEq) return false;

    return true;
}

bool tof_record()
{
    int index = gv_xi::xi_trk.track->bTofPidTraitsIndex();
    //std::cout << "index = " << index << endl;
    if(index < 0) return false;

    gv_xi::xi_trk.tofflag = (gv::picoDst->btofPidTraits(index))->btofMatchFlag();
    gv_xi::xi_trk.tof = (gv::picoDst->btofPidTraits(index))->btof();
    gv_xi::xi_trk.tofpos = (gv::picoDst->btofPidTraits(index))->btofHitPos();

    gv_xi::xi_trk.tofpathlen = -999.;
    if(gv_xi::xi_trk.tofflag > 0) gv_xi::xi_trk.tofpathlen = gv_xi::xi_trk.helix.pathLength(gv_xi::xi_trk.tofpos) - gv_xi::xi_trk.helix.pathLength(gv_xi::xi_calc.xxi);

    gv_xi::xi_trk.v0pathlen = gv_xi::xi_trk.helixv0.pathLength(gv_xi::xi_calc.xv0) - gv_xi::xi_trk.helixv0.pathLength(gv_xi::xi_calc.xxi);
    gv_xi::xi_trk.xipathlen = gv_xi::xi_trk.helixxi.pathLength(gv_xi::xi_calc.xxi) - gv_xi::xi_trk.helixxi.pathLength(gv_xi::xi_calc.pv);
    
    return true;
}

void fill_mXiDst(int nXi, int i)
{
    //save the xi information (including cut variables) into StXiDst.
    gv_xi::mXiDst.v0mass[nXi] = gv::p_lambda[i].mass;
    gv_xi::mXiDst.v0x[nXi]    = gv::p_lambda[i].x;
    gv_xi::mXiDst.v0y[nXi]    = gv::p_lambda[i].y;
    gv_xi::mXiDst.v0z[nXi]    = gv::p_lambda[i].z;
    gv_xi::mXiDst.v0px[nXi]    = gv::p_lambda[i].px;
    gv_xi::mXiDst.v0py[nXi]    = gv::p_lambda[i].py;
    gv_xi::mXiDst.v0pz[nXi]    = gv::p_lambda[i].pz;
    gv_xi::mXiDst.v0dca[nXi]    = gv::p_lambda[i].dcaglobal;
    gv_xi::mXiDst.v0pathlen[nXi]  = gv_xi::xi_trk.v0pathlen;

    gv_xi::mXiDst.dau1id[nXi]    = gv::p_dau[2 * i].trk_id;
    gv_xi::mXiDst.dau1dca[nXi]    = gv::p_dau[2 * i].dcaglobal;
    gv_xi::mXiDst.dau1px[nXi]    = gv::p_dau[2 * i].px;
    gv_xi::mXiDst.dau1py[nXi]    = gv::p_dau[2 * i].py;
    gv_xi::mXiDst.dau1pz[nXi]    = gv::p_dau[2 * i].pz;

    gv_xi::mXiDst.dau2id[nXi]    = gv::p_dau[2 * i + 1].trk_id;
    gv_xi::mXiDst.dau2dca[nXi]    = gv::p_dau[2 * i + 1].dcaglobal;
    gv_xi::mXiDst.dau2px[nXi]    = gv::p_dau[2 * i + 1].px;
    gv_xi::mXiDst.dau2py[nXi]    = gv::p_dau[2 * i + 1].py;
    gv_xi::mXiDst.dau2pz[nXi]    = gv::p_dau[2 * i + 1].pz;

    gv_xi::mXiDst.bachid[nXi] = gv_xi::xi_trk.track->id();
    gv_xi::mXiDst.bachdca[nXi]  = gv_xi::xi_calc.dca1;
    gv_xi::mXiDst.bachdca2d[nXi]  = gv_xi::xi_trk.helix.geometricSignedDistance(gv_xi::xi_calc.pv.X(), gv_xi::xi_calc.pv.Y());
    gv_xi::mXiDst.bachnhits[nXi] = gv_xi::xi_trk.track->nHits();
    gv_xi::mXiDst.bachdedx[nXi] = gv_xi::xi_trk.track->dEdx();
    gv_xi::mXiDst.bachnsigma[nXi] = (gv_xi::curr.mXiType == kXi || gv_xi::curr.mXiType == kAntiXi) ? gv_xi::xi_trk.track->nSigmaPion() : gv_xi::xi_trk.track->nSigmaKaon();
    gv_xi::mXiDst.bacheta[nXi]  = gv_xi::xi_calc.op1.PseudoRapidity();
    gv_xi::mXiDst.bachpt[nXi]   = gv_xi::xi_calc.op1.Perp();
    gv_xi::mXiDst.bachpx[nXi]   = gv_xi::xi_calc.op1.X();
    gv_xi::mXiDst.bachpy[nXi]   = gv_xi::xi_calc.op1.Y();
    gv_xi::mXiDst.bachpz[nXi]   = gv_xi::xi_calc.op1.Z();
    gv_xi::mXiDst.bachtpc[nXi]  = gv_xi::xi_trk.track->nHitsFit();
    gv_xi::mXiDst.bachssd[nXi]  = gv_xi::xi_trk.track->nHitsFit();
    gv_xi::mXiDst.bachsvt[nXi]  = gv_xi::xi_trk.track->nHitsFit();
    gv_xi::mXiDst.bachtofflag[nXi] = gv_xi::xi_trk.tofflag;
    gv_xi::mXiDst.bachtof[nXi] = gv_xi::xi_trk.tof;
    gv_xi::mXiDst.bachpathlen[nXi] = gv_xi::xi_trk.tofpathlen;


    gv_xi::mXiDst.dcav0tobach[nXi] = gv_xi::xi_calc.dca1tov0;

    gv_xi::mXiDst.ximass[nXi]   = gv_xi::xi_calc.ximass;
    gv_xi::mXiDst.xidca[nXi]    = gv_xi::xi_calc.dcaxitoPV;
    gv_xi::mXiDst.xidca2d[nXi]    = gv_xi::xi_trk.helixxi.geometricSignedDistance(gv_xi::xi_calc.pv.X(), gv_xi::xi_calc.pv.Y());
    gv_xi::mXiDst.xipt[nXi]   = gv_xi::xi_calc.pxi.Perp();
    gv_xi::mXiDst.xirapidity[nXi]    = log( (sqrt(gv_xi::xi_calc.ximass * gv_xi::xi_calc.ximass + gv_xi::xi_calc.pxi.Mag2()) + gv_xi::xi_calc.pxi.Z()) / sqrt(gv_xi::xi_calc.ximass * gv_xi::xi_calc.ximass + gv_xi::xi_calc.pxi.Perp2()));
    gv_xi::mXiDst.xieta[nXi]    = 0.5 * log( (gv_xi::xi_calc.pxi.Mag() + gv_xi::xi_calc.pxi.Z()) / (gv_xi::xi_calc.pxi.Mag() - gv_xi::xi_calc.pxi.Z()) );
    gv_xi::mXiDst.xix[nXi]      = gv_xi::xi_calc.xxi.X();
    gv_xi::mXiDst.xiy[nXi]      = gv_xi::xi_calc.xxi.Y();
    gv_xi::mXiDst.xiz[nXi]      = gv_xi::xi_calc.xxi.Z();
    gv_xi::mXiDst.xipx[nXi]     = gv_xi::xi_calc.pxi.X();
    gv_xi::mXiDst.xipy[nXi]     = gv_xi::xi_calc.pxi.Y();
    gv_xi::mXiDst.xipz[nXi]     = gv_xi::xi_calc.pxi.Z();
    gv_xi::mXiDst.xideclen[nXi] = gv_xi::xi_calc.xidecaylength;
    gv_xi::mXiDst.xisinth[nXi]  = gv_xi::xi_calc.sinth;
    gv_xi::mXiDst.xipathlen[nXi]  = gv_xi::xi_trk.xipathlen;
}
