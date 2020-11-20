/*
 * Author: Brian Chan
 * Date: February 27, 2019
 *
 * ProtonMaker tries to select primary protons for further analysis
 *
 *
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
#include "./particle.h"
#include "./event.h"

/// Analysis headers
#include "./StV0Type.h"
#include "./StDcaService.h"
#include "../structs.h"
#include "./ProtonMaker.h"
#include "./gv.h"

using namespace std;

void init_Proton(){
    cout << "init_Proton()" << endl;
    gv::all_hists.primary_p_perp = new TH1F( "primary_p_perp", "Transverse Momentum of Primary Protons", 100, 0, 10 ) ;
    gv::all_hists.primary_p_mag = new TH1F( "primary_p_mag", "Momentum of Primary Protons", 100, 0, 10 ) ;
    gv::all_hists.btofylocal_dist = new TH1F( "btofylocal_dist", "Distribution of BtofYLocal", 200, -10, 10 ) ;
    gv::all_hists.mass2_p = new TH1F( "mass2_p", "Distribution of Mass2", 100, 0, 5 ) ;
}

void make_Proton(){

    Int_t nTracks = gv::picoDst->numberOfTracks();

    for(Int_t iTrk = 0; iTrk < nTracks; iTrk++)
    {
        gv::track_forp = gv::picoDst->track(iTrk);
        if(!gv::track_forp) continue;
        
        if(!proton_cuts_decision()) continue;
        
        store_proton();
    }
}

bool proton_cuts_decision()
{
    if(!gv::track_forp->isPrimary()) return false;

    TVector3 primary_p = gv::track_forp->pMom();
    gv::all_hists.primary_p_perp->Fill(primary_p.Perp());
    gv::all_hists.primary_p_mag->Fill(primary_p.Mag());
    if((primary_p.Perp() < 0.2) || (primary_p.Mag() > 2.0) || (gv::track_forp->nHitsDedx() < 15)) return false;

    int index = gv::track_forp->bTofPidTraitsIndex();
    if(index < 0) return false;
    int tofflag = (gv::picoDst->btofPidTraits(index))->btofMatchFlag();
    double tof = (gv::picoDst->btofPidTraits(index))->btof();
    double BtofYLocal = (gv::picoDst->btofPidTraits(index))->btofYLocal();
    gv::all_hists.btofylocal_dist->Fill(BtofYLocal);
    if((tofflag < 1) || (tof <= 0) || (BtofYLocal < -1.8) || (BtofYLocal > 1.8)) return false;

    double nsigmaproton = gv::track_forp->nSigmaProton();
    TVector3 dca = gv::track_forp->gDCA(gv::picoEvent->primaryVertex());
    if((nsigmaproton < -2) || (nsigmaproton > 2) || (dca.Mag() > 1)) return false;

    StPicoPhysicalHelix helix = gv::track_forp->helix(gv::details->Magn);
    TVector3 pp = helix.momentum(gv::details->Magn * kilogauss);
    double beta = (gv::picoDst->btofPidTraits(index))->btofBeta();
    float mass2 = pp.Mag2() * (1.0 / beta / beta - 1.0);
    gv::all_hists.mass2_p->Fill(mass2);
    if((mass2 < 0.8) || ((mass2 > 1.0))) return false;

    return true;
}

void store_proton()
{
    TVector3 primary_p = gv::track_forp->pMom();
    TVector3 dca = gv::track_forp->gDCA(gv::picoEvent->primaryVertex());

    gv::p_proton[gv::details->n_Proton].px = primary_p.x();
    gv::p_proton[gv::details->n_Proton].py = primary_p.y();
    gv::p_proton[gv::details->n_Proton].pz = primary_p.z();
    gv::p_proton[gv::details->n_Proton].dcaglobal = dca.Mag();

    gv::p_proton[gv::details->n_Proton].Charge = gv::track_forp->charge();
    
    gv::p_proton[gv::details->n_Proton].trk_id = gv::track_forp->id();

    gv::details->n_Proton++;
}
