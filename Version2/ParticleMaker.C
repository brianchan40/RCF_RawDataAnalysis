/*
 * Author: Brian Chan
 * Date: February 27, 2019
 *
 * ParticleMaker tries to reconstruct Particles for further analysis
 * Current Particles:
 * Lambda
 * Ks
 * Proton
 * Xi
 * Omega
 *
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
//#include "TVector3.hh"
#include <assert.h>
#include <TTree.h>
#include <vector>
#include <TLorentzVector.h>
#include <iterator>
#include <sstream>
#include <fstream>

/// PicoDst headers
#include "./StRoot/StPicoDstReader.h"
#include "./StRoot/StPicoDst.h"
#include "./StRoot/StPicoEvent.h"
#include "./StRoot/StPicoTrack.h"
#include "./StRoot/StPicoBTofHit.h"
#include "./StRoot/StPicoBTowHit.h"
#include "./StRoot/StPicoHelix.h"
#include "./StRoot/StPicoPhysicalHelix.h"
#include "./StRoot/StPicoBTofPidTraits.h"
#include "./StRoot/StRefMultCorr/StRefMultCorr.h"
#include "./StRoot/StRefMultCorr/CentralityMaker.h"
#include "./StRoot/particle.h"
#include "./StRoot/event.h"

/// Analysis headers
/*#include "../StV0Type.h"
#include "../StDcaService.h"*/
#include "./StRoot/StV0Type.h"
#include "./StRoot/StXiType.h"
#include "./StRoot/StDcaService.h"
#include "./structs.h"
#include "./StRoot/gv.h"

// Particle Makers
#include "./StRoot/LambdaMaker.h"
#include "./StRoot/ProtonMaker.h"
#include "./StRoot/Run_by_RunQA.h"
#include "./StRoot/XiMaker.h"

//#include "./Gamma_112_module.C"


/// Load libraries
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0) // 393215
R__LOAD_LIBRARY(. / StRoot / libStPicoDst)
#endif

#ifdef __MAKECINT__
#pragma link C++ class vector<particle>+;
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Global Variables
TFile *outputfile;
TFile *outputfile_xi;
TFile *treefile_xi;

event *details_tmp = new event();

std::vector<particle> *p_lambda_tmp = new vector<particle>;
std::vector<particle> *p_all_tmp = new vector<particle>;
std::vector<particle> *p_dau_tmp = new vector<particle>;
std::vector<particle> *p_proton_tmp = new vector<particle>;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Supporting Functions
void init_output(TString JobIdName);
void write_hist();
void init_tree();
void write_hist_run_by_run();
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void ParticleMaker(const int cen = 1, const int opt_weight = 1, const Char_t *inFile = "./st_physics_12160020_raw_1010001.picoDst.root", const TString JobIdName = "1234")
{
    std::cout << "Hi! Lets do some physics, Master!" << std::endl;

    Char_t *InputFileList = "./new.list";
    std::ofstream fout(InputFileList); //create a file to write
    std::ifstream fin(inFile);

    string line;
    int file_count = 0;
    while(getline(fin, line)) //loop wiill run till end of file
    {
        //cout << line << endl;
        std::istringstream iss(line);
        std::string token;
        while(std::getline(iss, token, ' '))   // but we can specify a different one
        {
            if(file_count % 2 == 0)
            {
                fout << token << "\n";     //writing data to file
                cout << token << endl;
            }
            file_count++;
        }
    }

    fin.close();
    fout.close();


    StPicoDstReader *picoReader = new StPicoDstReader(InputFileList);
    //StPicoDstReader *picoReader = new StPicoDstReader(inFile);
    picoReader->Init();

    std::cout << "Now I know what to read, Master!" << std::endl;

    if( !picoReader->chain() ) std::cout << "No chain has been found." << std::endl;

    TTree *picoTree = picoReader->tree();

    Long64_t eventsInTree = picoTree->GetEntries();
    std::cout << "eventsInTree: "  << eventsInTree << std::endl;
    Long64_t events2read = picoReader->chain()->GetEntries();

    std::cout << "Number of events to read: " << events2read << std::endl;


    ////////////////////////////////////////////// Particle Reconstruction //////////////////////////////////////////////
    /////// Initialization ////////
    init_output(JobIdName); // Create output file
    init_tree();
    init_Lambda(); // from LambdaMaker.cpp, initializes LambdaMaker and ProtonMaker
    //std::cout << "init_Lambda()" << std::endl;
    init_Proton();
    //init_xi(); // from XiMaker.cpp, initializes XiMaker
    //init_Gamma112(cen);

    /////// ANALYSIS /////////

    /// Loop over events
    for(Long64_t iEvent = 0; iEvent < events2read; iEvent++)
    {
        if(iEvent % 100 == 0) std::cout << "Working on event #[" << (iEvent + 1) << "/" << events2read << "]" << std::endl;

        Bool_t readEvent = picoReader->readPicoEvent(iEvent);
        if( !readEvent )
        {
            std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
            break;
        }
        //std::cout << "gv::check = " << gv::check << endl;
        gv::check = 2;
        //std::cout << "gv::check = " << gv::check << endl;
        /// Retrieve picoDst
        gv::picoDst = picoReader->picoDst();
        //std::cout << "Read PicoDst" << endl;

        /// Retrieve event information
        gv::picoEvent = gv::picoDst->event();
        //std::cout << "Read PicoEvent" << endl;
        if( !gv::picoEvent )
        {
            std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
            break;
        }

        //std::cout << "1" << std::endl;
        if(!make_Lambda(cen)) continue;
        //std::cout << "Made Lambda" << endl; 
        if(gv::details->nLambda == 0)
        {
            if(p_lambda_tmp->size() != 0) p_lambda_tmp->clear();
            if(p_all_tmp->size() != 0) p_all_tmp->clear();
            if(p_dau_tmp->size() != 0) p_dau_tmp->clear();
            if(p_proton_tmp->size() != 0) p_proton_tmp->clear();
            continue;
        }

        //std::cout << "gv::check = " << gv::check << endl;
        //std::cout << "2" << std::endl;
        make_Proton();

        //std::cout << "gv::check = " << gv::check << endl;
        //std::cout << "3" << std::endl;
        //std::cout << "p_all[0].px = " << p_all[0].Px() << endl;
        details_tmp = gv::details;
        //details = gv::details;

        for(int pi = 0; pi < gv::details->num_trk; pi++)
        {
            p_all_tmp->push_back(gv::p_all[pi]);
        }
        for(int pl = 0; pl < (gv::details->nLambda / 2); pl++)
        {
            p_lambda_tmp->push_back(gv::p_lambda[pl]);
        }
        for(int pd = 0; pd < (gv::details->nLambda); pd++)
        {
            p_dau_tmp->push_back(gv::p_dau[pd]);
        }
        for(int pp = 0; pp < gv::details->n_Proton; pp++)
        {
            p_proton_tmp->push_back(gv::p_proton[pp]);
        }

        //Gamma_112_module(cen, opt_weight, JobIdName);

        //std::cout << details_tmp->cent << endl;

        //make_Xi();
        gv::corr_tree->Fill();
        //gv::mXiTree->Fill();

        p_lambda_tmp->clear();
        p_all_tmp->clear();
        p_dau_tmp->clear();
        p_proton_tmp->clear();
    }

    cout << "before write_hist()" << endl;
    write_hist();
    //finish_Gamma112(cen, opt_weight, JobIdName);

    treefile_xi->cd();
    gv::corr_tree->Write();
    //gv::mXiTree->CloneTree()->Write();
    treefile_xi->Close();
    ////////////////////////////////////////////// Particle Reconstruction //////////////////////////////////////////////

    ////////////////////////////////////////////// Run by Run QA //////////////////////////////////////////////
    /////// Initialization ////////
    /*init_output(JobIdName); // Create output file
    init_hist_run_by_run();

    /////// ANALYSIS /////////

    /// Loop over events
    for(Long64_t iEvent = 0; iEvent < events2read; iEvent++)
    {
        if(iEvent % 100 == 0) std::cout << "Working on event #[" << (iEvent + 1) << "/" << events2read << "]" << std::endl;

        Bool_t readEvent = picoReader->readPicoEvent(iEvent);
        if( !readEvent )
        {
            std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
            break;
        }
        
        /// Retrieve picoDst
        gv::picoDst = picoReader->picoDst();

        /// Retrieve event information
        gv::picoEvent = gv::picoDst->event();
        if( !gv::picoEvent )
        {
            std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
            break;
        }

        //std::cout << "1" << std::endl;
        if(!Run_by_Run()) continue;        
    }

    cout << "before write_hist_run_by_run()" << endl;
    write_hist_run_by_run();
    
    treefile_xi->cd();
    vz_dist->Write();
    treefile_xi->Close();*/
    ////////////////////////////////////////////// Run by Run QA //////////////////////////////////////////////


    picoReader->Finish();

    std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!" << std::endl;
}

//creating output files
void init_output(TString JobIdName)
{
    // Create output file for histograms
    TString Name = "sched";
    Name.Append(JobIdName);
    Name.Append("_plam.root") ;
    outputfile = new TFile(Name, "recreate");

    std::cout << Name << std::endl;

    // TString Name3 = "sched";
    // Name3.Append(JobIdName);
    // Name3.Append("_xi.root") ;
    // outputfile_xi = new TFile(Name3, "recreate");

    // std::cout << Name3 << std::endl;

    // Create output file for tree;
    TString Name4 = "sched";
    Name4.Append(JobIdName);
    Name4.Append("_tree.root");
    treefile_xi = new TFile(Name4, "recreate");
}

void write_hist()
{
    /// Filling the File
    outputfile->cd();
    cout << "output file in" << endl;
    gv::all_hists.hVertexZ->Write();
    gv::all_hists.hNRefMult->Write();
    gv::all_hists.hSelectNRefMultM0->Write();
    gv::all_hists.hSelectNRefMultM1->Write();
    gv::all_hists.hSelectNRefMultM2->Write();
    gv::all_hists.hSelectNRefMultFtpcE->Write();
    cout << "here 1" << endl;
    gv::all_hists.primary_vertex->Write();
    gv::all_hists.hPtRaw->Write();
    gv::all_hists.charge_dist->Write();
    gv::all_hists.hPtDiff->Write();
    gv::all_hists.hOrDiff->Write();
    gv::all_hists.hPDiff->Write();
    gv::all_hists.nHitsFitDist->Write();
    gv::all_hists.mom_track_gpt->Write();
    gv::all_hists.DCA_Dist->Write();
    cout << "here 2" << endl;
    gv::all_hists.hInvMass->Write();
    cout << "here 3" << endl;
    gv::all_hists.lambda_pt->Write();
    gv::all_hists.v0decaylengthQA->Write();
    gv::all_hists.dcav0toPVQA->Write();
    gv::all_hists.entireV0Mass->Write();
    gv::all_hists.dca1to2_Dist->Write();
    cout << "here 4" << endl;
    gv::all_hists.transversevertexmag->Write();
    gv::all_hists.RapidityQA->Write();
    gv::all_hists.pseudoRapidityQA1->Write();
    gv::all_hists.pseudoRapidityQA2->Write();
    gv::all_hists.lambdarapidity->Write();
    cout << "here 5" << endl;
    gv::all_hists.cent_check_9->Write();
    gv::all_hists.cent_check_16->Write();
    gv::all_hists.trigger_IDs->Write();
    gv::all_hists.hNSigmaPion->Write();
    gv::all_hists.hNSigmaProton->Write();
    gv::all_hists.hNSigmaKaon->Write();
    cout << "here 6" << endl;
    gv::all_hists.primary_p_perp->Write();
    cout << "here 7" << endl;
    gv::all_hists.primary_p_mag->Write();
    cout << "here 8" << endl;
    gv::all_hists.btofylocal_dist->Write();
    cout << "here 9" << endl;
    gv::all_hists.mass2_p->Write();

    cout << "before the 9?" << endl;
    for(int cen_id = 0; cen_id < 9; cen_id++)
    {
        for(int j = 0; j < 17; j++)
        {
            //cout << "cen_id = " << cen_id << endl;
            gv::all_hists.V0Mass_cent[cen_id][j]->Write();
        }
    }
    outputfile->Close();

    // outputfile_xi->cd();
    // gv::all_hists_xi.hVertexZ->Write();
    // gv::all_hists_xi.hPt->Write();
    // gv::all_hists_xi.hNPrimVertex->Write();
    // gv::all_hists_xi.hNRefMult->Write();
    // gv::all_hists_xi.hPtDiff->Write();
    // gv::all_hists_xi.hOrDiff->Write();
    // gv::all_hists_xi.hPDiff->Write();
    // gv::all_hists_xi.hDcaDiff->Write();
    // gv::all_hists_xi.hInvMass->Write();
    // gv::all_hists_xi.hInvMassSTD->Write();
    // outputfile_xi->Close();

}

void write_hist_run_by_run(){
    outputfile->cd();
    vz_dist->Write();
    vxy_dist->Write();
    refmult_dist->Write();
    mean_pt->Write();
    mean_eta->Write();
    mean_phi->Write();
    mean_dca->Write();
    mean_nHitsFit->Write();
    mean_nHitsDedx->Write();
    outputfile->Close();
}

void init_tree()
{
    gv::corr_tree->Branch("details", &details_tmp);
    // gv::corr_tree->Branch("p_lambda", p_lambda_tmp, "px[5000]/F:py[5000]/F:pz[5000]/F:x[5000]/F:y[5000]/F:z[5000]/F:Charge[5000]/F:TOFflag[5000]/I:dcaglobal[5000]/F:prim[5000]/I:nSigmaProton[5000]/F:isTofTrack[5000]/I:mass[5000]/D:trk_id[5000]/I:ToF[5000]/D:BTOFYLocal[5000]/D:Mass2[5000]/F");
    // gv::corr_tree->Branch("p_proton", p_proton_tmp, "px[5000]/F:py[5000]/F:pz[5000]/F:x[5000]/F:y[5000]/F:z[5000]/F:Charge[5000]/F:TOFflag[5000]/I:dcaglobal[5000]/F:prim[5000]/I:nSigmaProton[5000]/F:isTofTrack[5000]/I:mass[5000]/D:trk_id[5000]/I:ToF[5000]/D:BTOFYLocal[5000]/D:Mass2[5000]/F");
    // gv::corr_tree->Branch("p_all", p_all_tmp, "px[5000]/F:py[5000]/F:pz[5000]/F:x[5000]/F:y[5000]/F:z[5000]/F:Charge[5000]/F:TOFflag[5000]/I:dcaglobal[5000]/F:prim[5000]/I:nSigmaProton[5000]/F:isTofTrack[5000]/I:mass[5000]/D:trk_id[5000]/I:ToF[5000]/D:BTOFYLocal[5000]/D:Mass2[5000]/F");
    // gv::corr_tree->Branch("p_dau", p_dau_tmp, "px[5000]/F:py[5000]/F:pz[5000]/F:x[5000]/F:y[5000]/F:z[5000]/F:Charge[5000]/F:TOFflag[5000]/I:dcaglobal[5000]/F:prim[5000]/I:nSigmaProton[5000]/F:isTofTrack[5000]/I:mass[5000]/D:trk_id[5000]/I:ToF[5000]/D:BTOFYLocal[5000]/D:Mass2[5000]/F");
    gv::corr_tree->Branch("p_lambda", &p_lambda_tmp);
    gv::corr_tree->Branch("p_proton", &p_proton_tmp);
    gv::corr_tree->Branch("p_all", &p_all_tmp);
    gv::corr_tree->Branch("p_dau", &p_dau_tmp);
}
