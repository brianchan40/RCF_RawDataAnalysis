using namespace std;

/// C++ headers

#include "stdio.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include <string>
#include <cstring>
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

/// PicoDst headers
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoDstReader.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoDst.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoEvent.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoTrack.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoBTofHit.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoBTowHit.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoHelix.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoPhysicalHelix.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoBTofPidTraits.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StRefMultCorr/StRefMultCorr.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StRefMultCorr/CentralityMaker.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/particle.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/event.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoEmcTrigger.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StPicoTrackCovMatrix.h"

/// Analysis headers
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StV0Type.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StDcaService.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/structs.h"
#include "./namespaces/gv_gamma.h"
#include "./v2_module.h"

//void __attribute__((constructor)) LoadLib();
/////////////////////////////main program starts here/////////////////////////////////////
void v2_module(int cen = 1, int opt_weight = 1, const Char_t *inFile = "./st_physics_12160020_raw_1010001.picoDst.root", const TString JobIDName = "1234")      //main_function
{
    //TChain *corr_tree = new TChain ("tree");
    TChain *corr_tree = new TChain ("tr");
    std::ifstream fin(inFile);
    string line1;
    char line[200];
    while(getline(fin, line1)) //loop wiill run till end of file
    {
        //cout << "line1 = " << line1 << endl;
        strcpy(line, line1.c_str());
        corr_tree->AddFile(line);
    }
    fin.close();

    // corr_tree->SetBranchAddress("event", &event);
    // corr_tree->SetBranchAddress("refmult3", &refmult3);
    // corr_tree->SetBranchAddress("pid", &pid);
    // corr_tree->SetBranchAddress("px", &px);
    // corr_tree->SetBranchAddress("py", &py);
    // corr_tree->SetBranchAddress("pz", &pz);
    // corr_tree->SetBranchAddress("refmult2", &refmult2);
    // corr_tree->SetBranchAddress("imp", &imp);
    // corr_tree->SetBranchAddress("refmult", &refmult);
    // corr_tree->SetBranchAddress("npp", &npp);
    // corr_tree->SetBranchAddress("npt", &npt);
    // corr_tree->SetBranchAddress("nesp", &nesp);
    // corr_tree->SetBranchAddress("ninesp", &ninesp);
    // corr_tree->SetBranchAddress("nest", &nest);
    // corr_tree->SetBranchAddress("ninest", &ninest);

    // cout << "before assigment" << endl;
    corr_tree->SetBranchAddress("Event", &event);
    // corr_tree->SetBranchAddress("refmult3", &refmult3);
    corr_tree->SetBranchAddress("PID", pid);
    corr_tree->SetBranchAddress("Px", px);
    corr_tree->SetBranchAddress("Py", py);
    corr_tree->SetBranchAddress("Pz", pz);
    // corr_tree->SetBranchAddress("refmult2", &refmult2);
    corr_tree->SetBranchAddress("Imp", &imp);
    corr_tree->SetBranchAddress("Mult", &refmult);
    corr_tree->SetBranchAddress("Npartp", &npp);
    corr_tree->SetBranchAddress("Npartt", &npt);
    corr_tree->SetBranchAddress("Nesp", &nesp);
    corr_tree->SetBranchAddress("Ninesp", &ninesp);
    corr_tree->SetBranchAddress("Nest", &nest);
    corr_tree->SetBranchAddress("Ninest", &ninest);

    Int_t nentries = corr_tree->GetEntries();
    // cout << "nentries = " << nentries << endl;

    //loop through events
    for(int i = 0; i < nentries; i++)
    {
        if((i + 1) % 1000 == 0) cout << "Processing entry == " << i + 1 << " == out of " << nentries << ".\n";

        corr_tree->GetEntry(i);
        //cout << "entries?" << endl;
        if((npp + npt) < 3) continue;
        // if(imp > 20) continue;

        Fill_EventQAPlots();

        // for(int j = 0; j < px->size(); j++)
        for(int j = 0; j < refmult; j++)
        {
            // if(abs((int)pid->at(j)) == 42) continue;
            if(abs((int)pid[j]) == 42) continue;
            
            Get_BasicTrackInfo(j);
            
            Fill_TrackQAPlots(j);
            
        }

    }

    Write_QAPlots(JobIDName);
}

void Fill_EventQAPlots()
{
    
    imp_dist->Fill(imp);
    
    refmult_dist->Fill(refmult);
    // refmult2_dist->Fill(refmult2);
    
    // refmult3_dist->Fill(refmult3);
    
    // pmult_dist->Fill(px->size());
    imp_vs_refmult->Fill(imp, refmult);
    // imp_vs_refmult2->Fill(imp, refmult2);
    // imp_vs_refmult3->Fill(imp, refmult3);
    // imp_vs_pmult->Fill(imp, px->size());
    // pmult_vs_refmult->Fill(px->size(), refmult);
    // pmult_vs_refmult2->Fill(px->size(), refmult2);
    npp_dist->Fill(npp);
    npt_dist->Fill(npt);
    nesp_dist->Fill(nesp);
    ninesp_dist->Fill(ninesp);
    nest_dist->Fill(nest);
    ninest_dist->Fill(ninest);
    diff_npp_dist->Fill(npp-(nesp+ninesp));
    diff_npt_dist->Fill(npt-(nest+ninest));
}

void Get_BasicTrackInfo(int j)
{
    // gv_gamma::trk_mom->SetXYZ(px->at(j), py->at(j), pz->at(j));
    gv_gamma::trk_mom->SetXYZ(px[j], py[j], pz[j]);
    Pt    = gv_gamma::trk_mom->Pt();
    if(Pt != 0) Eta   = gv_gamma::trk_mom->Eta();
    Phi   = gv_gamma::trk_mom->Phi();
    Mom_Mag = gv_gamma::trk_mom->Mag();

    // p_info = db->GetParticle((int)pid->at(j));
    p_info = db->GetParticle((int)pid[j]);

    Charge = p_info->Charge() / 3;
}

void Fill_TrackQAPlots(int j)
{
    // px_dist->Fill(px->at(j));
    // py_dist->Fill(py->at(j));
    // pz_dist->Fill(pz->at(j));
    px_dist->Fill(px[j]);
    py_dist->Fill(py[j]);
    pz_dist->Fill(pz[j]);
    pt_dist->Fill(Pt);
    p_dist->Fill(Mom_Mag);
    eta_dist->Fill(Eta);
    phi_dist->Fill(Phi);
    charge_dist->Fill(Charge);
    pid_dist->Fill(pid[j]);
    // pid_dist->Fill(pid->at(j));

    // if((int)pid->at(j) == 3122) nLambda_dist->Fill(imp);
    // if((int)pid->at(j) == -3122) nAntiLambda_dist->Fill(imp);
    if((int)pid[j] == 3122) nLambda_dist->Fill(imp);
    if((int)pid[j] == -3122) nAntiLambda_dist->Fill(imp);
}

void Write_QAPlots(TString JobIDName)
{

    TString filename = JobIDName;
    filename.Append(".root");
    filename.Prepend("sched");
    TFile f(filename, "recreate");
    f.cd();

    imp_dist->Write();
    refmult_dist->Write();
    // refmult2_dist->Write();
    // refmult3_dist->Write();
    pmult_dist->Write();
    imp_vs_refmult->Write();
    // imp_vs_refmult2->Write();
    // imp_vs_refmult3->Write();
    imp_vs_pmult->Write();
    pmult_vs_refmult->Write();
    // pmult_vs_refmult2->Write();
    npp_dist->Write();
    npt_dist->Write();
    nesp_dist->Write();
    ninesp_dist->Write();
    nest_dist->Write();
    ninest_dist->Write();
    diff_npp_dist->Write();
    diff_npt_dist->Write();

    px_dist->Write();
    py_dist->Write();
    pz_dist->Write();
    pt_dist->Write();
    p_dist->Write();
    eta_dist->Write();
    phi_dist->Write();
    charge_dist->Write();
    pid_dist->Write();

    nLambda_dist->Write();
    nAntiLambda_dist->Write();

    f.Close();
}
