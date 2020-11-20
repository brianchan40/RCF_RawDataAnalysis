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
#include <TChain.h>
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

/// PicoDst headers
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/particle.h"
#include "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/event.h"

void simple_data_extract()
{
    TFile *outputfile = new TFile("plam_massonly.root", "recreate");
    int n_files = 28922;
    //int n_files = 10;

    TH1F *V0Mass_cent[9];
    for(int j = 0; j < 9; j++)
    {
        V0Mass_cent[j] = new TH1F(TString::Format("V0Mass_%d", j), TString::Format("V0 Inv. Mass for Cent-%d", j), 200, 1.115684 - 0.07, 1.115684 + 0.07);
    }

    std::vector<particle> *p_lambda = new vector<particle>;
    event *details = new event();

    TChain *corr_tree;
    char fname[200];

    int n_list = ceil((double)n_files / 100.0);

    for(int list_id = 0; list_id < n_list; list_id++)
    {
        int curr_n_files = 0;

        if(list_id != (n_list-1)) curr_n_files = list_id*100 + 100;
        else curr_n_files = n_files;

        cout << "Processing list == " << (list_id*100) << " to " << curr_n_files << " == out of " << n_files << ".\n";
        corr_tree = new TChain("corr_tree");

        for(int i = (list_id*100); i < curr_n_files; i++)
        {
            //if((i + 1) % 100 == 0) cout << "Processing entry == " << i + 1 << " == out of " << n_files << ".\n";

            sprintf(fname, "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/Results_lam_18/sched265153F42AB16B85EBA22DFB61D45A06_%d_tree.root", i);
            corr_tree->AddFile(fname);
        }

        //cout << "HERE?" << endl;

        corr_tree->SetBranchAddress("p_lambda", &p_lambda);
        corr_tree->SetBranchAddress("details", &details);

        Int_t nentries = corr_tree->GetEntries();
        cout << "nentries = " << nentries << endl;

        for(int l = 0; l < nentries; l++)
        {
            corr_tree->GetEntry(l);
            if((l + 1) % 1000 == 0) cout << "Processing entry == " << l + 1 << " == out of " << nentries << ".\n";

            int Centrality = details->cent;
            int n_Lambdas = (details->nLambda) / 2;

            //cout << "n_Lambdas = " << n_Lambdas << endl;

            for(int k = 0; k < n_Lambdas; k++)
            {
                V0Mass_cent[Centrality]->Fill(p_lambda->at(k).mass);
            }
        }

        p_lambda->clear();
        corr_tree->Clear();
    }

    //p_lambda->clear();

    
    outputfile->cd();
    for(int m = 0; m < 9; m++)
    {
        V0Mass_cent[m]->Write();
    }
    outputfile->Close();

}