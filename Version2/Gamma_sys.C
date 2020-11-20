using namespace std;

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
#include "TLorentzVector.h"

/// PicoDst headers
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StPicoEvent/StPicoDstReader.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StPicoEvent/StPicoDst.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StPicoEvent/StPicoEvent.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StPicoEvent/StPicoTrack.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StPicoEvent/StPicoBTofHit.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StPicoEvent/StPicoBTowHit.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StEpdUtil/StEpdEpFinder.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StRefMultCorr/StRefMultCorr.h"
#include "/star/data01/pwg/rexwg/27GeV_2018/StRoot/StRefMultCorr/CentralityMaker.h"

//class StRefMultCorr;
//class CentralityMaker;

const float PI = TMath::Pi();
const float MM = 2 / PI;
const int Phibin = 80;
const int order = 5;
const float mpi = 0.13957;

const int opt_sys = 0;
const int opt_HighOrder = 0;
const int opt_useEPD = 0;  //0 is TPC EP, 1 is EPD, 2 is BBC, 3 is ZDC
const int opt_boost = 0;
const int opt_ESE = 0; //0 means default, 1 will select 0.4 < q < 1.6
const int nHar = 3; //2nd order EP

const float Eta_EP_Cut = 0.55;  //particles in the event plane
const float Eta_ESE_Cut = 0.5;  //particles for event shape engineering

const float pt_asso_up = 2.0;//1;
const float pt_asso_lo = 0.15;//0.15;
const float EtaCut_asso = 1;
const float DcaCut_asso = 2;         //2
const int nFitHits_asso = 15;   //15, number of fit hits

const float pt_trig_up = 2;
const float pt_trig_lo = .2;
const float EtaCut_trig = 1;
const float DcaCut_trig = 2;
const int nFitHits_trig = 15;
const int nFitRatio_trig_lo = 0.52;  //0.52
const int nFitRatio_trig_up = 1.05;  //1.05

const Float_t Vz_cut = 70;      //30, Vz
const Float_t Vr_cut = 2;   //2, Vr
const float Vz_diff = 4;        //5, difference between Vz{TPC} and Vz{VPD}
const int run_sta = 71000;
const int run_end = 168000;
const int cenDef[9] = {6, 13, 25, 44, 72, 113, 169, 245, 295};  //27 GeV 2018

//bad runs
const int Nrun_MB = 20;
const int bad_Ref_day3[Nrun_MB] = {132063, 133018, 137003, 137050, 158013 //RefMult
                                   , 143013, 143014, 143015, 143016, 143017 //TOF
                                   , 137051, 157013 //NPA
                                   , 131050, 141019
                                   , 137052, 138010 //TPC sin2
                                   , 144019, 140016, 142034 //EPD east and west, cos2 and sin2
                                   , 147007 //BBC
                                  };
//new efficiency
const float PP0[9] = {8.75317e-01, 8.73691e-01, 8.69273e-01, 8.66705e-01, 8.61336e-01, 8.48794e-01, 8.33457e-01, 8.16553e-01, 8.01148e-01};
const float PP1[9] = {1.29625e-01, 1.29080e-01, 1.31148e-01, 1.30600e-01, 1.31081e-01, 1.33922e-01, 1.34175e-01, 1.34916e-01, 1.35806e-01};
const float PP2[9] = {6.88369e+00, 8.00026e+00, 6.25220e+00, 7.18563e+00, 7.06380e+00, 6.77852e+00, 6.91458e+00, 6.74690e+00, 7.10116e+00};

//defining histograms
struct WeightHists
{
    TProfile2D *pTPCmeanPhi_1 ;
    TProfile2D *pTPCmeanPhi_2 ;
    TProfile2D *pTPCmeanPhi_3 ;
    TProfile2D *pTPCmeanPhiAsso_1 ;
    TProfile2D *pTPCmeanPhiAsso_2 ;
    TProfile2D *pTPCmeanPhiAsso_3 ;
};
struct WeightHists WeightHist[10];
struct WeightHists ReadWeightHist[10];

struct PhiHists
{
    TH2D *Hist_Phi_1;
    TH2D *Hist_Phi_2;
    TH2D *Hist_Phi_3;
    TH2D *Hist_Phi_new_1;
    TH2D *Hist_Phi_new_2;
    TH2D *Hist_Phi_new_3;
};
struct PhiHists PhiHist[10];

TH2D *Ref_TOF = new TH2D("Ref_TOF", "Ref_TOF", 500, 0.5, 500.5, 5000, 0.5, 5000.5);
TProfile *Ref_Day3 = new TProfile("Ref_Day3", "RefMult vs Run", run_end - run_sta, run_sta, run_end, 0, 999);
TProfile *TOF_Day3 = new TProfile("TOF_Day3", "TOFMult vs Run", run_end - run_sta, run_sta, run_end, 0, 5000);
TProfile *NPT_Day3 = new TProfile("NPT_Day3", "NPTracks vs Run", run_end - run_sta, run_sta, run_end, 0, 5000);
TProfile *NPA_Day3 = new TProfile("NPA_Day3", "NPAsso vs Run", run_end - run_sta, run_sta, run_end, 0, 5000);
TProfile *TPC_Day3_cos2 = new TProfile("TPC_Day3_cos2", "cos(2*psi) vs Run", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *TPC_Day3_sin2 = new TProfile("TPC_Day3_sin2", "sin(2*psi) vs Run", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *BBCe_Day3_cos2 = new TProfile("BBCe_Day3_cos2", "BBCe_Day3_cos2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *BBCe_Day3_sin2 = new TProfile("BBCe_Day3_sin2", "BBCe_Day3_sin2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *BBCw_Day3_cos2 = new TProfile("BBCw_Day3_cos2", "BBCw_Day3_cos2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *BBCw_Day3_sin2 = new TProfile("BBCw_Day3_sin2", "BBCw_Day3_sin2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDe_Day3_cos2 = new TProfile("EPDe_Day3_cos2", "EPDe_Day3_cos2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDe_Day3_sin2 = new TProfile("EPDe_Day3_sin2", "EPDe_Day3_sin2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDw_Day3_cos2 = new TProfile("EPDw_Day3_cos2", "EPDw_Day3_cos2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPDw_Day3_sin2 = new TProfile("EPDw_Day3_sin2", "EPDw_Day3_sin2", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *ZDCe_Day3_cos = new TProfile("ZDCe_Day3_cos", "ZDCe_Day3_cos", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *ZDCe_Day3_sin = new TProfile("ZDCe_Day3_sin", "ZDCe_Day3_sin", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *ZDCw_Day3_cos = new TProfile("ZDCw_Day3_cos", "ZDCw_Day3_cos", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *ZDCw_Day3_sin = new TProfile("ZDCw_Day3_sin", "ZDCw_Day3_sin", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *BBC_cor_Day3 = new TProfile("BBC_cor_Day3", "BBC_cor_Day3", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *EPD_cor_Day3 = new TProfile("EPD_cor_Day3", "EPD_cor_Day3", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *ZDC_cor_Day3 = new TProfile("ZDC_cor_Day3", "ZDC_cor_Day3", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *v2_TPC_Day3  = new TProfile("v2_TPC_Day3", "v2_TPC_Day3", run_end - run_sta, run_sta, run_end, -100, 100);
TProfile *v2_EPD_Day3  = new TProfile("v2_EPD_Day3", "v2_EPD_Day3", run_end - run_sta, run_sta, run_end, -100, 100);
TProfile *v2_BBC_Day3  = new TProfile("v2_BBC_Day3", "v2_BBC_Day3", run_end - run_sta, run_sta, run_end, -100, 100);
TProfile *v2_ZDC_Day3  = new TProfile("v2_ZDC_Day3", "v2_ZDC_Day3", run_end - run_sta, run_sta, run_end, -100, 100);

TH1D *hTally = new TH1D("hTally", "hTally", 11, -0.5, 10.5);
TH1D *hTall  = new TH1D("hTall ", "hTall ", 12, -0.5, 11.5);
TH1D *hZDCcoin = new TH1D("hZDCcoin", "hZDCcoin", 1000, 0, 100000);
TH1D *hTrigger = new TH1D("hTrigger", "hTrigger", 200, 0.5, 200.5);
TH1D *hCentrality = new TH1D("hCentrality", "hCentrality", 11, -1.5, 9.5);
TH2D *hVertexXY = new TH2D("hVertexXY", "hVertexXY", 300, -3, 3, 300, -3, 3);
TH1D *hVertexZ = new TH1D("hVertexZ", "hVertexZ", 100, -100, 100);
TH1D *hVzDiff = new TH1D("hVzDiff", "hVzDiff", 120, -30, 30);
TH2D *hMult_Vz = new TH2D("hMult_Vz", "hMult_Vz", 1000, -0.5, 999.5, 100, -100, 100);
TH2D *hMult_Vz_new = new TH2D("hMult_Vz_new", "hMult_Vz_new", 1000, -0.5, 999.5, 100, -100, 100);
TH1D *BBC1 = new TH1D("BBC1", "BBC1", 6, 0.5, 6.5);
TH1D *BBC2 = new TH1D("BBC2", "BBC2", 10, 0.5, 10.5);
TH1D *BBC7 = new TH1D("BBC7", "BBC7", 10, 0.5, 10.5);
TH1D *BBC8 = new TH1D("BBC8", "BBC8", 6, 0.5, 6.5);
TH1D *ZDC_e_h = new TH1D("ZDC_e_h", "ZDC_e_h", 8, 0.5, 8.5);
TH1D *ZDC_e_v = new TH1D("ZDC_e_v", "ZDC_e_v", 8, 0.5, 8.5);
TH1D *ZDC_w_h = new TH1D("ZDC_w_h", "ZDC_w_h", 8, 0.5, 8.5);
TH1D *ZDC_w_v = new TH1D("ZDC_w_v", "ZDC_w_v", 8, 0.5, 8.5);
TProfile2D *pZDCcenter = new TProfile2D("pZDCcenter", "pZDCcenter",
                                        4, 0.5, 4.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, 0, 20);
TProfile2D *pZDCcenter_new = new TProfile2D("pZDCcenter_new", "pZDCcenter_new",
        4, 0.5, 4.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -20, 20);

TH1D *Hist_positive = new TH1D("Hist_positive", "Hist_positive", 500, -0.5, 499.5);
TH1D *Hist_negative = new TH1D("Hist_negative", "Hist_negative", 500, -0.5, 499.5);
TH1D *Hist_Ch       = new TH1D("Hist_Ch", "Hist_Ch", 1000, -0.5, 999.5);
TH1D *Hist_netCh    = new TH1D("Hist_netCh", "Hist_netCh", 999, -499.5, 499.5);
TH1D *Hist_netChAsym = new TH1D("Hist_netChAsym", "Hist_netChAsym", 600, -1.5 + 0.0025, 1.5 + 0.0025);
TH2D *Hist_netChAsym_bin    = new TH2D("Hist_netChAsym_bin", "Hist_netChAsym_bin", 5, 0.5, 5.5, 600, -1.5 + 0.0025, 1.5 + 0.0025);
TProfile *p_netChAsym_RefMult = new TProfile("p_netChAsym_RefMult", "p_netChAsym_RefMult", 300, -1.5 + 0.0025, 1.5 + 0.0025, 0., 900);
TProfile *p_netChAsym_cos     = new TProfile("p_netChAsym_cos", "p_netChAsym_cos", 300, -1.5 + 0.0025, 1.5 + 0.0025, -1, 1);
TProfile *Hist_cos_Ach = new TProfile("Hist_cos_Ach", "Hist_cos_Ach", 5, 0.5, 5.5, -1, 1);
TProfile *p_v2_Ach = new TProfile("p_v2_Ach", "p_v2_Ach", 300, -1.5 + 0.0025, 1.5 + 0.0025, -100, 100);
TH2D *Hist_pt_pos_Ach = new TH2D("Hist_pt_pos_Ach", "Hist_pt_pos_Ach", 5, 0.5, 5.5, 300, 0, 15);
TH2D *Hist_pt_neg_Ach = new TH2D("Hist_pt_neg_Ach", "Hist_pt_neg_Ach", 5, 0.5, 5.5, 300, 0, 15);
TProfile2D *p_v2_pt_pos_Ach = new TProfile2D("p_v2_pt_pos_Ach", "p_v2_pt_pos_Ach", 5, 0.5, 5.5, 300, 0, 15, -100, 100);
TProfile2D *p_v2_pt_neg_Ach = new TProfile2D("p_v2_pt_neg_Ach", "p_v2_pt_neg_Ach", 5, 0.5, 5.5, 300, 0, 15, -100, 100);

TH2D *Hist_TPC_EP_east = new TH2D("Hist_TPC_EP_east", "Hist_TPC_EP_east", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_west = new TH2D("Hist_TPC_EP_west", "Hist_TPC_EP_west", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_full = new TH2D("Hist_TPC_EP_full", "Hist_TPC_EP_full", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_for = new TH2D("Hist_TPC_EP_for", "Hist_TPC_EP_for", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_bac = new TH2D("Hist_TPC_EP_bac", "Hist_TPC_EP_bac", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_BBC_EP_east = new TH2D("Hist_BBC_EP_east", "Hist_BBC_EP_east", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_BBC_EP_west = new TH2D("Hist_BBC_EP_west", "Hist_BBC_EP_west", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP_east = new TH2D("Hist_EPD_EP_east", "Hist_EPD_EP_east", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP_west = new TH2D("Hist_EPD_EP_west", "Hist_EPD_EP_west", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_ZDC_EP_east = new TH2D("Hist_ZDC_EP_east", "Hist_ZDC_EP_east", 72, -PI, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_ZDC_EP_west = new TH2D("Hist_ZDC_EP_west", "Hist_ZDC_EP_west", 72, -PI, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_east_flat = new TH2D("Hist_TPC_EP_east_flat", "Hist_TPC_EP_east_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_west_flat = new TH2D("Hist_TPC_EP_west_flat", "Hist_TPC_EP_west_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_full_flat = new TH2D("Hist_TPC_EP_full_flat", "Hist_TPC_EP_full_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_for_flat = new TH2D("Hist_TPC_EP_for_flat", "Hist_TPC_EP_for_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_bac_flat = new TH2D("Hist_TPC_EP_bac_flat", "Hist_TPC_EP_bac_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_BBC_EP_east_flat = new TH2D("Hist_BBC_EP_east_flat", "Hist_BBC_EP_east_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_BBC_EP_west_flat = new TH2D("Hist_BBC_EP_west_flat", "Hist_BBC_EP_west_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP_east_flat = new TH2D("Hist_EPD_EP_east_flat", "Hist_EPD_EP_east_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_EPD_EP_west_flat = new TH2D("Hist_EPD_EP_west_flat", "Hist_EPD_EP_west_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_ZDC_EP_east_flat = new TH2D("Hist_ZDC_EP_east_flat", "Hist_ZDC_EP_east_flat", 72, -PI, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_ZDC_EP_west_flat = new TH2D("Hist_ZDC_EP_west_flat", "Hist_ZDC_EP_west_flat", 72, -PI, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH1D *Hist_TPC_EP_full_m1 = new TH1D("Hist_TPC_EP_full_m1", "Hist_TPC_EP_full_m1", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m2 = new TH1D("Hist_TPC_EP_full_m2", "Hist_TPC_EP_full_m2", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m1_flat = new TH1D("Hist_TPC_EP_full_m1_flat", "Hist_TPC_EP_full_m1_flat", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m2_flat = new TH1D("Hist_TPC_EP_full_m2_flat", "Hist_TPC_EP_full_m2_flat", 36, 0, PI);
TProfile2D *pTPC_EP_east = new TProfile2D("pTPC_EP_east", "pTPC_EP_east", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_west = new TProfile2D("pTPC_EP_west", "pTPC_EP_west", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_full = new TProfile2D("pTPC_EP_full", "pTPC_EP_full", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_for = new TProfile2D("pTPC_EP_for", "pTPC_EP_for", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_bac = new TProfile2D("pTPC_EP_bac", "pTPC_EP_bac", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pBBC_EP_east = new TProfile2D("pBBC_EP_east", "pBBC_EP_east", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pBBC_EP_west = new TProfile2D("pBBC_EP_west", "pBBC_EP_west", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pEPD_EP_east = new TProfile2D("pEPD_EP_east", "pEPD_EP_east", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pEPD_EP_west = new TProfile2D("pEPD_EP_west", "pEPD_EP_west", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pZDC_EP_east = new TProfile2D("pZDC_EP_east", "pZDC_EP_east", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pZDC_EP_west = new TProfile2D("pZDC_EP_west", "pZDC_EP_west", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");

TH2D *Hist_BBC_vs_TPC = new TH2D("Hist_BBC_vs_TPC", "Hist_BBC_vs_TPC", 30, 0, PI, 30, 0, PI);
TH2D *Hist_EPD_vs_TPC = new TH2D("Hist_EPD_vs_TPC", "Hist_EPD_vs_TPC", 30, 0, PI, 30, 0, PI);

TProfile *Hist_cos = new TProfile("Hist_cos", "Hist_cos", 4, 0.5, 4.5, -1, 1, "");
TProfile *Hist_cos_BBC = new TProfile("Hist_cos_BBC", "Hist_cos_BBC", 4, 0.5, 4.5, -1, 1, "");
TProfile *Hist_cos_EPD = new TProfile("Hist_cos_EPD", "Hist_cos_EPD", 4, 0.5, 4.5, -1, 1, "");
TProfile *Hist_cos_ZDC = new TProfile("Hist_cos_ZDC", "Hist_cos_ZDC", 4, 0.5, 4.5, -1, 1, "");

TH1D *Hist_DCA = new TH1D("Hist_DCA", "Hist_DCA", 100, 0, 10);
TH2D *hEtaPtDist = new TH2D("EtaPtDist", "EtaPtDist", 26, -1.3, 1.3, 300, 0, 15);
TH2D *hEtaPhiDist = new TH2D("hEtaPhiDist", "hEtaPhiDist", 26, -1.3, 1.3, Phibin, -PI, PI);
TH2D *hPhiPtDist = new TH2D("PhiPtDist", "PhiPtDist", Phibin, -PI, PI, 300, 0, 15);
TH1D *Hist_Pt = new TH1D("Hist_Pt", "Hist_Pt", 300, 0, 15);
TH1D *Hist_Pt_TOF = new TH1D("Hist_Pt_TOF", "Hist_Pt_TOF", 300, 0, 15);
TH1D *rc;
TH2D *wt = new TH2D("Order1etaWeight", "Order1etaWeight", 100, 1.5, 6.5, 9, 0, 9);
TH2D *wt2 = new TH2D("Order2etaWeight", "Order2etaWeight", 100, 1.5, 6.5, 9, 0, 9);
TH2D *Hist_Phi = new TH2D("Hist_Phi", "Hist_Phi", Phibin, -PI, PI, 4, 0.5, 4.5);
TH1D *hDpt   = new TH1D("hDpt", "hDpt", 200, 0, 2);

TProfile *Hist_v2_pt_obs1 = new TProfile("Hist_v2_pt_obs1", "Hist_v2_pt_obs1", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs2 = new TProfile("Hist_v2_pt_obs2", "Hist_v2_pt_obs2", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_eta_obs1 = new TProfile("Hist_v2_eta_obs1", "Hist_v2_eta_obs1", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_obs2 = new TProfile("Hist_v2_eta_obs2", "Hist_v2_eta_obs2", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_obs3 = new TProfile("Hist_v2_eta_obs3", "Hist_v2_eta_obs3", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_pt_EPD_obs = new TProfile("Hist_v2_pt_EPD_obs", "Hist_v2_pt_EPD_obs", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_BBC_obs = new TProfile("Hist_v2_pt_BBC_obs", "Hist_v2_pt_BBC_obs", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_ZDC_obs = new TProfile("Hist_v2_pt_ZDC_obs", "Hist_v2_pt_ZDC_obs", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_eta_EPD_obs = new TProfile("Hist_v2_eta_EPD_obs", "Hist_v2_eta_EPD_obs", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_EPDe_obs = new TProfile("Hist_v2_eta_EPDe_obs", "Hist_v2_eta_EPDe_obs", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_EPDw_obs = new TProfile("Hist_v2_eta_EPDw_obs", "Hist_v2_eta_EPDw_obs", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_BBC_obs = new TProfile("Hist_v2_eta_BBC_obs", "Hist_v2_eta_BBC_obs", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_BBCe_obs = new TProfile("Hist_v2_eta_BBCe_obs", "Hist_v2_eta_BBCe_obs", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_BBCw_obs = new TProfile("Hist_v2_eta_BBCw_obs", "Hist_v2_eta_BBCw_obs", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_ZDC_obs = new TProfile("Hist_v2_eta_ZDC_obs", "Hist_v2_eta_ZDC_obs", 40, -1, 1, -100, 100, "");

TProfile *pTemp_origin_obs = new TProfile("pTemp_origin_obs", "pTemp_origin_obs", 4, 0.5, 4.5);
TProfile *pTemp_shuffle_obs = new TProfile("pTemp_shuffle_obs", "pTemp_shuffle_obs", 4, 0.5, 4.5);
TH1D *Hist_dS_sin_obs = new TH1D("Hist_dS_sin_obs", "Hist_dS_sin_obs", 200, -1, 1);
TH1D *Hist_dS_cos_obs = new TH1D("Hist_dS_cos_obs", "Hist_dS_cos_obs", 200, -1, 1);
TH1D *Hist_dS_sin_shuffle_obs = new TH1D("Hist_dS_sin_shuffle_obs", "Hist_dS_sin_shuffle_obs", 200, -1, 1);
TH1D *Hist_dS_cos_shuffle_obs = new TH1D("Hist_dS_cos_shuffle_obs", "Hist_dS_cos_shuffle_obs", 200, -1, 1);

TH1D *Hist_dQ_in1  = new TH1D("Hist_dQ_in1", "Hist_dQ_in1", 201, -100.5, 100.5);
TH1D *Hist_dQ_out1 = new TH1D("Hist_dQ_out1", "Hist_dQ_out1", 201, -100.5, 100.5);
TH1D *Hist_dQ_in2  = new TH1D("Hist_dQ_in2", "Hist_dQ_in2", 201, -100.5, 100.5);
TH1D *Hist_dQ_out2 = new TH1D("Hist_dQ_out2", "Hist_dQ_out2", 201, -100.5, 100.5);
TProfile *pParity_int_MSC_same_in1 = new TProfile("pParity_int_MSC_same_in1", "pParity_int_MSC_same_in1", 201, -100.5, 100.5, -100, 100, "");
TProfile *pParity_int_MSC_oppo_in1 = new TProfile("pParity_int_MSC_oppo_in1", "pParity_int_MSC_oppo_in1", 201, -100.5, 100.5, -100, 100, "");
TProfile *pParity_int_MSC_same_out1 = new TProfile("pParity_int_MSC_same_out1", "pParity_int_MSC_same_out1", 201, -100.5, 100.5, -100, 100, "");
TProfile *pParity_int_MSC_oppo_out1 = new TProfile("pParity_int_MSC_oppo_out1", "pParity_int_MSC_oppo_out1", 201, -100.5, 100.5, -100, 100, "");
TProfile *pParity_int_MSC_same_in2 = new TProfile("pParity_int_MSC_same_in2", "pParity_int_MSC_same_in2", 201, -100.5, 100.5, -100, 100, "");
TProfile *pParity_int_MSC_oppo_in2 = new TProfile("pParity_int_MSC_oppo_in2", "pParity_int_MSC_oppo_in2", 201, -100.5, 100.5, -100, 100, "");
TProfile *pParity_int_MSC_same_out2 = new TProfile("pParity_int_MSC_same_out2", "pParity_int_MSC_same_out2", 201, -100.5, 100.5, -100, 100, "");
TProfile *pParity_int_MSC_oppo_out2 = new TProfile("pParity_int_MSC_oppo_out2", "pParity_int_MSC_oppo_out2", 201, -100.5, 100.5, -100, 100, "");
TProfile *pParity_int_MSC1_obs1 = new TProfile("pParity_int_MSC1_obs1", "pParity_int_MSC1_obs1", 4, 0.5, 4.5, -200, 200, "");
TProfile *pParity_int_MSC1_obs2 = new TProfile("pParity_int_MSC1_obs2", "pParity_int_MSC1_obs2", 4, 0.5, 4.5, -200, 200, "");
TProfile *pParity_int_MSC2_obs1 = new TProfile("pParity_int_MSC2_obs1", "pParity_int_MSC2_obs1", 4, 0.5, 4.5, -200, 200, "");
TProfile *pParity_int_MSC2_obs2 = new TProfile("pParity_int_MSC2_obs2", "pParity_int_MSC2_obs2", 4, 0.5, 4.5, -200, 200, "");

TProfile *pParity_int_obs1 = new TProfile("Parity_int_obs1", "Parity_int_obs1", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_obs3 = new TProfile("Parity_int_obs3", "Parity_int_obs3", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_ss_obs1 = new TProfile("Parity_int_ss_obs1", "Parity_int_ss_obs1", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_ss_obs3 = new TProfile("Parity_int_ss_obs3", "Parity_int_ss_obs3", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_int_ss_same_run = new TProfile("Parity_int_ss_same_run", "Parity_int_ss_same_run", (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -100, 100, "");
TProfile *pParity_int_ss_oppo_run = new TProfile("Parity_int_ss_oppo_run", "Parity_int_ss_oppo_run", (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -100, 100, "");
TProfile2D *pParity_eta_ss_obs1 = new TProfile2D("Parity_eta_ss_obs1", "Parity_eta_ss_obs1", 12, 0.5, 12.5, 20, -1, 1, -100, 100, "");
TProfile2D *pParity_eta_ss_obs3 = new TProfile2D("Parity_eta_ss_obs3", "Parity_eta_ss_obs3", 12, 0.5, 12.5, 20, -1, 1, -100, 100, "");
TProfile2D *pParity_Deta_ss_obs1 = new TProfile2D("Parity_Deta_ss_obs1", "Parity_Deta_ss_obs1", 12, 0.5, 12.5, 40, 0, 2, -100, 100, "");
TProfile2D *pParity_Deta_ss_obs3 = new TProfile2D("Parity_Deta_ss_obs3", "Parity_Deta_ss_obs3", 12, 0.5, 12.5, 40, 0, 2, -100, 100, "");
TProfile2D *pParity_pt_ss_obs1  = new TProfile2D("Parity_pt_ss_obs1", "Parity_pt_ss_obs1", 12, 0.5, 12.5, 20, 0, 2.0, -100, 100, "");
TProfile2D *pParity_pt_ss_obs3  = new TProfile2D("Parity_pt_ss_obs3", "Parity_pt_ss_obs3", 12, 0.5, 12.5, 20, 0, 2.0, -100, 100, "");
TProfile2D *pParity_Dpt_ss_obs1 = new TProfile2D("Parity_Dpt_ss_obs1", "Parity_Dpt_ss_obs1", 12, 0.5, 12.5, 200, 0, 2.0, -100, 100, "");
TProfile2D *pParity_Dpt_ss_obs3 = new TProfile2D("Parity_Dpt_ss_obs3", "Parity_Dpt_ss_obs3", 12, 0.5, 12.5, 200, 0, 2.0, -100, 100, "");
TProfile *pParity_noHBT_ss_obs1   = new TProfile("Parity_noHBT_ss_obs1", "Parity_noHBT_ss_obs1", 4, 0.5, 4.5, -100, 100, "");
TProfile *pParity_noHBT_ss_obs3   = new TProfile("Parity_noHBT_ss_obs3", "Parity_noHBT_ss_obs3", 4, 0.5, 4.5, -100, 100, "");
TProfile2D *pParity_Deta_highDpt_ss_obs1 = new TProfile2D("pParity_Deta_highDpt_ss_obs1", "pParity_Deta_highDpt_ss_obs1", 4, 0.5, 4.5, 40, 0, 2, -100, 100, "");
TProfile2D *pParity_Deta_highDpt_ss_obs3 = new TProfile2D("pParity_Deta_highDpt_ss_obs3", "pParity_Deta_highDpt_ss_obs3", 4, 0.5, 4.5, 40, 0, 2, -100, 100, "");

TProfile *pDelta_int_ss_obs1 = new TProfile("Delta_int_ss_obs1", "Delta_int_ss_obs1", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_obs3 = new TProfile("Delta_int_ss_obs3", "Delta_int_ss_obs3", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_int_ss_same_run = new TProfile("Delta_int_ss_same_run", "Delta_int_ss_same_run", (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -100, 100, "");
TProfile *pDelta_int_ss_oppo_run = new TProfile("Delta_int_ss_oppo_run", "Delta_int_ss_oppo_run", (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -100, 100, "");
TProfile2D *pDelta_eta_ss_obs1 = new TProfile2D("Delta_eta_ss_obs1", "Delta_eta_ss_obs1", 12, 0.5, 12.5, 20, -1, 1, -100, 100, "");
TProfile2D *pDelta_eta_ss_obs3 = new TProfile2D("Delta_eta_ss_obs3", "Delta_eta_ss_obs3", 12, 0.5, 12.5, 20, -1, 1, -100, 100, "");
TProfile2D *pDelta_Deta_ss_obs1 = new TProfile2D("Delta_Deta_ss_obs1", "Delta_Deta_ss_obs1", 12, 0.5, 12.5, 40, 0, 2, -100, 100, "");
TProfile2D *pDelta_Deta_ss_obs3 = new TProfile2D("Delta_Deta_ss_obs3", "Delta_Deta_ss_obs3", 12, 0.5, 12.5, 40, 0, 2, -100, 100, "");
TProfile2D *pDelta_pt_ss_obs1  = new TProfile2D("Delta_pt_ss_obs1", "Delta_pt_ss_obs1", 12, 0.5, 12.5, 20, 0, 2.0, -100, 100, "");
TProfile2D *pDelta_pt_ss_obs3  = new TProfile2D("Delta_pt_ss_obs3", "Delta_pt_ss_obs3", 12, 0.5, 12.5, 20, 0, 2.0, -100, 100, "");
TProfile2D *pDelta_Dpt_ss_obs1 = new TProfile2D("Delta_Dpt_ss_obs1", "Delta_Dpt_ss_obs1", 12, 0.5, 12.5, 200, 0, 2.0, -100, 100, "");
TProfile2D *pDelta_Dpt_ss_obs3 = new TProfile2D("Delta_Dpt_ss_obs3", "Delta_Dpt_ss_obs3", 12, 0.5, 12.5, 200, 0, 2.0, -100, 100, "");
TProfile *pDelta_noHBT_ss_obs1   = new TProfile("Delta_noHBT_ss_obs1", "Delta_noHBT_ss_obs1", 4, 0.5, 4.5, -100, 100, "");
TProfile *pDelta_noHBT_ss_obs3   = new TProfile("Delta_noHBT_ss_obs3", "Delta_noHBT_ss_obs3", 4, 0.5, 4.5, -100, 100, "");
TProfile2D *pDelta_Deta_highDpt_ss_obs1 = new TProfile2D("pDelta_Deta_highDpt_ss_obs1", "pDelta_Deta_highDpt_ss_obs1", 4, 0.5, 4.5, 40, 0, 2, -100, 100, "");
TProfile2D *pDelta_Deta_highDpt_ss_obs3 = new TProfile2D("pDelta_Deta_highDpt_ss_obs3", "pDelta_Deta_highDpt_ss_obs3", 4, 0.5, 4.5, 40, 0, 2, -100, 100, "");

TProfile *pTemp_v2 = new TProfile("pTemp_v2", "pTemp_v2", 10, 0.5, 10.5, -100, 100, "");
TProfile *pTemp_parity_e = new TProfile("pTemp_parity_e", "pTemp_parity_e", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_parity_w = new TProfile("pTemp_parity_w", "pTemp_parity_w", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_delta = new TProfile("pTemp_delta", "pTemp_delta", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_v2_noHBT = new TProfile("pTemp_v2_noHBT", "pTemp_v2_noHBT", 6, 0.5, 6.5, -100, 100, "");
TProfile *pTemp_parity_e_noHBT = new TProfile("pTemp_parity_e_noHBT", "pTemp_parity_e_noHBT", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_parity_w_noHBT = new TProfile("pTemp_parity_w_noHBT", "pTemp_parity_w_noHBT", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_delta_noHBT = new TProfile("pTemp_delta_noHBT", "pTemp_delta_noHBT", 8, 0.5, 8.5, -100, 100, "");

TH1D *Hist_Q  = new TH1D("Hist_Q", "Hist_Q", 250, 0, 5);
TH1D *Hist_Q2 = new TH1D("Hist_Q2", "Hist_Q2", 250, 0, 25);
TH2D *Hist_RefMult_Q2 = new TH2D("Hist_RefMult_Q2", "Hist_RefMult_Q2", 250, 0, 25, 500, 0, 500);
TProfile *p_RefMult_Q2 = new TProfile("p_RefMult_Q2", "p_RefMult_Q2", 250, 0, 25, 0, 1000, "");
TProfile *p_cos_Q2 = new TProfile("p_cos_Q2", "p_cos_Q2", 250, 0, 25, -1, 1, "");
TProfile *p_v2e_Q2_obs1 = new TProfile("p_v2e_Q2_obs1", "p_v2e_Q2_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_Q2_obs2 = new TProfile("p_v2e_Q2_obs2", "p_v2e_Q2_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_Q2_obs1 = new TProfile("p_v2w_Q2_obs1", "p_v2w_Q2_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_Q2_obs2 = new TProfile("p_v2w_Q2_obs2", "p_v2w_Q2_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_Q2_P_obs1 = new TProfile("p_v2e_Q2_P_obs1", "p_v2e_Q2_P_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_Q2_P_obs2 = new TProfile("p_v2e_Q2_P_obs2", "p_v2e_Q2_P_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_Q2_P_obs1 = new TProfile("p_v2w_Q2_P_obs1", "p_v2w_Q2_P_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_Q2_P_obs2 = new TProfile("p_v2w_Q2_P_obs2", "p_v2w_Q2_P_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_Q2_N_obs1 = new TProfile("p_v2e_Q2_N_obs1", "p_v2e_Q2_N_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_Q2_N_obs2 = new TProfile("p_v2e_Q2_N_obs2", "p_v2e_Q2_N_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_Q2_N_obs1 = new TProfile("p_v2w_Q2_N_obs1", "p_v2w_Q2_N_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_Q2_N_obs2 = new TProfile("p_v2w_Q2_N_obs2", "p_v2w_Q2_N_obs2", 250, 0, 25, -100, 100, "");
TProfile2D *pParity_e_Q2_obs1 = new TProfile2D("Parity_e_Q2_obs1", "Parity_e_Q2_obs1", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_e_Q2_obs2 = new TProfile2D("Parity_e_Q2_obs2", "Parity_e_Q2_obs2", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_w_Q2_obs1 = new TProfile2D("Parity_w_Q2_obs1", "Parity_w_Q2_obs1", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_w_Q2_obs2 = new TProfile2D("Parity_w_Q2_obs2", "Parity_w_Q2_obs2", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pDelta_Q2_obs1 = new TProfile2D("pDelta_Q2_obs1", "pDelta_Q2_obs1", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pDelta_Q2_obs2 = new TProfile2D("pDelta_Q2_obs2", "pDelta_Q2_obs2", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_obs1 = new TProfile("p_v2e_noHBT_Q2_obs1", "p_v2e_noHBT_Q2_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_obs2 = new TProfile("p_v2e_noHBT_Q2_obs2", "p_v2e_noHBT_Q2_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_obs1 = new TProfile("p_v2w_noHBT_Q2_obs1", "p_v2w_noHBT_Q2_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_obs2 = new TProfile("p_v2w_noHBT_Q2_obs2", "p_v2w_noHBT_Q2_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_P_obs1 = new TProfile("p_v2e_noHBT_Q2_P_obs1", "p_v2e_noHBT_Q2_P_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_P_obs2 = new TProfile("p_v2e_noHBT_Q2_P_obs2", "p_v2e_noHBT_Q2_P_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_P_obs1 = new TProfile("p_v2w_noHBT_Q2_P_obs1", "p_v2w_noHBT_Q2_P_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_P_obs2 = new TProfile("p_v2w_noHBT_Q2_P_obs2", "p_v2w_noHBT_Q2_P_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_N_obs1 = new TProfile("p_v2e_noHBT_Q2_N_obs1", "p_v2e_noHBT_Q2_N_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2e_noHBT_Q2_N_obs2 = new TProfile("p_v2e_noHBT_Q2_N_obs2", "p_v2e_noHBT_Q2_N_obs2", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_N_obs1 = new TProfile("p_v2w_noHBT_Q2_N_obs1", "p_v2w_noHBT_Q2_N_obs1", 250, 0, 25, -100, 100, "");
TProfile *p_v2w_noHBT_Q2_N_obs2 = new TProfile("p_v2w_noHBT_Q2_N_obs2", "p_v2w_noHBT_Q2_N_obs2", 250, 0, 25, -100, 100, "");
TProfile2D *pParity_e_noHBT_Q2_obs1 = new TProfile2D("Parity_e_noHBT_Q2_obs1", "Parity_e_noHBT_Q2_obs1", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_e_noHBT_Q2_obs2 = new TProfile2D("Parity_e_noHBT_Q2_obs2", "Parity_e_noHBT_Q2_obs2", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_w_noHBT_Q2_obs1 = new TProfile2D("Parity_w_noHBT_Q2_obs1", "Parity_w_noHBT_Q2_obs1", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pParity_w_noHBT_Q2_obs2 = new TProfile2D("Parity_w_noHBT_Q2_obs2", "Parity_w_noHBT_Q2_obs2", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pDelta_noHBT_Q2_obs1 = new TProfile2D("pDelta_noHBT_Q2_obs1", "pDelta_noHBT_Q2_obs1", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");
TProfile2D *pDelta_noHBT_Q2_obs2 = new TProfile2D("pDelta_noHBT_Q2_obs2", "pDelta_noHBT_Q2_obs2", 4, 0.5, 4.5, 250, 0, 25, -100, 100, "");

//defining variables
float pVx, pVy, pVz, VPDvz, BBCco, ZDCcoin, net_Nch_Asym, mQx, mQy;                   //run, event info
int   Run, Day, Day2, Day3, Trigger, RefMult, TOFMult, Ntofmatch, Centrality, NPTracks, Fcount, Scount, POIcount;   //
int   Charge, Charge2, ChargeAsso;
float ndEdx, nSigma_p, nSigma_pi, DCAGlobal, Eta, Theta, Phi, Pt, eff, TOFflag;             //track info
float ndEdx2, nSigma_p2, nSigma_pi2, DCAGlobal2, Eta2, Theta2, Phi2, Pt2, eff2, TOFflag2;   //2nd track info
float DCAGlobalAsso, EtaAsso, PhiAsso, PtAsso, TOFflagAsso;
TVector3 pV;

float Eweight = 1;
int Weight_Read = 0;
TH1D *TOF_eff = 0;
float PhiMean_sin[order] = {0}, PhiMean_cos[order] = {0}, PhiMeanAsso_sin[order] = {0}, PhiMeanAsso_cos[order] = {0};
vector<float> PhiAsso_new, Phi_new;
vector<int> iCharge, match;
TProfile2D *Read_TPC_EP_full = 0, *Read_TPC_EP_east = 0, *Read_TPC_EP_west = 0, *Read_TPC_EP_for = 0, *Read_TPC_EP_bac = 0;
float PsiMean_F[2 * order] = {0}, PsiMean_E[2 * order] = {0}, PsiMean_W[2 * order] = {0}, PsiMean_f[2 * order] = {0}, PsiMean_b[2 * order] = {0};
TProfile2D *Read_BBC_EP_east = 0, *Read_BBC_EP_west = 0, *Read_EPD_EP_east = 0, *Read_EPD_EP_west = 0;
float Psi_BBC_E[2 * order] = {0}, Psi_BBC_W[2 * order] = {0}, Psi_EPD_E[2 * order] = {0}, Psi_EPD_W[2 * order] = {0};
float MeanNetChargeAsym, RMSNetChargeAsym;
float BBC_gain_east[16] = {0}, BBC_gain_west[16] = {0};
float ZDC_gain_east_v[7] = {0}, ZDC_gain_east_h[8] = {0}, ZDC_gain_west_v[7] = {0}, ZDC_gain_west_h[8] = {0};
TProfile2D *Read_ZDCcenter = 0, *Read_ZDC_EP_east = 0, *Read_ZDC_EP_west = 0;
float zdcsmd_x0_e = 0, zdcsmd_y0_e = 0, zdcsmd_x0_w = 0, zdcsmd_y0_w = 0;
float Psi_ZDC_E[2 * order] = {0}, Psi_ZDC_W[2 * order] = {0};

float TPC_EP_full = 0, TPC_EP_east = 0, TPC_EP_west = 0, TPC_EP_for = 0, TPC_EP_bac = 0;
float TPC_EP_full_new = 0, TPC_EP_east_new = 0, TPC_EP_west_new = 0, TPC_EP_for_new = 0, TPC_EP_bac_new = 0;
float BBC_EP_east = 0, BBC_EP_west = 0, EPD_EP_east = 0, EPD_EP_west = 0, ZDC_EP_east = 0, ZDC_EP_west = 0;
float BBC_EP_east_new = 0, BBC_EP_west_new = 0, EPD_EP_east_new = 0, EPD_EP_west_new = 0, ZDC_EP_east_new = 0, ZDC_EP_west_new = 0;
float Q2 = 0;
int N_L_P1 = 0, N_L_N1 = 0, N_R_P1 = 0, N_R_N1 = 0, N_T_P1 = 0, N_T_N1 = 0, N_B_P1 = 0, N_B_N1 = 0;
int N_L_P2 = 0, N_L_N2 = 0, N_R_P2 = 0, N_R_N2 = 0, N_T_P2 = 0, N_T_N2 = 0, N_B_P2 = 0, N_B_N2 = 0;
float v2 = 0, v2e = 0, v2w = 0, v2_sub = 0, v2_EPDe = 0, v2_EPDw = 0, v2_BBCe = 0, v2_BBCw = 0, v2_ZDC = 0;
float Rsin = 0, Rcos = 0;
int iCh = 0;
float correlator0 = 0, correlator3 = 0, correlator4 = 0, correlator4e = 0, correlator4w = 0;

//StRefMultCorr* refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;

//sub-functions
void DefineHistogram();                     //define histograms
void WriteHistogram(int c, int o);              //into result ROOT file
void WriteWeight(char *OutFileName);                //into weight file
int  ReadWeight(char *InFileName);              //input from weight file
bool IsGoodEvent(int c);                    //select good events and fill event level histograms
bool IsGoodBBC(StPicoEvent *e);                     //cuts on BBC ADCs
bool IsGoodZDC(StPicoEvent *e);                     //cuts on ZDC ADCs
bool IsGoodAsso(StPicoTrack *p);            //cuts on EP particles
bool IsGoodPOI(StPicoTrack *p);             //cuts on particles of interest
bool IsGoodTrack(StPicoTrack *p);               //cuts on tracks
bool IsGoodPionNoTOF(StPicoTrack *p);               //cuts on low pT pions without TOF
bool IsGoodPion(StPicoDst *d, StPicoTrack *p, int opt);     //cuts on pions
bool IsGoodKaon(StPicoDst *d, StPicoTrack *p, int opt);         //cuts on kaons
bool IsGoodProton(StPicoDst *d, StPicoTrack *p, int opt);       //cuts on protons
bool CountCharge(StPicoDst *d);                 //count good tracks, charge asymmetry
void MakeTPC_EP(StPicoDst *d, int *iTr);            //reconstruct TPC EPs
void MakeBBC_EP(StPicoEvent *ev);               //reconstruct BBC EPs
void MakeZDC_EP(StPicoEvent *ev);               //reconstruct ZDC EPs
void FillEP_resolution();                   //Fill profiles for EP resolution
void FillPhiPOI();                          //particles of interest, shift parameters to make phi distribution flat
void FillPhiAsso();                         //particles for EP, shift parameters to make phi distribution flat
void ShiftPhiAsso(int tr);                      //flatten the phi distribution
void ShiftPhiPOI(int tr);                   //flatten the phi distribution
void ShiftPsi();                        //flatten EPs
void FillCMW();                         //v2 for CMW
void Fillv2();                          //v2 with different EPs
void FillRcorr();                       //R correlator
void WrapUpRcorr();
void FillGamma(int ord);                    //correlations
void CountMSC(int *iTr, float Phi_ne);              //MSC correlator
void WrapUpMSC();                       //MSC
void WrapUpESE();                                   //Event Shape Engineering
float GetPhiInBBC(int e_w, int bbcN);               //input est_wst(0,1),BBC PMT#
float GetXYInZDC(int e_w, int v_h, int zdcN, int opt_raw = 0);  //input ver_hor(0,1),ZDCSMD slat#

//void __attribute__((constructor)) LoadLib();
/////////////////////////////main program starts here/////////////////////////////////////
void Gamma_sys(int cen = 1, int opt_weight = 1, const Char_t *inFile = "test.list")     //main_function
{
    delete gRandom;
    gRandom = new TRandom(0);

    DefineHistogram();
    char fname_old[200], fname_new[200];
    sprintf(fname_new, "cen%d.weight_1%d%d_new.root", cen, nHar - 1, nHar);
    sprintf(fname_old, "cen%d.weight_1%d%d.root", cen, nHar - 1, nHar);
    Weight_Read = ReadWeight(fname_old);

    StPicoDstReader *picoReader = new StPicoDstReader(inFile);
    picoReader->Init();
    if( !picoReader->chain() )
    {
        std::cout << "No chain has been found." << std::endl;
    }
    Int_t nentries = picoReader->chain()->GetEntries();

    StRefMultCorr *refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;

    // This is a way if you want to spead up IO
    std::cout << "Explicit read status for some branches" << std::endl;
    picoReader->SetStatus("*", 0);
    picoReader->SetStatus("Event", 1);
    picoReader->SetStatus("Track", 1);
    picoReader->SetStatus("BTofHit", 1);
    picoReader->SetStatus("BTofPidTraits", 1);
    // Prepare EPD
    TClonesArray *mEpdHits = new TClonesArray("StPicoEpdHit");
    unsigned int found;
    picoReader->chain()->SetBranchStatus("EpdHit*", 1, &found);
    cout << "StPicoDstMaker::SetStatus " << found << " to EpdHit" << endl;
    picoReader->chain()->SetBranchAddress("EpdHit", &mEpdHits);
    StEpdEpFinder *mEpFinder = new StEpdEpFinder(9, fname_new, fname_old);
    mEpFinder->SetnMipThreshold(0.3);       // recommended by EPD group
    mEpFinder->SetMaxTileWeight(2.0);       // recommended by EPD group
    mEpFinder->SetEpdHitFormat(2);          // 2=pico
    mEpFinder->SetEtaWeights(1, wt);        // eta weight for 1st-order EP
    //        mEpFinder->SetEtaWeights(2,wt2);  // eta weight for 2nd-order EP, select different eta range

    //loop through events
    for(int i = 0; i < nentries; i++)
    {

        if((i + 1) % 1000 == 0) cout << "Processing entry == " << i + 1 << " == out of " << nentries << ".\n";
        Bool_t readEvent = picoReader->readPicoEvent(i);
        if( !readEvent )
        {
            cout << "Something went wrong, Master! Nothing to analyze..." << endl;
            break;
        }
        // Retrieve picoDst
        StPicoDst *dst = picoReader->picoDst();
        // Retrieve event information
        StPicoEvent *event = dst->event();
        if( !event )
        {
            cout << "Something went wrong, Master! Event is hiding from me..." << endl;
            break;
        }

        hTally->Fill(0);
        if(!event->isTrigger(610001) && !event->isTrigger(610011) && !event->isTrigger(610021)
                && !event->isTrigger(610031) && !event->isTrigger(610041) && !event->isTrigger(610051) ) continue;

        Run = event->runId();
        pV  = event->primaryVertex();
        pVz = pV.Z();
        pVx = pV.X();
        pVy = pV.Y();
        VPDvz   = event->vzVpd();
        RefMult = event->refMult();
        TOFMult = event->btofTrayMultiplicity();
        Ntofmatch = event->nBTOFMatch();
        NPTracks = dst->numberOfTracks();
        BBCco   = event->BBCx();
        ZDCcoin = event->ZDCx();
        Day     = (int)((Run % 1000000) / 1000);
        Day2    = (int)((Run % 1000000) / 10);
        Day3    = (int)((Run % 1000000) / 1);

        refmultCorrUtil->init(Run);
        if ( refmultCorrUtil->isBadRun(Run) ) continue;
        hTally->Fill(1);
        refmultCorrUtil->initEvent(RefMult, pVz, ZDCcoin) ;
        if(!refmultCorrUtil->passnTofMatchRefmultCut(1.*RefMult, 1.*Ntofmatch)) continue;
        hTally->Fill(2);
        Centrality  = 1 + refmultCorrUtil->getCentralityBin9() ;
        Eweight = refmultCorrUtil->getWeight();
        if((i + 1) % 1000 == 0) cout << "centraility = " << Centrality << " Eweight = " << Eweight << endl;

        if(!IsGoodEvent(cen)) continue;
        if(!CountCharge(dst)) continue;;
        //shuffle tracks for random EPs
        int iTrack[Fcount];
        Scount = Fcount / 2 - 1;
        for(int q = 0; q < Fcount; q++) iTrack[q] = q;
        random_shuffle(iTrack, iTrack + Fcount);

        //TPC EP reconstruction
        MakeTPC_EP(dst, iTrack);
        random_shuffle(iCharge.begin(), iCharge.end()); //shuffle charge for POI
        if(TPC_EP_for < -9 || TPC_EP_bac < -9) continue;
        if(opt_ESE == 1 && (sqrt(Q2) < 0.4 || sqrt(Q2) > 1.6)) continue;
        //BBC EP
        if(!IsGoodBBC(event)) continue;
        MakeBBC_EP(event);

        //EPD EP
        StEpdEpInfo result = mEpFinder->Results(mEpdHits, pV, (Centrality > 0) ? Centrality - 1 : 0);
        EPD_EP_east = result.EastPhiWeightedPsi(nHar);
        EPD_EP_west = result.WestPhiWeightedPsi(nHar);
        if(EPD_EP_east == EPD_EP_west) continue;
        EPDe_Day3_cos2->Fill(Day3, cos(nHar * EPD_EP_east));
        EPDe_Day3_sin2->Fill(Day3, sin(nHar * EPD_EP_east));
        EPDw_Day3_cos2->Fill(Day3, cos(nHar * EPD_EP_west));
        EPDw_Day3_sin2->Fill(Day3, sin(nHar * EPD_EP_west));

        //ZDC EP
        if(!IsGoodZDC(event)) continue;
        MakeZDC_EP(event);
        //flatten EPs
        ShiftPsi();
        FillEP_resolution();
        //if(i>80000) cout<<"5th"<<endl;
        //Store the flattened phi for POI
        Phi_new.resize(NPTracks);
        for(int trki = 0; trki < NPTracks; trki++)
        {
            StPicoTrack *picoTrack = dst->track(trki);
            if(!IsGoodTrack(picoTrack)) continue;
            Pt        = picoTrack->pMom().Pt();
            Eta       = picoTrack->pMom().Eta();
            Charge    = picoTrack->charge();
            Phi       = picoTrack->pMom().Phi();
            DCAGlobal = picoTrack->gDCA(pV).Mag();

            if(!IsGoodPOI(picoTrack)) continue;
            FillPhiPOI();   //Charge is needed here
            ShiftPhiPOI(trki);
        }
        //if(i>80000) cout<<"6th"<<endl;
        //////////Real analysis begins here//////////////////////////////////////////////////////////////////////////////////////
        N_L_P1 = 0, N_L_N1 = 0, N_R_P1 = 0, N_R_N1 = 0, N_T_P1 = 0, N_T_N1 = 0, N_B_P1 = 0, N_B_N1 = 0;
        N_L_P2 = 0, N_L_N2 = 0, N_R_P2 = 0, N_R_N2 = 0, N_T_P2 = 0, N_T_N2 = 0, N_B_P2 = 0, N_B_N2 = 0;
        Fcount = 0;
        if(opt_useEPD == 1 || opt_useEPD == 2) Eweight *= 2;
        //loop for the real analysis
        for(int trki = 0; trki < NPTracks; trki++)
        {
            StPicoTrack *picoTrack = dst->track(trki);
            if(!IsGoodTrack(picoTrack)) continue;
            if(!IsGoodPOI(picoTrack)) continue;
            Pt    = picoTrack->pMom().Pt();
            Eta   = picoTrack->pMom().Eta();
            Theta     = 2.*atan(exp(-Eta));
            Charge    = picoTrack->charge();
            Phi   = picoTrack->pMom().Phi();
            ndEdx     = picoTrack->nHitsDedx();
            DCAGlobal = picoTrack->gDCA(pV).Mag();
            eff = (cen > 0) ? PP0[cen - 1] * exp(-pow(PP1[cen - 1] / Pt, PP2[cen - 1])) : PP0[0] * exp(-pow(PP1[0] / Pt, PP2[0]));
            float eff_tof = 1;
            //                      if(TOF_eff && TOF_eff->GetEntries()) eff_tof = TOF_eff->GetBinContent(TOF_eff->FindBin(Pt));
            eff *= eff_tof;

            hEtaPtDist->Fill(Eta, Pt, Eweight);
            hEtaPhiDist->Fill(Eta, Phi, Eweight);
            hPhiPtDist->Fill(Phi, Pt, Eweight);
            if(DCAGlobal < 1 && fabs(Eta) < 0.9)    // && Pt>0.2 && Pt*cosh(Eta)<2) {
            {
                Hist_Pt->Fill(Pt, Eweight);
                if(picoTrack->isTofTrack() && dst->btofPidTraits(picoTrack->bTofPidTraitsIndex())->btof() > 0) Hist_Pt_TOF->Fill(Pt, Eweight);
            }

            //          if(!IsGoodPion(dst,picoTrack,1)) continue;
            v2_sub = (Eta > 0) ? cos(nHar * Phi_new[trki] - nHar * TPC_EP_bac_new) * 100 : cos(nHar * Phi_new[trki] - nHar * TPC_EP_for_new) * 100;
            if(opt_useEPD == 1) v2_sub = (cos(nHar * Phi_new[trki] - nHar * EPD_EP_east_new) + cos(nHar * Phi_new[trki] - nHar * EPD_EP_west_new)) * 50;
            if(opt_useEPD == 2) v2_sub = (cos(nHar * Phi_new[trki] - nHar * BBC_EP_east_new) + cos(nHar * Phi_new[trki] - nHar * BBC_EP_west_new)) * 50;
            FillCMW();

            float mQx_i = mQx, mQy_i = mQy;
            if(IsGoodAsso(picoTrack))
            {
                mQx_i -= Pt * cos(PhiAsso_new[trki] * nHar);
                mQy_i -= Pt * sin(PhiAsso_new[trki] * nHar);
            }
            TVector2 mQ_i(mQx_i, mQy_i);
            float psi_F = mQ_i.Phi() / nHar;
            Hist_TPC_EP_full_m1->Fill(psi_F);
            float psi_F_new = psi_F;
            for(int jj = 0; jj < order; jj++)
                psi_F_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F) / nHar / (jj + 1);
            Hist_TPC_EP_full_m1_flat->Fill(psi_F_new);

            v2  = cos(nHar * Phi_new[trki] - nHar * psi_F_new) * 100;
            v2e = cos(nHar * Phi_new[trki] - nHar * TPC_EP_for_new) * 100;
            v2w = cos(nHar * Phi_new[trki] - nHar * TPC_EP_bac_new) * 100;
            v2_EPDe = cos(nHar * Phi_new[trki] - nHar * EPD_EP_east_new) * 100;
            v2_EPDw = cos(nHar * Phi_new[trki] - nHar * EPD_EP_west_new) * 100;
            v2_BBCe = cos(nHar * Phi_new[trki] - nHar * BBC_EP_east_new) * 100;
            v2_BBCw = cos(nHar * Phi_new[trki] - nHar * BBC_EP_west_new) * 100;
            v2_ZDC  = cos(nHar * Phi_new[trki] - ZDC_EP_east_new - ZDC_EP_west_new) * 100;

            iCh = iCharge[match[trki]];
            Rsin = sin(0.5 * nHar * (Phi_new[trki] - psi_F_new));
            Rcos = sin(0.5 * nHar * (Phi_new[trki] - psi_F_new - PI / nHar));
            if(opt_useEPD == 1)
            {
                Rsin = 0.5 * (sin(0.5 * nHar * (Phi_new[trki] - EPD_EP_east_new)) + sin(0.5 * nHar * (Phi_new[trki] - EPD_EP_west_new)));
                Rcos = 0.5 * (sin(0.5 * nHar * (Phi_new[trki] - EPD_EP_east_new - PI / nHar)) + sin(0.5 * nHar * (Phi_new[trki] - EPD_EP_west_new - PI / nHar)));
            }
            if(opt_useEPD == 2)
            {
                Rsin = 0.5 * (sin(0.5 * nHar * (Phi_new[trki] - BBC_EP_east_new)) + sin(0.5 * nHar * (Phi_new[trki] - BBC_EP_west_new)));
                Rcos = 0.5 * (sin(0.5 * nHar * (Phi_new[trki] - BBC_EP_east_new - PI / nHar)) + sin(0.5 * nHar * (Phi_new[trki] - BBC_EP_west_new - PI / nHar)));
            }

            Fillv2();
            FillRcorr();
            CountMSC(iTrack, Phi_new[trki]);

            Fcount++;
            if(opt_weight == 1) continue;
            for(int trkj = trki + 1; trkj < NPTracks; trkj++)
            {
                StPicoTrack *picoTrack2 = dst->track(trkj);
                if(!IsGoodTrack(picoTrack2)) continue;
                if(!IsGoodPOI(picoTrack2)) continue;
                Pt2        = picoTrack2->pMom().Pt();
                Eta2       = picoTrack2->pMom().Eta();
                Charge2    = picoTrack2->charge();
                Phi2       = picoTrack2->pMom().Phi();
                DCAGlobal2 = picoTrack2->gDCA(pV).Mag();
                eff2 = (cen > 0) ? PP0[cen - 1] * exp(-pow(PP1[cen - 1] / Pt2, PP2[cen - 1])) : PP0[0] * exp(-pow(PP1[0] / Pt2, PP2[0]));
                float eff_tof2 = 1;
                //                                if(TOF_eff && TOF_eff->GetEntries()) eff_tof2 = TOF_eff->GetBinContent(TOF_eff->FindBin(Pt2));
                eff2 *= eff_tof2;

                //              if(!IsGoodPion(dst,picoTrack2,1)) continue;
                hDpt->Fill(fabs(Pt - Pt2), Eweight);

                float mQx_j = mQx_i, mQy_j = mQy_i;
                if(IsGoodAsso(picoTrack2))
                {
                    mQx_j -= Pt2 * cos(PhiAsso_new[trkj] * nHar);
                    mQy_j -= Pt2 * sin(PhiAsso_new[trkj] * nHar);
                }
                TVector2 mQ_j(mQx_j, mQy_j);
                float psi_F2 = mQ_j.Phi() / nHar;
                Hist_TPC_EP_full_m2->Fill(psi_F2);
                float psi_F2_new = psi_F2;
                for(int jj = 0; jj < order; jj++)
                    psi_F2_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * psi_F2) / nHar / (jj + 1);
                if(psi_F2_new < 0) psi_F2_new += PI;
                if(psi_F2_new > PI) psi_F2_new -= PI;
                Hist_TPC_EP_full_m2_flat->Fill(psi_F2_new);

                float Phi1_new = Phi_new[trki];
                float Phi2_new = Phi_new[trkj];
                TLorentzVector lA;
                lA.SetPxPyPzE(Pt * cos(Phi1_new), Pt * sin(Phi1_new), Pt * sinh(Eta), sqrt(mpi * mpi + Pt * Pt * cosh(Eta)*cosh(Eta)));
                TLorentzVector lB;
                lB.SetPxPyPzE(Pt2 * cos(Phi2_new), Pt2 * sin(Phi2_new), Pt2 * sinh(Eta2), sqrt(mpi * mpi + Pt2 * Pt2 * cosh(Eta2)*cosh(Eta2)));
                TVector3 boost = -(lA + lB).BoostVector();
                lA.Boost(boost);
                lB.Boost(boost);
                if(opt_boost)
                {
                    Phi1_new = lA.Phi();
                    Phi2_new = lB.Phi();
                }

                //correlations
                if(opt_useEPD == 1)
                {
                    TPC_EP_for_new = EPD_EP_west_new;
                    TPC_EP_bac_new = EPD_EP_east_new;
                }
                if(opt_useEPD == 2)
                {
                    TPC_EP_for_new = BBC_EP_west_new;
                    TPC_EP_bac_new = BBC_EP_east_new;
                }
                correlator3 = cos(Phi1_new - Phi2_new) * 100;
                correlator0 = cos(Phi + (nHar - 1) * Phi2 - nHar * psi_F2_new) * 100;
                correlator4e = cos(Phi1_new + (nHar - 1) * Phi2_new - nHar * TPC_EP_for_new) * 100;
                correlator4w = cos(Phi1_new + (nHar - 1) * Phi2_new - nHar * TPC_EP_bac_new) * 100;
                correlator4 = cos(Phi1_new + (nHar - 1) * Phi2_new - nHar * psi_F2_new) * 100;
                if(opt_HighOrder == 1)
                {
                    correlator0 = cos(Phi - (nHar + 1) * Phi2 + nHar * psi_F2_new) * 100;
                    correlator4e = cos(Phi1_new - (nHar + 1) * Phi2_new + nHar * TPC_EP_for_new) * 100;
                    correlator4w = cos(Phi1_new - (nHar + 1) * Phi2_new + nHar * TPC_EP_bac_new) * 100;
                    correlator4 = cos(Phi1_new - (nHar + 1) * Phi2_new + nHar * psi_F2_new) * 100;
                }
                if(opt_useEPD == 1 || opt_useEPD == 2) correlator4 = 0.5 * (correlator4e+correlator4w);

                if(Charge > 0 && Charge2 > 0) FillGamma(1);
                if(Charge < 0 && Charge2 < 0) FillGamma(2);
                if(Charge * Charge2 > 0) FillGamma(3);
                if(Charge * Charge2 < 0) FillGamma(4);

                if(fabs(Pt - Pt2) > 0.15 && fabs(Eta - Eta2) > 0.15)
                {
                    Hist_v2_eta_obs3->Fill(Eta, v2, Eweight / eff);
                    if((opt_useEPD == 0 && fabs(Eta) < Eta_ESE_Cut && fabs(Eta2) < Eta_ESE_Cut) || opt_useEPD > 0)
                    {
                        pTemp_v2_noHBT->Fill(1, v2e, 1. / eff);
                        pTemp_v2_noHBT->Fill(2, v2w, 1. / eff);
                        if(Charge > 0)
                        {
                            pTemp_v2_noHBT->Fill(3, v2e, 1. / eff);
                            pTemp_v2_noHBT->Fill(4, v2w, 1. / eff);
                        }
                        if(Charge < 0)
                        {
                            pTemp_v2_noHBT->Fill(5, v2e, 1. / eff);
                            pTemp_v2_noHBT->Fill(6, v2w, 1. / eff);
                        }
                    }
                }
            } // 2nd track

        }  //1st Track

        v2_TPC_Day3->Fill(Day3, pTemp_v2->GetBinContent(7));
        v2_EPD_Day3->Fill(Day3, pTemp_v2->GetBinContent(8));
        v2_BBC_Day3->Fill(Day3, pTemp_v2->GetBinContent(9));
        v2_ZDC_Day3->Fill(Day3, pTemp_v2->GetBinContent(10));

        WrapUpRcorr();
        WrapUpMSC();
        WrapUpESE();

        PhiAsso_new.clear();
        Phi_new.clear();
    } // Event

    rc = (TH1D *)Hist_Pt_TOF->Clone();
    rc->SetName("rc");
    rc->Divide(Hist_Pt);
    WriteHistogram(cen, opt_weight);
    if(opt_weight == 1)
    {
        mEpFinder->Finish();
        WriteWeight(fname_new);
    }
    return;
}
//////////////////////////////////
bool IsGoodEvent(int c)
{
    hZDCcoin->Fill(ZDCcoin);
    hVertexXY->Fill(pVx, pVy);
    if((pVx) * (pVx) + (pVy) * (pVy) > Vr_cut) return false;
    hTally->Fill(3);
    hVzDiff->Fill(pVz - VPDvz);
    if(TMath::Abs(pVz - VPDvz) > Vz_diff) return false;
    hMult_Vz->Fill(RefMult, pVz);
    hVertexZ->Fill(pVz);
    hTally->Fill(4);
    if(TMath::Abs(pVz) > Vz_cut) return false;
    //  if(pVz<0) return false;
    hMult_Vz_new->Fill(RefMult, pVz, Eweight);
    Ref_TOF->Fill(RefMult, TOFMult);
    Ref_Day3->Fill(Day3, RefMult);
    TOF_Day3->Fill(Day3, TOFMult);
    NPT_Day3->Fill(Day3, NPTracks);
    hTally->Fill(5);
    int Bad = 0;
    for(int jj = 0; jj < Nrun_MB; jj++) if(Day3 == bad_Ref_day3[jj])
        {
            Bad = 1;
            break;
        }
    if(Bad) return false;  //bad run
    hTally->Fill(6);

    hCentrality->Fill(Centrality, Eweight);
    if(c && Centrality != c) return false;
    hTally->Fill(7);

    return true;
}
//////////////////////////////////////////////
bool CountCharge(StPicoDst *d)
{
    int Ntof = 0, Npos = 0, Nneg = 0;
    Fcount = 0, POIcount = 0;
    for(int trk = 0; trk < NPTracks; trk++)
    {
        StPicoTrack *picoTrack = d->track(trk);
        if(!IsGoodTrack(picoTrack)) continue;
        EtaAsso   = picoTrack->pMom().Eta();
        PtAsso    = picoTrack->pMom().Pt();
        ChargeAsso = picoTrack->charge();
        DCAGlobalAsso = picoTrack->gDCA(pV).Mag();
        nSigma_p  = picoTrack->nSigmaProton();

        if(ChargeAsso > 0 && fabs(EtaAsso) < 1 && PtAsso > 0.15 && DCAGlobalAsso < 1 && !(fabs(nSigma_p) < 3 && PtAsso < 0.4)) Npos++;
        if(ChargeAsso < 0 && fabs(EtaAsso) < 1 && PtAsso > 0.15 && DCAGlobalAsso < 1 && !(fabs(nSigma_p) < 3 && PtAsso < 0.4)) Nneg++;
        if(picoTrack->isTofTrack()) Ntof++;
        Hist_DCA->Fill(DCAGlobalAsso);

        if(IsGoodAsso(picoTrack)) Fcount++;
        if(IsGoodPOI(picoTrack))  POIcount++;
    }
    if(Ntof < 2) return false; //at least 2 tracks match TOF
    NPA_Day3->Fill(Day3, Fcount);
    hTally->Fill(9);
    int net_Nch = Npos - Nneg;
    net_Nch_Asym = -99;
    if((Npos + Nneg) > 0) net_Nch_Asym = (Npos - Nneg) / float(Npos + Nneg);
    Hist_positive->Fill(Npos);
    Hist_negative->Fill(Nneg);
    Hist_Ch->Fill(Npos + Nneg);
    Hist_netCh->Fill(net_Nch);
    if (net_Nch_Asym > -99)
    {
        Hist_netChAsym->Fill(net_Nch_Asym);
        p_netChAsym_RefMult->Fill(net_Nch_Asym, RefMult);
        if(net_Nch_Asym < (MeanNetChargeAsym - RMSNetChargeAsym)) Hist_netChAsym_bin->Fill(1, net_Nch_Asym);
        else if(net_Nch_Asym < (MeanNetChargeAsym - (0.3 * RMSNetChargeAsym))) Hist_netChAsym_bin->Fill(2, net_Nch_Asym);
        else if(net_Nch_Asym < (MeanNetChargeAsym + (0.3 * RMSNetChargeAsym))) Hist_netChAsym_bin->Fill(3, net_Nch_Asym);
        else if(net_Nch_Asym < (MeanNetChargeAsym + RMSNetChargeAsym)) Hist_netChAsym_bin->Fill(4, net_Nch_Asym);
        else Hist_netChAsym_bin->Fill(5, net_Nch_Asym);

    }
    return true;
}
/////////////////////////////////////////////
void MakeTPC_EP(StPicoDst *d, int *iTr)
{
    iCharge.resize(POIcount);
    match.resize(NPTracks);
    PhiAsso_new.resize(NPTracks);

    TVector2 mQ, mQ1, mQ2, mQ3, mQ4;
    mQx = 0., mQy = 0.;
    float mQQx = 0., mQQy = 0., mQx1 = 0., mQy1 = 0., mQx2 = 0., mQy2 = 0., mQx3 = 0., mQy3 = 0., mQx4 = 0., mQy4 = 0.;
    Fcount = 0, POIcount = 0;
    int Qcount = 0;

    for(int trk = 0; trk < NPTracks; trk++)
    {
        StPicoTrack *picoTrack = d->track(trk);
        if(!IsGoodTrack(picoTrack)) continue;
        EtaAsso   = picoTrack->pMom().Eta();
        PtAsso    = picoTrack->pMom().Pt();
        PhiAsso   = picoTrack->pMom().Phi();
        DCAGlobalAsso = picoTrack->gDCA(pV).Mag();
        ChargeAsso = picoTrack->charge();

        //to shuffle charge of POI
        match[trk] = 0;
        if(IsGoodPOI(picoTrack))
        {
            match[trk] = POIcount;
            iCharge[POIcount] = ChargeAsso;
            POIcount++;
        }
        if(!IsGoodAsso(picoTrack)) continue;

        FillPhiAsso();  //ChargeAsso is needed here
        ShiftPhiAsso(trk);

        mQx += PtAsso * cos(PhiAsso_new[trk] * nHar);
        mQy += PtAsso * sin(PhiAsso_new[trk] * nHar);
        if(fabs(EtaAsso) < Eta_ESE_Cut && IsGoodPOI(picoTrack))
        {
            mQQx += cos(PhiAsso_new[trk] * nHar);
            mQQy += sin(PhiAsso_new[trk] * nHar);
            Qcount++;
        }
        if(iTr[Fcount] > Scount)
        {
            mQx1 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy1 += PtAsso * sin(PhiAsso_new[trk] * nHar);
        }
        else
        {
            mQx2 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy2 += PtAsso * sin(PhiAsso_new[trk] * nHar);
        }
        if(EtaAsso > Eta_EP_Cut)
        {
            mQx3 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy3 += PtAsso * sin(PhiAsso_new[trk] * nHar);
        }
        if(EtaAsso < -Eta_EP_Cut)
        {
            mQx4 += PtAsso * cos(PhiAsso_new[trk] * nHar);
            mQy4 += PtAsso * sin(PhiAsso_new[trk] * nHar);
        }
        Fcount++;
    }
    mQ.Set(mQx, mQy);
    mQ1.Set(mQx1, mQy1);
    mQ2.Set(mQx2, mQy2);
    mQ3.Set(mQx3, mQy3);
    mQ4.Set(mQx4, mQy4);
    TPC_EP_full = mQ.Phi() / nHar;
    TPC_EP_east = mQ1.Phi() / nHar;
    TPC_EP_west = mQ2.Phi() / nHar;
    TPC_EP_for  = mQ3.Phi() / nHar;
    TPC_EP_bac  = mQ4.Phi() / nHar;
    if(mQx3 == 0 || mQy3 == 0) TPC_EP_for = -10;
    if(mQx4 == 0 || mQy4 == 0) TPC_EP_bac = -10;
    Q2 = (mQQx * mQQx + mQQy * mQQy) / float(Qcount);
    TPC_Day3_cos2->Fill(Day3, cos(nHar * TPC_EP_full));
    TPC_Day3_sin2->Fill(Day3, sin(nHar * TPC_EP_full));
}
///////////////////////////////////////////
bool IsGoodBBC(StPicoEvent *e)
{
    float BBC_sum_east = 0, BBC_sum_west = 0;
    for(int j = 0; j < 16; j++) BBC_sum_east += e->bbcAdcEast(j);
    for(int j = 0; j < 16; j++) BBC_sum_west += e->bbcAdcWest(j);
    if(BBC_sum_east < 150 || BBC_sum_west < 150) return false;
    int saturate = 0;
    for(int j = 0; j < 16; j++) if(e->bbcAdcEast(j) > 4500 || e->bbcAdcWest(j) > 4500)
        {
            saturate = 1;
            break;
        }
    if(saturate == 1) return false;
    hTally->Fill(10);

    return true;
}
//////////////////////////////////////////////
void MakeBBC_EP(StPicoEvent *ev)
{
    for(int j = 0; j < 6; j++)
    {
        BBC1->Fill(j + 1, ev->bbcAdcEast(j) / 1000.);
        BBC8->Fill(j + 1, ev->bbcAdcWest(j) / 1000.);
    }
    for(int j = 0; j < 10; j++)
    {
        BBC2->Fill(j + 1, ev->bbcAdcEast(j + 6) / 1000.);
        BBC7->Fill(j + 1, ev->bbcAdcWest(j + 6) / 1000.);
    }
    TVector2 mBe, mBw;
    float mBe_x = 0, mBe_y = 0, mBw_x = 0, mBw_y = 0;
    for(int trk = 0; trk < 16; trk++)
    {
        float sig_bbc = ev->bbcAdcEast(trk);
        float phi_bbc = GetPhiInBBC(0, trk + 1);
        mBe_x += cos(nHar * phi_bbc) * sig_bbc * BBC_gain_east[trk];
        mBe_y += sin(nHar * phi_bbc) * sig_bbc * BBC_gain_east[trk];
    }
    for(int trk = 0; trk < 16; trk++)
    {
        float sig_bbc = ev->bbcAdcWest(trk);
        float phi_bbc = GetPhiInBBC(1, trk + 1);
        mBw_x += cos(nHar * phi_bbc) * sig_bbc * BBC_gain_west[trk];
        mBw_y += sin(nHar * phi_bbc) * sig_bbc * BBC_gain_west[trk];
    }
    if(mBe_x == 0 || mBe_y == 0 || mBw_x == 0 || mBw_y == 0) return;
    mBe.Set(mBe_x, mBe_y);
    mBw.Set(mBw_x, mBw_y);
    BBC_EP_east = mBe.Phi() / nHar;
    BBC_EP_west = mBw.Phi() / nHar;
    BBCe_Day3_cos2->Fill(Day3, cos(nHar * BBC_EP_east));
    BBCe_Day3_sin2->Fill(Day3, sin(nHar * BBC_EP_east));
    BBCw_Day3_cos2->Fill(Day3, cos(nHar * BBC_EP_west));
    BBCw_Day3_sin2->Fill(Day3, sin(nHar * BBC_EP_west));
}
//////////////////////////////////////
bool IsGoodZDC(StPicoEvent *e)
{

    return true;
}
/////////////////////////////////////////////
void MakeZDC_EP(StPicoEvent *ev)
{

    for(int j = 0; j < 8; j++)
    {
        ZDC_e_h->Fill(j + 1, ev->ZdcSmdEastHorizontal(j) / 1000.);
        ZDC_e_v->Fill(j + 1, ev->ZdcSmdEastVertical(j) / 1000.);
        ZDC_w_h->Fill(j + 1, ev->ZdcSmdWestHorizontal(j) / 1000.);
        ZDC_w_v->Fill(j + 1, ev->ZdcSmdWestVertical(j) / 1000.);
    }

    float eXsum = 0., eYsum = 0., eXWgt = 0., eYWgt = 0.;
    float wXsum = 0., wYsum = 0., wXWgt = 0., wYWgt = 0.;
    for(int v = 0; v <= 6; v++)
    {
        eXsum += GetXYInZDC(0, 0, v + 1, 1) * ZDC_gain_east_v[v] * ev->ZdcSmdEastVertical(v);
        eXWgt += ZDC_gain_east_v[v] * ev->ZdcSmdEastVertical(v);
        wXsum += GetXYInZDC(1, 0, v + 1, 1) * ZDC_gain_west_v[v] * ev->ZdcSmdWestVertical(v);
        wXWgt += ZDC_gain_west_v[v] * ev->ZdcSmdWestVertical(v);
    }
    for(int v = 0; v <= 7; v++)
    {
        eYsum += GetXYInZDC(0, 1, v + 1, 1) * ZDC_gain_east_h[v] * ev->ZdcSmdEastHorizontal(v);
        eYWgt += ZDC_gain_east_h[v] * ev->ZdcSmdEastHorizontal(v);
        wYsum += GetXYInZDC(1, 1, v + 1, 1) * ZDC_gain_west_h[v] * ev->ZdcSmdWestHorizontal(v);
        wYWgt += ZDC_gain_west_h[v] * ev->ZdcSmdWestHorizontal(v);
    }
    if(eXWgt) pZDCcenter->Fill(1, Day2, eXsum / eXWgt);
    if(eYWgt) pZDCcenter->Fill(2, Day2, eYsum / eYWgt);
    if(wXWgt) pZDCcenter->Fill(3, Day2, wXsum / wXWgt);
    if(wYWgt) pZDCcenter->Fill(4, Day2, wYsum / wYWgt);

    if(Weight_Read && Read_ZDCcenter->ProjectionY("center1", 1, 1)->GetBinContent(Day2 - run_sta / 10 + 1) > 0)
    {
        zdcsmd_x0_e = Read_ZDCcenter->ProjectionY("center1", 1, 1)->GetBinContent(Day2 - run_sta / 10 + 1);
        zdcsmd_y0_e = Read_ZDCcenter->ProjectionY("center2", 2, 2)->GetBinContent(Day2 - run_sta / 10 + 1);
        zdcsmd_x0_w = Read_ZDCcenter->ProjectionY("center3", 3, 3)->GetBinContent(Day2 - run_sta / 10 + 1);
        zdcsmd_y0_w = Read_ZDCcenter->ProjectionY("center4", 4, 4)->GetBinContent(Day2 - run_sta / 10 + 1);
    }
    eXsum = 0., eYsum = 0., eXWgt = 0., eYWgt = 0.;
    wXsum = 0., wYsum = 0., wXWgt = 0., wYWgt = 0.;
    for(int v = 0; v <= 6; v++)
    {
        eXsum += GetXYInZDC(0, 0, v + 1, 0) * ZDC_gain_east_v[v] * ev->ZdcSmdEastVertical(v);
        eXWgt += ZDC_gain_east_v[v] * ev->ZdcSmdEastVertical(v);
        wXsum += GetXYInZDC(1, 0, v + 1, 0) * ZDC_gain_west_v[v] * ev->ZdcSmdWestVertical(v);
        wXWgt += ZDC_gain_west_v[v] * ev->ZdcSmdWestVertical(v);
    }
    for(int v = 0; v <= 7; v++)
    {
        eYsum += GetXYInZDC(0, 1, v + 1, 0) * ZDC_gain_east_h[v] * ev->ZdcSmdEastHorizontal(v);
        eYWgt += ZDC_gain_east_h[v] * ev->ZdcSmdEastHorizontal(v);
        wYsum += GetXYInZDC(1, 1, v + 1, 0) * ZDC_gain_west_h[v] * ev->ZdcSmdWestHorizontal(v);
        wYWgt += ZDC_gain_west_h[v] * ev->ZdcSmdWestHorizontal(v);
    }
    //these 4 profiles should all give zero
    if(eXWgt) pZDCcenter_new->Fill(1, Day2, eXsum / eXWgt);
    if(eYWgt) pZDCcenter_new->Fill(2, Day2, eYsum / eYWgt);
    if(wXWgt) pZDCcenter_new->Fill(3, Day2, wXsum / wXWgt);
    if(wYWgt) pZDCcenter_new->Fill(4, Day2, wYsum / wYWgt);
    ZDC_EP_east = atan2(eYsum / eYWgt, eXsum / eXWgt);
    ZDC_EP_west = atan2(wYsum / wYWgt, wXsum / wXWgt);
    ZDCe_Day3_cos->Fill(Day3, cos(ZDC_EP_east));
    ZDCe_Day3_sin->Fill(Day3, sin(ZDC_EP_east));
    ZDCw_Day3_cos->Fill(Day3, cos(ZDC_EP_west));
    ZDCw_Day3_sin->Fill(Day3, sin(ZDC_EP_west));
}
/////////////////////////////////////////////////
void FillEP_resolution()
{
    float cos_ew = cos(nHar * TPC_EP_east_new - nHar * TPC_EP_west_new);
    float cos_fb = cos(nHar * TPC_EP_for_new - nHar * TPC_EP_bac_new);

    Hist_cos->Fill(1, cos_fb, Eweight);
    Hist_cos->Fill(2, cos_ew, Eweight);
    Hist_cos->Fill(3, cos(nHar * TPC_EP_for_new - nHar * TPC_EP_full_new));
    Hist_cos->Fill(4, cos(nHar * TPC_EP_bac_new - nHar * TPC_EP_full_new));
    Hist_cos_BBC->Fill(1, cos(nHar * BBC_EP_east_new - nHar * BBC_EP_west_new), Eweight);
    Hist_cos_BBC->Fill(2, cos(nHar * BBC_EP_east_new - nHar * TPC_EP_full_new), Eweight);
    Hist_cos_BBC->Fill(3, cos(nHar * BBC_EP_west_new - nHar * TPC_EP_full_new), Eweight);
    Hist_cos_EPD->Fill(1, cos(nHar * EPD_EP_east_new - nHar * EPD_EP_west_new), Eweight);
    Hist_cos_EPD->Fill(2, cos(nHar * EPD_EP_east_new - nHar * TPC_EP_full_new), Eweight);
    Hist_cos_EPD->Fill(3, cos(nHar * EPD_EP_west_new - nHar * TPC_EP_full_new), Eweight);
    Hist_cos_ZDC->Fill(1, cos(ZDC_EP_east_new - ZDC_EP_west_new), Eweight);
    BBC_cor_Day3->Fill(Day3, cos(nHar * BBC_EP_east_new - nHar * BBC_EP_west_new));
    EPD_cor_Day3->Fill(Day3, cos(nHar * EPD_EP_east_new - nHar * EPD_EP_west_new));
    ZDC_cor_Day3->Fill(Day3, cos(ZDC_EP_east_new - ZDC_EP_west_new));

    Hist_Q->Fill(sqrt(Q2));
    Hist_Q2->Fill(Q2);
    if(opt_useEPD == 1) cos_fb = cos(nHar * EPD_EP_east_new - nHar * EPD_EP_west_new);
    if(opt_useEPD == 2) cos_fb = cos(nHar * BBC_EP_east_new - nHar * BBC_EP_west_new);
    p_cos_Q2->Fill(Q2, cos_fb, Eweight);
    p_RefMult_Q2->Fill(Q2, Fcount);
    Hist_RefMult_Q2->Fill(Q2, Fcount);

    if(net_Nch_Asym > -99)
    {
        p_netChAsym_cos->Fill(net_Nch_Asym, cos_fb);
        if(net_Nch_Asym < (MeanNetChargeAsym - RMSNetChargeAsym)) Hist_cos_Ach->Fill(1, cos_fb, Eweight);
        else if(net_Nch_Asym < (MeanNetChargeAsym - (0.3 * RMSNetChargeAsym))) Hist_cos_Ach->Fill(2, cos_fb, Eweight);
        else if(net_Nch_Asym < (MeanNetChargeAsym + (0.3 * RMSNetChargeAsym))) Hist_cos_Ach->Fill(3, cos_fb, Eweight);
        else if(net_Nch_Asym < (MeanNetChargeAsym + RMSNetChargeAsym)) Hist_cos_Ach->Fill(4, cos_fb, Eweight);
        else Hist_cos_Ach->Fill(5, cos_fb, Eweight);
    }
}
/////////////////////////////////////////////////
bool IsGoodAsso(StPicoTrack *p)
{
    if(p->nHitsFit() < nFitHits_asso) return false;
    if(p->pMom().Pt() > pt_asso_up || p->pMom().Pt() < pt_asso_lo) return false;
    if(p->gDCA(pV).Mag() > DcaCut_asso) return false;
    if(fabs(p->pMom().Eta()) >= EtaCut_asso) return false;
    return true;
}
/////////////////////////////////////////////////////
bool IsGoodPOI(StPicoTrack *p)
{
    int nHF = p->nHitsFit();
    if(nHF < nFitHits_trig) return false;
    int nHM = p->nHitsMax();
    if(nHF / float(nHM) < nFitRatio_trig_lo || nHF / float(nHM) >= nFitRatio_trig_up) return false;
    if(p->pMom().Pt() > pt_trig_up || p->pMom().Pt() < pt_trig_lo) return false;
    if(p->gDCA(pV).Mag() > DcaCut_trig) return false;
    if(fabs(p->pMom().Eta()) >= EtaCut_trig) return false;
    return true;
}
/////////////////////////////////////////////////
bool IsGoodTrack(StPicoTrack *p)
{
    if(!p->isPrimary()) return false;
    return true;
}
/////////////////////////////////////////////////
bool IsGoodPionNoTOF(StPicoTrack *p)
{
    if(p->gDCA(pV).Mag() > 1) return false;
    float nSig_i = p->nSigmaPion();
    float ndEdx_i = p->nHitsDedx();
    if(ndEdx_i < 10 || nSig_i > 2 || nSig_i < -2) return false;
    return true;
}
//////////////////////////////////////////////////
bool IsGoodPion(StPicoDst *d, StPicoTrack *p, int opt)
{
    //hTall->Fill(1);
    if(p->gDCA(pV).Mag() > 1)   return false;
    //hTall->Fill(2);
    float eta_i = p->pMom().Eta();
    if(fabs(eta_i) > 0.9) return false;
    //hTall->Fill(3);
    float pt_i =  p->pMom().Pt();
    float p_i = pt_i * cosh(eta_i);
    if(pt_i < 0.2 || p_i > 1.6) return false;
    //hTall->Fill(4);
    float nSig_i = p->nSigmaPion();
    float ndEdx_i = p->nHitsDedx();
    if(ndEdx_i < 15 || nSig_i > 2 || nSig_i < -2) return false;
    //hTall->Fill(5);
    if(opt == 0) return true;
    if(!(p->isTofTrack())) return false;
    //hTall->Fill(6);
    StPicoBTofPidTraits *trait = d->btofPidTraits( p->bTofPidTraitsIndex() );
    if(!trait) return false;
    //hTall->Fill(7);
    if(trait->btof() <= 0) return false;
    //hTall->Fill(8);
    if(fabs(trait->btofYLocal()) > 1.8) return false;
    //hTall->Fill(9);
    float beta_i = trait->btofBeta();
    if(beta_i == 0) return false;
    //hTall->Fill(10);
    float mass2_i = p_i * p_i * (1.0 / beta_i / beta_i - 1.0);
    if(mass2_i < -0.01 || mass2_i > 0.1) return false;
    //hTall->Fill(11);
    return true;
}
///////////////////////////////////////////////////
bool IsGoodKaon(StPicoDst *d, StPicoTrack *p, int opt)
{
    if(p->gDCA(pV).Mag() > 1) return false;
    float eta_i = p->pMom().Eta();
    if(fabs(eta_i) > 0.9) return false;
    float pt_i =  p->pMom().Pt();
    float p_i = pt_i * cosh(eta_i);
    if(pt_i < 0.2 || p_i > 1.6) return false;

    float nSig_i = p->nSigmaKaon();
    float ndEdx_i = p->nHitsDedx();
    if(ndEdx_i < 15 || nSig_i > 2 || nSig_i < -2) return false;

    if(opt == 0) return true;
    if(!(p->isTofTrack())) return false;
    StPicoBTofPidTraits *trait = d->btofPidTraits( p->bTofPidTraitsIndex() );
    if(!trait) return false;
    if(trait->btof() <= 0) return false;
    if(fabs(trait->btofYLocal()) > 1.8) return false;
    float beta_i = trait->btofBeta();
    if(beta_i == 0) return false;
    float mass2_i = p_i * p_i * (1.0 / beta_i / beta_i - 1.0);
    if(mass2_i < 0.2 || mass2_i > 0.35) return false;
    return true;
}
////////////////////////////////////////////////////
bool IsGoodProton(StPicoDst *d, StPicoTrack *p, int opt)
{
    if(p->gDCA(pV).Mag() > 1) return false;
    float eta_i = p->pMom().Eta();
    if(fabs(eta_i) > 0.9) return false;
    float pt_i =  p->pMom().Pt();
    float p_i = pt_i * cosh(eta_i);
    if(pt_i < 0.4 || p_i > 2) return false;

    float nSig_i = p->nSigmaProton();
    float ndEdx_i = p->nHitsDedx();
    if(ndEdx_i < 15 || nSig_i > 2 || nSig_i < -2) return false;

    if(opt == 0) return true;
    if(!(p->isTofTrack())) return false;
    StPicoBTofPidTraits *trait = d->btofPidTraits( p->bTofPidTraitsIndex() );
    if(!trait) return false;
    if(trait->btof() <= 0) return false;
    if(fabs(trait->btofYLocal()) > 1.8) return false;
    float beta_i = trait->btofBeta();
    if(beta_i == 0) return false;
    float mass2_i = p_i * p_i * (1.0 / beta_i / beta_i - 1.0);
    if(mass2_i < 0.8 || mass2_i > 1) return false;
    return true;
}
////////////////////////////////////
void ShiftPsi()
{
    Hist_TPC_EP_full->Fill(TPC_EP_full, Day);
    Hist_TPC_EP_east->Fill(TPC_EP_east, Day);
    Hist_TPC_EP_west->Fill(TPC_EP_west, Day);
    Hist_TPC_EP_for->Fill(TPC_EP_for, Day);
    Hist_TPC_EP_bac->Fill(TPC_EP_bac, Day);
    Hist_BBC_EP_east->Fill(BBC_EP_east, Day);
    Hist_BBC_EP_west->Fill(BBC_EP_west, Day);
    Hist_EPD_EP_east->Fill(EPD_EP_east, Day);
    Hist_EPD_EP_west->Fill(EPD_EP_west, Day);
    Hist_ZDC_EP_east->Fill(ZDC_EP_east, Day);
    Hist_ZDC_EP_west->Fill(ZDC_EP_west, Day);

    for(int kk = 0; kk < order; kk++)
    {
        pTPC_EP_full->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*TPC_EP_full), Eweight);
        pTPC_EP_full->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*TPC_EP_full), Eweight);
        pTPC_EP_east->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*TPC_EP_east), Eweight);
        pTPC_EP_east->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*TPC_EP_east), Eweight);
        pTPC_EP_west->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*TPC_EP_west), Eweight);
        pTPC_EP_west->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*TPC_EP_west), Eweight);
        pTPC_EP_for->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*TPC_EP_for), Eweight);
        pTPC_EP_for->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*TPC_EP_for), Eweight);
        pTPC_EP_bac->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*TPC_EP_bac), Eweight);
        pTPC_EP_bac->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*TPC_EP_bac), Eweight);
        pBBC_EP_east->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*BBC_EP_east), Eweight);
        pBBC_EP_east->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*BBC_EP_east), Eweight);
        pBBC_EP_west->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*BBC_EP_west), Eweight);
        pBBC_EP_west->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*BBC_EP_west), Eweight);
        pEPD_EP_east->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*EPD_EP_east), Eweight);
        pEPD_EP_east->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*EPD_EP_east), Eweight);
        pEPD_EP_west->Fill(1 + 2 * kk, Day2, cos(nHar * (kk + 1)*EPD_EP_west), Eweight);
        pEPD_EP_west->Fill(2 + 2 * kk, Day2, sin(nHar * (kk + 1)*EPD_EP_west), Eweight);
        pZDC_EP_east->Fill(1 + 2 * kk, Day2, cos((kk + 1)*ZDC_EP_east), Eweight);
        pZDC_EP_east->Fill(2 + 2 * kk, Day2, sin((kk + 1)*ZDC_EP_east), Eweight);
        pZDC_EP_west->Fill(1 + 2 * kk, Day2, cos((kk + 1)*ZDC_EP_west), Eweight);
        pZDC_EP_west->Fill(2 + 2 * kk, Day2, sin((kk + 1)*ZDC_EP_west), Eweight);
    }

    if(Weight_Read && Read_TPC_EP_full->GetEntries())
    {
        for(int k = 0; k < 2 * order; k++)
        {
            PsiMean_F[k] = Read_TPC_EP_full->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_E[k] = Read_TPC_EP_east->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_W[k] = Read_TPC_EP_west->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_f[k] = Read_TPC_EP_for->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            PsiMean_b[k] = Read_TPC_EP_bac->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            Psi_BBC_E[k] = Read_BBC_EP_east->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            Psi_BBC_W[k] = Read_BBC_EP_west->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            Psi_EPD_E[k] = Read_EPD_EP_east->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            Psi_EPD_W[k] = Read_EPD_EP_west->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            Psi_ZDC_E[k] = Read_ZDC_EP_east->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
            Psi_ZDC_W[k] = Read_ZDC_EP_west->GetBinContent(k + 1, Day2 - run_sta / 10 + 1);
        }
    }

    TPC_EP_full_new = TPC_EP_full, TPC_EP_east_new = TPC_EP_east, TPC_EP_west_new = TPC_EP_west;
    TPC_EP_for_new = TPC_EP_for, TPC_EP_bac_new = TPC_EP_bac;
    BBC_EP_east_new = BBC_EP_east, BBC_EP_west_new = BBC_EP_west;
    EPD_EP_east_new = EPD_EP_east, EPD_EP_west_new = EPD_EP_west;
    ZDC_EP_east_new = ZDC_EP_east, ZDC_EP_west_new = ZDC_EP_west;
    for(int jj = 0; jj < order; jj++)
    {
        TPC_EP_full_new += -2 * PsiMean_F[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_full) / nHar / (jj + 1) + 2 * PsiMean_F[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_full) / nHar / (jj + 1);
        TPC_EP_east_new += -2 * PsiMean_E[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_east) / nHar / (jj + 1) + 2 * PsiMean_E[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_east) / nHar / (jj + 1);
        TPC_EP_west_new += -2 * PsiMean_W[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_west) / nHar / (jj + 1) + 2 * PsiMean_W[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_west) / nHar / (jj + 1);
        TPC_EP_for_new  += -2 * PsiMean_f[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_for) / nHar / (jj + 1)  + 2 * PsiMean_f[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_for) / nHar / (jj + 1);
        TPC_EP_bac_new  += -2 * PsiMean_b[1 + 2 * jj] * cos(nHar * (jj + 1) * TPC_EP_bac) / nHar / (jj + 1)  + 2 * PsiMean_b[0 + 2 * jj] * sin(nHar * (jj + 1) * TPC_EP_bac) / nHar / (jj + 1);
        BBC_EP_east_new += -2 * Psi_BBC_E[1 + 2 * jj] * cos(nHar * (jj + 1) * BBC_EP_east) / nHar / (jj + 1) + 2 * Psi_BBC_E[0 + 2 * jj] * sin(nHar * (jj + 1) * BBC_EP_east) / nHar / (jj + 1);
        BBC_EP_west_new += -2 * Psi_BBC_W[1 + 2 * jj] * cos(nHar * (jj + 1) * BBC_EP_west) / nHar / (jj + 1) + 2 * Psi_BBC_W[0 + 2 * jj] * sin(nHar * (jj + 1) * BBC_EP_west) / nHar / (jj + 1);
        EPD_EP_east_new += -2 * Psi_EPD_E[1 + 2 * jj] * cos(nHar * (jj + 1) * EPD_EP_east) / nHar / (jj + 1) + 2 * Psi_EPD_E[0 + 2 * jj] * sin(nHar * (jj + 1) * EPD_EP_east) / nHar / (jj + 1);
        EPD_EP_west_new += -2 * Psi_EPD_W[1 + 2 * jj] * cos(nHar * (jj + 1) * EPD_EP_west) / nHar / (jj + 1) + 2 * Psi_EPD_W[0 + 2 * jj] * sin(nHar * (jj + 1) * EPD_EP_west) / nHar / (jj + 1);
        ZDC_EP_east_new += -2 * Psi_ZDC_E[1 + 2 * jj] * cos((jj + 1) * ZDC_EP_east) / (jj + 1) + 2 * Psi_ZDC_E[0 + 2 * jj] * sin((jj + 1) * ZDC_EP_east) / (jj + 1);
        ZDC_EP_west_new += -2 * Psi_ZDC_W[1 + 2 * jj] * cos((jj + 1) * ZDC_EP_west) / (jj + 1) + 2 * Psi_ZDC_W[0 + 2 * jj] * sin((jj + 1) * ZDC_EP_west) / (jj + 1);
    }
    if(TPC_EP_full_new > 2.*PI / nHar) TPC_EP_full_new -= 2.*PI / nHar;
    if(TPC_EP_full_new < 0) TPC_EP_full_new += 2.*PI / nHar;
    if(TPC_EP_east_new > 2.*PI / nHar) TPC_EP_east_new -= 2.*PI / nHar;
    if(TPC_EP_east_new < 0) TPC_EP_east_new += 2.*PI / nHar;
    if(TPC_EP_west_new > 2.*PI / nHar) TPC_EP_west_new -= 2.*PI / nHar;
    if(TPC_EP_west_new < 0) TPC_EP_west_new += 2.*PI / nHar;
    if(TPC_EP_for_new > 2.*PI / nHar)  TPC_EP_for_new  -= 2.*PI / nHar;
    if(TPC_EP_for_new < 0)  TPC_EP_for_new  += 2.*PI / nHar;
    if(TPC_EP_bac_new > 2.*PI / nHar)  TPC_EP_bac_new  -= 2.*PI / nHar;
    if(TPC_EP_bac_new < 0)  TPC_EP_bac_new  += 2.*PI / nHar;
    if(BBC_EP_east_new > 2.*PI / nHar) BBC_EP_east_new -= 2.*PI / nHar;
    if(BBC_EP_east_new < 0) BBC_EP_east_new += 2.*PI / nHar;
    if(BBC_EP_west_new > 2.*PI / nHar) BBC_EP_west_new -= 2.*PI / nHar;
    if(BBC_EP_west_new < 0) BBC_EP_west_new += 2.*PI / nHar;
    if(EPD_EP_east_new > 2.*PI / nHar) EPD_EP_east_new -= 2.*PI / nHar;
    if(EPD_EP_east_new < 0) EPD_EP_east_new += 2.*PI / nHar;
    if(EPD_EP_west_new > 2.*PI / nHar) EPD_EP_west_new -= 2.*PI / nHar;
    if(EPD_EP_west_new < 0) EPD_EP_west_new += 2.*PI / nHar;
    if(ZDC_EP_east_new > PI) ZDC_EP_east_new -= 2 * PI;
    if(ZDC_EP_east_new < -PI) ZDC_EP_east_new += 2 * PI;
    if(ZDC_EP_west_new > PI) ZDC_EP_west_new -= 2 * PI;
    if(ZDC_EP_west_new < -PI) ZDC_EP_west_new += 2 * PI;
    Hist_TPC_EP_east_flat->Fill(TPC_EP_east_new, Day);
    Hist_TPC_EP_west_flat->Fill(TPC_EP_west_new, Day);
    Hist_TPC_EP_for_flat->Fill(TPC_EP_for_new, Day);
    Hist_TPC_EP_bac_flat->Fill(TPC_EP_bac_new, Day);
    Hist_TPC_EP_full_flat->Fill(TPC_EP_full_new, Day);
    Hist_BBC_EP_east_flat->Fill(BBC_EP_east_new, Day);
    Hist_BBC_EP_west_flat->Fill(BBC_EP_west_new, Day);
    Hist_EPD_EP_east_flat->Fill(EPD_EP_east_new, Day);
    Hist_EPD_EP_west_flat->Fill(EPD_EP_west_new, Day);
    Hist_ZDC_EP_east_flat->Fill(ZDC_EP_east_new, Day);
    Hist_ZDC_EP_west_flat->Fill(ZDC_EP_west_new, Day);

    float BBC_EP_full = atan2(sin(nHar * BBC_EP_east_new) + sin(nHar * BBC_EP_west_new), cos(nHar * BBC_EP_east_new) + cos(nHar * BBC_EP_west_new)) / nHar;
    if(BBC_EP_full < 0) BBC_EP_full += 2.*PI / nHar;
    float EPD_EP_full = atan2(sin(nHar * EPD_EP_east_new) + sin(nHar * EPD_EP_west_new), cos(nHar * EPD_EP_east_new) + cos(nHar * EPD_EP_west_new)) / nHar;
    if(EPD_EP_full < 0) EPD_EP_full += 2.*PI / nHar;
    Hist_BBC_vs_TPC->Fill(TPC_EP_full_new, BBC_EP_full);
    Hist_EPD_vs_TPC->Fill(TPC_EP_full_new, EPD_EP_full);
}
//////////////////////////////////
void ShiftPhiPOI(int tr)
{
    int index = 0, index2 = 0;
    if(pVz > 0) index2 = (Charge > 0) ? 1 : 2;
    else index2 = (Charge > 0) ? 3 : 4;
    int Eta_index = int((Eta + 1) * 5);
    if(Weight_Read && (ReadWeightHist[Eta_index].pTPCmeanPhi_1->GetEntries()))
    {
        for(int kk = 0; kk < order; kk++)
        {
            if(pVz > 0) index = (Charge > 0) ? 1 + 8 * kk  : 3 + 8 * kk;
            else      index = (Charge > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;

            if(Pt < 0.5)
            {
                PhiMean_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_1->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                PhiMean_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_1->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
            }
            else if(Pt < 1)
            {
                PhiMean_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_2->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                PhiMean_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_2->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
            }
            else
            {
                PhiMean_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_3->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                PhiMean_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhi_3->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
            }
        }
    }

    Phi_new[tr] = Phi;   //store the shifted angles
    for(int jj = 0; jj < order; jj++) Phi_new[tr] += -2 * PhiMean_sin[jj] * cos(jj * Phi + Phi) / (jj + 1) + 2 * PhiMean_cos[jj] * sin(jj * Phi + Phi) / (jj + 1);
    if(Phi_new[tr] > PI) Phi_new[tr] -= 2 * PI;
    if(Phi_new[tr] < -PI) Phi_new[tr] += 2 * PI;

    if(Pt < 0.5) PhiHist[Eta_index].Hist_Phi_new_1->Fill(Phi_new[tr], index2, Eweight);
    else if(Pt < 1) PhiHist[Eta_index].Hist_Phi_new_2->Fill(Phi_new[tr], index2, Eweight);
    else PhiHist[Eta_index].Hist_Phi_new_3->Fill(Phi_new[tr], index2, Eweight);
}
//////////////////////////////////////////
void ShiftPhiAsso(int tr)
{
    int index = 0;
    int Eta_index = int((EtaAsso + 1) * 5);
    if(Weight_Read && (ReadWeightHist[Eta_index].pTPCmeanPhiAsso_1->GetEntries()))
    {
        for(int kk = 0; kk < order; kk++)
        {
            if(pVz > 0) index = (ChargeAsso > 0) ? 1 + 8 * kk  : 3 + 8 * kk;
            else      index = (ChargeAsso > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;

            if(PtAsso < 0.5)
            {
                PhiMeanAsso_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_1->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                PhiMeanAsso_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_1->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
            }
            else if(PtAsso < 1)
            {
                PhiMeanAsso_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_2->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                PhiMeanAsso_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_2->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
            }
            else
            {
                PhiMeanAsso_cos[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_3->GetBinContent(index,  Day2 - run_sta / 10 + 1);
                PhiMeanAsso_sin[kk] = ReadWeightHist[Eta_index].pTPCmeanPhiAsso_3->GetBinContent(index + 1, Day2 - run_sta / 10 + 1);
            }
        }
    }

    PhiAsso_new[tr] = PhiAsso;
    for(int jj = 0; jj < order; jj++) PhiAsso_new[tr] += -2 * PhiMeanAsso_sin[jj] * cos(jj * PhiAsso + PhiAsso) / (jj + 1) + 2 * PhiMeanAsso_cos[jj] * sin(jj * PhiAsso + PhiAsso) / (jj + 1);
}
/////////////////////////////////////
void FillPhiAsso()    // shift parameters for Particles of EP
{
    int index = 0;
    int Eta_index = int((EtaAsso + 1) * 5);
    for(int kk = 0; kk < order; kk++)
    {
        if(pVz > 0) index = (ChargeAsso > 0) ? 1 + 8 * kk  : 3 + 8 * kk;
        else      index = (ChargeAsso > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;

        if(PtAsso < 0.5)
        {
            WeightHist[Eta_index].pTPCmeanPhiAsso_1->Fill(index,  Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
            WeightHist[Eta_index].pTPCmeanPhiAsso_1->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
        }
        else if(PtAsso < 1)
        {
            WeightHist[Eta_index].pTPCmeanPhiAsso_2->Fill(index,  Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
            WeightHist[Eta_index].pTPCmeanPhiAsso_2->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
        }
        else
        {
            WeightHist[Eta_index].pTPCmeanPhiAsso_3->Fill(index,  Day2, cos(kk * PhiAsso + PhiAsso), Eweight);
            WeightHist[Eta_index].pTPCmeanPhiAsso_3->Fill(index + 1, Day2, sin(kk * PhiAsso + PhiAsso), Eweight);
        }
    }
}
///////////////////////////////////
void FillPhiPOI()   // shift parameters for Particles of Interest
{
    int index = 0, index2 = 0;
    if(pVz > 0) index2 = (Charge > 0) ? 1 : 2;
    else index2 = (Charge > 0) ? 3 : 4;
    int Eta_index = int((Eta + 1) * 5);
    for(int kk = 0; kk < order; kk++)
    {
        if(pVz > 0) index = (Charge > 0) ? 1 + 8 * kk  : 3 + 8 * kk;
        else      index = (Charge > 0) ? 1 + 4 + 8 * kk : 3 + 4 + 8 * kk;
        if(Pt < 0.5)
        {
            WeightHist[Eta_index].pTPCmeanPhi_1->Fill(index,  Day2, cos(kk * Phi + Phi), Eweight);
            WeightHist[Eta_index].pTPCmeanPhi_1->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
            if(kk == 0) PhiHist[Eta_index].Hist_Phi_1->Fill(Phi, index2, Eweight);
        }
        else if(Pt < 1)
        {
            WeightHist[Eta_index].pTPCmeanPhi_2->Fill(index,  Day2, cos(kk * Phi + Phi), Eweight);
            WeightHist[Eta_index].pTPCmeanPhi_2->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
            if(kk == 0) PhiHist[Eta_index].Hist_Phi_2->Fill(Phi, index2, Eweight);
        }
        else
        {
            WeightHist[Eta_index].pTPCmeanPhi_3->Fill(index,  Day2, cos(kk * Phi + Phi), Eweight);
            WeightHist[Eta_index].pTPCmeanPhi_3->Fill(index + 1, Day2, sin(kk * Phi + Phi), Eweight);
            if(kk == 0) PhiHist[Eta_index].Hist_Phi_3->Fill(Phi, index2, Eweight);
        }
    }
}
/////////////////////////////////
void WriteHistogram(int c, int o)
{
    char fname_out[200];
    sprintf(fname_out, "cen%d.gamma1%d%d_EP%d_Boost%d.root", c, nHar - 1, nHar, opt_useEPD, opt_boost);
    if(opt_HighOrder == 1) sprintf(fname_out, "cen%d.gamma1%d%d_EP%d_Boost%d.root", c, nHar + 1, nHar, opt_useEPD, opt_boost);
    if(o == 1) sprintf(fname_out, "cen%d.v%d.root", c, nHar);
    TFile *fout = new TFile(fname_out, "RECREATE");

    //        gROOT->GetList()->ls();
    TList *list = gROOT->GetList();//GetListOfKeys();
    TIter next(list);
    TKey *key;
    TObject *obj;
    while ((key = (TKey *)next()))
    {
        TString tempStr(key->GetName());
        if (tempStr.Contains("Temp")) continue;
        if(o == 1 && (tempStr.Contains("Parity") || tempStr.Contains("Delta"))) continue;
        obj = gROOT->Get(key->GetName());
        if (!obj) continue;
        if(obj->IsA() == TDirectory::Class())
        {
            delete obj;
            obj = NULL;
            continue;
        }
        obj->Write();
    }
    fout->Write();
    fout->Close();
}
//////////////////////////////////////////////////
void WriteWeight(char *OutFileName)
{
    TFile *fWgtNew = new TFile(OutFileName, "UPDATE");
    for(int i = 0; i < 10; i++)
    {
        WeightHist[i].pTPCmeanPhi_1->Write();
        WeightHist[i].pTPCmeanPhi_2->Write();
        WeightHist[i].pTPCmeanPhi_3->Write();
        WeightHist[i].pTPCmeanPhiAsso_1->Write();
        WeightHist[i].pTPCmeanPhiAsso_2->Write();
        WeightHist[i].pTPCmeanPhiAsso_3->Write();
    }
    BBC1->Write();
    BBC2->Write();
    BBC7->Write();
    BBC8->Write();
    ZDC_e_h->Write();
    ZDC_e_v->Write();
    ZDC_w_h->Write();
    ZDC_w_v->Write();
    pZDCcenter->Write();
    Hist_netChAsym->Write();
    pTPC_EP_east->Write();
    pTPC_EP_west->Write();
    pTPC_EP_for->Write();
    pTPC_EP_bac->Write();
    pTPC_EP_full->Write();
    pBBC_EP_east->Write();
    pBBC_EP_west->Write();
    pEPD_EP_east->Write();
    pEPD_EP_west->Write();
    pZDC_EP_east->Write();
    pZDC_EP_west->Write();
    //  rc->Write();
    fWgtNew->Close();
}
/////////////////////////////////////////////////////
int ReadWeight(char *InFileName)
{
    TFile *fWgt = new TFile(InFileName, "READ");
    if(!fWgt->IsOpen()) return 0;
    if(fWgt->IsOpen())
    {
        TString *histTitle;
        for(int i = 0; i < 10; i++)
        {
            histTitle = new TString("TPCmeanPhi_1_");
            *histTitle += i + 1;
            ReadWeightHist[i].pTPCmeanPhi_1 = (TProfile2D *)fWgt->Get(histTitle->Data());
            delete histTitle;

            histTitle = new TString("TPCmeanPhi_2_");
            *histTitle += i + 1;
            ReadWeightHist[i].pTPCmeanPhi_2 = (TProfile2D *)fWgt->Get(histTitle->Data());
            delete histTitle;

            histTitle = new TString("TPCmeanPhi_3_");
            *histTitle += i + 1;
            ReadWeightHist[i].pTPCmeanPhi_3 = (TProfile2D *)fWgt->Get(histTitle->Data());
            delete histTitle;

            histTitle = new TString("TPCmeanPhiAsso_1_");
            *histTitle += i + 1;
            ReadWeightHist[i].pTPCmeanPhiAsso_1 = (TProfile2D *)fWgt->Get(histTitle->Data());
            delete histTitle;

            histTitle = new TString("TPCmeanPhiAsso_2_");
            *histTitle += i + 1;
            ReadWeightHist[i].pTPCmeanPhiAsso_2 = (TProfile2D *)fWgt->Get(histTitle->Data());
            delete histTitle;

            histTitle = new TString("TPCmeanPhiAsso_3_");
            *histTitle += i + 1;
            ReadWeightHist[i].pTPCmeanPhiAsso_3 = (TProfile2D *)fWgt->Get(histTitle->Data());
            delete histTitle;
        }

        TOF_eff = (TH1D *)fWgt->Get("rc");
        if(TOF_eff && TOF_eff->GetEntries())
        {
            float cont = TOF_eff->GetBinContent(20);
            TOF_eff->Scale(1.25 / cont);
        }
        Read_TPC_EP_full = (TProfile2D *)fWgt->Get("pTPC_EP_full");
        Read_TPC_EP_east = (TProfile2D *)fWgt->Get("pTPC_EP_east");
        Read_TPC_EP_west = (TProfile2D *)fWgt->Get("pTPC_EP_west");
        Read_TPC_EP_for = (TProfile2D *)fWgt->Get("pTPC_EP_for");
        Read_TPC_EP_bac = (TProfile2D *)fWgt->Get("pTPC_EP_bac");
        Read_BBC_EP_east = (TProfile2D *)fWgt->Get("pBBC_EP_east");
        Read_BBC_EP_west = (TProfile2D *)fWgt->Get("pBBC_EP_west");
        Read_EPD_EP_east = (TProfile2D *)fWgt->Get("pEPD_EP_east");
        Read_EPD_EP_west = (TProfile2D *)fWgt->Get("pEPD_EP_west");
        Read_ZDC_EP_east = (TProfile2D *)fWgt->Get("pZDC_EP_east");
        Read_ZDC_EP_west = (TProfile2D *)fWgt->Get("pZDC_EP_west");
        cout << "Loaded: TPC/BBC/EPD/ZDC EP corrections" << endl;
        TH1D *Read_netChAsym   = (TH1D *)fWgt->Get("Hist_netChAsym");
        MeanNetChargeAsym = Read_netChAsym->GetMean();
        RMSNetChargeAsym = Read_netChAsym->GetRMS();
        cout << "Loaded: charge asymmetry " << endl;
        TH1D *BBC_east_1 = (TH1D *)fWgt->Get("BBC1;1");
        TH1D *BBC_east_2 = (TH1D *)fWgt->Get("BBC2;1");
        TH1D *BBC_west_2 = (TH1D *)fWgt->Get("BBC7;1");
        TH1D *BBC_west_1 = (TH1D *)fWgt->Get("BBC8;1");
        float east_mean1 = (BBC_east_1->GetSum()) / 6.;
        float east_mean2 = (BBC_east_2->GetSum()) / 10.;
        float west_mean2 = (BBC_west_2->GetSum()) / 10.;
        float west_mean1 = (BBC_west_1->GetSum()) / 6.;

        float content;
        for(int i = 0; i < 6; i++)
        {
            content = BBC_east_1->GetBinContent(i + 1);
            BBC_gain_east[i] = (content > 0) ? east_mean1 / content : 1;
            content = BBC_west_1->GetBinContent(i + 1);
            BBC_gain_west[i] = (content > 0) ? west_mean1 / content : 1;
        }
        for(int i = 0; i < 10; i++)
        {
            content = BBC_east_2->GetBinContent(i + 1);
            BBC_gain_east[i + 6] = (content > 0) ? 0.2 * east_mean2 / content : 1;
            content = BBC_west_2->GetBinContent(i + 1);
            BBC_gain_west[i + 6] = (content > 0) ? 0.2 * west_mean2 / content : 1;
        }
        cout << "Loaded: BBC gains" << endl;
        float lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
        float cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};
        float par1[9] = {474.651, 474.651, 474.651, 474.651, 474.651, 3.27243e+02, 1.72351, 1.72351, 1.72351};
        float par2[9] = {3.55515, 3.55515, 3.55515, 3.55515, 3.55515, 3.56938, -7.69075, -7.69075, -7.69075};
        float par3[9] = {1.80162, 1.80162, 1.80162, 1.80162, 1.80162, 1.67113, 4.72770, 4.72770, 4.72770};

        for (int iy = 1; iy <= 9; iy++)
        {
            for (int ix = 1; ix < 101; ix++)
            {
                double eta = wt->GetXaxis()->GetBinCenter(ix);
                wt->SetBinContent(ix, iy, lin[iy - 1]*eta + cub[iy - 1]*pow(eta, 3));
                //              wt2->SetBinContent(ix,iy, sqrt(1-1/par1[iy-1]/par1[iy-1]/cosh(eta)/cosh(eta))
                //                              /(1+exp((abs(eta)-par2[iy-1])/par3[iy-1])));
                wt2->SetBinContent(ix, iy, (fabs(eta) > 3.8) ? 1 : 0);
            }
        }
        cout << "Loaded: EPD 1st-order EP weights" << endl;
        TH1D *Read_ZDC_e_h = (TH1D *)fWgt->Get("ZDC_e_h;1");
        TH1D *Read_ZDC_e_v = (TH1D *)fWgt->Get("ZDC_e_v;1");
        TH1D *Read_ZDC_w_h = (TH1D *)fWgt->Get("ZDC_w_h;1");
        TH1D *Read_ZDC_w_v = (TH1D *)fWgt->Get("ZDC_w_v;1");
        east_mean1 = (Read_ZDC_e_v->GetSum()) / 7.;
        east_mean2 = (Read_ZDC_e_h->GetSum()) / 8.;
        west_mean1 = (Read_ZDC_w_v->GetSum()) / 7.;
        west_mean2 = (Read_ZDC_w_h->GetSum()) / 8.;
        for(int i = 0; i < 7; i++)
        {
            content = Read_ZDC_e_v->GetBinContent(i + 1);
            ZDC_gain_east_v[i] = (content > 0) ? east_mean1 / content : 1;
            content = Read_ZDC_w_v->GetBinContent(i + 1);
            ZDC_gain_west_v[i] = (content > 0) ? west_mean1 / content : 1;
        }
        for(int i = 0; i < 8; i++)
        {
            content = Read_ZDC_e_h->GetBinContent(i + 1);
            ZDC_gain_east_h[i] = (content > 0) ? east_mean2 / content : 1;
            content = Read_ZDC_w_h->GetBinContent(i + 1);
            ZDC_gain_west_h[i] = (content > 0) ? west_mean2 / content : 1;
        }
        cout << "Loaded: ZDCSMD gains" << endl;
        Read_ZDCcenter = (TProfile2D *)fWgt->Get("pZDCcenter");
        cout << "Loaded: ZDCSMD beam center" << endl;

    }
    return 1;
}
///////////////////////////////////////////////////
void FillCMW()
{
    p_v2_Ach->Fill(net_Nch_Asym, v2_sub);
    if(net_Nch_Asym <= -99) return;
    if(net_Nch_Asym < (MeanNetChargeAsym - RMSNetChargeAsym))
    {
        if(Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(1, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(1, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(1, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(1, Pt);
        }
    }
    else if(net_Nch_Asym < (MeanNetChargeAsym - (0.3 * RMSNetChargeAsym)))
    {
        if(Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(2, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(2, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(2, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(2, Pt);
        }
    }
    else if(net_Nch_Asym < (MeanNetChargeAsym + (0.3 * RMSNetChargeAsym)))
    {
        if(Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(3, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(3, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(3, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(3, Pt);
        }
    }
    else if(net_Nch_Asym < (MeanNetChargeAsym + RMSNetChargeAsym))
    {
        if(Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(4, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(4, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(4, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(4, Pt);
        }
    }
    else
    {
        if(Charge > 0)
        {
            p_v2_pt_pos_Ach->Fill(5, Pt, v2_sub);
            Hist_pt_pos_Ach->Fill(5, Pt);
        }
        else
        {
            p_v2_pt_neg_Ach->Fill(5, Pt, v2_sub);
            Hist_pt_neg_Ach->Fill(5, Pt);
        }
    }
}
/////////////////////////////////////////////////////
void Fillv2()
{
    Hist_v2_pt_obs1->Fill(Pt, v2);
    Hist_v2_pt_obs2->Fill(Pt, v2, Eweight);
    Hist_v2_eta_obs1->Fill(Eta, v2, 1. / eff);
    Hist_v2_eta_obs2->Fill(Eta, v2, Eweight / eff);
    Hist_v2_pt_EPD_obs->Fill(Pt, 0.5 * (v2_EPDe + v2_EPDw), Eweight);
    Hist_v2_pt_BBC_obs->Fill(Pt, 0.5 * (v2_BBCe + v2_BBCw), Eweight);
    Hist_v2_pt_ZDC_obs->Fill(Pt, v2_ZDC, Eweight);
    Hist_v2_eta_EPD_obs->Fill(Eta, 0.5 * (v2_EPDe + v2_EPDw), Eweight / eff);
    Hist_v2_eta_EPDe_obs->Fill(Eta, v2_EPDe, Eweight / eff);
    Hist_v2_eta_EPDw_obs->Fill(Eta, v2_EPDw, Eweight / eff);
    Hist_v2_eta_BBC_obs->Fill(Eta, 0.5 * (v2_BBCe + v2_BBCw), Eweight / eff);
    Hist_v2_eta_BBCe_obs->Fill(Eta, v2_BBCe, Eweight / eff);
    Hist_v2_eta_BBCw_obs->Fill(Eta, v2_BBCw, Eweight / eff);
    Hist_v2_eta_ZDC_obs->Fill(Eta, v2_ZDC, Eweight / eff);
    if(fabs(Eta) < 0.5)
    {
        pTemp_v2->Fill(1, v2e, 1. / eff);
        pTemp_v2->Fill(2, v2w, 1. / eff);
        if(Charge > 0)
        {
            pTemp_v2->Fill(3, v2e, 1. / eff);
            pTemp_v2->Fill(4, v2w, 1. / eff);
        }
        if(Charge < 0)
        {
            pTemp_v2->Fill(5, v2e, 1. / eff);
            pTemp_v2->Fill(6, v2w, 1. / eff);
        }
    }
    pTemp_v2->Fill(7, v2, 1. / eff);
    pTemp_v2->Fill(8, 0.5 * (v2_EPDe + v2_EPDw), 1. / eff);
    pTemp_v2->Fill(9, 0.5 * (v2_BBCe + v2_BBCw), 1. / eff);
    pTemp_v2->Fill(10, v2_ZDC, 1. / eff);
}
//////////////////////////////////////////////////
void FillRcorr()
{
    if(Charge > 0)
    {
        pTemp_origin_obs->Fill(1, Rsin);
        pTemp_origin_obs->Fill(3, Rcos);
    }
    else
    {
        pTemp_origin_obs->Fill(2, Rsin);
        pTemp_origin_obs->Fill(4, Rcos);
    }
    if(iCh > 0)
    {
        pTemp_shuffle_obs->Fill(1, Rsin);
        pTemp_shuffle_obs->Fill(3, Rcos);
    }
    else
    {
        pTemp_shuffle_obs->Fill(2, Rsin);
        pTemp_shuffle_obs->Fill(4, Rcos);
    }

}
void WrapUpRcorr()
{
    float ds_sin_obs = pTemp_origin_obs->GetBinContent(1) - pTemp_origin_obs->GetBinContent(2);
    float ds_cos_obs = pTemp_origin_obs->GetBinContent(3) - pTemp_origin_obs->GetBinContent(4);
    Hist_dS_sin_obs->Fill(ds_sin_obs, Eweight);
    Hist_dS_cos_obs->Fill(ds_cos_obs, Eweight);
    float ds_sin_shuffle_obs = pTemp_shuffle_obs->GetBinContent(1) - pTemp_shuffle_obs->GetBinContent(2);
    float ds_cos_shuffle_obs = pTemp_shuffle_obs->GetBinContent(3) - pTemp_shuffle_obs->GetBinContent(4);
    Hist_dS_sin_shuffle_obs->Fill(ds_sin_shuffle_obs, Eweight);
    Hist_dS_cos_shuffle_obs->Fill(ds_cos_shuffle_obs, Eweight);

    pTemp_origin_obs->Reset();
    pTemp_shuffle_obs->Reset();
    iCharge.clear();
    match.clear();
}
//////////////////////////////////////////////////
void FillGamma(int ord)
{
    if(ord == 4)
    {
        pParity_int_ss_oppo_run->Fill(Day2, correlator4, Eweight / eff / eff2);
        pDelta_int_ss_oppo_run->Fill(Day2, correlator3, Eweight / eff / eff2);
    }
    if(ord == 3)
    {
        pParity_int_ss_same_run->Fill(Day2, correlator4, Eweight / eff / eff2);
        pDelta_int_ss_same_run->Fill(Day2, correlator3, Eweight / eff / eff2);
    }
    pParity_int_obs1->Fill(ord, correlator0);
    pParity_int_obs3->Fill(ord, correlator0, Eweight / eff / eff2);
    pParity_int_ss_obs1->Fill(ord, correlator4);
    pParity_int_ss_obs3->Fill(ord, correlator4, Eweight / eff / eff2);
    pDelta_int_ss_obs1->Fill(ord, correlator3);
    pDelta_int_ss_obs3->Fill(ord, correlator3, Eweight / eff / eff2);
    pParity_eta_ss_obs1->Fill(ord, 0.5 * (Eta + Eta2), correlator4);
    pParity_eta_ss_obs3->Fill(ord, 0.5 * (Eta + Eta2), correlator4, Eweight / eff / eff2);
    pParity_Deta_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator4);
    pParity_Deta_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator4, Eweight / eff / eff2);
    pParity_pt_ss_obs1->Fill(ord, 0.5 * (Pt + Pt2), correlator4);
    pParity_pt_ss_obs3->Fill(ord, 0.5 * (Pt + Pt2), correlator4, Eweight / eff / eff2);
    pParity_Dpt_ss_obs1->Fill(ord, fabs(Pt - Pt2), correlator4);
    pParity_Dpt_ss_obs3->Fill(ord, fabs(Pt - Pt2), correlator4, Eweight / eff / eff2);
    pDelta_eta_ss_obs1->Fill(ord, 0.5 * (Eta + Eta2), correlator3);
    pDelta_eta_ss_obs3->Fill(ord, 0.5 * (Eta + Eta2), correlator3, Eweight / eff / eff2);
    pDelta_Deta_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator3);
    pDelta_Deta_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator3, Eweight / eff / eff2);
    pDelta_pt_ss_obs1->Fill(ord, 0.5 * (Pt + Pt2), correlator3);
    pDelta_pt_ss_obs3->Fill(ord, 0.5 * (Pt + Pt2), correlator3, Eweight / eff / eff2);
    pDelta_Dpt_ss_obs1->Fill(ord, fabs(Pt - Pt2), correlator3);
    pDelta_Dpt_ss_obs3->Fill(ord, fabs(Pt - Pt2), correlator3, Eweight / eff / eff2);
    if((opt_useEPD == 0 && fabs(Eta) < 0.5 && fabs(Eta2) < 0.5) || opt_useEPD > 0)
    {
        pTemp_parity_e->Fill(ord, correlator4e, 1. / eff / eff2);
        pTemp_parity_w->Fill(ord, correlator4w, 1. / eff / eff2);
        pTemp_delta->Fill(ord, correlator3, 1. / eff / eff2);
        if(fabs(Pt - Pt2) > 0.15 && fabs(Eta - Eta2) > 0.15)
        {
            pTemp_parity_e_noHBT->Fill(ord, correlator4e, 1. / eff / eff2);
            pTemp_parity_w_noHBT->Fill(ord, correlator4w, 1. / eff / eff2);
            pTemp_delta_noHBT->Fill(ord, correlator3, 1. / eff / eff2);
        }
    }
    if(fabs(Pt - Pt2) > 0.15)
    {
        pParity_Deta_highDpt_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator4);
        pParity_Deta_highDpt_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator4, Eweight / eff / eff2);
        pDelta_Deta_highDpt_ss_obs1->Fill(ord, fabs(Eta - Eta2), correlator4);
        pDelta_Deta_highDpt_ss_obs3->Fill(ord, fabs(Eta - Eta2), correlator4, Eweight / eff / eff2);
    }
    if(fabs(Pt - Pt2) > 0.15 && fabs(Eta - Eta2) > 0.15)
    {
        pParity_noHBT_ss_obs1->Fill(ord, correlator4);
        pParity_noHBT_ss_obs3->Fill(ord, correlator4, Eweight / eff / eff2);
        pDelta_noHBT_ss_obs1->Fill(ord, correlator3);
        pDelta_noHBT_ss_obs3->Fill(ord, correlator3, Eweight / eff / eff2);
    }
}
//////////////////////////////////////////////////
void CountMSC(int *iTr, float Phi_ne)
{
    int jud = (iTr[Fcount] > Scount) ? 1 : 0;
    if(jud && Charge > 0)
    {
        if(sin(Phi_ne - TPC_EP_west_new) > 0) N_T_P1++;
        if(sin(Phi_ne - TPC_EP_west_new) < 0) N_B_P1++;
        if(cos(Phi_ne - TPC_EP_west_new) > 0) N_R_P1++;
        if(cos(Phi_ne - TPC_EP_west_new) < 0) N_L_P1++;
    }
    if(jud && Charge < 0)
    {
        if(sin(Phi_ne - TPC_EP_west_new) > 0) N_T_N1++;
        if(sin(Phi_ne - TPC_EP_west_new) < 0) N_B_N1++;
        if(cos(Phi_ne - TPC_EP_west_new) > 0) N_R_N1++;
        if(cos(Phi_ne - TPC_EP_west_new) < 0) N_L_N1++;
    }
    if(!jud && Charge > 0)
    {
        if(sin(Phi_ne - TPC_EP_east_new) > 0) N_T_P2++;
        if(sin(Phi_ne - TPC_EP_east_new) < 0) N_B_P2++;
        if(cos(Phi_ne - TPC_EP_east_new) > 0) N_R_P2++;
        if(cos(Phi_ne - TPC_EP_east_new) < 0) N_L_P2++;
    }
    if(!jud && Charge < 0)
    {
        if(sin(Phi_ne - TPC_EP_east_new) > 0) N_T_N2++;
        if(sin(Phi_ne - TPC_EP_east_new) < 0) N_B_N2++;
        if(cos(Phi_ne - TPC_EP_east_new) > 0) N_R_N2++;
        if(cos(Phi_ne - TPC_EP_east_new) < 0) N_L_N2++;
    }
}
///////////////////////////////////////////////////
void WrapUpMSC()
{
    float PP_in1  = 100.*(N_L_P1 * (N_L_P1 - 1) + N_R_P1 * (N_R_P1 - 1) - 2 * N_L_P1 * N_R_P1) / (N_L_P1 + N_R_P1) / (N_L_P1 + N_R_P1 - 1);
    float PP_out1 = 100.*(N_T_P1 * (N_T_P1 - 1) + N_B_P1 * (N_B_P1 - 1) - 2 * N_T_P1 * N_B_P1) / (N_T_P1 + N_B_P1) / (N_T_P1 + N_B_P1 - 1);
    float NN_in1  = 100.*(N_L_N1 * (N_L_N1 - 1) + N_R_N1 * (N_R_N1 - 1) - 2 * N_L_N1 * N_R_N1) / (N_L_N1 + N_R_N1) / (N_L_N1 + N_R_N1 - 1);
    float NN_out1 = 100.*(N_T_N1 * (N_T_N1 - 1) + N_B_N1 * (N_B_N1 - 1) - 2 * N_T_N1 * N_B_N1) / (N_T_N1 + N_B_N1) / (N_T_N1 + N_B_N1 - 1);
    float PN_in1  = 100.*(N_L_P1 * N_L_N1 + N_R_P1 * N_R_N1 - N_L_P1 * N_R_N1 - N_L_N1 * N_R_P1) / (N_L_P1 + N_R_P1) / (N_L_N1 + N_R_N1);
    float PN_out1 = 100.*(N_T_P1 * N_T_N1 + N_B_P1 * N_B_N1 - N_T_P1 * N_B_N1 - N_T_N1 * N_B_P1) / (N_T_P1 + N_B_P1) / (N_T_N1 + N_B_N1);
    float NP_in1  = 100.*(N_L_N1 * N_L_P1 + N_R_N1 * N_R_P1 - N_L_N1 * N_R_P1 - N_L_P1 * N_R_N1) / (N_L_N1 + N_R_N1) / (N_L_P1 + N_R_P1);
    float NP_out1 = 100.*(N_T_N1 * N_T_P1 + N_B_N1 * N_B_P1 - N_T_N1 * N_B_P1 - N_T_P1 * N_B_N1) / (N_T_N1 + N_B_N1) / (N_T_P1 + N_B_P1);
    float PP_in2  = 100.*(N_L_P2 * (N_L_P2 - 1) + N_R_P2 * (N_R_P2 - 1) - 2 * N_L_P2 * N_R_P2) / (N_L_P2 + N_R_P2) / (N_L_P2 + N_R_P2 - 1);
    float PP_out2 = 100.*(N_T_P2 * (N_T_P2 - 1) + N_B_P2 * (N_B_P2 - 1) - 2 * N_T_P2 * N_B_P2) / (N_T_P2 + N_B_P2) / (N_T_P2 + N_B_P2 - 1);
    float NN_in2  = 100.*(N_L_N2 * (N_L_N2 - 1) + N_R_N2 * (N_R_N2 - 1) - 2 * N_L_N2 * N_R_N2) / (N_L_N2 + N_R_N2) / (N_L_N2 + N_R_N2 - 1);
    float NN_out2 = 100.*(N_T_N2 * (N_T_N2 - 1) + N_B_N2 * (N_B_N2 - 1) - 2 * N_T_N2 * N_B_N2) / (N_T_N2 + N_B_N2) / (N_T_N2 + N_B_N2 - 1);
    float PN_in2  = 100.*(N_L_P2 * N_L_N2 + N_R_P2 * N_R_N2 - N_L_P2 * N_R_N2 - N_L_N2 * N_R_P2) / (N_L_P2 + N_R_P2) / (N_L_N2 + N_R_N2);
    float PN_out2 = 100.*(N_T_P2 * N_T_N2 + N_B_P2 * N_B_N2 - N_T_P2 * N_B_N2 - N_T_N2 * N_B_P2) / (N_T_P2 + N_B_P2) / (N_T_N2 + N_B_N2);
    float NP_in2  = 100.*(N_L_N2 * N_L_P2 + N_R_N2 * N_R_P2 - N_L_N2 * N_R_P2 - N_L_P2 * N_R_N2) / (N_L_N2 + N_R_N2) / (N_L_P2 + N_R_P2);
    float NP_out2 = 100.*(N_T_N2 * N_T_P2 + N_B_N2 * N_B_P2 - N_T_N2 * N_B_P2 - N_T_P2 * N_B_N2) / (N_T_N2 + N_B_N2) / (N_T_P2 + N_B_P2);

    int Nin1  = (N_L_P1 - N_L_N1) - (N_R_P1 - N_R_N1);
    int Nout1 = (N_T_P1 - N_T_N1) - (N_B_P1 - N_B_N1);
    int Nin2  = (N_L_P2 - N_L_N2) - (N_R_P2 - N_R_N2);
    int Nout2 = (N_T_P2 - N_T_N2) - (N_B_P2 - N_B_N2);
    Hist_dQ_in1->Fill(Nin1, Eweight);
    Hist_dQ_out1->Fill(Nout1, Eweight);
    Hist_dQ_in2->Fill(Nin2, Eweight);
    Hist_dQ_out2->Fill(Nout2, Eweight);

    if(!isnan(PP_in1)) pParity_int_MSC_same_in1->Fill(Nin1, MM * PP_in1, Eweight);
    if(!isnan(NN_in1)) pParity_int_MSC_same_in1->Fill(Nin1, MM * NN_in1, Eweight);
    if(!isnan(PN_in1)) pParity_int_MSC_oppo_in1->Fill(Nin1, MM * PN_in1, Eweight);
    if(!isnan(NP_in1)) pParity_int_MSC_oppo_in1->Fill(Nin1, MM * NP_in1, Eweight);
    if(!isnan(PP_out1)) pParity_int_MSC_same_out1->Fill(Nout1, MM * PP_out1, Eweight);
    if(!isnan(NN_out1)) pParity_int_MSC_same_out1->Fill(Nout1, MM * NN_out1, Eweight);
    if(!isnan(PN_out1)) pParity_int_MSC_oppo_out1->Fill(Nout1, MM * PN_out1, Eweight);
    if(!isnan(NP_out1)) pParity_int_MSC_oppo_out1->Fill(Nout1, MM * NP_out1, Eweight);
    if(!isnan(PP_in2)) pParity_int_MSC_same_in2->Fill(Nin2, MM * PP_in2, Eweight);
    if(!isnan(NN_in2)) pParity_int_MSC_same_in2->Fill(Nin2, MM * NN_in2, Eweight);
    if(!isnan(PN_in2)) pParity_int_MSC_oppo_in2->Fill(Nin2, MM * PN_in2, Eweight);
    if(!isnan(NP_in2)) pParity_int_MSC_oppo_in2->Fill(Nin2, MM * NP_in2, Eweight);
    if(!isnan(PP_out2)) pParity_int_MSC_same_out2->Fill(Nout2, MM * PP_out2, Eweight);
    if(!isnan(NN_out2)) pParity_int_MSC_same_out2->Fill(Nout2, MM * NN_out2, Eweight);
    if(!isnan(PN_out2)) pParity_int_MSC_oppo_out2->Fill(Nout2, MM * PN_out2, Eweight);
    if(!isnan(NP_out2)) pParity_int_MSC_oppo_out2->Fill(Nout2, MM * NP_out2, Eweight);

    if(!isnan(PP_in1) && !isnan(PP_out1)) pParity_int_MSC1_obs1->Fill(1, MM * (PP_in1 - PP_out1));
    if(!isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs1->Fill(2, MM * (NN_in1 - NN_out1));
    if(!isnan(PP_in1) && !isnan(PP_out1) && !isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs1->Fill(3, MM * (PP_in1 - PP_out1));
    if(!isnan(PP_in1) && !isnan(PP_out1) && !isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs1->Fill(3, MM * (NN_in1 - NN_out1));
    if(!isnan(PN_in1) && !isnan(PN_out1) && !isnan(NP_in1) && !isnan(NP_out1)) pParity_int_MSC1_obs1->Fill(4, MM * (PN_in1 - PN_out1));
    if(!isnan(PN_in1) && !isnan(PN_out1) && !isnan(NP_in1) && !isnan(NP_out1)) pParity_int_MSC1_obs1->Fill(4, MM * (NP_in1 - NP_out1));
    if(!isnan(PP_in1) && !isnan(PP_out1)) pParity_int_MSC1_obs2->Fill(1, MM * (PP_in1 - PP_out1), Eweight);
    if(!isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs2->Fill(2, MM * (NN_in1 - NN_out1), Eweight);
    if(!isnan(PP_in1) && !isnan(PP_out1) && !isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs2->Fill(3, MM * (PP_in1 - PP_out1), Eweight);
    if(!isnan(PP_in1) && !isnan(PP_out1) && !isnan(NN_in1) && !isnan(NN_out1)) pParity_int_MSC1_obs2->Fill(3, MM * (NN_in1 - NN_out1), Eweight);
    if(!isnan(PN_in1) && !isnan(PN_out1) && !isnan(NP_in1) && !isnan(NP_out1)) pParity_int_MSC1_obs2->Fill(4, MM * (PN_in1 - PN_out1), Eweight);
    if(!isnan(PN_in1) && !isnan(PN_out1) && !isnan(NP_in1) && !isnan(NP_out1)) pParity_int_MSC1_obs2->Fill(4, MM * (NP_in1 - NP_out1), Eweight);
    if(!isnan(PP_in2) && !isnan(PP_out2)) pParity_int_MSC2_obs1->Fill(1, MM * (PP_in2 - PP_out2));
    if(!isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs1->Fill(2, MM * (NN_in2 - NN_out2));
    if(!isnan(PP_in2) && !isnan(PP_out2) && !isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs1->Fill(3, MM * (PP_in2 - PP_out2));
    if(!isnan(PP_in2) && !isnan(PP_out2) && !isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs1->Fill(3, MM * (NN_in2 - NN_out2));
    if(!isnan(PN_in2) && !isnan(PN_out2) && !isnan(NP_in2) && !isnan(NP_out2)) pParity_int_MSC2_obs1->Fill(4, MM * (PN_in2 - PN_out2));
    if(!isnan(PN_in2) && !isnan(PN_out2) && !isnan(NP_in2) && !isnan(NP_out2)) pParity_int_MSC2_obs1->Fill(4, MM * (NP_in2 - NP_out2));
    if(!isnan(PP_in2) && !isnan(PP_out2)) pParity_int_MSC2_obs2->Fill(1, MM * (PP_in2 - PP_out2), Eweight);
    if(!isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs2->Fill(2, MM * (NN_in2 - NN_out2), Eweight);
    if(!isnan(PP_in2) && !isnan(PP_out2) && !isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs2->Fill(3, MM * (PP_in2 - PP_out2), Eweight);
    if(!isnan(PP_in2) && !isnan(PP_out2) && !isnan(NN_in2) && !isnan(NN_out2)) pParity_int_MSC2_obs2->Fill(3, MM * (NN_in2 - NN_out2), Eweight);
    if(!isnan(PN_in2) && !isnan(PN_out2) && !isnan(NP_in2) && !isnan(NP_out2)) pParity_int_MSC2_obs2->Fill(4, MM * (PN_in2 - PN_out2), Eweight);
    if(!isnan(PN_in2) && !isnan(PN_out2) && !isnan(NP_in2) && !isnan(NP_out2)) pParity_int_MSC2_obs2->Fill(4, MM * (NP_in2 - NP_out2), Eweight);
}
//////////////////////////////////////////////////////
void WrapUpESE()
{
    float Temp_v2e = pTemp_v2->GetBinContent(1);
    float Temp_v2w = pTemp_v2->GetBinContent(2);
    float Temp_v2e_P = pTemp_v2->GetBinContent(3);
    float Temp_v2w_P = pTemp_v2->GetBinContent(4);
    float Temp_v2e_N = pTemp_v2->GetBinContent(5);
    float Temp_v2w_N = pTemp_v2->GetBinContent(6);
    float Temp_parity1e = pTemp_parity_e->GetBinContent(1);
    float Temp_parity2e = pTemp_parity_e->GetBinContent(2);
    float Temp_parity3e = pTemp_parity_e->GetBinContent(3);
    float Temp_parity4e = pTemp_parity_e->GetBinContent(4);
    float Temp_parity1w = pTemp_parity_w->GetBinContent(1);
    float Temp_parity2w = pTemp_parity_w->GetBinContent(2);
    float Temp_parity3w = pTemp_parity_w->GetBinContent(3);
    float Temp_parity4w = pTemp_parity_w->GetBinContent(4);
    float Temp_delta1 = pTemp_delta->GetBinContent(1);
    float Temp_delta2 = pTemp_delta->GetBinContent(2);
    float Temp_delta3 = pTemp_delta->GetBinContent(3);
    float Temp_delta4 = pTemp_delta->GetBinContent(4);

    if(Temp_v2e != 0)
    {
        p_v2e_Q2_obs1->Fill(Q2, Temp_v2e);
        p_v2e_Q2_obs2->Fill(Q2, Temp_v2e, Eweight);
    }
    if(Temp_v2w != 0)
    {
        p_v2w_Q2_obs1->Fill(Q2, Temp_v2w);
        p_v2w_Q2_obs2->Fill(Q2, Temp_v2w, Eweight);
    }
    if(Temp_v2e_P != 0)
    {
        p_v2e_Q2_P_obs1->Fill(Q2, Temp_v2e_P);
        p_v2e_Q2_P_obs2->Fill(Q2, Temp_v2e_P, Eweight);
    }
    if(Temp_v2w_P != 0)
    {
        p_v2w_Q2_P_obs1->Fill(Q2, Temp_v2w_P);
        p_v2w_Q2_P_obs2->Fill(Q2, Temp_v2w_P, Eweight);
    }
    if(Temp_v2e_N != 0)
    {
        p_v2e_Q2_N_obs1->Fill(Q2, Temp_v2e_N);
        p_v2e_Q2_N_obs2->Fill(Q2, Temp_v2e_N, Eweight);
    }
    if(Temp_v2w_N != 0)
    {
        p_v2w_Q2_N_obs1->Fill(Q2, Temp_v2w_N);
        p_v2w_Q2_N_obs2->Fill(Q2, Temp_v2w_N, Eweight);
    }
    if(Temp_parity1e != 0)
    {
        pParity_e_Q2_obs1->Fill(1, Q2, Temp_parity1e);
        pParity_e_Q2_obs2->Fill(1, Q2, Temp_parity1e, Eweight);
    }
    if(Temp_parity2e != 0)
    {
        pParity_e_Q2_obs1->Fill(2, Q2, Temp_parity2e);
        pParity_e_Q2_obs2->Fill(2, Q2, Temp_parity2e, Eweight);
    }
    if(Temp_parity3e != 0)
    {
        pParity_e_Q2_obs1->Fill(3, Q2, Temp_parity3e);
        pParity_e_Q2_obs2->Fill(3, Q2, Temp_parity3e, Eweight);
    }
    if(Temp_parity4e != 0)
    {
        pParity_e_Q2_obs1->Fill(4, Q2, Temp_parity4e);
        pParity_e_Q2_obs2->Fill(4, Q2, Temp_parity4e, Eweight);
    }
    if(Temp_parity1w != 0)
    {
        pParity_w_Q2_obs1->Fill(1, Q2, Temp_parity1w);
        pParity_w_Q2_obs2->Fill(1, Q2, Temp_parity1w, Eweight);
    }
    if(Temp_parity2w != 0)
    {
        pParity_w_Q2_obs1->Fill(2, Q2, Temp_parity2w);
        pParity_w_Q2_obs2->Fill(2, Q2, Temp_parity2w, Eweight);
    }
    if(Temp_parity3w != 0)
    {
        pParity_w_Q2_obs1->Fill(3, Q2, Temp_parity3w);
        pParity_w_Q2_obs2->Fill(3, Q2, Temp_parity3w, Eweight);
    }
    if(Temp_parity4w != 0)
    {
        pParity_w_Q2_obs1->Fill(4, Q2, Temp_parity4w);
        pParity_w_Q2_obs2->Fill(4, Q2, Temp_parity4w, Eweight);
    }
    if(Temp_delta1 != 0)
    {
        pDelta_Q2_obs1->Fill(1, Q2, Temp_delta1);
        pDelta_Q2_obs2->Fill(1, Q2, Temp_delta1, Eweight);
    }
    if(Temp_delta2 != 0)
    {
        pDelta_Q2_obs1->Fill(2, Q2, Temp_delta2);
        pDelta_Q2_obs2->Fill(2, Q2, Temp_delta2, Eweight);
    }
    if(Temp_delta3 != 0)
    {
        pDelta_Q2_obs1->Fill(3, Q2, Temp_delta3);
        pDelta_Q2_obs2->Fill(3, Q2, Temp_delta3, Eweight);
    }
    if(Temp_delta4 != 0)
    {
        pDelta_Q2_obs1->Fill(4, Q2, Temp_delta4);
        pDelta_Q2_obs2->Fill(4, Q2, Temp_delta4, Eweight);
    }

    Temp_v2e = pTemp_v2_noHBT->GetBinContent(1);
    Temp_v2w = pTemp_v2_noHBT->GetBinContent(2);
    Temp_v2e_P = pTemp_v2_noHBT->GetBinContent(3);
    Temp_v2w_P = pTemp_v2_noHBT->GetBinContent(4);
    Temp_v2e_N = pTemp_v2_noHBT->GetBinContent(5);
    Temp_v2w_N = pTemp_v2_noHBT->GetBinContent(6);
    Temp_parity1e = pTemp_parity_e_noHBT->GetBinContent(1);
    Temp_parity2e = pTemp_parity_e_noHBT->GetBinContent(2);
    Temp_parity3e = pTemp_parity_e_noHBT->GetBinContent(3);
    Temp_parity4e = pTemp_parity_e_noHBT->GetBinContent(4);
    Temp_parity1w = pTemp_parity_w_noHBT->GetBinContent(1);
    Temp_parity2w = pTemp_parity_w_noHBT->GetBinContent(2);
    Temp_parity3w = pTemp_parity_w_noHBT->GetBinContent(3);
    Temp_parity4w = pTemp_parity_w_noHBT->GetBinContent(4);
    Temp_delta1 = pTemp_delta_noHBT->GetBinContent(1);
    Temp_delta2 = pTemp_delta_noHBT->GetBinContent(2);
    Temp_delta3 = pTemp_delta_noHBT->GetBinContent(3);
    Temp_delta4 = pTemp_delta_noHBT->GetBinContent(4);

    if(Temp_v2e != 0)
    {
        p_v2e_noHBT_Q2_obs1->Fill(Q2, Temp_v2e);
        p_v2e_noHBT_Q2_obs2->Fill(Q2, Temp_v2e, Eweight);
    }
    if(Temp_v2w != 0)
    {
        p_v2w_noHBT_Q2_obs1->Fill(Q2, Temp_v2w);
        p_v2w_noHBT_Q2_obs2->Fill(Q2, Temp_v2w, Eweight);
    }
    if(Temp_v2e_P != 0)
    {
        p_v2e_noHBT_Q2_P_obs1->Fill(Q2, Temp_v2e_P);
        p_v2e_noHBT_Q2_P_obs2->Fill(Q2, Temp_v2e_P, Eweight);
    }
    if(Temp_v2w_P != 0)
    {
        p_v2w_noHBT_Q2_P_obs1->Fill(Q2, Temp_v2w_P);
        p_v2w_noHBT_Q2_P_obs2->Fill(Q2, Temp_v2w_P, Eweight);
    }
    if(Temp_v2e_N != 0)
    {
        p_v2e_noHBT_Q2_N_obs1->Fill(Q2, Temp_v2e_N);
        p_v2e_noHBT_Q2_N_obs2->Fill(Q2, Temp_v2e_N, Eweight);
    }
    if(Temp_v2w_N != 0)
    {
        p_v2w_noHBT_Q2_N_obs1->Fill(Q2, Temp_v2w_N);
        p_v2w_noHBT_Q2_N_obs2->Fill(Q2, Temp_v2w_N, Eweight);
    }
    if(Temp_parity1e != 0)
    {
        pParity_e_noHBT_Q2_obs1->Fill(1, Q2, Temp_parity1e);
        pParity_e_noHBT_Q2_obs2->Fill(1, Q2, Temp_parity1e, Eweight);
    }
    if(Temp_parity2e != 0)
    {
        pParity_e_noHBT_Q2_obs1->Fill(2, Q2, Temp_parity2e);
        pParity_e_noHBT_Q2_obs2->Fill(2, Q2, Temp_parity2e, Eweight);
    }
    if(Temp_parity3e != 0)
    {
        pParity_e_noHBT_Q2_obs1->Fill(3, Q2, Temp_parity3e);
        pParity_e_noHBT_Q2_obs2->Fill(3, Q2, Temp_parity3e, Eweight);
    }
    if(Temp_parity4e != 0)
    {
        pParity_e_noHBT_Q2_obs1->Fill(4, Q2, Temp_parity4e);
        pParity_e_noHBT_Q2_obs2->Fill(4, Q2, Temp_parity4e, Eweight);
    }
    if(Temp_parity1w != 0)
    {
        pParity_w_noHBT_Q2_obs1->Fill(1, Q2, Temp_parity1w);
        pParity_w_noHBT_Q2_obs2->Fill(1, Q2, Temp_parity1w, Eweight);
    }
    if(Temp_parity2w != 0)
    {
        pParity_w_noHBT_Q2_obs1->Fill(2, Q2, Temp_parity2w);
        pParity_w_noHBT_Q2_obs2->Fill(2, Q2, Temp_parity2w, Eweight);
    }
    if(Temp_parity3w != 0)
    {
        pParity_w_noHBT_Q2_obs1->Fill(3, Q2, Temp_parity3w);
        pParity_w_noHBT_Q2_obs2->Fill(3, Q2, Temp_parity3w, Eweight);
    }
    if(Temp_parity4w != 0)
    {
        pParity_w_noHBT_Q2_obs1->Fill(4, Q2, Temp_parity4w);
        pParity_w_noHBT_Q2_obs2->Fill(4, Q2, Temp_parity4w, Eweight);
    }
    if(Temp_delta1 != 0)
    {
        pDelta_noHBT_Q2_obs1->Fill(1, Q2, Temp_delta1);
        pDelta_noHBT_Q2_obs2->Fill(1, Q2, Temp_delta1, Eweight);
    }
    if(Temp_delta2 != 0)
    {
        pDelta_noHBT_Q2_obs1->Fill(2, Q2, Temp_delta2);
        pDelta_noHBT_Q2_obs2->Fill(2, Q2, Temp_delta2, Eweight);
    }
    if(Temp_delta3 != 0)
    {
        pDelta_noHBT_Q2_obs1->Fill(3, Q2, Temp_delta3);
        pDelta_noHBT_Q2_obs2->Fill(3, Q2, Temp_delta3, Eweight);
    }
    if(Temp_delta4 != 0)
    {
        pDelta_noHBT_Q2_obs1->Fill(4, Q2, Temp_delta4);
        pDelta_noHBT_Q2_obs2->Fill(4, Q2, Temp_delta4, Eweight);
    }

    pTemp_v2->Reset();
    pTemp_v2_noHBT->Reset();
    pTemp_parity_e->Reset();
    pTemp_parity_e_noHBT->Reset();
    pTemp_parity_w->Reset();
    pTemp_parity_w_noHBT->Reset();
    pTemp_delta->Reset();
    pTemp_delta_noHBT->Reset();
}
/////////////////////////////
float GetPhiInBBC(int e_w, int bbcN)   //bbcN=1 to 24
{
    const float phi_div = PI / 6;
    float bbc_phi = phi_div;
    switch(bbcN)
    {
    case 1:
        bbc_phi = 3 * phi_div;
        break;
    case 2:
        bbc_phi = phi_div;
        break;
    case 3:
        bbc_phi = -1 * phi_div;
        break;
    case 4:
        bbc_phi = -3 * phi_div;
        break;
    case 5:
        bbc_phi = -5 * phi_div;
        break;
    case 6:
        bbc_phi = 5 * phi_div;
        break;
    case 7:
        bbc_phi = (gRandom->Rndm() > 0.5) ? 2 * phi_div : 4 * phi_div;
        break;
    case 8:
        bbc_phi = 3 * phi_div;
        break;
    case 9:
        bbc_phi = phi_div;
        break;
    case 10:
        bbc_phi = 0.;
        break;
    case 11:
        bbc_phi = -phi_div;
        break;
    case 12:
        bbc_phi = (gRandom->Rndm() > 0.5) ? -2 * phi_div : -4 * phi_div;
        break;
    case 13:
        bbc_phi = -3 * phi_div;
        break;
    case 14:
        bbc_phi = -5 * phi_div;
        break;
    case 15:
        bbc_phi = PI;
        break;
    case 16:
        bbc_phi = 5 * phi_div;
        break;
    case 17:
        bbc_phi = 3 * phi_div;
        break;
    case 18:
        bbc_phi = 0.;
        break;
    case 19:
        bbc_phi = -3 * phi_div;
        break;
    case 20:
        bbc_phi = PI;
        break;
    case 21:
        bbc_phi = 3 * phi_div;
        break;
    case 22:
        bbc_phi = 0.;
        break;
    case 23:
        bbc_phi = -3 * phi_div;
        break;
    case 24:
        bbc_phi = PI;
        break;
    }
    if(e_w == 0)
    {
        if(bbc_phi > -0.001)
        {
            bbc_phi = PI - bbc_phi;
        }
        else
        {
            bbc_phi = -PI - bbc_phi;
        }
    }
    return bbc_phi;
}
///////////////////////////////////////////
float GetXYInZDC(int e_w, int v_h, int zdcN, int opt_raw )   //zdcN=1 to 7 or 8
{
    //get position of each slat;strip starts from 1
    float zdcsmd_x[7] = {0.5, 2, 3.5, 5, 6.5, 8, 9.5};
    float zdcsmd_y[8] = {1.25, 3.25, 5.25, 7.25, 9.25, 11.25, 13.25, 15.25};
    float zdcsmd_xx = zdcsmd_x0_e, zdcsmd_yy = zdcsmd_y0_e;
    if(1 == e_w)
    {
        zdcsmd_xx = zdcsmd_x0_w;
        zdcsmd_yy = zdcsmd_y0_w;
    }
    if(1 == opt_raw)
    {
        if(0 == v_h) return zdcsmd_x[zdcN - 1];
        else if(1 == v_h) return zdcsmd_y[zdcN - 1];
    }
    else
    {
        if(0 == v_h) return zdcsmd_x[zdcN - 1] - zdcsmd_xx;
        else if(1 == v_h) return zdcsmd_y[zdcN - 1] - zdcsmd_yy;
    }
    return 0;
}
/////////////////////////////////////////////
void DefineHistogram()
{
    TString *histTitle;
    for(int i = 0; i < 10; i++)
    {
        histTitle = new TString("TPCmeanPhi_1_");
        *histTitle += i + 1;
        WeightHist[i].pTPCmeanPhi_1 = new TProfile2D(histTitle->Data(), histTitle->Data(), 8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
        delete histTitle;

        histTitle = new TString("TPCmeanPhi_2_");
        *histTitle += i + 1;
        WeightHist[i].pTPCmeanPhi_2 = new TProfile2D(histTitle->Data(), histTitle->Data(), 8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
        delete histTitle;

        histTitle = new TString("TPCmeanPhi_3_");
        *histTitle += i + 1;
        WeightHist[i].pTPCmeanPhi_3 = new TProfile2D(histTitle->Data(), histTitle->Data(), 8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
        delete histTitle;

        histTitle = new TString("TPCmeanPhiAsso_1_");
        *histTitle += i + 1;
        WeightHist[i].pTPCmeanPhiAsso_1 = new TProfile2D(histTitle->Data(), histTitle->Data(), 8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
        delete histTitle;

        histTitle = new TString("TPCmeanPhiAsso_2_");
        *histTitle += i + 1;
        WeightHist[i].pTPCmeanPhiAsso_2 = new TProfile2D(histTitle->Data(), histTitle->Data(), 8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
        delete histTitle;

        histTitle = new TString("TPCmeanPhiAsso_3_");
        *histTitle += i + 1;
        WeightHist[i].pTPCmeanPhiAsso_3 = new TProfile2D(histTitle->Data(), histTitle->Data(), 8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
        delete histTitle;

        histTitle = new TString("Hist_Phi_1_");
        *histTitle += i + 1;
        PhiHist[i].Hist_Phi_1 = new TH2D(histTitle->Data(), histTitle->Data(), Phibin, -PI, PI, 4, 0.5, 4.5);
        delete histTitle;

        histTitle = new TString("Hist_Phi_2_");
        *histTitle += i + 1;
        PhiHist[i].Hist_Phi_2 = new TH2D(histTitle->Data(), histTitle->Data(), Phibin, -PI, PI, 4, 0.5, 4.5);
        delete histTitle;

        histTitle = new TString("Hist_Phi_3_");
        *histTitle += i + 1;
        PhiHist[i].Hist_Phi_3 = new TH2D(histTitle->Data(), histTitle->Data(), Phibin, -PI, PI, 4, 0.5, 4.5);
        delete histTitle;

        histTitle = new TString("Hist_Phi_new_1_");
        *histTitle += i + 1;
        PhiHist[i].Hist_Phi_new_1 = new TH2D(histTitle->Data(), histTitle->Data(), Phibin, -PI, PI, 4, 0.5, 4.5);
        delete histTitle;

        histTitle = new TString("Hist_Phi_new_2_");
        *histTitle += i + 1;
        PhiHist[i].Hist_Phi_new_2 = new TH2D(histTitle->Data(), histTitle->Data(), Phibin, -PI, PI, 4, 0.5, 4.5);
        delete histTitle;

        histTitle = new TString("Hist_Phi_new_3_");
        *histTitle += i + 1;
        PhiHist[i].Hist_Phi_new_3 = new TH2D(histTitle->Data(), histTitle->Data(), Phibin, -PI, PI, 4, 0.5, 4.5);
        delete histTitle;
    }

}
