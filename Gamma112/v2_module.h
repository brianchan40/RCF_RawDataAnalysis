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

const float PI = TMath::Pi();
const float MM = 2 / PI;
const int opt_useBBC = 1;
const int Phibin = 80;
const int order = 10;
const int nHar = 2; //2nd order EP
const int run_sta = 71000;
const int run_end = 168000;
const float pt_trig_up = 2;
const float pt_trig_lo = .2;
const float pt_asso_up = 2.0;//1;
const float pt_asso_lo = 0.15;//0.15;
const float EtaCut = 1.3;
const float DcaCut = 2;
const float Eta_EP_Cut = 0;  //particles in the event plane
const Float_t Vz_cut = 40;    //40 for 39 Gev, 75 for 11 GeV and 70 for 7GeV
//const int cenDef[9] = {10, 22, 43, 76, 125, 193, 281, 396, 466}; // 200 GeV run11
const int cenDef[9] = {19, 31, 46, 67, 93, 128, 172, 230, 267}; //isobar

//bad runs
const int Nrun_MB = 22;
const int bad_day3_MB[Nrun_MB] = {170107, 171012, 171017, 171039, 174023, 174024, 174025, 174026, 174027, 174028, 176030, 178052, 178053, 178055, //TOF
                                  171081, 172136, 179021, 184011, 184012, 184013, 185004, 185008 //Cos
                                 };

//new efficiency
const float PP0[9] = {8.63097e-01, 8.69215e-01, 8.74097e-01, 8.69024e-01, 8.55400e-01, 8.28717e-01, 7.85676e-01, 7.37110e-01, 6.93514e-01};
const float PP1[9] = {1.30629e-01, 1.36947e-01, 1.28525e-01, 1.30408e-01, 1.32583e-01, 1.34462e-01, 1.38044e-01, 1.41042e-01, 1.44047e-01};
const float PP2[9] = {8.17412e+00, 5.23658e+00, 6.77459e+00, 7.23513e+00, 6.23855e+00, 6.88205e+00, 6.37277e+00, 7.25640e+00, 6.58271e+00};

//defining histograms
TH2D *Ref_TOF = new TH2D("Ref_TOF", "Ref_TOF", 500, 0.5, 500.5, 5000, 0.5, 5000.5);
TProfile *Ref_Day3 = new TProfile("Ref_Day3", "RefMult vs Run", run_end - run_sta, run_sta, run_end, 0, 999);
TProfile *TOF_Day3 = new TProfile("TOF_Day3", "TOFMult vs Run", run_end - run_sta, run_sta, run_end, 0, 5000);
TProfile *NPT_Day3 = new TProfile("NPT_Day3", "NPTracks vs Run", run_end - run_sta, run_sta, run_end, 0, 5000);
TProfile *NPA_Day3 = new TProfile("NPA_Day3", "NPAsso vs Run", run_end - run_sta, run_sta, run_end, 0, 5000);
TProfile *TPC_Day3_cos2 = new TProfile("TPC_Day3_cos2", "cos(2*psi) vs Run", run_end - run_sta, run_sta, run_end, -1, 1);
TProfile *TPC_Day3_sin2 = new TProfile("TPC_Day3_sin2", "sin(2*psi) vs Run", run_end - run_sta, run_sta, run_end, -1, 1);

TH1D *hTally = new TH1D("hTally", "hTally", 10, 0.5, 10.5);
TH1D *hTall  = new TH1D("hTall ", "hTall ", 11, 0.5, 11.5);
TH1D *hZDCcoin = new TH1D("hZDCcoin", "hZDCcoin", 1000, 0, 100000);
TH1D *hTrigger = new TH1D("hTrigger", "hTrigger", 200, 0.5, 200.5);
TH1D *hCentrality = new TH1D("hCentrality", "hCentrality", 11, -1.5, 9.5);
TH2D *hVertexXY = new TH2D("hVertexXY", "hVertexXY", 300, -3, 3, 300, -3, 3);
TH1D *hVertexZ = new TH1D("hVertexZ", "hVertexZ", 100, -100, 100);
TH1D *hVzDiff = new TH1D("hVzDiff", "hVzDiff", 120, -30, 30);
TH2D *hMult_Vz = new TH2D("hMult_Vz", "hMult_Vz", 1000, -0.5, 999.5, 100, -100, 100);
TH2D *hMult_Vz_new = new TH2D("hMult_Vz_new", "hMult_Vz_new", 1000, -0.5, 999.5, 100, -100, 100);

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

TProfile2D *pTPCmeanPhi_FF_1 = new TProfile2D("TPCmeanPhi_FF_1", "TPCmeanPhi_FF_1",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_1_p = new TProfile2D("TPCmeanPhi_FF_1_p", "TPCmeanPhi_FF_1_p",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_1 = new TProfile2D("TPCmeanPhi_RF_1", "TPCmeanPhi_RF_1",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_1_p = new TProfile2D("TPCmeanPhi_RF_1_p", "TPCmeanPhi_RF_1_p",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_FF_1 = new TProfile2D("TPCmeanPhiAsso_FF_1", "TPCmeanPhiAsso_FF_1",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_RF_1 = new TProfile2D("TPCmeanPhiAsso_RF_1", "TPCmeanPhiAsso_RF_1",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_2 = new TProfile2D("TPCmeanPhi_FF_2", "TPCmeanPhi_FF_2",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_2_p = new TProfile2D("TPCmeanPhi_FF_2_p", "TPCmeanPhi_FF_2_p",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_2 = new TProfile2D("TPCmeanPhi_RF_2", "TPCmeanPhi_RF_2",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_2_p = new TProfile2D("TPCmeanPhi_RF_2_p", "TPCmeanPhi_RF_2_p",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_FF_2 = new TProfile2D("TPCmeanPhiAsso_FF_2", "TPCmeanPhiAsso_FF_2",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_RF_2 = new TProfile2D("TPCmeanPhiAsso_RF_2", "TPCmeanPhiAsso_RF_2",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_3 = new TProfile2D("TPCmeanPhi_FF_3", "TPCmeanPhi_FF_3",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_FF_3_p = new TProfile2D("TPCmeanPhi_FF_3_p", "TPCmeanPhi_FF_3_p",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_3 = new TProfile2D("TPCmeanPhi_RF_3", "TPCmeanPhi_RF_3",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhi_RF_3_p = new TProfile2D("TPCmeanPhi_RF_3_p", "TPCmeanPhi_RF_3_p",
                                              8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_FF_3 = new TProfile2D("TPCmeanPhiAsso_FF_3", "TPCmeanPhiAsso_FF_3",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPCmeanPhiAsso_RF_3 = new TProfile2D("TPCmeanPhiAsso_RF_3", "TPCmeanPhiAsso_RF_3",
        8 * order, 0.5, 8 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TH2D *Hist_TPC_EP_east = new TH2D("Hist_TPC_EP_east", "Hist_TPC_EP_east", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_west = new TH2D("Hist_TPC_EP_west", "Hist_TPC_EP_west", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_full = new TH2D("Hist_TPC_EP_full", "Hist_TPC_EP_full", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_for = new TH2D("Hist_TPC_EP_for", "Hist_TPC_EP_for", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_bac = new TH2D("Hist_TPC_EP_bac", "Hist_TPC_EP_bac", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_east_flat = new TH2D("Hist_TPC_EP_east_flat", "Hist_TPC_EP_east_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_west_flat = new TH2D("Hist_TPC_EP_west_flat", "Hist_TPC_EP_west_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_full_flat = new TH2D("Hist_TPC_EP_full_flat", "Hist_TPC_EP_full_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_for_flat = new TH2D("Hist_TPC_EP_for_flat", "Hist_TPC_EP_for_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH2D *Hist_TPC_EP_bac_flat = new TH2D("Hist_TPC_EP_bac_flat", "Hist_TPC_EP_bac_flat", 36, 0, PI, (run_end - run_sta) / 1000, run_sta / 1000, run_end / 1000);
TH1D *Hist_TPC_EP_full_m1 = new TH1D("Hist_TPC_EP_full_m1", "Hist_TPC_EP_full_m1", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m2 = new TH1D("Hist_TPC_EP_full_m2", "Hist_TPC_EP_full_m2", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m1_flat = new TH1D("Hist_TPC_EP_full_m1_flat", "Hist_TPC_EP_full_m1_flat", 36, 0, PI);
TH1D *Hist_TPC_EP_full_m2_flat = new TH1D("Hist_TPC_EP_full_m2_flat", "Hist_TPC_EP_full_m2_flat", 36, 0, PI);
TProfile2D *pTPC_EP_east = new TProfile2D("pTPC_EP_east", "pTPC_EP_east", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_west = new TProfile2D("pTPC_EP_west", "pTPC_EP_west", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_full = new TProfile2D("pTPC_EP_full", "pTPC_EP_full", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_for = new TProfile2D("pTPC_EP_for", "pTPC_EP_for", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
TProfile2D *pTPC_EP_bac = new TProfile2D("pTPC_EP_bac", "pTPC_EP_bac", 2 * order, 0.5, 2 * order + 0.5, (run_end - run_sta) / 10, run_sta / 10, run_end / 10, -1, 1, "");
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
TH2D *Hist_Phi = new TH2D("Hist_Phi", "Hist_Phi", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_1 = new TH2D("Hist_Phi_FF_1", "Hist_Phi_FF_1", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_1_p = new TH2D("Hist_Phi_FF_1_p", "Hist_Phi_FF_1_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_1 = new TH2D("Hist_Phi_RF_1", "Hist_Phi_RF_1", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_1_p = new TH2D("Hist_Phi_RF_1_p", "Hist_Phi_RF_1_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_1 = new TH2D("Hist_Phi_FF_new_1", "Hist_Phi_FF_new_1", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_1 = new TH2D("Hist_Phi_RF_new_1", "Hist_Phi_RF_new_1", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_1_p = new TH2D("Hist_Phi_FF_new_1_p", "Hist_Phi_FF_new_1_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_1_p = new TH2D("Hist_Phi_RF_new_1_p", "Hist_Phi_RF_new_1_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_2 = new TH2D("Hist_Phi_FF_2", "Hist_Phi_FF_2", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_2_p = new TH2D("Hist_Phi_FF_2_p", "Hist_Phi_FF_2_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_2 = new TH2D("Hist_Phi_RF_2", "Hist_Phi_RF_2", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_2_p = new TH2D("Hist_Phi_RF_2_p", "Hist_Phi_RF_2_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_2 = new TH2D("Hist_Phi_FF_new_2", "Hist_Phi_FF_new_2", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_2 = new TH2D("Hist_Phi_RF_new_2", "Hist_Phi_RF_new_2", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_2_p = new TH2D("Hist_Phi_FF_new_2_p", "Hist_Phi_FF_new_2_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_2_p = new TH2D("Hist_Phi_RF_new_2_p", "Hist_Phi_RF_new_2_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_3 = new TH2D("Hist_Phi_FF_3", "Hist_Phi_FF_3", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_3_p = new TH2D("Hist_Phi_FF_3_p", "Hist_Phi_FF_3_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_3 = new TH2D("Hist_Phi_RF_3", "Hist_Phi_RF_3", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_3_p = new TH2D("Hist_Phi_RF_3_p", "Hist_Phi_RF_3_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_3 = new TH2D("Hist_Phi_FF_new_3", "Hist_Phi_FF_new_3", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_3 = new TH2D("Hist_Phi_RF_new_3", "Hist_Phi_RF_new_3", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_FF_new_3_p = new TH2D("Hist_Phi_FF_new_3_p", "Hist_Phi_FF_new_3_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH2D *Hist_Phi_RF_new_3_p = new TH2D("Hist_Phi_RF_new_3_p", "Hist_Phi_RF_new_3_p", Phibin, -PI, PI, 4, 0.5, 4.5);
TH1D *hDpt   = new TH1D("hDpt", "hDpt", 200, 0, 2);

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

TProfile *Hist_v2_pt_obs1 = new TProfile("Hist_v2_pt_obs1", "Hist_v2_pt_obs1", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs1_p = new TProfile("Hist_v2_pt_obs1_p", "Hist_v2_pt_obs1_p", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_Ach = new TProfile("Hist_v2_Ach", "V2 against Charge Assymetry", 300, -1.5 + 0.0025, 1.5 + 0.0025, -100, 100);
TProfile *Hist_v2_Ach_p = new TProfile("Hist_v2_Ach_p", "V2 against Charge Assymetry protons", 300, -1.5 + 0.0025, 1.5 + 0.0025, -100, 100);
TProfile *Hist_v2_pt_obs2 = new TProfile("Hist_v2_pt_obs2", "Hist_v2_pt_obs2", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs2_caysm[5];
TProfile *Hist_v2_pt_obs2_p = new TProfile("Hist_v2_pt_obs2_p", "Hist_v2_pt_obs2_p", 300, 0, 15, -100, 100, "");
TProfile *Hist_v2_pt_obs2_p_caysm[5];
TProfile *Hist_v2_eta_obs1 = new TProfile("Hist_v2_eta_obs1", "Hist_v2_eta_obs1", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_obs2 = new TProfile("Hist_v2_eta_obs2", "Hist_v2_eta_obs2", 40, -1, 1, -100, 100, "");
TProfile *Hist_v2_eta_obs3 = new TProfile("Hist_v2_eta_obs3", "Hist_v2_eta_obs3", 40, -1, 1, -100, 100, "");

TProfile *pTemp_v2 = new TProfile("pTemp_v2", "pTemp_v2", 6, 0.5, 6.5, -100, 100, "");
TProfile *pTemp_parity_e = new TProfile("pTemp_parity_e", "pTemp_parity_e", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_parity_w = new TProfile("pTemp_parity_w", "pTemp_parity_w", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_delta = new TProfile("pTemp_delta", "pTemp_delta", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_v2_noHBT = new TProfile("pTemp_v2_noHBT", "pTemp_v2_noHBT", 6, 0.5, 6.5, -100, 100, "");
TProfile *pTemp_parity_e_noHBT = new TProfile("pTemp_parity_e_noHBT", "pTemp_parity_e_noHBT", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_parity_w_noHBT = new TProfile("pTemp_parity_w_noHBT", "pTemp_parity_w_noHBT", 8, 0.5, 8.5, -100, 100, "");
TProfile *pTemp_delta_noHBT = new TProfile("pTemp_delta_noHBT", "pTemp_delta_noHBT", 8, 0.5, 8.5, -100, 100, "");

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
TH1D *lam_pt[9], *p_pt[9];

//defining variables
float pVx, pVy, pVz, VPDvz, BBCco, ZDCcoin, net_Nch_Asym, mQx, mQy;                   //run, event info
int   Run, Day, Day2, Day3, Trigger, RefMult, TOFMult, Centrality, NPTracks, n_Lambdas, n_protons, Fcount, Scount;   //
int   Charge, Charge2, ChargeAsso;
int net_charge_asym_bin = -1;
float ndEdx, nSigma_p, nSigma_pi, DCAGlobal, Eta, Theta, Phi, Pt, eff, TOFflag;             //track info
float ndEdx2, nSigma_p2, nSigma_pi2, DCAGlobal2, Eta2, Theta2, Phi2, Pt2, eff2, TOFflag2;   //2nd track info
float DCAGlobalAsso, EtaAsso, PhiAsso, PtAsso, TOFflagAsso;
TVector3 pV;

//struct StRefMultCorr refmultCorrUtil  = StRefMultCorr("refmult") ;
float Eweight = 1;
int Weight_Read = 0;
TH1D *TOF_eff = 0;
TProfile2D *TPCmean_FF_1 = 0, *TPCmean_RF_1 = 0, *TPCmean_FF_1_p = 0, *TPCmean_RF_1_p = 0, *TPCmeanAsso_FF_1 = 0, *TPCmeanAsso_RF_1 = 0;
TProfile2D *TPCmean_FF_2 = 0, *TPCmean_RF_2 = 0, *TPCmean_FF_2_p = 0, *TPCmean_RF_2_p = 0, *TPCmeanAsso_FF_2 = 0, *TPCmeanAsso_RF_2 = 0;
TProfile2D *TPCmean_FF_3 = 0, *TPCmean_RF_3 = 0, *TPCmean_FF_3_p = 0, *TPCmean_RF_3_p = 0, *TPCmeanAsso_FF_3 = 0, *TPCmeanAsso_RF_3 = 0;
float PhiMean_sin[order] = {0}, PhiMean_cos[order] = {0}, PhiMean_sin_p[order] = {0}, PhiMean_cos_p[order] = {0}, PhiMeanAsso_sin[order] = {0}, PhiMeanAsso_cos[order] = {0};
vector<float> PhiAsso_new, Phi_new, Phi_pnew;
TProfile2D *Read_TPC_EP_full = 0, *Read_TPC_EP_east = 0, *Read_TPC_EP_west = 0, *Read_TPC_EP_for = 0, *Read_TPC_EP_bac = 0;
float PsiMean_F[2 * order] = {0}, PsiMean_E[2 * order] = {0}, PsiMean_W[2 * order] = {0}, PsiMean_f[2 * order] = {0}, PsiMean_b[2 * order] = {0};
TProfile2D *Read_BBC_EP_east = 0, *Read_BBC_EP_west = 0, *Read_EPD_EP_east = 0, *Read_EPD_EP_west = 0;
float Psi_BBC_E[2 * order] = {0}, Psi_BBC_W[2 * order] = {0}, Psi_EPD_E[2 * order] = {0}, Psi_EPD_W[2 * order] = {0};
float MeanNetChargeAsym, RMSNetChargeAsym, StdDevNetChargeAsym;

float TPC_EP_full = 0, TPC_EP_east = 0, TPC_EP_west = 0, TPC_EP_for = 0, TPC_EP_bac = 0;
float TPC_EP_full_new = 0, TPC_EP_east_new = 0, TPC_EP_west_new = 0, TPC_EP_for_new = 0, TPC_EP_bac_new = 0;
float Q2 = 0;
float v2 = 0, v2e = 0, v2w = 0, v2_sub = 0, v2p = 0, v2pe = 0, v2pw = 0;
float correlator0 = 0, correlator3 = 0, correlator4 = 0, correlator4e = 0, correlator4w = 0;

bool debug_1, debug_2;

//Structure 2 - particles
//particle *p_lambda = new particle[5000];
//particle *p_proton = new particle[5000];
//particle *p_all = new particle[5000];
//particle *p_dau = new particle[5000];

std::vector<particle> *p_lambda = new vector<particle>;
std::vector<particle> *p_all = new vector<particle>;
std::vector<particle> *p_dau = new vector<particle>;
std::vector<particle> *p_proton = new vector<particle>;

//Structure 3 - Event details
event *details = new event();

//sub-functions
void WriteHistogram(int c, int o, TString JobIDName);                              //into result ROOT file
void WriteWeight(TString OutFileName);                            //into weight file
int  ReadWeight(char *InFileName);                              //input from weight file
bool IsGoodEvent(int c);                                        //select good events and fill event level histograms
//bool IsGoodBBC(StPicoEvent *e);                                 //cuts on BBC ADCs
//bool IsGoodZDC(StPicoEvent *e);                                 //cuts on ZDC ADCs
bool IsGoodAsso(float p, float e, float d);                     //cuts on EP particles
bool IsGoodPOI(float p, float e, float d);                      //cuts on particles of interest
//bool IsGoodTrack(StPicoTrack *p);                               //cuts on tracks
bool IsGoodPion(StPicoDst *d, StPicoTrack *p, int opt);         //cuts on pions
bool IsGoodKaon(StPicoDst *d, StPicoTrack *p, int opt);         //cuts on kaons
bool IsGoodProton(StPicoDst *d, StPicoTrack *p, int opt);       //cuts on protons
bool CountCharge();                                 //count good tracks, charge asymmetry
void MakeTPC_EP();                        //reconstruct TPC EPs
//void MakeBBC_EP(StPicoEvent *ev);                               //reconstruct BBC EPs
//void MakeZDC_EP(StPicoEvent *ev);                               //reconstruct ZDC EPs
void FillEP_resolution();                                       //Fill profiles for EP resolution
void FillPhiPOI();                                              //particles of interest, shift parameters to make phi distribution flat
void FillPhiPOI_p();                                              //particles of interest, shift parameters to make phi distribution flat
void FillPhiAsso();                                             //particles for EP, shift parameters to make phi distribution flat
void ShiftPhiAsso(int tr);                                      //flatten the phi distribution
void ShiftPhiPOI(int tr);                                       //flatten the phi distribution
void ShiftPhiPOI_p(int tr);                                       //flatten the phi distribution
void ShiftPsi();                                                //flatten EPs
void FillCMW();                                                 //v2 for CMW
void FillGamma(int ord);                                        //correlations
//void CountMSC(int *iTr, float Phi_ne);                          //MSC correlator
//void WrapUpMSC();                                               //MSC
//void WrapUpESE();                                               //Event Shape Engineering
//float GetPhiInBBC(int e_w, int bbcN);                           //input est_wst(0,1),BBC PMT#
//float GetXYInZDC(int e_w, int v_h, int zdcN, int opt_raw = 0);  //input ver_hor(0,1),ZDCSMD slat#
void set_debug();
