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

// Data From Tree
//std::vector<particle> *p_all = new vector<particle>;
//std::vector<particle> *particle1 = new vector<particle>;
//std::vector<particle> *particle2 = new vector<particle>;
int event = 0, pmult = 0, refmult2 = 0, refmult3 = 0, refmult = 0, npp = 0, npt = 0, nesp = 0, ninesp = 0, nest = 0, ninest = 0;
float imp = 0;
int pid[10000] = {0};
float px[10000] = {0.};
float py[10000] = {0.};
float pz[10000] = {0.};
//std::vector<int> *pid = new vector<int>;
//std::vector<float> *px = new vector<float>;
//std::vector<float> *py = new vector<float>;
//std::vector<float> *pz = new vector<float>;

// Track Data for Analysis
float Eta = 0, Pt = 0, Charge = 0, Phi = 0, Mom_Mag = 0;
TDatabasePDG *db = new TDatabasePDG();
TParticlePDG *p_info = new TParticlePDG();


// Event QA Plots
TH1F *imp_dist = new TH1F("imp_dist", "Distribution of Impact Parameter", 250000, 0, 25);
TH1F *nLambda_dist = new TH1F("nLambda_dist", "Distribution of Lambdas", 250, 0, 25);
TH1F *nAntiLambda_dist = new TH1F("nAntiLambda_dist", "Distribution of AntiLambdas", 250, 0, 25);
TH1F *refmult_dist = new TH1F("refmult_dist", "Distribution of Refmult", 450, 0, 450);
TH1F *refmult2_dist = new TH1F("refmult2_dist", "Distribution of Refmult2", 450, 0, 450);
TH1F *refmult3_dist = new TH1F("refmult3_dist", "Distribution of Refmult3", 450, 0, 450);
TH1F *pmult_dist = new TH1F("pmult_dist", "Distribution of pmult", 450, 0, 450);
TH2D *imp_vs_refmult2 = new TH2D("imp_vs_refmult2", "Impact Parameter vs. Refmult2", 250, 0, 25, 450, 0, 450);
TH2D *imp_vs_refmult3 = new TH2D("imp_vs_refmult3", "Impact Parameter vs. Refmult3", 250, 0, 25, 450, 0, 450);
TH2D *imp_vs_pmult = new TH2D("imp_vs_pmult", "Impact Parameter vs. pmult", 250, 0, 25, 450, 0, 450);
TH2D *imp_vs_refmult = new TH2D("imp_vs_refmult", "Impact Parameter vs. Refmult", 250, 0, 25, 450, 0, 450);
TH2D *pmult_vs_refmult = new TH2D("pmult_vs_refmult", "pmult vs. Refmult", 450, 0, 450, 450, 0, 450);
TH2D *pmult_vs_refmult2 = new TH2D("pmult_vs_refmult2", "pmult vs. Refmult2", 450, 0, 450, 450, 0, 450);
TH1F *npp_dist = new TH1F("npp_dist", "Distribution of npp", 450, 0, 450);
TH1F *npt_dist = new TH1F("npt_dist", "Distribution of npt", 450, 0, 450);
TH1F *nesp_dist = new TH1F("nesp_dist", "Distribution of nesp", 450, 0, 450);
TH1F *ninesp_dist = new TH1F("ninesp_dist", "Distribution of ninesp", 450, 0, 450);
TH1F *nest_dist = new TH1F("nest_dist", "Distribution of nest", 450, 0, 450);
TH1F *ninest_dist = new TH1F("ninest_dist", "Distribution of ninest", 450, 0, 450);
TH1F *diff_npp_dist = new TH1F("diff_npp_dist", "Difference between npp and (nesp + ninesp)", 400, -200, 200);
TH1F *diff_npt_dist = new TH1F("diff_npt_dist", "Difference between npt and (nest + ninest)", 400, -200, 200);

// Track QA Plots
TH1F *px_dist = new TH1F("px_dist", "Distribution of Px", 100, -5, 5);
TH1F *py_dist = new TH1F("py_dist", "Distribution of Py", 100, -5, 5);
TH1F *pz_dist = new TH1F("pz_dist", "Distribution of Pz", 100, -5, 5);
TH1F *pt_dist = new TH1F("pt_dist", "Distribution of Pt", 100, -5, 5);
TH1F *p_dist = new TH1F("p_dist", "Distribution of |P|", 100, -5, 5);
TH1F *eta_dist = new TH1F("eta_dist", "Distribution of Eta", 100, -TMath::Pi(), TMath::Pi());
TH1F *phi_dist = new TH1F("phi_dist", "Distribution of Phi", 100, -TMath::Pi(), TMath::Pi());
TH1F *charge_dist = new TH1F("charge_dist", "Distribution of Charge", 10, -5, 5);
TH1I *pid_dist = new TH1I("pid_dist", "Distribution of PID", 10000, -5000, 5000);




void Fill_EventQAPlots();
void Get_BasicTrackInfo(int j);
void Fill_TrackQAPlots(int j);
void Write_QAPlots(TString JobIDName);
