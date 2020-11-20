#ifndef structs_h
#define structs_h

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
//#include "./StRoot/StRefMultCorr.h"
//#include "./StRoot/CentralityMaker.h"
#include "./StRoot/StRefMultCorr/StRefMultCorr.h"
#include "./StRoot/StRefMultCorr/CentralityMaker.h"

/// Analysis headers
#include "./StRoot/StV0Type.h"
#include "./StRoot/StXiType.h"
#include "./StRoot/StDcaService.h"
#include "./StRoot/StV0Type.h"

struct condition
{
    StV0Type  mV0Type;
    StXiType mXiType;
    bool  mRotate;
    bool  mSameSignPlus;
    bool  mSameSignMinus;
    bool  mDcaAlgoLong;
    bool  mDumpNull;
    
    //choosing what to debug
    bool debug_Event;
    bool debug_priTrack;
    bool debug_gTrack;
    bool debug_gTrack2;
    bool debug_reconstruct;
};

struct event_cuts
{
    //statistic information
    UInt_t        mEventsProcessed;                 //  Number of Events read and processed
    
    //some diagnosing variables
    Float_t   mTestVZ;
    UInt_t    mTestNTrack;
    
    //section of parameters for v0 analysis
    float  cutAbsVertexZLeEq;
};

struct v0_const
{
    //section of v0type related constants
    Float_t   mMass1;
    Float_t mMass2;
    Float_t mMassV0;
    Int_t mCharge1;
    Int_t mCharge2;
    Float_t	mMassBachelor;
    Int_t 	mChargeBachelor;
    Float_t	mMassXi;
    Int_t	mChargeXi;
};

struct track_cuts
{
    //section of parameters for v0 analysis
    int    cutNHitsGr;
    float  cutPtGrEq;
    float  cutPtGrEq1;
    float  cutPtGrEq2;
    
    float  cutAbsNSigma1Le;
    float  cutAbsNSigma2Le;
    float  cutDca1GrEq;
    float  cutDca2GrEq;
    float  cutDca1LeEq;
    float  cutDca2LeEq;
    
    float  cutDca1to2LeEq;
    float  cutV0MassWidthLeEq;
    float  cutDauPtArmLeEq;
    float  cutAbsDausPtShoulderDiffLeEq;
    float  cutDau1DecAngGr;
    float  cutDau2DecAngGr;
    float  cutV0rdotpGr;
    float  cutDcaV0Le;
    float  cutV0DecLenGrEq;
    float  cutDau1Dau2Ang3DLe;
    float  cutDau1Dau2DipAngDiffLe;
    
    double nsigma;
    
    float  cutV0DcaGrEq;
    float  cutXiMassWidthLeEq;
    float  cutXirdotpGr;
    float  cutXircrosspLeEq;
    float  cutDcaXiLe;
    float  cutXiDecLenGrEq;
};

struct histograms_created
{
    /// Initializing Histograms
    TH1F   *hNPrimVertex;
    TH1F   *hVertexZ;
    TH1F   *hNRefMult;
    TH1F   *hSelectNRefMultM0;
    TH1F   *hSelectNRefMultM1;
    TH1F   *hSelectNRefMultM2;
    TH1F         *hSelectNRefMultM3;
    TH1F         *hSelectNRefMultM4;
    TH1F         *hSelectNRefMultM5;
    TH1F         *hSelectNRefMultM6;
    TH1F         *hSelectNRefMultM7;
    TH1F         *hSelectNRefMultM8;
    TH1F         *hSelectNRefMultM9;
    
    TH1F   *hSelectNRefMultFtpcE;
    TH1F   *hSelectBbcEAdcSum;
    TH1F   *hEta;
    TH1F   *hPhi;
    TH1F   *hPhiLowPt;
    TH1F   *hPhiHighPt;
    TH2F   *hDedxP;
    TH1F   *hPtRaw;
    TH1F   *hEtaRaw;
    TH1F   *hPhiRaw;
    TH1F   *hNSigmaPion;
    TH1F   *hNSigmaProton;
    TH1F   *hNSigmaKaon;
    TH1F   *hPtDiff;
    TH1F   *hOrDiff;
    TH1F   *hPDiff;
    TH1F   *hDcaDiff;
    
    TH1F   *hVertexZDiff;
    
    TH1D       *charge_dist;
    
    
    TH1F *hInvMass;
    TH1F *mom_track_gpt;
    TH1F *nHitsFitDist;
    TH1D *nsigmadist;
    TH1F *DCA_Dist;
    TH1F *dca1to2_Dist;
    
    TH1F *lambda_pt;
    
    //TH1F *V0Mass_pt_cent[14][7];
    TH1F* V0Mass_cent[9][17];
    
    TH3F *primary_vertex;
    
    TH1D *v0decaylengthQA;
    TH1D *dcav0toPVQA;
    TH1D *entireV0Mass;
    TH1D *transversevertexmag;
    TH1D *pseudoRapidityQA1;
    TH1D *pseudoRapidityQA2;
    TH1D *RapidityQA;
    TH1D *lambdarapidity;
    
    TH1I *cent_check_9;
    TH1I *cent_check_16;
    
    TH1I *trigger_IDs;
    
    TH1F *primary_p_perp;
    TH1F *primary_p_mag;
    TH1F *btofylocal_dist;
    TH1F *mass2_p;
};

struct v0_dau_vectors
{
    std::vector<StPicoTrack *> mDauVec;
    std::vector<double> mDauDcaVec;
        
    int PrimaryTrackID[5000];
    double PrimaryTrackPx[5000], PrimaryTrackPy[5000], PrimaryTrackPz[5000];
    int nPrimary;
    
    double mag_field;
};

struct TRK{
    StPicoTrack *track;
    StPicoPhysicalHelix helix;
    StPicoPhysicalHelix helixv0, helixxi;
    TVector3 p;
    double pt;
    double dca;
    int tofflag;
    double tof, tofpathlen, v0pathlen, xipathlen;
    TVector3 tofpos;
};

struct crucial_v0calc{
    TVector3 xv0, op1, op2;
    double dca1to2;
    TVector3 ox1, ox2;
    TLorentzVector Lambda_fvec;
};

struct crucial_xicalc{
    TVector3 pv; //Event details used in calculations
    TVector3 xxi, op1, op2, xv0, pv0, pxi, xxitoPV, pro_cross, dca;
    double dca1;
    TVector3 p1;
    double dca1tov0;
    double oe1, oev0, ximass, ang1, rdotp, pathlength, dcaxitoPV, xidecaylength, sinth;
};

struct histograms_xi{
    TH1F* hVertexZ;
    TH1F* hPt;
    TH1F* hNPrimVertex;
    TH1F* hNRefMult;
    TH1F* hPtDiff;
    TH1F* hOrDiff;
    TH1F* hPDiff;
    TH1F* hDcaDiff;
    TH1F* hInvMass;
    TH1F* hInvMassSTD;
};

struct lambda_info{
    double mass;
    double dcaV0toPv;
    
    Float_t x, y, z;
};

#endif
