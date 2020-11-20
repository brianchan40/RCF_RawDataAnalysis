#include "event.h"

event::event(Int_t in_cent, Int_t in_num_trk, Int_t in_nLambda, Int_t in_Run, Int_t in_TOFMult, Int_t in_RefMult, Int_t in_n_Proton, Double_t in_VPDvz, Float_t in_PVtxz, Float_t in_PVtxx, Float_t in_PVtxy, Float_t in_Eweight, Double_t in_Magn, Int_t in_EventID)
{
    cent = in_cent;
    num_trk = in_num_trk;
    nLambda = in_nLambda;
    Run = in_Run;
    TOFMult = in_TOFMult;
    RefMult = in_RefMult;
    n_Proton = in_n_Proton;
    VPDvz = in_VPDvz;
    PVtxz = in_PVtxz;
    PVtxx = in_PVtxx;
    PVtxy = in_PVtxy;
    Eweight = in_Eweight;
    
    Magn = in_Magn;
    
    EventID = in_EventID;
}
