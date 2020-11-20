#ifndef event_hh
#define event_hh

class event
{
public:
    event()
    {
        cent = -999;
        num_trk = -999;
        nLambda = -999;
        Run = -999;
        TOFMult = -999;
        RefMult = -999;
        n_Proton = -999;
        VPDvz = -3.14999;
        PVtxz = -3.14999;
        PVtxx = -3.14999;
        PVtxy = -3.14999;
        Eweight = -3.14999;
        EventID = -999;

        Magn = -3.14999;
    }
    event(Int_t in_cent, Int_t in_num_trk, Int_t in_nLambda, Int_t in_Run, Int_t in_TOFMult, Int_t in_RefMult, Int_t in_n_Proton, Double_t in_VPDvz, Float_t in_PVtxz, Float_t in_PVtxx, Float_t in_PVtxy, Float_t in_Eweight, Double_t in_Magn, Int_t in_EventID);

    virtual       ~event() { }

    Int_t cent, num_trk, nLambda, Run, TOFMult, RefMult, n_Proton, EventID;
    Double_t VPDvz;
    Float_t PVtxz, PVtxx, PVtxy, Eweight;
    
    Double_t Magn;
};



#endif
