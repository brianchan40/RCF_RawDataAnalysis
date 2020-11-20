#ifndef particle_hh
#define particle_hh

class particle
{
public:
    particle() {
        px = -3.14999;
        py = -3.14999;
        pz = -3.14999;
        x = -3.14999;
        y = -3.14999;
        z = -3.14999;
        Charge = -3.14999;
        TOFflag = -999;
        dcaglobal = -3.14999;
        prim = -999;
        nSigmaProton = -3.14999;
        isTofTrack = -999;
        mass = -3.14999;
        trk_id = -999;
        ToF = -3.14999;
        BTOFYLocal = -3.14999;
        Mass2 = -3.14999;
    }
    particle(Float_t in_px, Float_t in_py, Float_t in_pz, Float_t in_x, Float_t in_y, Float_t in_z, Float_t in_charge, Int_t in_TOFflag, Float_t in_dcaglobal, Int_t in_prim, float in_nSigmaProton, Int_t in_isTofTrack, double in_mass, int in_trk_id, double in_ToF, double in_BTOFYLocal, float in_Mass2);

    virtual       ~particle() { }

    Float_t px, py, pz;
    Float_t x, y, z;
    Float_t Charge;
    Int_t TOFflag;
    Float_t dcaglobal;
    Int_t prim;
    float nSigmaProton;
    Int_t isTofTrack;
    double mass;
    int trk_id;
    double ToF;
    double BTOFYLocal;
    float Mass2;
};



#endif
