#include "particle.h"

particle::particle(Float_t in_px, Float_t in_py, Float_t in_pz, Float_t in_x, Float_t in_y, Float_t in_z, Float_t in_charge, Int_t in_TOFflag, Float_t in_dcaglobal, Int_t in_prim, float in_nSigmaProton, Int_t in_isTofTrack, double in_mass, int in_trk_id, double in_ToF, double in_BTOFYLocal, float in_Mass2)
{
    px = in_px;
    py = in_py;
    pz = in_pz;
    x = in_x;
    y = in_y;
    z = in_z;
    Charge = in_charge;
    TOFflag = in_TOFflag;
    dcaglobal = in_dcaglobal;
    prim = in_prim;
    nSigmaProton = in_nSigmaProton;
    isTofTrack = in_isTofTrack;
    mass = in_mass;
    trk_id = in_trk_id;
    ToF = in_ToF;
    BTOFYLocal = in_BTOFYLocal;
    Mass2 = in_Mass2;
}
