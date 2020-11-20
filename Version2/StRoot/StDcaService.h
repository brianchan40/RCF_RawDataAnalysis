#ifndef StDcaService_def
#define StDcaService_def

#include <TROOT.h>
#include <TMath.h>
#include "StPicoHelix.h"
#include "StPicoPhysicalHelix.h"

//this struct is for Hui LONG's dcaPToPi function ONLY. 
//it is just the 'Helix' in his code. 
//some useless items are commented out.
struct StTrackHelix {
  //float  pid;
  //short  id;
  double Xc;
  double Yc;
  double Zc;
  double pt;
  //double X;
  //double Y;
  //double Z;
  //double Px;
  //double Py;
  double Pz;
  double r;
  float  theta;
  float  Vz;
  int    h;
  //int    h2;
  int    Flag;
  //int    q;
  //int    nhits;
  //float  dca;
  //float  dedx;
};

extern bool kStHelixDca;
extern bool kMinimize;
extern float kShiftConnect;
extern float kShiftContain;

//wrapper to Hui LONG's two helix dca code
double closestDistance(const StPicoPhysicalHelix& helix1, const StPicoPhysicalHelix& helix2, double magnet, const TVector3& pv, TVector3& xv0, TVector3& op1, TVector3& op2); 
void   dcaPToPi(float *fi_p,float *fi_pi, StTrackHelix* proton,StTrackHelix* pion,float *d_root,float *v0,float alfa); 			  

//wrapper to Hui LONG's dca to primary vertex code (OBSOLETE)
double getDcaToPV(const StPicoPhysicalHelix& helix, const TVector3& pv);
//wrapper to Hui LONG's helix to straight line dca code (for Xi or Omega reconstruction)(OBSOLETE)
double closestDistance(const TVector3& xv0, const TVector3& pv0, const StPicoPhysicalHelix& helix, double magnet, TVector3& xxi, TVector3& op );

#endif
