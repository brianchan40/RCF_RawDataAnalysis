#ifndef StXiDst_def
#define StXiDst_def

#define MAX_NUM_Xi  1000
struct StXiDst{ 
   //event information, copy StV0Dst
   int runnumber;
   int evtnumber;
   int trgmode;
   int nrefmult; 
   int nrefmultftpceast;
   int bbcadcsumeast;
   float zdcadc0;
   double bbccirate;
   float primvertexX;
   float primvertexY;
   float primvertexZ;
   float magn;

   int nxi;

   //v0 information, copy StV0Dst
   float v0mass[MAX_NUM_Xi];
   //float v0pt[MAX_NUM_Xi];
   //float v0rapidity[MAX_NUM_Xi];
   //float v0eta[MAX_NUM_Xi];
   float v0x[MAX_NUM_Xi];
   float v0y[MAX_NUM_Xi];
   float v0z[MAX_NUM_Xi];
   float v0px[MAX_NUM_Xi];
   float v0py[MAX_NUM_Xi];
   float v0pz[MAX_NUM_Xi];
   //float v0declen[MAX_NUM_Xi];
   float v0dca[MAX_NUM_Xi];
   //float v0dca2d[MAX_NUM_Xi];
   float v0pathlen[MAX_NUM_Xi];

   int   dau1id[MAX_NUM_Xi];
   float dau1dca[MAX_NUM_Xi];
   //float dau1dca2d[MAX_NUM_Xi];
   //int dau1nhits[MAX_NUM_Xi];
   //float dau1dedx[MAX_NUM_Xi];
   //float dau1nsigma[MAX_NUM_Xi];
   //float dau1eta[MAX_NUM_Xi];
   //float dau1pt[MAX_NUM_Xi];
   float dau1px[MAX_NUM_Xi];
   float dau1py[MAX_NUM_Xi];
   float dau1pz[MAX_NUM_Xi];
   /*int dau1tpc[MAX_NUM_Xi];
   int dau1ssd[MAX_NUM_Xi];
   int dau1svt[MAX_NUM_Xi];
   int dau1tofflag[MAX_NUM_Xi];
   float dau1tof[MAX_NUM_Xi];
   float dau1pathlen[MAX_NUM_Xi];*/

   int   dau2id[MAX_NUM_Xi];
   float dau2dca[MAX_NUM_Xi];
   //float dau2dca2d[MAX_NUM_Xi];
   //int dau2nhits[MAX_NUM_Xi];
   //float dau2dedx[MAX_NUM_Xi];
   //float dau2nsigma[MAX_NUM_Xi];
   //float dau2eta[MAX_NUM_Xi];
   //float dau2pt[MAX_NUM_Xi];
   float dau2px[MAX_NUM_Xi];
   float dau2py[MAX_NUM_Xi];
   float dau2pz[MAX_NUM_Xi];
   /*int dau2tpc[MAX_NUM_Xi];
   int dau2ssd[MAX_NUM_Xi];
   int dau2svt[MAX_NUM_Xi];
   int dau2tofflag[MAX_NUM_Xi];
   float dau2tof[MAX_NUM_Xi];
   float dau2pathlen[MAX_NUM_Xi];*/

   //float dca1to2[MAX_NUM_Xi];

   //bachelor information
   int   bachid[MAX_NUM_Xi];
   float bachdca[MAX_NUM_Xi];
   float bachdca2d[MAX_NUM_Xi];
   int bachnhits[MAX_NUM_Xi];
   float bachdedx[MAX_NUM_Xi];
   float bachnsigma[MAX_NUM_Xi];
   float bacheta[MAX_NUM_Xi];
   float bachpt[MAX_NUM_Xi];
   float bachpx[MAX_NUM_Xi];
   float bachpy[MAX_NUM_Xi];
   float bachpz[MAX_NUM_Xi];
   int bachtpc[MAX_NUM_Xi];
   int bachssd[MAX_NUM_Xi];
   int bachsvt[MAX_NUM_Xi];
   int bachtofflag[MAX_NUM_Xi];
   float bachtof[MAX_NUM_Xi];
   float bachpathlen[MAX_NUM_Xi];
   
   float dcav0tobach[MAX_NUM_Xi];
   //float stdcav0tobach[MAX_NUM_Xi];

   //xi information
   float ximass[MAX_NUM_Xi];
   float xipt[MAX_NUM_Xi];
   float xirapidity[MAX_NUM_Xi];
   float xieta[MAX_NUM_Xi];
   float xix[MAX_NUM_Xi];
   float xiy[MAX_NUM_Xi];
   float xiz[MAX_NUM_Xi];
   float xipx[MAX_NUM_Xi];
   float xipy[MAX_NUM_Xi];
   float xipz[MAX_NUM_Xi];
   float xideclen[MAX_NUM_Xi];
   float xidca[MAX_NUM_Xi];
   float xidca2d[MAX_NUM_Xi];
   float xisinth[MAX_NUM_Xi];
   float xipathlen[MAX_NUM_Xi];
   
};

#endif
