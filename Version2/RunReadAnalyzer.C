#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <./structs.h>

class StRefMultCorr;
class CentralityMaker;

//_________________
void RunReadAnalyzer(const int cen = 1, const int opt_weight = 1, const Char_t *inFileName = "./schedB954995E01B88BE24718672311AF0A4E_800.list", const TString JobIDName = "1234") {
  gROOT->ProcessLine("#define _VANILLA_ROOT_");
  //gROOT->ProcessLine("#include <vector>");
  //gSystem->Load("./StRoot/StRefMultCorr/StRefMultCorr.so");
  //gSystem->Load("./StRoot/StRefMultCorr/CentralityMaker_cxx.so");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("./StRoot/event_cxx.so");
  gSystem->Load("./StRoot/particle_cxx.so");
  gSystem->Load("./StRoot/gv_cpp.so");
  gSystem->Load("./StRoot/gv_lam_cpp.so");
  gSystem->Load("./StRoot/gv_xi_cpp.so");
  gSystem->Load("./StRoot/Run_by_RunQA_cpp.so");
  gSystem->Load("./StRoot/LambdaMaker_cpp.so");
  gSystem->Load("./StRoot/ProtonMaker_cpp.so");
  gSystem->Load("./StRoot/XiMaker_cpp.so");
  gSystem->Load("./StRoot/libStPicoDst.so");
  gSystem->Load("./namespaces/gv_gamma_cpp.so");
  TString str;
  str = ".x ParticleMaker.C+(";
  str += cen;
  str += ",";
  str += opt_weight;
  str += ", \"";
  str += inFileName;
  str += "\", \"";
  str += JobIDName;
  str += "\")";
  gROOT->ProcessLine( str.Data() );
  gROOT->ProcessLine(".!rm -f ParticleMaker_C* ");
}
