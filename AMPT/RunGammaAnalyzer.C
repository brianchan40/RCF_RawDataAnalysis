/**
 *  Author: Grigory Nigmatkulov
 *  Date:   July 5, 2018
 *
 *  Description:
 *  This macros takes inFileName argument with a picoDst.root file
 *  or with a list of files (name.lis or name.list). It sets _VANILLA_ROOT_
 *  (necessary for standalone and can be skipped on RACF), loads pre compiled
 *  libStPicoDst.so (from StPicoEvent), compiles and executes a text
 *  PicoDstAnalyzer.C macro with passing inFileName to it, and
 *  cleans up the directory from the compilation products at the end.
 *
 *  Some details:
 *    inFileName - is a name of name.picoDst.root file or a name
 *                 of a name.lis(t) files that contains a list of
 *                 name1.picoDst.root files.
 *    NOTE: inFileName should contain either /absolutePath/inFileName
 *          or /relative2currentDir/inFileName
 *  It is assumed that PicoDstAnalyzer.C is placed in the same
 *  directory where the RunAnalyzer.C is stored.
 **/

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include </star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/structs.h>

//_________________
void RunGammaAnalyzer(const int cen = 1, const int opt_weight = 1, const Char_t *inFile = "/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/file.list", const TString JobIDName = "1234", const TString lam_type = "lam") {
  // Next line is not needed if you are not running in a standalone mode
  gROOT->ProcessLine("#define _VANILLA_ROOT_");
  //gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StRefMultCorr/StRefMultCorr_cxx.so");
  //gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/StRefMultCorr/CentralityMaker_cxx.so");
  gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/particle_cxx.so");
  gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/event_cxx.so");
  gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/AMPT/namespaces/gv_gamma_cpp.so");
  gSystem->Load("/star/u/brian40/PicoDst/PicoDst_NoMaker/Version2/StRoot/libStPicoDst.so");
  TString str;
  str = ".x Gamma_112_module.C+(";
  str += cen;
  str += ",";
  str += opt_weight;
  str += ", \"";
  str += inFile;
  str += "\", \"";
  str += JobIDName;
  str += "\", \"";
  str += lam_type;
  str += "\")";
  cout<<str.Data()<<endl;
  gROOT->ProcessLine( str.Data() );
  // Next line should be commented if you run in a batch mode
  gROOT->ProcessLine(".!rm -f Gamma_112_module_C* ");
}
