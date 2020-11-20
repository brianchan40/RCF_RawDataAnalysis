#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFitResultPtr.h"

using namespace std;


void readoutnrms(char const * inputfile, char const * histname, char const * regionfile){

  char * infile0 = new char[400];
  sprintf(infile0,"%s",inputfile);
  TFile *f = new TFile(infile0);
  sprintf(infile0,"%s",histname);
  char * regionfile0 = new char[400];
  sprintf(regionfile0,"%s",regionfile);

  TProfile *myHist = (TProfile *)f->Get(infile0);
  //  myHist->Sumw2(); // it's no need to use sumw2 here

  //////////////////////
  int goodrunno=1;
  int goodrunno2=1;
  double content=0.;
  double entries=0.;
  double error=0.;
  double rms=0.;
  double factor=1.;
  for(int i=1;i<=myHist->GetNbinsX();i++){
    if(myHist->GetBinContent(i)!=0){
      goodrunno++;
      content += myHist->GetBinContent(i)*myHist->GetBinEntries(i);
      entries += myHist->GetBinEntries(i);
      error += pow(myHist->GetBinError(i),2)*pow(myHist->GetBinEntries(i),2);
    }
  }

  content = (entries>0 ? content/entries : 0);
  error = (entries>0 ? error/pow(entries,2) : 0);

  //Check this link for the weighted St Dev
  //https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
  
  for(int i=1;i<=myHist->GetNbinsX();i++){
    if(myHist->GetBinContent(i)!=0){
      rms += pow(content-myHist->GetBinContent(i),2)*myHist->GetBinEntries(i);
      //      rms += pow(myHist->GetMean(2)-myHist->GetBinContent(i),2);
      //      rms += pow(myHist->GetBinContent(i),2);
    }
  }
  factor = (double(goodrunno)-1)/double(goodrunno);
  rms = (entries>0 ? rms/entries/factor : 0);  
  //  rms = (entries>0 ? rms/goodrunno : 0);
  //    cout<< "Runs " << goodrunno << "; Mean " << content << "; Error " << sqrt(error) << "; RMS " << myHist->GetStdDev(2) << "; rms "<< sqrt(rms) <<endl;

    double correction_factor=1.;
    double SanityIndex=0.;
    // if we use 27GeV as the standard, we can correct the 5-weighted-error into a different n value, in case the quantity is so bad and get so many jumps. 288.881/5=57.776107 is the ratios for 27 GeV.
    // 57.776107 ~ 1-weighted-error; 
    correction_factor =  sqrt(rms)/sqrt(error)/57.776107;
    SanityIndex =  sqrt(rms)/sqrt(error);
    //    correction_factor =  sqrt(rms)/sqrt(error)/50;
    correction_factor = (correction_factor>5. ? correction_factor : 5. );

    cout << int(correction_factor) <<  ' ' << SanityIndex << endl;
  

  f->Close();  


}  






