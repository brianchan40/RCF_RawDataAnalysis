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

using namespace std;


void readoutroot(char const * inputfile, char const * histname){

  char * infile0 = new char[400];
  sprintf(infile0,"%s",inputfile);
  TFile *f = new TFile(infile0);
  sprintf(infile0,"%s",histname);

  TProfile *myHist = (TProfile *)f->Get(infile0);
    //myHist->Rebin($REBINX);
  fprintf(stdout,"## %s %s \n",myHist->GetXaxis()->GetTitle(),myHist->GetYaxis()->GetTitle());
    
  ///  cout<<" I am going to read this file "<<infile<<endl;
  ///char * his3 = new char[400];
  //////////////////////
  int goodrunno=1;
  int goodrunno2=1;

  for(int i=1;i<=myHist->GetNbinsX();i++){if(myHist->GetBinContent(i)!=0)goodrunno++;}

  fprintf(stdout,"# %d %d %g %g %g %g %g %g %g \n",0,goodrunno,1.*(myHist->GetMinimum()),1.*(myHist->GetMaximum()),myHist->GetEntries(),myHist->GetMean(1),myHist->GetMean(2),myHist->GetRMS(2),myHist->GetMeanError(2)*sqrt(goodrunno));

  ////// Here we prepare the RMS by fill a new TH1D
  double xmean=myHist->GetMean(2);
  double xsigma=myHist->GetMeanError(2)*sqrt(goodrunno);
  double xlowlimit=xmean-20*xsigma;
  double xhighlimit=xmean+20*xsigma;

  for(int i=1;i<myHist->GetNbinsX();i++){if(myHist->GetBinContent(i)!=0){ if(xlowlimit>=myHist->GetBinContent(i)){xlowlimit=myHist->GetBinContent(i);}}};
  for(int i=1;i<myHist->GetNbinsX();i++){if(myHist->GetBinContent(i)!=0){ if(xhighlimit<=myHist->GetBinContent(i)){xhighlimit=myHist->GetBinContent(i);}}};
  
  TH1D *h1 = new TH1D("h1","Check RMS",100,xlowlimit,xhighlimit);
  for(int i=1;i<myHist->GetNbinsX();i++){if((myHist->GetBinContent(i)!=0)&&(myHist->GetBinError(i)!=0)){h1->Fill(myHist->GetBinContent(i), 1./(myHist->GetBinError(i))/(myHist->GetBinError(i)));}};
  ////

  for(int i=1;i<myHist->GetNbinsX();i++){if(myHist->GetBinContent(i)!=0){fprintf(stdout,"%ld %g %g %g %g %d %g \n",long(myHist->GetBinCenter(i)),myHist->GetBinContent(i),myHist->GetBinError(i),myHist->GetMean(2),h1->GetRMS(),goodrunno2,myHist->GetMeanError(2)*sqrt(goodrunno)); goodrunno2++;}};

  
  ///  for(int i=1;i<myHist->GetNbinsX();i++){if(myHist->GetBinContent(i)!=0){fprintf(stdout,"%ld %g %g %g %g %d %g \n",long(myHist->GetBinCenter(i)),myHist->GetBinContent(i),myHist->GetBinError(i),myHist->GetMean(2),myHist->GetRMS(2),goodrunno2,myHist->GetMeanError(2)*sqrt(goodrunno)); goodrunno2++;}};


  f->Close();  


}  






