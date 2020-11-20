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


void readoutdata(char const * inputfile, char const * histname, char const * regionfile){

  char * infile0 = new char[400];
  sprintf(infile0,"%s",inputfile);
  TFile *f = new TFile(infile0);
  sprintf(infile0,"%s",histname);
  char * regionfile0 = new char[400];
  sprintf(regionfile0,"%s",regionfile);

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

  int cutids[100];
  int j=0;  
  ifstream myfile;
  myfile.open (regionfile0);
  if (myfile.is_open())
    {
      while(!myfile.eof())
	{
	  string number;
	  int data;
	  getline (myfile,number);
	  data = atoi(number.c_str());
	  if (data!=0)
	    {
	      cutids[j]=data;
	      //   cout << data << '\n';
	      j++;
	    }
	}
      myfile.close();
    }
  

  double xmean=myHist->GetMean(2);
  double xsigma=myHist->GetMeanError(2)*sqrt(goodrunno);
  double xlowlimit=xmean-20*xsigma;
  double xhighlimit=xmean+20*xsigma;

  for(int i=1;i<myHist->GetNbinsX();i++){if(myHist->GetBinContent(i)!=0){ if(xlowlimit>=myHist->GetBinContent(i)){xlowlimit=myHist->GetBinContent(i);}}};
  for(int i=1;i<myHist->GetNbinsX();i++){if(myHist->GetBinContent(i)!=0){ if(xhighlimit<=myHist->GetBinContent(i)){xhighlimit=myHist->GetBinContent(i);}}};
  
  TH1D *h1 = new TH1D("h1","Check RMS",300,xlowlimit,xhighlimit);
  double MM, RR, meanh1, rmsh1;
  // int ii=1;
  for(int jj=0;jj<(j-1);jj++){

    for(int i=1;i<myHist->GetNbinsX();i++){if(myHist->GetBinCenter(i)>=cutids[jj]&&myHist->GetBinCenter(i)<cutids[jj+1]){if((myHist->GetBinContent(i)!=0)&&(myHist->GetBinError(i)!=0)){h1->Fill(myHist->GetBinContent(i), 1./(myHist->GetBinError(i))/(myHist->GetBinError(i)));}}};

    /// Use fit to get the RMS may different from GetRMS()
    MM = h1->GetMean();
    RR = h1->GetRMS();
    /// Need to make sure the RR is not to small
    //    if (RR < 0.001){RR = 0.001;}
    //    TF1 *fun = new TF1("fun","[0]*exp(-(x-[1])*(x-[1])/2./[2]/[2])",0,10000);
    //    fun->SetParameters(h1->Integral(),MM,RR);
    //    TFitResultPtr r = h1->Fit("fun","QNSE","",MM-10*RR,MM+10*RR);
    //    Double_t ndf = r->Ndf();
    //    Double_t prob = r->Prob();
    // Normally the fitting will be good, but if the histgrams has only 1 entry
    // it will be have some issues
    //    if (ndf==0 || prob==0){
      meanh1 = MM;
      rmsh1 = RR;
      //    }else{
      //      meanh1 = fun->GetParameter(1);
      //      rmsh1 = fabs(fun->GetParameter(2));
      //    }
    /// Fitting End

    ///old code, use h1->GetRMS()
    //   for(int i=1;i<myHist->GetNbinsX();i++){if(myHist->GetBinCenter(i)>=cutids[jj]&&myHist->GetBinCenter(i)<cutids[jj+1]){if(myHist->GetBinContent(i)!=0){fprintf(stdout,"%ld %g %g %g %g %d %g \n",long(myHist->GetBinCenter(i)),myHist->GetBinContent(i),myHist->GetBinError(i),myHist->GetMean(2),h1->GetRMS(),goodrunno2,myHist->GetMeanError(2)*sqrt(goodrunno)); goodrunno2++;}}};

   ///new code, use fitting function to get the better RMS, and mean
   for(int i=1;i<myHist->GetNbinsX();i++){if(myHist->GetBinCenter(i)>=cutids[jj]&&myHist->GetBinCenter(i)<cutids[jj+1]){if(myHist->GetBinContent(i)!=0){fprintf(stdout,"%ld %g %g %g %g %d %g \n",long(myHist->GetBinCenter(i)),myHist->GetBinContent(i),myHist->GetBinError(i),meanh1,rmsh1,goodrunno2,myHist->GetMeanError(2)*sqrt(goodrunno)); goodrunno2++;}}};
      h1->Reset();    
  }
  

  f->Close();  


}  






