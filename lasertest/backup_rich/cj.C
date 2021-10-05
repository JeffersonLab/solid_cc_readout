#include <iostream> 
#include <fstream>
#include <cmath> 
#include "math.h" 
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TPaveText.h"
#include "TText.h"
#include "TSystem.h"
#include "TArc.h"
#include "TString.h"
#include <vector>
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"

using namespace std;


void cj(string dir)
{
gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
// gStyle->SetPadRightMargin(0.32);

const int n=192;

TH1F *hh[n];
 	for(int j=0;j<n;j++){	
 	  char hstname[100];
		sprintf(hstname,"h%03i",j);
		hh[j]=new TH1F(hstname,Form("pixel%03i;injected;adc counts",j),81,0,4050);
	}	

int m = 81;
//int m = 81*64;

for(int i=1;i<m;i++){

  char filename[100];
	sprintf(filename,"%s/run_%06i.bin.hist.root",dir.c_str(),i);

  TFile *file=new TFile(filename);
 
 	for(int j=0;j<n;j++){	
    if (j !=i/81 && j !=i/81+64 && j !=i/81+128) continue; 
		//cout << j << endl;
	  char hstname[100];
		sprintf(hstname,"hspe%03i",j);

  //cout << hstname << endl;
  h=(TH1F*) file->Get(hstname);
  hh[j]->Fill(i*50,h->GetMean(1));
	}
  
 file->Close();
}

TCanvas *c[3];
for(int i=0;i<3;i++){
c[i] = new TCanvas(Form("c%i",i),Form("c%i",i),2000,1100);
c[i]->Divide(8,8);
}

 	for(int j=0;j<n;j++){	
	i=j/64;
	c[i]->cd(j%64+1);
	hh[j]->SetMaximum(3000);
	hh[j]->SetMinimum(0);
	hh[j]->Draw();
	}

}



