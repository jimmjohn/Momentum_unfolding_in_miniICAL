#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TTree.h>
#include <TCanvas.h>
#include "TLegend.h"

#include "TVector.h"
#include <vector>
#include <TF1.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TStyle.h>
#include "TPostScript.h"
#include <TPad.h>
#include <TLine.h>
#include <TRandom.h>

#include "TUnfold.h"
#include "TUnfoldDensity.h"


using namespace std;


double probCut = 0.015;
//bool muMinus = false;
bool mlApplied = false;

bool var_p = true;

//int version = 108;

double chisquare_p = 0;
double chisquare_phi = 0;

//theta
int nthetabins_gen = 40;
int nthetabins_reco = 80;
double theta_low = 0.;
double theta_high = 90.;

//phi
int nphibins_gen = 40;
int nphibins_reco = 80;
double phi_low = -180.;
double phi_high = 180.;

//Momentum
int npbins_gen_init=8;
int npbins_reco = 16;
double p_low = 0.8;
double p_high = 3.0;


double binwidth_gen = (p_high-p_low)/npbins_gen_init;
double binwidth_reco = (p_high-p_low)/npbins_reco;

int extra_bins = (p_low - 0.5)/binwidth_gen;
int npbins_gen = npbins_gen_init;//+2*extra_bins;

void DivideHistogramByBinWidth(TH1D *histogram) {
  for (int i = 1; i <= histogram->GetNbinsX(); ++i) {
    double binContent = histogram->GetBinContent(i);
    double binWidth = histogram->GetBinWidth(i);
    histogram->SetBinContent(i, binContent / binWidth);
    histogram->SetBinError(i, histogram->GetBinError(i) / binWidth);
  }
}

int main(int argc, char** argv){
//void unfolding(){

  //RSA read to rootfile containing response matrix from all model
  bool muMinus = (atoi(argv[3])==1) ? true : false;
  
  int version = atoi(argv[1]);
  string saveDir = "fileOut";
  //string filename = "Mc_v108_2018_Data_Magnetic.root";
  string filename = argv[2];
  cout<<filename<<" "<<version<<endl;
  int myndof = atoi(argv[4]);

  
  
  
  

  char name[200];

  double p_low_temp = p_low;

  if(muMinus) {p_low = -1*p_high; p_high = -1*p_low_temp;}

  double pival = acos(-1.);
  const int numThetaRanges = 4;
  //double thetaRanges[numThetaRanges][2] = {{0, 17}, {17, 26}, {26, 34}, {34, 44}, {44,90}};
  double thetaRanges[numThetaRanges][2] = {{0, 17}, {17, 26}, {26, 34}, {34, 50}};

  //double thetaRanges[numThetaRanges][2] = {{0,90}};

  //double reco_bin_sch[23] = {0,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3,3.1};
  //double gen_bin_sch[13] = {0,1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3,3.1};

  //double reco_bin_sch[35] = {0, 1, 1.06149, 1.12298, 1.18446, 1.24595, 1.30744, 1.36893, 1.43041, 1.4919, 1.55339, 1.61488, 1.67636, 1.73785, 1.79934, 1.86083, 1.92231, 1.9838, 2.04529, 2.10678, 2.16826, 2.22975, 2.29124, 2.35273, 2.41421, 2.4757, 2.53719, 2.59868, 2.66016, 2.72165, 2.78314, 2.84463, 2.90611, 2.9676, 3.1};
  //double gen_bin_sch[17] = {0, 1, 1.14142, 1.28284, 1.42426, 1.56569, 1.70711, 1.84853, 1.98995, 2.13137, 2.27279, 2.41421, 2.55563, 2.69706, 2.83848, 2.9799, 3.1};

  //double reco_bin_sch[20] = {0, 1, 1.11314, 1.22627, 1.33941, 1.45255, 1.56569, 1.67882, 1.79196, 1.9051, 2.01823, 2.13137, 2.24451, 2.35765, 2.47078, 2.58392, 2.69706, 2.81019, 2.92333, 3.1};
  //double gen_bin_sch[10] = {0, 1, 1.28284, 1.56569, 1.84853, 2.13137, 2.41421, 2.69706, 2.9799, 3.1};

  //Variable bin width based on resolution MuPlus
  //double reco_bin_sch_muplus[18] = {0, 1, 1.11314, 1.22819, 1.3452, 1.4642, 1.58521, 1.70828, 1.83344, 1.96072, 2.09016, 2.2218, 2.35567, 2.49181, 2.63027, 2.77107, 2.91427, 3.1};
  //double gen_bin_sch_muplus[9] = {0, 1, 1.28284, 1.57769, 1.88504, 2.20543, 2.53941, 2.88757, 3.1};
  // double reco_bin_sch_muplus[47] = {0.8, 0.82, 0.85, 0.87, 0.9, 0.92, 0.95, 0.98, 1.01, 1.04, 1.07, 1.1, 1.13, 1.16, 1.2, 1.23, 1.27, 1.3, 1.34, 1.38, 1.42, 1.46, 1.51, 1.55, 1.59, 1.64, 1.69, 1.74, 1.79, 1.84, 1.89, 1.95, 2.01, 2.06, 2.13, 2.19, 2.25, 2.32, 2.38, 2.45, 2.52, 2.6, 2.67, 2.75, 2.83, 2.92, 3};
  // double gen_bin_sch_muplus[23] = {0.8, 0.85, 0.9, 0.96, 1.02, 1.08, 1.15, 1.22, 1.29, 1.37, 1.46, 1.55, 1.65, 1.75, 1.86, 1.97, 2.09, 2.22, 2.36, 2.51, 2.66, 2.83, 3};
  double reco_bin_sch_muplus[17] = {0.8, 0.87, 0.94, 1.02, 1.11, 1.21, 1.31, 1.43, 1.55, 1.68, 1.83, 1.98, 2.16, 2.34, 2.54, 2.76, 3};
  double gen_bin_sch_muplus[9] = {0.8, 0.94, 1.11, 1.31, 1.55, 1.83, 2.16, 2.54, 3};

  //Start from 0.6 - 3.5
  //0.7-0.85, 0.85-1.0

  //suppose 0.8-3.0 its 20 bins total bins =22, number of points =23
  //suppose 0.8-3.0 its 8 bins total bins =10 , number of points =11


  //Variable bin width based on resolution MuMinus
  //double reco_bin_sch_muminus[47] = {-3, -2.92, -2.83, -2.75, -2.67, -2.6, -2.52, -2.45, -2.38, -2.32, -2.25, -2.19, -2.13, -2.06, -2.01, -1.95, -1.89, -1.84, -1.79, -1.74, -1.69, -1.64, -1.59, -1.55, -1.51, -1.46, -1.42, -1.38, -1.34, -1.3, -1.27, -1.23, -1.2, -1.16, -1.13, -1.1, -1.07, -1.04, -1.01, -0.98, -0.95, -0.92, -0.9, -0.87, -0.85, -0.82, -0.8};
  //double gen_bin_sch_muminus[23] = {-3, -2.83, -2.66, -2.51, -2.36, -2.22, -2.09, -1.97, -1.86, -1.75, -1.65, -1.55, -1.46, -1.37, -1.29, -1.22, -1.15, -1.08, -1.02, -0.96, -0.9, -0.85, -0.8,};
  double reco_bin_sch_muminus[17] = {-3, -2.76, -2.54, -2.34, -2.16, -1.98, -1.83, -1.68, -1.55, -1.43, -1.31, -1.21, -1.11, -1.02, -0.94, -0.87, -0.8};
  double gen_bin_sch_muminus[9] = {-3, -2.54, -2.16, -1.83, -1.55, -1.31, -1.11, -0.94, -0.8};

  double reco_bin_sch[17];
  double gen_bin_sch[9];

  if(!muMinus) {
    std::copy(std::begin(reco_bin_sch_muplus), std::end(reco_bin_sch_muplus), std::begin(reco_bin_sch));
    std::copy(std::begin(gen_bin_sch_muplus), std::end(gen_bin_sch_muplus), std::begin(gen_bin_sch));
  } else {
    std::copy(std::begin(reco_bin_sch_muminus), std::end(reco_bin_sch_muminus), std::begin(reco_bin_sch));
    std::copy(std::begin(gen_bin_sch_muminus), std::end(gen_bin_sch_muminus), std::begin(gen_bin_sch));
  }


  //Momentum

  TH1D *hist_p_reco[numThetaRanges];
  for (int j = 0; j < numThetaRanges; ++j) {
    if(!muMinus) {sprintf(name,"h_p_mc_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_p_reco[j] = new TH1D(name,name,npbins_reco,reco_bin_sch); hist_p_reco[j]->Sumw2();}
    if(muMinus)  {sprintf(name,"h_p_mc_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_p_reco[j] = new TH1D(name,name,npbins_reco,reco_bin_sch); hist_p_reco[j]->Sumw2();}
  }


  //TH1D *hist_p_data = new TH1D(name,"",npbins_reco,reco_bin_sch); hist_p_data->Sumw2();
  TH1D *hist_p_data[numThetaRanges];
  for (int j = 0; j < numThetaRanges; ++j) {
    if(!muMinus) {sprintf(name,"h_p_data_reco_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_p_data[j] = new TH1D(name,name,npbins_reco,reco_bin_sch); hist_p_data[j]->Sumw2();}
    if(muMinus)  {sprintf(name,"h_p_data_reco_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      hist_p_data[j] = new TH1D(name,name,npbins_reco,reco_bin_sch); hist_p_data[j]->Sumw2();}
  }


  
  int nbins_gen, nbins_reco;
  double lowedge, upedge;

  double thgen, phgen, pgen;
  double threco, phreco, preco;
  double chi2reco;
  int ndfreco;

  //In tree
  Float_t momin[1], thein[1], phiin[1];
  Float_t trkmm[1];  
  Float_t momrf[1], therf[1], phirf[1];
  float chisquare[1];
  UInt_t ndof[1];
  UInt_t fitfailed;
  Float_t learnedVal, learnedErr;
  Float_t trk_para_before_ext_roof[10];
  Float_t trk_para_after_ext_roof[10];

  Float_t momrfd[1], therfd[1], phirfd[1];
  Float_t trkmmd[1];
  float chisquared[1];
  UInt_t ndofd[1];
  UInt_t fitfailedd;
  Float_t learnedVald, learnedErrd;
  Float_t trk_para_before_ext_roofd[10];
  Float_t trk_para_after_ext_roofd[10];

  //  sprintf(name,"%s",filename.c_str());
  sprintf(name,"/var/nfscondor/rajshah/paper2/mom_unfolding_fiducial/%s",filename.c_str());

  TFile *file1 = new TFile(name,"read");
  file1->cd();
  TTree *T1;
  T1 = (TTree*)file1->Get("T3");
  T1->SetBranchAddress("thein",thein);
  T1->SetBranchAddress("phiin",phiin);
  T1->SetBranchAddress("momin",momin);
  T1->SetBranchAddress("therf",therf);
  T1->SetBranchAddress("phirf",phirf);
  T1->SetBranchAddress("trkmm",trkmm);
  T1->SetBranchAddress("momrf",momrf);
  T1->SetBranchAddress("chisquare",chisquare);
  T1->SetBranchAddress("ndof",ndof);
  T1->SetBranchAddress("fitfailed",&fitfailed);
  T1->SetBranchAddress("learnedVal",&learnedVal);
  T1->SetBranchAddress("learnedErr",&learnedErr);
  T1->SetBranchAddress("trk_para_before_ext_roof",trk_para_before_ext_roof);
  T1->SetBranchAddress("trk_para_after_ext_roof",trk_para_after_ext_roof);


  TTree *T2;
  T2 = (TTree*)file1->Get("T5");
  T2->SetBranchAddress("therf",therfd);
  T2->SetBranchAddress("phirf",phirfd);
  T2->SetBranchAddress("trkmm",trkmmd);
  T2->SetBranchAddress("momrf",momrfd);
  T2->SetBranchAddress("chisquare",chisquared);
  T2->SetBranchAddress("ndof",ndofd);
  T2->SetBranchAddress("fitfailed",&fitfailedd);
  T2->SetBranchAddress("learnedVal",&learnedVald);
  T2->SetBranchAddress("learnedErr",&learnedErrd);
  T2->SetBranchAddress("trk_para_before_ext_roof",trk_para_before_ext_roofd);
  T2->SetBranchAddress("trk_para_after_ext_roof",trk_para_after_ext_roofd);

  int totalevents= 0;

  
	
  

  //for(int ientry=0; ientry<T1->GetEntries(); ientry++) {
  for(int ientry=0; ientry<2620000; ientry++) {//RSA    
    T1->GetEntry(ientry);
    //if(abs(momin[0])>3) {continue;}   To check to see it has any effect on assymettry

    if(!mlApplied && !muMinus) {learnedVal=1.;}
    if(!mlApplied && muMinus)  {learnedVal=0.;}

    thgen = thein[0]*180./pival; phgen = phiin[0]*180./pival; pgen = momin[0];
    threco = therf[0]*180./pival; phreco = phirf[0]*180./pival;
    preco = momrf[0];
    chi2reco = chisquare[0]; ndfreco = ndof[0];

    double prob = TMath::Prob(3.*chi2reco,ndfreco);

    bool selected =  (fitfailed==1 && ndfreco>=myndof && chi2reco/ndfreco<2)? true: false;
    //bool selected =  (fitfailed==1 && ndfreco>=5 && prob > probCut && abs(trkmm[0])>0.03 && abs(trk_para_before_ext_roof[2])<15 && abs(trk_para_before_ext_roof[3])<15 )? true: false;

    if(fitfailed==1 || fitfailed !=1){
      //if(totalevents>1000000) break;

      //////////////////MuPlus/////////////////////////////////////////////////////////////////////////

      if(selected && learnedVal>0.9 && !muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    hist_p_reco[j]->Fill(preco);
	  }
	}
	totalevents++;   
      } else if(!muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    hist_p_reco[j]->Fill(-100);
	  }
	}
      }

      //////////////MuMinus----------------------------------------------------------
      if(selected && learnedVal<0.1 && muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    hist_p_reco[j]->Fill(preco);
	  }
	}
	totalevents++;
      } else if(muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    hist_p_reco[j]->Fill(100);
	  }
	}
      }

      
    }
  }

  cout<<"MC TotalEvents Reco good"<<totalevents<<endl;

  totalevents=0;

  for(int ientry=0; ientry<T2->GetEntries(); ientry++){

    T2->GetEntry(ientry);

    if(!mlApplied && !muMinus) {learnedVald=1.;}
    if(!mlApplied && muMinus)  {learnedVald=0.;}

    threco = therfd[0]*180./pival; phreco = phirfd[0]*180./pival;
    //if(momrfd[0]>0) {preco = -0.014765 + 1.60720*momrfd[0];}//pol1 fit in range 0.6 to 2.0
    ////if(momrfd[0]< 0){preco = -0.0398542 + 1.46176*momrfd[0];}//pol1 fit in range -2.0 to -0.6
    ////if(momrfd[0]< 0){preco = -0.0162093 + 1.59788*momrfd[0];}//pol1 fit in range -2.0 to -0.6
    //if(momrfd[0]< 0){preco = -0.185643 + 1.21464*momrfd[0];}//pol1 fit in range -2.0 to -0.6
    ////if(momrfd[0]< 0){preco = momrfd[0]-0.2;}
    preco = momrfd[0];

    chi2reco = chisquared[0]; ndfreco = ndofd[0];

    double probd = TMath::Prob(4.*chi2reco, ndfreco);

    bool selected = (fitfailedd==1 && ndfreco>=myndof && chi2reco/ndfreco<2)? true: false;
    //bool selected =  (fitfailedd==1 && ndfreco>=5 && probd > probCut && abs(trkmmd[0])>0.03 && abs(trk_para_before_ext_roofd[2])<15 && abs(trk_para_before_ext_roofd[3])<15 )? true: false;

    if(fitfailedd==1 || fitfailedd!=1){
      //if(totalevents>1000000) break;
      if(selected && learnedVald>0.9 && !muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    hist_p_data[j]->Fill(preco);
	  }
	}
	totalevents++;
      } else if(!muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    hist_p_data[j]->Fill(-100.);
	  }
	}
      }

      if(selected && learnedVald<0.1 && muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    hist_p_data[j]->Fill(preco);
	  }
	}
	totalevents++;
      } else if(muMinus) {
	for (int j = 0; j < numThetaRanges; ++j) {
	  if (threco > thetaRanges[j][0] && threco < thetaRanges[j][1]) {
	    hist_p_data[j]->Fill(100.);
	  }
	}
      }


    }
  }

  cout<<"Data TotalEvents Reco good"<<totalevents<<endl;


  for (int j = 0; j < numThetaRanges; ++j) {
    hist_p_data[j]->Scale(0.8*2620000./T2->GetEntries());
  }

  //Get Fakes and effi, RM, gen distribution from response matrix rootfile
  TFile* infileRM;
  if(!muMinus){
    infileRM = TFile::Open(TString::Format("response_matrix_new/fileOut_allmodels_chisq/momentum_response_matrix_v%d_muplus_ndof%d.root", 100, myndof),"READ");
  }
  else{
    infileRM = TFile::Open(TString::Format("response_matrix_new/fileOut_allmodels_chisq/momentum_response_matrix_v%d_muminus_ndof%d.root", 100, myndof),"READ");
  }
  
  
  TH1D *hist_p_gen[numThetaRanges];
  TH2D *mat_p_rm[numThetaRanges];
  TH1D* hFakeRate[numThetaRanges];
  TH1D* hIneffi[numThetaRanges];
  TH1D* hEffi[numThetaRanges];

  infileRM->cd();
  for(int ij=0;ij<numThetaRanges;ij++){
    cout<<ij<<endl;
    if(!muMinus){
      sprintf(name,"FakeEvents_muplus_theta_%d_%d", (int)thetaRanges[ij][0], (int)thetaRanges[ij][1]);  
      hFakeRate[ij]= (TH1D*)infileRM->Get(name);
      sprintf(name,"InEffi_Events_muplus_theta_%d_%d", (int)thetaRanges[ij][0], (int)thetaRanges[ij][1]);
      hIneffi[ij] = (TH1D*)infileRM->Get(name);
      sprintf(name,"Effi_Events_muplus_theta_%d_%d", (int)thetaRanges[ij][0], (int)thetaRanges[ij][1]); 
      hEffi[ij] = (TH1D*)infileRM->Get(name);
    }
    else{
      sprintf(name,"FakeEvents_muminus_theta_%d_%d", (int)thetaRanges[ij][0], (int)thetaRanges[ij][1]); 
      hFakeRate[ij]= (TH1D*)infileRM->Get(name);
      sprintf(name,"InEffi_Events_muminus_theta_%d_%d", (int)thetaRanges[ij][0], (int)thetaRanges[ij][1]);
      hIneffi[ij] = (TH1D*)infileRM->Get(name);
      sprintf(name,"Effi_Events_muminus_theta_%d_%d", (int)thetaRanges[ij][0], (int)thetaRanges[ij][1]);
      hEffi[ij] = (TH1D*)infileRM->Get(name);
    }
  }

  
  cout<<"Clone response, gen from root file"<<endl;
   for (int j = 0; j < numThetaRanges; ++j) {
     cout<<j<<endl;
     if(!muMinus) {
       sprintf(name,"h_p_gen_muplus_v%d_theta_%d_%d", 100, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
       hist_p_gen[j] = (TH1D*)infileRM->Get(name);
     }
     if(muMinus)  {
       sprintf(name,"h_p_gen_muminus_v%d_theta_%d_%d", 100, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
        hist_p_gen[j] = (TH1D*)infileRM->Get(name);
     }

     cout<<"name "<<j<<" "<<name<<endl;

     if(!muMinus) {sprintf(name,"h_p_gen_reco_muplus_v%d_theta_%d_%d", 100, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
       
       mat_p_rm[j] = (TH2D*)infileRM->Get(name);
     }
     if(muMinus)  {sprintf(name,"h_p_gen_reco_muminus_v%d_theta_%d_%d", 100, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
       mat_p_rm[j] = (TH2D*)infileRM->Get(name);}

     cout<<"name: "<<name<<endl;

     //Scale Only when using all models:   
     mat_p_rm[j]->Scale(1./8.);
     hist_p_gen[j]->Scale(1./8.);

   }
   

   cout<<"Create new rootfile to save unfolded distribution"<<endl;
   
  TFile *fileout ;
  if(!muMinus)sprintf(name,"%s/TUnfold_momentum_muplus_v%d.root", saveDir.c_str(), version);
  if(muMinus)sprintf(name,"%s/TUnfold_momentum_muminus_v%d.root", saveDir.c_str(), version);
  fileout = new TFile(name,"recreate") ;
  cout<<"fileout created"<<endl;

  //Momentum
  TH1D *hist_unf_reco[numThetaRanges];
  for (int j = 0; j < numThetaRanges; ++j) {
	if(!muMinus) {sprintf(name,"h_unfolded_mc_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
	hist_unf_reco[j] = new TH1D(name,name,npbins_gen,gen_bin_sch); hist_unf_reco[j]->Sumw2();}
	if(muMinus)  {sprintf(name,"h_unfolded_mc_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
	hist_unf_reco[j] = new TH1D(name,name,npbins_gen,gen_bin_sch); hist_unf_reco[j]->Sumw2();}
  }

  TH1D *hist_unf_data[numThetaRanges];
  for (int j = 0; j < numThetaRanges; ++j) {
	if(!muMinus) {sprintf(name,"h_unfolded_data_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
	hist_unf_data[j] = new TH1D(name,name,npbins_gen,gen_bin_sch); hist_unf_data[j]->Sumw2();}
	if(muMinus)  {sprintf(name,"h_unfolded_data_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
	hist_unf_data[j] = new TH1D(name,name,npbins_gen,gen_bin_sch); hist_unf_data[j]->Sumw2();}
  }

  
    if(!muMinus){
      sprintf(name,"Lcurve_mc_reco_v%d_muplus",version); 
    }
    else{
      sprintf(name,"Lcurve_mc_reco_v%d_muminus",version); 
    }

    cout<<"check1"<<endl;
    TCanvas* canvasm = new TCanvas(name,name,800,1200);                
    canvasm->Divide(3,5,1.e-6,1.e-6);
    gStyle->SetLabelSize(0.06,"XYZ");
    gStyle->SetTitleSize(0.06,"XYZ");

    if(!muMinus){
      sprintf(name,"Lcurve_data_reco_v%d_muplus",version); 
    }
    else{
      sprintf(name,"Lcurve_data_reco_v%d_muminus",version); 
    }


    TCanvas* canvasd = new TCanvas(name,name,800,1200);
    canvasd->Divide(3,5,1.e-6,1.e-6);
    gStyle->SetLabelSize(0.06,"XYZ");
    gStyle->SetTitleSize(0.06,"XYZ");

    cout<<"check2"<<endl;

    //RSA
    TGraph *LCurve[numThetaRanges];
    TGraph *BestLCurve[numThetaRanges];
  
    
    for (int j = 0; j < numThetaRanges; ++j) {
      cout<<"check3 "<<j<<endl;

      // regularize curvature                                                                                
  TUnfold::ERegMode regMode=TUnfold::kRegModeDerivative;
  // preserve the area                                                                                   
  TUnfold::EConstraint constraintMode=TUnfold::kEConstraintArea;
  //TUnfold::EConstraint constraintMode=TUnfold::kEConstraintNone;                                       
  // bin content is divided by the bin width                                                             
  TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;

  TUnfoldDensity unfoldBbBreco(mat_p_rm[j],TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
  TUnfoldDensity unfoldBbBdata(mat_p_rm[j],TUnfold::kHistMapOutputVert, regMode, constraintMode, densityFlags);
    unfoldBbBreco.SetInput(hist_p_reco[j],1.0);
    unfoldBbBdata.SetInput(hist_p_data[j],1.0);	

    
    //Regulatisation
    
    Double_t estimatedPeakPosition=1.6;//1.2
    if(muMinus) estimatedPeakPosition = -1.6;
    Int_t nPeek=1;
    Int_t iPeek=(Int_t)(npbins_gen*(estimatedPeakPosition-p_low)/(p_high-p_low)
			// offset 1.5
			// accounts for start bin 1
			// and rounding errors +0.5
			+1.5);
    iPeek = hist_p_gen[j]->GetMaximumBin();
    //if(!muminus) iPeek=6;
    //else iPeek = 4;
    // regularize output bins 1..iPeek-nPeek
    ////unfoldBbB.RegularizeBins(iPeek+10,1,npbins_gen,regMode);

    //Assuming we are putting this condition in gen distribution avoid peak bins
    if(!muMinus){unfoldBbBreco.RegularizeBins(3,1,9, regMode);}
    else{unfoldBbBreco.RegularizeBins(3,1,9, regMode);}

    //RSA
    //unfoldBbBreco.RegularizeBins(1,1,npbins_gen,regMode);
    //unfoldBbBreco.RegularizeBins(1,1,iPeek-nPeek,regMode);
    //unfoldBbBreco.RegularizeBins(iPeek+nPeek,1,npbins_gen-(iPeek+nPeek),regMode);


    //RSA
    if(!muMinus){unfoldBbBdata.RegularizeBins(3,1,9, regMode);}
    else{unfoldBbBdata.RegularizeBins(3,1,9, regMode);}
    
    //iPeek = hist_p_data[j]->GetMaximumBin();
    //unfoldBbBdata.RegularizeBins(1,1,npbins_gen,regMode);
    //unfoldBbBdata.RegularizeBins(1,1,iPeek-nPeek,regMode);
    //unfoldBbBdata.RegularizeBins(iPeek+nPeek,1,npbins_gen-(iPeek+nPeek),regMode);

    
    // unfoldBbBreco.DoUnfold(pow(10,-1));
    // unfoldBbBdata.DoUnfold(pow(10,-1));

    // hist_unf_reco[j] = (TH1D*)unfoldBbBreco.GetOutput("Unfolded");
    // hist_unf_data[j] = (TH1D*)unfoldBbBdata.GetOutput("Unfolded");
    
    Double_t tauMin=pow(10,-1);//-1
    Double_t tauMax=pow(10,3);//3

    //Double_t tauMin=0.0;//-1
    //Double_t tauMax=0.0;//3
    
    // if(muMinus){
    //   tauMin= pow(10,-1);
    //   tauMax = pow(10,0.);
    // }
    // else{
    //   tauMin= pow(10,-1.);
    //   tauMax = pow(10,0.);      
    // }
    
    Int_t nScan=50;
    Int_t iBest;
    TSpline *logTauX,*logTauY;
    TSpline *logTauCurvature;
    TGraph *lCurvem;
    TGraph *lCurved;

    //naming
    char title1[300];
    char title2[300];
    if(!muMinus){
      sprintf(title1,"h_unfolded_mc_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      sprintf(title2,"h_unfolded_data_muplus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
    }
    if(muMinus){
      sprintf(title1,"h_unfolded_mc_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
      sprintf(title2,"h_unfolded_data_muminus_v%d_theta_%d_%d",version, (int)thetaRanges[j][0], (int)thetaRanges[j][1]);
    }

    iBest = unfoldBbBreco.ScanLcurve(nScan,tauMin,tauMax,&lCurvem,&logTauX,&logTauY, &logTauCurvature);
    std::cout<<"chi**2_reco="<<unfoldBbBreco.GetChi2A()<<"+"<<unfoldBbBreco.GetChi2L()<<" / "<<unfoldBbBreco.GetNdf()<<"\n";
    Double_t t[1],x[1],y[1];
    logTauX->GetKnot(iBest,t[0],x[0]);
    logTauY->GetKnot(iBest,t[0],y[0]);
    TGraph *bestLcurvem=new TGraph(1,x,y);
    TGraph *bestLogTauLogChi2m=new TGraph(1,t,x);
    
    cout<<"tau = "<<t[0] <<"\tx = "<<x[0]<<"\ty = "<<y[0]<<endl;
    hist_unf_reco[j] = (TH1D*)unfoldBbBreco.GetOutput(title1);
    cout<<"check1"<<endl;

    
    cout<<"check2"<<endl;
    
    //plot RSA
    canvasm->cd(3*j+1);
    logTauX->Draw();
    bestLogTauLogChi2m->SetMarkerColor(kRed);
    bestLogTauLogChi2m->Draw("*");
    // show the L curve
    canvasm->cd(3*j+2);
    lCurvem->Draw("AL");
    bestLcurvem->SetMarkerColor(kRed);
    bestLcurvem->Draw("*");   
    //logtaucurvature
    canvasm->cd(3*j+3);
    logTauCurvature->Draw();

    
    if(j==3){
      sprintf(name,"%s/%s.png",saveDir.c_str(),canvasm->GetName());
      canvasm->SaveAs(name);
      fileout->cd();
      canvasm->Write();
    }

    //    delete bestLcurvem;
    //    delete bestLogTauLogChi2m;

    cout<<"check3"<<endl;

    
    iBest = unfoldBbBdata.ScanLcurve(nScan,tauMin,tauMax,&lCurved,&logTauX,&logTauY,&logTauCurvature);
    std::cout<<"chi**2_reco="<<unfoldBbBdata.GetChi2A()<<"+"<<unfoldBbBdata.GetChi2L()<<" / "<<unfoldBbBdata.GetNdf()<<"\n";
    logTauX->GetKnot(iBest,t[0],x[0]);
    logTauY->GetKnot(iBest,t[0],y[0]);
    TGraph *bestLcurved=new TGraph(1,x,y);
    TGraph *bestLogTauLogChi2d=new TGraph(1,t,x);

    
    cout<<"tau = "<<t[0] <<"\tx = "<<x[0]<<"\ty = "<<y[0]<<endl;
    hist_unf_data[j] = (TH1D*)unfoldBbBdata.GetOutput(title2);

    canvasd->cd(3*j+1);
    logTauX->Draw();
    bestLogTauLogChi2d->SetMarkerColor(kRed);
    bestLogTauLogChi2d->Draw("*");


    // show the L curve
    canvasd->cd(3*j+2);
    
    lCurved->Draw("AL");

    bestLcurved->SetMarkerColor(kRed);
    bestLcurved->Draw("*");

    LCurve[j] = (TGraph*)lCurved->Clone();
    BestLCurve[j] = (TGraph*)bestLcurved->Clone();

    canvasd->cd(3*j+3);
    logTauCurvature->Draw();

    if(j==4){sprintf(name,"%s/%s.png",saveDir.c_str(),canvasd->GetName());
    canvasd->SaveAs(name);
    fileout->cd();
    canvasd->Write();

    }
    //    delete bestLcurved;
    //    delete bestLogTauLogChi2d;

    cout<<"check13"<<endl;

    
	
    hist_p_data[j]->SetLineColor(kBlack);
    hist_p_reco[j]->SetLineColor(kBlue);
    hist_p_gen[j]->SetLineColor(kRed);
    hist_unf_reco[j]->SetLineColor(kMagenta);
    hist_unf_data[j]->SetLineColor(kGreen);

    hist_p_reco[j]->SetStats(0);
    hist_p_reco[j]->SetStats(0);
    hist_p_gen[j]->SetStats(1);
    hist_unf_reco[j]->SetStats(0);
    hist_unf_data[j]->SetStats(0);
    hist_p_gen[j]->SetLineStyle(kDashed);
  
    DivideHistogramByBinWidth(hist_p_data[j]);
    DivideHistogramByBinWidth(hist_p_reco[j]);
    DivideHistogramByBinWidth(hist_p_gen[j]);
    DivideHistogramByBinWidth(hist_unf_reco[j]);
    DivideHistogramByBinWidth(hist_unf_data[j]);

    // unfoldBbBreco.Delete();
    // unfoldBbBdata.Delete();
    
    }
  //Momentum
  if(!muMinus){sprintf(name,"Tunfold_output_mc_reco_momentum_muplus_v%d",version);}
  else{
    sprintf(name,"Tunfold_output_mc_reco_momentum_muminus_v%d",version);
  }
  TCanvas *c1 = new TCanvas(name,name,800,1200);
  c1->cd();
  c1->Divide(2,3);

  for (int j = 0; j < numThetaRanges; ++j) {
    TLegend *leg1 = new TLegend(0.4,0.2,0.6,0.4);
    leg1->SetBorderSize(0);
    c1->cd(j+1);

    //gPad->SetLogx(1);
    gPad->SetLogy(1);

    hist_p_gen[j]->GetYaxis()->SetRangeUser(1000,1000000);
    //hist_p_gen[j]->SetTitle("MC Reco Unfolding Without ML");

    hist_p_gen[j]->Draw("hist");  leg1->AddEntry(hist_p_gen[j],"GEN","l");
    hist_p_data[j]->Draw("hist:same");  leg1->AddEntry(hist_p_data[j],"Data Reco","l");
    hist_p_reco[j]->Draw("hist:sames"); leg1->AddEntry(hist_p_reco[j],"MC Reco","l");
    hist_unf_reco[j]->Draw("sames"); leg1->AddEntry(hist_unf_reco[j],"Unfolded","l");

    leg1->Draw();
		
  }

  sprintf(name,"%s/%s.png", saveDir.c_str(), c1->GetName());
  c1->SaveAs(name);
  
  if(!muMinus){
    sprintf(name,"Tunfold_output_data_reco_momentum_muplus_v%d",version);
  }else{
    sprintf(name,"Tunfold_output_data_reco_momentum_muminus_v%d",version);
  }
  
  TCanvas *c2 = new TCanvas(name,name,800,1200);
  c2->cd();
  c2->Divide(2,3);
  for (int j = 0; j < numThetaRanges; ++j) {
    TLegend *leg2 = new TLegend(0.4,0.2,0.6,0.4);
    leg2->SetBorderSize(0);
    c2->cd(j+1);
    gPad->SetLogy(1);

    hist_p_gen[j]->GetYaxis()->SetRangeUser(1000,1000000);
    //hist_p_gen[j]->SetTitle("Data Reco Unfolding Without ML");

    hist_p_gen[j]->Draw("hist");  leg2->AddEntry(hist_p_gen[j],"GEN","l");
    hist_p_data[j]->Draw("hist:same");  leg2->AddEntry(hist_p_data[j],"Data Reco","l");
    hist_p_reco[j]->Draw("hist:sames"); leg2->AddEntry(hist_p_reco[j],"MC Reco","l");
    hist_unf_data[j]->Draw("sames"); leg2->AddEntry(hist_unf_data[j],"Unfolded","l");

    leg2->Draw();
  }

  sprintf(name,"%s/%s.png",saveDir.c_str(),c2->GetName());
  c2->SaveAs(name);



  for (int j = 0; j < numThetaRanges; ++j) {
    hist_p_gen[j]->Write();
    hist_p_reco[j]->Write();
    hist_p_data[j]->Write();
    hist_unf_reco[j]->Write();
    hist_unf_data[j]->Write();
    mat_p_rm[j]->Write();

    hFakeRate[j]->Write();
    hIneffi[j]->Write();
    hEffi[j]->Write();
    
    LCurve[j]->Write();
    BestLCurve[j]->Write();

  }
  fileout->Write();
  fileout->Close();


    }


/*
- 0.8 Scalig 
(base) [rajshah@rpc2 dir_corr]$ root -l Mc_v107_2018_Data_Magnetic.root 
root [0] 
Attaching file Mc_v107_2018_Data_Magnetic.root as _file0...
(TFile *) 0x369dc40
root [1] T3->GetEntries()
(long long) 2680000
root [2] .q
(base) [rajshah@rpc2 dir_corr]$ root -l Mc_v108_2018_Data_Magnetic.root 
root [0] 
Attaching file Mc_v108_2018_Data_Magnetic.root as _file0...
(TFile *) 0x2db5630
root [1] T3->GetEntries()
(long long) 2740000
root [2] .q
(base) [rajshah@rpc2 dir_corr]$ root -l Mc_v110_2018_Data_Magnetic.root 
root [0] 
Attaching file Mc_v110_2018_Data_Magnetic.root as _file0...
(TFile *) 0x3016b40
root [1] T3->GetEntries()
(long long) 2620000
root [2] .q
(base) [rajshah@rpc2 dir_corr]$ root -l Mc_v111_2018_Data_Magnetic.root 
root [0] 
Attaching file Mc_v111_2018_Data_Magnetic.root as _file0...
(TFile *) 0x3c9c840
root [1] T3->GetEntries()
(long long) 2710000
root [2] .q
(base) [rajshah@rpc2 dir_corr]$ root -l Mc_v112_2018_Data_Magnetic.root 
root [0] 
Attaching file Mc_v112_2018_Data_Magnetic.root as _file0...
(TFile *) 0x39d9420
root [1] T3->GetEntries()
(long long) 2800000
root [2] .q
(base) [rajshah@rpc2 dir_corr]$ root -l Mc_v113_2018_Data_Magnetic.root 
root [0] 
Attaching file Mc_v113_2018_Data_Magnetic.root as _file0...
(TFile *) 0x287c400
root [1] T3->GetEntries()
(long long) 2780000
root [2] .q
(base) [rajshah@rpc2 dir_corr]$ root -l Mc_v114_2018_Data_Magnetic.root 
root [0] 
Attaching file Mc_v114_2018_Data_Magnetic.root as _file0...
(TFile *) 0x1eaf1e0
root [1] T3->GetEntries()
(long long) 2870000
root [2] .q



 */
