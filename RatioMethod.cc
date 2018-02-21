#include <iostream>
#include <fstream>
#include <iomanip>
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom3.h"
#include "TVirtualFFT.h"
#include "TFitter.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include <TMatrixDSym.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

// Function which does Ratio Method and produces R(t), D(t), S(t):

//TH1F* GetRatioMethodHistFromPseudoExp(TH1F* PseudoExp, double Total_Time, double BinWidth, int BinNum, double binshift TH1F* RatDist_hist, TH1F* Ddist_hist, TH1F* Sdist_hist){

vector<TH1F*> GetRatioMethodHistFromPseudoExp(TH1F* PseudoExpA, TH1F* PseudoExpB, TH1F* PseudoExpC, TH1F* PseudoExpD, double Total_Time, double BinWidth, int BinNum, int choice = 0){
  
  // For now shifting by 1 bin width:
  // Float_t Toff    = BinWidth;           //cross check
  
  Float_t Toff    = 2100.0;   // Toff = pi/omega -> omega = 2pi/4200 -> Toff would be 2100 ns
  
  Double_t Ddist[int(Total_Time/BinWidth)];
  Double_t Sdist[int(Total_Time/BinWidth)];
  
  int binshift = Toff / BinWidth;

  //cout << "bin shift= " << binshift << " width " << BinWidth << " BinNum " << BinNum << " total time " << Total_Time << endl;
  
  // Book histograms for R(t), D(t), S(t)
  
  TH1F* RatDist_hist = new TH1F("" , "" , BinNum - (2 * binshift) , Toff - (BinWidth/2.0) , Total_Time - Toff - (BinWidth/2.0));
  TH1F* Ddist_hist   = new TH1F("" , "" , BinNum - (2 * binshift) , Toff - (BinWidth/2.0) , Total_Time - Toff - (BinWidth/2.0));
  TH1F* Sdist_hist   = new TH1F("" , "" , BinNum - (2 * binshift) , Toff - (BinWidth/2.0) , Total_Time - Toff - (BinWidth/2.0));

  // without half bin shift
  //TH1F* RatDist_hist = new TH1F("" , "" , BinNum - (2 * binshift) , Toff , Total_Time - (Toff) );
  //TH1F* Ddist_hist   = new TH1F("" , "" , BinNum - (2 * binshift) , Toff , Total_Time - (Toff) );
  //TH1F* Sdist_hist   = new TH1F("" , "" , BinNum - (2 * binshift) , Toff , Total_Time - (Toff) );
  
  // Loop over bins of pseudo exp
  
  //for (int ibin(BinWidth/2.0 + binshift); ibin <= (BinNum  - binshift); ibin++ ){
  for (int ibin(0 + binshift); ibin <= (BinNum  - binshift); ibin++ ){
    
    int ibin_up = ibin + binshift;
    int ibin_do = ibin - binshift;
    
    // Divide bin content by 4, apply shifts and add to correct array
    double d1 = PseudoExpA -> GetBinContent(ibin) ;
    double d2 = PseudoExpB -> GetBinContent(ibin) ;   
    double d3 = PseudoExpC -> GetBinContent(ibin_up) ;
    double d4 = PseudoExpD -> GetBinContent(ibin_do) ;

    double s1 = PseudoExpA -> GetBinContent(ibin) ;
    double s2 = PseudoExpB -> GetBinContent(ibin) ;   
    double s3 = PseudoExpC -> GetBinContent(ibin_up) ;
    double s4 = PseudoExpD -> GetBinContent(ibin_do) ;
    
    double Ddist = d3 + d4 - d1 - d2;
    double Sdist = s1 + s2 + s3 + s4;
        
    RatDist_hist -> SetBinContent(ibin-binshift, Ddist/Sdist);
    Ddist_hist   -> SetBinContent(ibin-binshift, Ddist);
    Sdist_hist   -> SetBinContent(ibin-binshift, Sdist);

    // Find error on R(t), D(t), S(t)

    double RatioValue = RatDist_hist -> GetBinContent(ibin);
    double RatioError = sqrt( (1 - pow(RatioValue,2)) / (s1+s2+s3+s4) );

    double DValue = Ddist_hist -> GetBinContent(ibin);
    double DError = sqrt(d3 + d4 + d1 + d2);
    
    double SValue = Sdist_hist -> GetBinContent(ibin);
    double SError = sqrt(s1 + s2 + s3 + s4);

    // Set error on R(t), D(t), S(t)
    
    RatDist_hist -> SetBinError(ibin-binshift, RatioError);
    //RatDist_hist -> SetBinError(ibin, RatioError);
    gStyle       -> SetEndErrorSize(2);
    RatDist_hist -> SetMarkerStyle(20);
    RatDist_hist -> SetMarkerSize(0.5);
    RatDist_hist -> SetMarkerColor(kBlack);
    RatDist_hist -> Draw("E1");

    Ddist_hist -> SetBinError(ibin-binshift , DError);
    gStyle     -> SetEndErrorSize(2);
    Ddist_hist -> SetMarkerStyle(20);
    Ddist_hist -> SetMarkerSize(0.5);
    Ddist_hist -> SetMarkerColor(kBlack);
    Ddist_hist -> Draw("E1");

    Sdist_hist -> SetBinError(ibin-binshift , SError);
    gStyle     -> SetEndErrorSize(2);
    Sdist_hist -> SetMarkerStyle(1);
    Sdist_hist -> SetMarkerSize(0.5);
    Sdist_hist -> SetMarkerColor(kBlack);
    Sdist_hist -> Draw("E1");
  }
  
  vector<TH1F*> tmp = {RatDist_hist, Ddist_hist, Sdist_hist};

  return tmp;

}

int main(){
  
  Double_t w = 1500;
  Double_t h = 1000;
  TCanvas* c2 = new TCanvas("c2", "c2", w, h);
  
  // Reading ROOT file and extracting histograms (choosen by user)
  
  int N_Pseudo_Exp;
  
  cout << "Please enter the number of histograms you want to use for the Ratio Method: " << endl;
  cin >>  N_Pseudo_Exp;
  
  TFile *file = TFile::Open("PseudoExp.root", "READ");
  
  if (file == 0) {
    
    // if we cannot open the file, print an error message and return immediatly
    printf("Error: cannot open PseudoExp.root");
    return 1;
    
  }

  TFile* file1 = new TFile("RatioMethodData.root", "RECREATE"); // Create a ROOT file containg Ratio Method data

  for (int i(0); i < (N_Pseudo_Exp); i++){                      // Loop over number of pseudo exp

    // Using stringstream to accommodate the number of pseudo experiements when naming and giving titles to plots
    
    string Pseudo_num;
    ostringstream conv;
    conv << i;
    string name1 = "PseudoA_" + conv.str();
    string name2 = "PseudoB_" + conv.str();
    string name3 = "PseudoC_" + conv.str();
    string name4 = "PseudoD_" + conv.str();

    string outname1 = "RatDist_Pseudo_" + conv.str();
    string outname2 = "DDist_Pseudo_" + conv.str();
    string outname3 = "SDist_Pseudo_" + conv.str();
    string outtitle1 = "R(t) for Pseudo Experiement " + conv.str();
    string outtitle2 = "D(t) for Pseudo Experiement " + conv.str();
    string outtitle3 = "S(t) for Pseudo Experiement " + conv.str();
    
    TH1F* h1 = (TH1F*)file -> Get(name1.c_str());   // Gets the pseudo wiggle plot
    TH1F* h2 = (TH1F*)file -> Get(name2.c_str());   // Gets the pseudo wiggle plot
    TH1F* h3 = (TH1F*)file -> Get(name3.c_str());   // Gets the pseudo wiggle plot
    TH1F* h4 = (TH1F*)file -> Get(name4.c_str());   // Gets the pseudo wiggle plot
    
    Double_t BinWidth = h1 -> GetBinWidth(1);      // Get bin width of pseudo wiggle plot
    
    Double_t BinNum = h1 -> GetNbinsX();           // Get bin number of pseudo wiggle plot
    
    Double_t Total_Time = BinWidth * BinNum;       // Find total time
    
    vector<TH1F*> tmp = GetRatioMethodHistFromPseudoExp(h1, h2, h3, h4, Total_Time , BinWidth , BinNum , 0);

    tmp.at(0) -> SetNameTitle(outname1.c_str(), outtitle1.c_str());
    tmp.at(0) -> GetXaxis() -> SetTitle("Time (ns)");
    tmp.at(0) -> GetYaxis() -> SetTitle("R(t)");
    

    tmp.at(1) -> SetNameTitle(outname2.c_str(), outtitle2.c_str());
    tmp.at(1) -> GetXaxis() -> SetTitle("Time (ns)");
    tmp.at(1) -> GetYaxis() -> SetTitle("D(t)");
    
    tmp.at(2) -> SetNameTitle(outname3.c_str(), outtitle3.c_str());
    tmp.at(2) -> GetXaxis() -> SetTitle("Time (ns)");
    tmp.at(2) -> GetYaxis() -> SetTitle("S(t)");

    tmp.at(0) -> Write();
    tmp.at(1) -> Write();
    tmp.at(2) -> Write();
    
  }

}
