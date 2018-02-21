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

// Function for producing wiggle plot:
// All units are in ns

Double_t Total_Time = 30.0 * 4200.0;

Double_t BinWidth = 150;

//int rebinFactor = 1;

Double_t myfunction(Double_t *x, Double_t *par){

  Float_t time     = x[0];
  double tau       = 2200;              // Lifetime of the muon at rest (ns)
  double gamma     = 29.3;              // Magic gamma
  double A         = 0.05;              // Amplitude
  double omega     = 2 * M_PI / 4200;   // Time of a single wiggle ~ 4200 ns (omega=2pi/t)
  double phase     = M_PI/2;
  double tau_gamma = tau * gamma;       //Making life easier

  //Double_t Npositrons =  N * exp(- time / tau_gamma ) * (1 + A * cos ((omega * time) + phase));
  //          N(t)    =  N(0) exp (-t/(tau * gamma )) * (1 + A * cos(omega * t + phase))
  //return Npositrons;

  double Ntot      = 2E11;               // Total number of events 
  double N         = Ntot *  ( (1 - exp( - (1 / tau_gamma) * BinWidth )) / (1 - exp( - (1 / tau_gamma) * Total_Time)) );
  
  Double_t Npositrons =  N * exp(- time / tau_gamma ) * (1 + A * cos ((omega * time) + phase));
  //          N(t)    =  N(0) exp (-t/(tau * gamma )) * (1 + A * cos(omega * t + phase))
  
  return Npositrons;
  
}

int main(){

  Double_t w = 1500;
  Double_t h = 1000;
  TCanvas* c2 = new TCanvas("c2", "c2", w, h);
  
  //  Double_t Total_Time = 30.0 * 4200.0;

  TF1 *f1 = new TF1("myfunc", myfunction, 0, Total_Time, 0); // Book histogram for ideal wiggle plot

  // Loop over full time of experiment
  
  Double_t par[] = {};
  
  //  Double_t BinWidth = 150;
  
  Double_t Nfunc[int(Total_Time / BinWidth)];
 
  for (int ibin = 0; ibin <= (Total_Time / BinWidth); ibin++){

    double time = BinWidth * ibin;
    Nfunc[ibin] = myfunction(&time, par); // Time distribution for decay positrons

    //double timeLow  = time - (BinWidth/2.0);
    //double timeHigh = time + (BinWidth/2.0);

    //vector<double> values;
    //double integral = 0.0;

    //// per ns
    //int ndiv = 1;
    //double fracTime = BinWidth / double(ndiv);
    //for (int i = 0; i <= ndiv; i++){
    //  double tmpTime = timeLow + double(i)*fracTime;
    //  values.push_back(tmpTime);
    //  integral += myfunction(&tmpTime, par);
    //  if (ibin == 0) {
    //	cout << "tmpTime: " << tmpTime << "\n";
    //  }
    //}
    //
    //integral = integral / double(values.size());//f1->Integral(timeLow, timeHigh);
    //
    //if (ibin == 0) {
    //  cout << "time: " << time << " Nfunc: " << Nfunc[ibin] << "\n";
    //  cout << "timeLow: " << timeLow << " timeHigh: " << timeHigh << ", integral: " << integral << "\n";
    //}

    // if we want to use average across bin...
    //Nfunc[ibin] = integral;
    
  }
    
  TRandom3 *r = new TRandom3(1); // Set seed to get same random numbers (0 means no seed)

  // Loop over number of psuedo experiments (choosen by the user)

  int N_Pseudo_Exp;
  
  cout << "Please enter the number of pseudo experiements: " << endl;
  cin >>  N_Pseudo_Exp;   

  TFile* file = new TFile("PseudoExp.root", "RECREATE"); // Create a ROOT file containing all pseudo experiments
  
  for (int i(0); i <= (N_Pseudo_Exp); i++){

    // Using stringstream to accommodate the number of pseudo experiements when naming and giving titles to plots
    
    string Pseudo_num;
    ostringstream conv;
    conv << i;
    string t1name  = "PseudoA_" + conv.str();
    string t1title = "PseudoA Wiggle Plot " + conv.str();
    string t2name  = "PseudoB_" + conv.str();
    string t2title = "PseudoB Wiggle Plot " + conv.str();
    string t3name  = "PseudoC_" + conv.str();
    string t3title = "PseudoC Wiggle Plot " + conv.str();
    string t4name  = "PseudoD_" + conv.str();
    string t4title = "PseudoD Wiggle Plot " + conv.str();
    string tTotname  = "PseudoTot_" + conv.str();
    string tTottitle = "PseudoTot Wiggle Plot " + conv.str();
    
    string tDiffname  = "Pseudo_diff_" + conv.str();
    string tDifftitle = "Wiggle Difference Plot " + conv.str();

    // Book a histogram for pseudo wiggle plot

    TH1F* t1 = new TH1F("t1", "", int(Total_Time / BinWidth), -(BinWidth/2), Total_Time - (BinWidth/2));
    TH1F* t2 = new TH1F("t2", "", int(Total_Time / BinWidth), -(BinWidth/2), Total_Time - (BinWidth/2));
    TH1F* t3 = new TH1F("t3", "", int(Total_Time / BinWidth), -(BinWidth/2), Total_Time - (BinWidth/2));
    TH1F* t4 = new TH1F("t4", "", int(Total_Time / BinWidth), -(BinWidth/2), Total_Time - (BinWidth/2));
    TH1F* tTot = new TH1F("tT0", "", int(Total_Time / BinWidth), -(BinWidth/2), Total_Time - (BinWidth/2));
    
    //TH1F* t1 = new TH1F("t1", "", int(Total_Time / BinWidth), -(BinWidth*rebinFactor /2), Total_Time - (BinWidth* rebinFactor/2));
    //TH1F* t2 = new TH1F("t2", "", int(Total_Time / BinWidth), -(BinWidth*rebinFactor /2), Total_Time - (BinWidth* rebinFactor/2));
    //TH1F* t3 = new TH1F("t3", "", int(Total_Time / BinWidth), -(BinWidth*rebinFactor /2), Total_Time - (BinWidth* rebinFactor/2));
    //TH1F* t4 = new TH1F("t4", "", int(Total_Time / BinWidth), -(BinWidth*rebinFactor /2), Total_Time - (BinWidth* rebinFactor/2));
    //TH1F* tTot = new TH1F("tT0", "", int(Total_Time / BinWidth), -(BinWidth*rebinFactor /2), Total_Time - (BinWidth * rebinFactor/2));
    
    t1 -> SetNameTitle(t1name.c_str(), t1title.c_str());
    t1 -> GetXaxis() -> SetTitle("Time (ns)");
    t1 -> GetYaxis() -> SetTitle("Number of Positrons N");

    t2 -> SetNameTitle(t2name.c_str(), t2title.c_str());
    t2 -> GetXaxis() -> SetTitle("Time (ns)");
    t2 -> GetYaxis() -> SetTitle("Number of Positrons N");

    t3 -> SetNameTitle(t3name.c_str(), t3title.c_str());
    t3 -> GetXaxis() -> SetTitle("Time (ns)");
    t3 -> GetYaxis() -> SetTitle("Number of Positrons N");

    t4 -> SetNameTitle(t4name.c_str(), t4title.c_str());
    t4 -> GetXaxis() -> SetTitle("Time (ns)");
    t4 -> GetYaxis() -> SetTitle("Number of Positrons N");

    tTot -> SetNameTitle(tTotname.c_str(), tTottitle.c_str());
    tTot -> GetXaxis() -> SetTitle("Time (ns)");
    tTot -> GetYaxis() -> SetTitle("Number of Positrons N");

    // Book a histogram for difference between ideal wiggle and pseudo wiggle plots
    
    TH1F* tDiff = new TH1F("", "", 100, -10, 10);
    
    tDiff -> SetNameTitle(tDiffname.c_str(), tDifftitle.c_str());
    tDiff -> GetXaxis() -> SetTitle("N/#sqrt{N}");
    tDiff -> GetYaxis() -> SetTitle("Number of Entries");    
    
    for (int ibin = 0; ibin <= (Total_Time / BinWidth); ibin++){

      Double_t PseudoWiggle1 = r -> Gaus(Nfunc[ibin]/4.0 , sqrt( Nfunc[ibin]/4.0 ));  // Gaussian smearing
      Double_t PseudoWiggle2 = r -> Gaus(Nfunc[ibin]/4.0 , sqrt( Nfunc[ibin]/4.0 ));  // Gaussian smearing
      Double_t PseudoWiggle3 = r -> Gaus(Nfunc[ibin]/4.0 , sqrt( Nfunc[ibin]/4.0 ));  // Gaussian smearing
      Double_t PseudoWiggle4 = r -> Gaus(Nfunc[ibin]/4.0 , sqrt( Nfunc[ibin]/4.0 ));  // Gaussian smearing
      
      Double_t PseudoWiggle = PseudoWiggle1 + PseudoWiggle2 + PseudoWiggle3 + PseudoWiggle4;
      t1 -> SetBinContent(ibin+1, PseudoWiggle1);
      t2 -> SetBinContent(ibin+1, PseudoWiggle2);
      t3 -> SetBinContent(ibin+1, PseudoWiggle3);
      t4 -> SetBinContent(ibin+1, PseudoWiggle4);
      
      tTot -> SetBinContent(ibin+1, PseudoWiggle);

      //if (ibin == 100 ) {
      //	cout << "p1: " << PseudoWiggle1 << " p2: " << PseudoWiggle2 << "p3: " << PseudoWiggle3 << " p4: " << PseudoWiggle4 << "\n";
      //	cout << "and p1: " << t1->GetBinContent(ibin+1) << " p2: " << t2->GetBinContent(ibin+1) << "p3: " << t3->GetBinContent(ibin+1) << " p4: " << t4->GetBinContent(ibin+1) << "\n";	
      //}
      
      Double_t diff =  (Nfunc[ibin] - PseudoWiggle) / (sqrt(Nfunc[ibin])); // Difference between ideal and pseudo wiggle
      tDiff -> Fill(diff);
      
    }

    //try rebinning
    //    t1->Rebin(rebinFactor);
    //    t2->Rebin(rebinFactor);
    //    t3->Rebin(rebinFactor);
    //    t4->Rebin(rebinFactor);
    //    tTot->Rebin(rebinFactor);
    //tDiff->Rebin(rebinFactor);
    
    // Write TH1F to ROOT file
    t1 -> Write();
    t2 -> Write();
    t3 -> Write();
    t4 -> Write();
    tTot -> Write();
    tDiff -> Write();

  }
  
  // Draw ideal wiggle plot
  
  f1 -> SetTitle("Ideal Wiggle Plot");
  f1 -> GetXaxis() -> SetTitle("Time (ns)");
  f1 -> GetYaxis() -> SetTitle("Number of Positrons N");
  //f1 -> GetYaxis() -> SetRangeUser(0, 2E11* 1.1);
  f1 -> Draw();
  c2 -> SaveAs("IdealFit.eps");
  
}
