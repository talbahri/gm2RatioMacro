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

Double_t myfunction(Double_t *x, Double_t *par){

  Float_t time     = x[0];
  double N         = par[0];
  double tau_gamma = par[1];
  double A         = par[2];
  double omega     = par[3];
  double phase     = par[4];

  Double_t Npositrons =  par[0] * exp( -x[0] / par[1] ) * (1 + par[2] * cos ((par[3] * x[0]) + par[4]));
  //N(t)    = N(0) exp (-t/(tau * gamma )) * (1 + A * cos(omega * t + phase)); 

  return Npositrons;
  
}

int main(){
  
  Double_t w = 1500;
  Double_t h = 1000;
  TCanvas* c2 = new TCanvas("c2", "c2", w, h);

  // Reading Root file and extracting histograms (choosen by user)

  int N_Pseudo_Exp;

  cout << "Please enter the number of histograms you want to use for the five parameter fit: " << endl;
  cin >>  N_Pseudo_Exp;

  // Setup correlation plots for fitted parameters
  // Setup fits for correlation plots

  //  TH2F* k1  = new TH2F("", "", 100, 9.9998e8, 1.00002e9, 100, 6.4458e4, 6.4462e4);       // N & TauGamma
  //TF1* CorrFit1 = NULL;
  
  //TH2F* k2  = new TH2F("", "", 100, 9.9998e8, 1.00002e9, 100, 0.049985, 0.050015);        // N & A
  //TF1* CorrFit2 = NULL;

  //TH2F* k3  = new TH2F("", "", 100, 9.9998e8, 1.00002e9, 100, 1.495985e-3, 1.49601e-3);    // N & Omega
  //TF1* CorrFit3 = NULL;

  //TH2F* k4  = new TH2F("", "", 100, 9.9998e8, 1.00002e9, 100, 1.5702, 1.5714);            // N & Phase
  //TF1* CorrFit4 = NULL;

  //TH2F* k5  = new TH2F("", "", 100, 6.4458e4, 6.4462e4, 100, 0.049985, 0.050015);          // TauGamma & A
  //TF1* CorrFit5 = NULL;

  //TH2F* k6  = new TH2F("", "", 100, 6.4458e4, 6.4462e4, 100, 1.495985e-3, 1.49601e-3);      // TauGamma & Omega
  //TF1* CorrFit6 = NULL;

  //TH2F* k7  = new TH2F("", "", 100, 6.4458e4, 6.4462e4, 100, 1.5702, 1.5714);              // TauGamma & Phase
  //TF1* CorrFit7 = NULL;

  //TH2F* k8  = new TH2F("", "", 100, 0.049985, 0.050015, 100, 1.495985e-3, 1.49601e-3);       // A & Omega
  //TF1* CorrFit8 = NULL;

  //TH2F* k9  = new TH2F("", "", 100, 0.049985, 0.050015, 100, 1.5702, 1.5714);               // A & Phase
  //TF1* CorrFit9 = NULL;

  //TH2F* k10 = new TH2F("", "", 100, 1.495985e-3, 1.49601e-3, 100, 1.5702, 1.5714);           // Omega & Phase
  //TF1* CorrFit10 = NULL; 

  // Setup histogram for Chi Square distribution of the 5-parameter fits

  //TH1F* cs  = new TH1F("", "", 100, 0, 2);        // Chi Square

  // Setup different histograms, 1 for each parameter

  //TH1F* g1 = new TH1F("", "", 100, -0.05, 0.05);     //N

  //TH1F* g2 = new TH1F("", "", 100, -0.05, 0.05);      //tau_gamma
 
  //TH1F* g3 = new TH1F("", "", 100, -0.05, 0.05);          //A

  //TH1F* g4 = new TH1F("", "", 100, -0.05, 0.05);      //omega

  //TH1F* g5 = new TH1F("", "", 100, -0.05, 0.05);          //phase
  
  // Setup different histograms, 1 for each fitted parameter error
  
  //TH1F* e1 = new TH1F("", "", 100, 2.830E3, 2.832E3);            //N

  //TH1F* e2 = new TH1F("", "", 100, 3.69E-1, 3.698E-1);                //tau_gamma

  //TH1F* e3 = new TH1F("", "", 100, 2.23E-6, 2.233E-6);  //A

  //TH1F* e4 = new TH1F("", "", 100, 2.495e-9, 2.52e-9);        //omega

  //TH1F* e5 = new TH1F("", "", 100, 7.955E-5, 7.975E-5);        //phase

  // for debugging
  //  vector<int> display_hist= {1,2,4};
  //vector<int> display_hist= {1,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000};

  TFile *file = TFile::Open("PseudoExp.root", "READ");
    
  if (file == 0) {
    // if we cannot open the file, print an error message and return immediatly
    printf("Error: cannot open PseudoExp.root");
    return 1;
  }

  TFile* file2 = new TFile("FiveParameterFit.root", "RECREATE");  // Create a ROOT file containg five parameter fits of the pseudo wiggle plots
  
  for (int i(0); i < (N_Pseudo_Exp); i++){                        // Loop over number of pseudo exp

    // Using stringstream to accommodate the number of pseudo experiements when naming and giving titles to plots
    
    string Pseudo_num;
    ostringstream conv;
    conv << i;
    string name = "PseudoTot_" + conv.str();
    string outname = name + ".eps";
      
    TH1F* h1 = (TH1F*)file -> Get(name.c_str());                 // Gets the pseudo wiggle plot
      
    Double_t BinWidth = h1 -> GetBinWidth(1);                    // Get bin width of pseudo wiggle plot
            
    Double_t BinNum = h1 -> GetNbinsX();                         // Get bin number of pseudo wiggle plot
      
    Double_t Total_Time = BinWidth * BinNum;                     // Find total time
      
    TF1 *fit = new TF1 ("fit", myfunction, 0, Total_Time, 5);    // Five parameter fit

    double tmp = 1.03;

    double Ntot      = 2E11;               // The number inputted into the generate wiggle code
    // 
    double N         = Ntot *  ( (1 - exp( - (1 / (2200 * 29.3)) * BinWidth)) / (1 - exp( - (1 / (2200 * 29.3)) * Total_Time)) );
    
    fit -> SetParameter(0, 0.99*N);
    fit -> SetParameter(1, 2200*29.3); // set to predicted value for now
    fit -> SetParameter(2, 0.99*0.05);
    fit -> SetParameter(3, tmp*(2*M_PI/4200));
    fit -> SetParameter(4, 0.986*M_PI/2);

    fit -> SetParName(0,"N");
    fit -> SetParName(1,"TauGamma");
    fit -> SetParName(2,"A");
    fit -> SetParName(3,"Omega");
    fit -> SetParName(4,"Phi");

    fit -> FixParameter(0, N);
    fit -> FixParameter(1, 2200*29.3); // set to predicted value for now
    fit -> FixParameter(2, 0.05);
    fit -> FixParameter(3, (2*M_PI/4200));
    fit -> FixParameter(4, M_PI/2);

    fit -> SetNpx(1000);
    fit -> SetNDF(5);
    
    fit -> SetLineColor(kMagenta);
    fit -> SetLineStyle(1);
    fit -> SetLineWidth(1);

    h1 -> Fit("fit");
    
    gStyle->SetOptFit(111);

    // Get fitter
    
    //TVirtualFitter* f = TVirtualFitter::GetFitter();

    // Get number of total parameters in the fit
    
    //const int n = f -> GetNumberTotalParameters();
    
    //if (0 == f->GetCovarianceMatrix()) {
    //std::cout << "Fit failed, we therefore have no covariance matrix." <<
    //std::endl;
    //return 1;
    //}
    
    // Get covariance matrix
    
    //TMatrixDSym cov(n, f->GetCovarianceMatrix());
    
    // Print correlation matrix (conversion from covariance matrix to correlation matrix is done on the fly)
    
    //std::cout << "Correlation Matrix of the fit is: " << std::endl;

    //std::cout << setw(15) << right << " N " << setw(20) << right << " TauGamma " << setw(10) << right << " A " << setw(15) << right << " Omega " << setw(15) << right << " Phi " << std::endl;

    //for (int i = 0; i < n; ++i) {
      
    //switch(i){
	
    //case 0:
    //cout << setw(10) << left << "N" ;
    //break;
    //case 1:
    //cout << setw(10) << left << "TauGamma" ;
    //break;
    //case 2:
    //cout << setw(10) << left << "A" ;
    //break;
    //case 3:
    //cout << setw(10) << left << "Omega" ;
    //break;
    //case 4:
    //cout << setw(10) << left << "Phi" ;
    //break;
	
    //}
      
    //for (int j = 0; j < n; ++j)
    //std::cout << " " <<
    // std::setw(12) << std::setprecision(6) << 
    //  cov(i, j) / std::sqrt(cov(i, i) * cov(j, j));
    //std::cout << std::endl;
      
    //}

    // Ideal values

    //double ideal_tau_gamma = 2200*29.3;
    //double ideal_A = 0.05;
    //double ideal_omega = 2 * M_PI / 4200;
    //double ideal_phase = M_PI/2;
    //double ideal_N = 1E9; //-1*(100000.0*ideal_tau_gamma/BinWidth) *(exp(-BinWidth/ideal_tau_gamma) - 1);

    // Fitted parameters
    
    //double fitted_N         = fit -> GetParameter(0);
    //double fitted_tau_gamma = fit -> GetParameter(1);
    //double fitted_A         = fit -> GetParameter(2);
    //double fitted_omega     = fit -> GetParameter(3);
    //double fitted_phase     = fit -> GetParameter(4);

    // Correlation plots and fits for fitted parameters

    //CorrFit1 = new TF1 ("CorrFit1", "pol1", 9.9998e8, 1.00002e9);   // N & TauGamma correlation fit
    //CorrFit1 -> SetParName(0, "C");
    //CorrFit1 -> SetParName(1, "Grad");

    //CorrFit2 = new TF1 ("CorrFit2", "pol1", 9.9998e8, 1.00002e9);   // N & A correlation fit
    //CorrFit2 -> SetParName(0, "C");
    //CorrFit2 -> SetParName(1, "Grad");

    //CorrFit3 = new TF1 ("CorrFit3", "pol1", 9.9998e8, 1.00002e9);   // N & Omega correlation fit
    //CorrFit3 -> SetParName(0, "C");
    //CorrFit3 -> SetParName(1, "Grad");

    //CorrFit4 = new TF1 ("CorrFit4", "pol1", 9.9998e8, 1.00002e9);   // N & Phi correlation fit
    //CorrFit4 -> SetParName(0, "C");
    //CorrFit4 -> SetParName(1, "Grad");

    //CorrFit5 = new TF1 ("CorrFit5", "pol1", 6.4458e4, 6.4462e4);   // TauGamma & A correlation fit
    //CorrFit5 -> SetParName(0, "C");
    //CorrFit5 -> SetParName(1, "Grad");

    //CorrFit6 = new TF1 ("CorrFit6", "pol1", 6.4458e4, 6.4462e4);   // TauGamma & Omega correlation fit
    //CorrFit6 -> SetParName(0, "C");
    //CorrFit6 -> SetParName(1, "Grad");

    //CorrFit7 = new TF1 ("CorrFit7", "pol1", 6.4458e4, 6.4462e4);   // TauGamma & Phi correlation fit
    //CorrFit7 -> SetParName(0, "C");
    //CorrFit7 -> SetParName(1, "Grad");

    //CorrFit8 = new TF1 ("CorrFit8", "pol1", 0.049985, 0.050015);   // A & Omega correlation fit
    //CorrFit8 -> SetParName(0, "C");
    //CorrFit8 -> SetParName(1, "Grad");

    //CorrFit9 = new TF1 ("CorrFit9", "pol1", 0.049985, 0.050015);   // A & Phi correlation fit
    //CorrFit9 -> SetParName(0, "C");
    //CorrFit9 -> SetParName(1, "Grad");

    //CorrFit10 = new TF1 ("CorrFit10", "pol1", 1.495985e-3, 1.49601e-3);   // Omega & Phi correlation fit
    //CorrFit10 -> SetParName(0, "C");
    //CorrFit10 -> SetParName(1, "Grad");

    //k1 -> Fill(fitted_N, fitted_tau_gamma);
    //k1 -> Fit("CorrFit1","R");
    //k1 -> SetTitle("Correlation between N_{Fitted} and #tau#gamma_{Fitted} (5 Parameter Fit)");
    //k1 -> GetXaxis() -> SetTitle("N_{Fitted}");
    //k1 -> GetYaxis() -> SetTitle("#tau#gamma_{Fitted}");
    //k1 -> SetMarkerStyle(20);
    //k1 -> SetMarkerSize(0.5);
    //k1 -> SetMarkerColor(kBlack);

    //k2 -> Fill(fitted_N, fitted_A);
    //k2 -> Fit("CorrFit2","R");
    //k2 -> SetTitle("Correlation between N_{Fitted} and A_{Fitted} (5 Parameter Fit)");
    //k2 -> GetXaxis() -> SetTitle("N_{Fitted}");
    //k2 -> GetYaxis() -> SetTitle("A_{Fitted}");
    //k2 -> SetMarkerStyle(20);
    //k2 -> SetMarkerSize(0.5);
    //k2 -> SetMarkerColor(kBlack);
    
    //k3 -> Fill(fitted_N, fitted_omega);
    //k3 -> Fit("CorrFit3","R");
    //k3 -> SetTitle("Correlation between N_{Fitted} and #omega_{Fitted} (5 Parameter Fit)");
    //k3 -> GetXaxis() -> SetTitle("N_{Fitted}");
    //k3 -> GetYaxis() -> SetTitle("#omega_{Fitted}");
    //k3 -> SetMarkerStyle(20);
    //k3 -> SetMarkerSize(0.5);
    //k3 -> SetMarkerColor(kBlack);

    //k4 -> Fill(fitted_N, fitted_phase);
    //k4 -> Fit("CorrFit4","R");
    //k4 -> SetTitle("Correlation between N_{Fitted} and #phi_{Fitted} (5 Parameter Fit)");
    //k4 -> GetXaxis() -> SetTitle("N_{Fitted}");
    //k4 -> GetYaxis() -> SetTitle("#phi_{Fitted}");
    //k4 -> SetMarkerStyle(20);
    //k4 -> SetMarkerSize(0.5);
    //k4 -> SetMarkerColor(kBlack);

    //k5 -> Fill(fitted_tau_gamma, fitted_A);
    //k5 -> Fit("CorrFit5","R");
    //k5 -> SetTitle("Correlation between #tau#gamma_{Fitted} and A_{Fitted} (5 Parameter Fit)");
    //k5 -> GetXaxis() -> SetTitle("#tau#gamma_{Fitted}");
    //k5 -> GetYaxis() -> SetTitle("A_{Fitted}");
    //k5 -> SetMarkerStyle(20);
    //k5 -> SetMarkerSize(0.5);
    //k5 -> SetMarkerColor(kBlack);

    //k6 -> Fill(fitted_tau_gamma, fitted_omega);
    //k6 -> Fit("CorrFit6","R");
    //k6 -> SetTitle("Correlation between #tau#gamma_{Fitted} and #omega_{Fitted} (5 Parameter Fit)");
    //k6 -> GetXaxis() -> SetTitle("#tau#gamma_{Fitted}");
    //k6 -> GetYaxis() -> SetTitle("#omega_{Fitted}");
    //k6 -> SetMarkerStyle(20);
    //k6 -> SetMarkerSize(0.5);
    //k6 -> SetMarkerColor(kBlack);

    //k7 -> Fill(fitted_tau_gamma, fitted_phase);
    //k7 -> Fit("CorrFit7","R");
    //k7 -> SetTitle("Correlation between #tau#gamma_{Fitted} and #phi_{Fitted} (5 Parameter Fit)");
    //k7 -> GetXaxis() -> SetTitle("#tau#gamma_{Fitted}");
    //k7 -> GetYaxis() -> SetTitle("#phi_{Fitted}");
    //k7 -> SetMarkerStyle(20);
    //k7 -> SetMarkerSize(0.5);
    //k7 -> SetMarkerColor(kBlack);

    //k8 -> Fill(fitted_A, fitted_omega);
    //k8 -> Fit("CorrFit8","R");
    //k8 -> SetTitle("Correlation between A_{Fitted} and #omega_{Fitted} (5 Parameter Fit)");
    //k8 -> GetXaxis() -> SetTitle("A_{Fitted}");
    //k8 -> GetYaxis() -> SetTitle("#omega_{Fitted}");
    //k8 -> SetMarkerStyle(20);
    //k8 -> SetMarkerSize(0.5);
    //k8 -> SetMarkerColor(kBlack);

    //k9 -> Fill(fitted_A, fitted_phase);
    //k9 -> Fit("CorrFit9","R");
    //k9 -> SetTitle("Correlation between A_{Fitted} and #phi_{Fitted} (5 Parameter Fit)");
    //k9 -> GetXaxis() -> SetTitle("A_{Fitted}");
    //k9 -> GetYaxis() -> SetTitle("#phi_{Fitted}");
    //k9 -> SetMarkerStyle(20);
    //k9 -> SetMarkerSize(0.5);
    //k9 -> SetMarkerColor(kBlack);

    //k10 -> Fill(fitted_omega, fitted_phase);
    //k10 -> Fit("CorrFit10","R");
    //k10 -> SetTitle("Correlation between #omega_{Fitted} and #phi_{Fitted} (5 Parameter Fit)");
    //k10 -> GetXaxis() -> SetTitle("#omega_{Fitted}");
    //k10 -> GetYaxis() -> SetTitle("#phi_{Fitted}");
    //k10 -> SetMarkerStyle(20);
    //k10 -> SetMarkerSize(0.5);
    //k10 -> SetMarkerColor(kBlack);

    // ChiSquare distribution histogram

    //double ChiSq = fit -> GetChisquare()/fit->GetNDF();
    //cs -> Fill(ChiSq);
    //cs -> SetTitle("#chi^{2} distribution (5 Parameter Fit)");
    //cs -> GetXaxis() -> SetTitle("#chi^{2}");
    //cs -> GetYaxis() -> SetTitle("Number of Entries");
   
    // Error on fitted parameters histograms

    //double ErrorPar0 = fit -> GetParError(0);
    //e1 -> Fill(ErrorPar0);
    //e1 -> SetTitle("Error on N_{Fitted} (5 Parameter Fit)");
    //e1 -> GetXaxis() -> SetTitle("Error");
    //e1 -> GetYaxis() -> SetTitle("Number of Entries");
    
    //double ErrorPar1 = fit -> GetParError(1);
    //e2 -> Fill(ErrorPar1);
    //e2 -> SetTitle("Error on #tau#gamma_{Fitted} (5 Parameter Fit)");
    //e2 -> GetXaxis() -> SetTitle("Error");
    //e2 -> GetYaxis() -> SetTitle("Number of Entries");

    //double ErrorPar2 = fit -> GetParError(2);
    //e3 -> Fill(ErrorPar2);
    //e3 -> SetTitle("Error on A_{Fitted} (5 Parameter Fit)");
    //e3 -> GetXaxis() -> SetTitle("Error");
    //e3 -> GetYaxis() -> SetTitle("Number of Entries");

    //double ErrorPar3 = fit -> GetParError(3);
    //e4 -> Fill(ErrorPar3);
    //e4 -> SetTitle("Error on #omega_{Fitted} (5 Parameter Fit)");
    //e4 -> GetXaxis() -> SetTitle("Error");
    //e4 -> GetYaxis() -> SetTitle("Number of Entries");

    //double ErrorPar4 = fit -> GetParError(4);
    //e5 -> Fill(ErrorPar4);
    //e5 -> SetTitle("Error on #phi_{Fitted} (5 Parameter Fit)");
    //e5 -> GetXaxis() -> SetTitle("Error");
    //e5 -> GetYaxis() -> SetTitle("Number of Entries");

    // Fractional difference histograms
    
    //double diff_N = (ideal_N - fitted_N)/(ideal_N);
    //g1 -> Fill(diff_N);
    //g1 -> SetTitle("Fractional difference in N (5 Parameter Fit)");
    //g1 -> GetXaxis() -> SetTitle("#frac{N - N'}{N}");
    //g1 -> GetXaxis() -> SetTitleOffset(1.2);
    //g1 -> GetYaxis() -> SetTitle("Number of Entries");
    
    //double diff_tau_gamma = (ideal_tau_gamma - fitted_tau_gamma)/(ideal_tau_gamma);
    //g2 -> Fill(diff_tau_gamma);
    //g2 -> SetTitle("Fractional difference in #tau#gamma (5 Parameter Fit)");
    //g2 -> GetXaxis() -> SetTitle("#frac{#tau#gamma - #tau#gamma'}{#tau#gamma}");
    //g2 -> GetXaxis() -> SetTitleOffset(1.2);
    //g2 -> GetYaxis() -> SetTitle("Number of Entries");
    
    //double diff_A = (ideal_A - fitted_A)/(ideal_A);
    //g3 -> Fill(diff_A);
    //g3 -> SetTitle("Fractional difference in A (5 Parameter Fit)");
    //g3 -> GetXaxis() -> SetTitle("#frac{A - A'}{A}");
    //g3 -> GetXaxis() -> SetTitleOffset(1.2);
    //g3 -> GetYaxis() -> SetTitle("Number of Entries");
    
    //double diff_omega = (ideal_omega - fitted_omega)/(ideal_omega);
    //g4 -> Fill(diff_omega);
    //g4 -> SetTitle("Fractional difference in #omega (5 Parameter Fit)");
    //g4 -> GetXaxis() -> SetTitle("#frac{#omega - #omega'}{#omega}");
    //g4 -> GetXaxis() -> SetTitleOffset(1.2);
    //g4 -> GetYaxis() -> SetTitle("Number of Entries");

    //double diff_phase = (ideal_phase - fitted_phase)/(ideal_phase);
    //g5 -> Fill(diff_phase);
    //g5 -> SetTitle("Fractional difference in #phi (5 Parameter Fit)");
    //g5 -> GetXaxis() -> SetTitle("#frac{#phi - #phi'}{#phi}");
    //g5 -> GetXaxis() -> SetTitleOffset(1.2);
    //g5 -> GetYaxis() -> SetTitle("Number of Entries");

    h1 -> Write();
    
    //  if (display_hist.size() > 0 && i == display_hist.at(0)){
    
    // h1 -> Draw();
    // c2 -> SaveAs(outname.c_str());
    // display_hist.erase(display_hist.begin());
      
    // }
    
  }

  // k1 -> Draw();
  //CorrFit1 -> Draw("SAME");
  //c2 -> SaveAs("correlation_N_tau_gamma_5PF.eps");
  
  //k2 -> Draw();
  //CorrFit2 -> Draw("SAME");
  //c2 -> SaveAs("correlation_N_A_5PF.eps");
  
  //k3 -> Draw();
  //CorrFit3 -> Draw("SAME");
  //c2 -> SaveAs("correlation_N_omega_5PF.eps");

  //k4 -> Draw();
  //CorrFit4 -> Draw("SAME");
  //c2 -> SaveAs("correlation_N_phase_5PF.eps");
  
  //k5 -> Draw();
  //CorrFit5 -> Draw("SAME");
  //c2 -> SaveAs("correlation_tau_gamma_A_5PF.eps");
  
  //k6 -> Draw();
  //CorrFit6 -> Draw("SAME");
  //c2 -> SaveAs("correlation_tau_gamma_omega_5PF.eps");
  
  //k7 -> Draw();
  //CorrFit7 -> Draw("SAME");
  //c2 -> SaveAs("correlation_tau_gamma_phase_5PF.eps");
  
  //k8 -> Draw();
  //CorrFit8 -> Draw("SAME");
  //c2 -> SaveAs("correlation_A_omega_5PF.eps");
  
  //k9 -> Draw();
  //CorrFit9 -> Draw("SAME");
  //c2 -> SaveAs("correlation_A_phase_5PF.eps");
  
  //k10 -> Draw();
  //CorrFit10 -> Draw("SAME");
  //c2 -> SaveAs("correlation_phase_omega_5PF.eps");

  ///

  //cs -> Draw();
  //c2 -> SaveAs("ChiSquare_5PF.eps");
  
  ///
  
  //e1 -> Draw();
  //c2 -> SaveAs("Error_N_5PF.eps");
  
  //e2 -> Draw();
  //c2 -> SaveAs("Error_TauGamma_5PF.eps");
  
  //e3 -> Draw();
  //c2 -> SaveAs("Error_A_5PF.eps");
  
  //e4 -> Draw();
  //c2 -> SaveAs("Error_Omega_5PF.eps");
  
  //e5 -> Draw();
  //c2 -> SaveAs("Error_Phase_5PF.eps");  

  ///
  
  //g1 -> Draw();
  //c2 -> SaveAs("Diff_N_5PF.eps");
  
  //g2 -> Draw();
  //c2 -> SaveAs("Diff_Tau_Gamma_5PF.eps");
  
  //g3 -> Draw();
  //c2 -> SaveAs("Diff_A_5PF.eps");
  
  //g4 -> Draw();
  //c2 -> SaveAs("Diff_Omega_5PF.eps");
  
  //g5 -> Draw();
  //c2 -> SaveAs("Diff_Phase_5PF.eps");
    
}
