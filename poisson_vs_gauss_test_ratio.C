#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TVirtualFitter.h"
#include <iostream>
#include <string>

using namespace std;

// Global pointer for the custom FCNs to access the data
TH1D* gHistA = nullptr;
TH1D* gHistB = nullptr;

// --- Default chisquare with Neyman bias ---
void DefaultFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistA->GetNbinsX(); i++) {
		double A_i = gHistA->GetBinContent(i);
    double B_i = gHistB->GetBinContent(i);
    if (A_i > 0 && B_i > 0) chi2 += pow(A_i - E*B_i, 2) / (A_i*(1+A_i/B_i)); // Standard Neyman chi^2
  }
  f = chi2;
}

// --- Pearson chisquare ---
void PearsonFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistA->GetNbinsX(); i++) {
		double A_i = gHistA->GetBinContent(i);
    double B_i = gHistB->GetBinContent(i);
    if (E > 0 && B_i > 0) chi2 += pow(A_i - E*B_i, 2) / (E*B_i*(1+A_i/B_i));  // Expected variance
  }
  f = chi2;
}

// --- Yates chisquare ---
void YatesFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistA->GetNbinsX(); i++) {
		double A_i = gHistA->GetBinContent(i);
    double B_i = gHistB->GetBinContent(i);
    double diff = abs(A_i - E*B_i) - 0.5;
	  if(diff < 0) diff = 0;
    if (A_i > 0 && B_i > 0) chi2 += diff * diff / (A_i*(1+A_i/B_i));  // Correction with observed variance
  }
  f = chi2;
}

// --- Yates & Pearson chisquare ---
void YatesPearsonFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistA->GetNbinsX(); i++) {
    if (E > 0) {
		  double A_i = gHistA->GetBinContent(i);
      double B_i = gHistB->GetBinContent(i);
      double diff = abs(A_i - E*B_i) - 0.5;
			if(diff<0) diff = 0;
      if (E > 0 && B_i > 0) chi2 += diff * diff / (E*B_i*(1+A_i/B_i));   // Continuity correction
    }
  }
  f = chi2;
}

// --- Nagy, Csanad, et al. (Gaussian to Poisson best fit) chisquare ---
void NagyCsanadFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistA->GetNbinsX(); i++) {
		double A_i = gHistA->GetBinContent(i);
    double B_i = gHistB->GetBinContent(i);
    if (B_i > 0) chi2 += pow(A_i - E*B_i, 2) / (A_i*(1+A_i/B_i)+0.5); // Modified variance for Gaussian to Poisson
  }
  f = chi2;
}

// --- Log Likelihood ---
void LikelihoodFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistA->GetNbinsX(); i++) {
		double A_i = gHistA->GetBinContent(i);
    double B_i = gHistB->GetBinContent(i);
    if (E > 0 && A_i > 0 && B_i > 0) chi2 += (E*B_i - A_i + A_i * log(A_i / (E*B_i))); // Poisson LL
  }
  f = chi2;
}

// --- The "Neater" Runner Function ---
void ExecuteFit(TH1D* hA, TH1D* hB, void (*fcn)(int&, double*, double&, double*, int), 
        string name, string formula, string saveName)
{  
  gHistA = hA;
  gHistB = hB;
  TH1D *hC = (TH1D*)hA->Clone(hA->GetName());
  for (int i = 1; i <= hC->GetNbinsX(); i++)
	{
		hC->SetBinContent(i,hA->GetBinContent(i)*1.0/hB->GetBinContent(i));
		hC->SetBinError(i,hC->GetBinContent(i)*sqrt(1.0/hA->GetBinContent(i)+1.0/hB->GetBinContent(i)));
	}
  TVirtualFitter *fitter = TVirtualFitter::Fitter(hC);
  fitter->SetFCN(fcn);
  
  fitter->SetParameter(0, "Constant", 1.0, 0.01, 0, 0);
  double args[2] = {5000, 0.01};
  fitter->ExecuteCommand("MIGRAD", args, 2);
  fitter->ExecuteCommand("HESSE", args, 2);

  double val = fitter->GetParameter(0);
  double err = fitter->GetParError(0);
  double chi2, edm, errdef;
  int nv, np;
  fitter->GetStats(chi2, edm, errdef, nv, np);
  int ndf = hC->GetNbinsX() - 1;

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  hC->SetTitle(Form("%s: %s;Bin Index;Value", name.c_str(), formula.c_str()));
  hC->SetMinimum(0.0);
  hC->SetMaximum(2.5);
  hC->Draw("E");

  
  TLatex l; l.SetNDC(); l.SetTextSize(0.03);
  l.DrawLatex(0.15, 0.85, Form("Method: %s", name.c_str()));
  l.DrawLatex(0.15, 0.81, Form("Fitted: %.3f #pm %.3f", val, err));
  double prob = TMath::Prob(chi2, ndf);
  if(prob>0.001)
    l.DrawLatex(0.15, 0.77, Form("#chi^{2}/ndf: %.2f/%d (P: %.1f%%)", chi2, ndf, prob*100));
  else
    l.DrawLatex(0.15, 0.77, Form("#chi^{2}/ndf: %.2f/%d (P: %.3e%%)", chi2, ndf, prob*100));
  
  c->Print(saveName.c_str());
  delete c;
  cerr << "TABLE:\t| " << name << " | " << (int)(hA->GetEntries()/1e6) << "M | " << Form("%.3f±%.3f",val,err) << " | " << Form("%.2f",chi2) << " | " << Form("%.4f",chi2/ndf) << " | " << Form("%.2f%%",prob*100) << " |" << endl;
}

// --- Main Test Function ---
void poisson_vs_gauss_test_ratio(int mean = 100)
{
  gStyle->SetOptStat(0);
  int nbins = 4e5;
		int hitsM = nbins*mean/1e6;

  // Data generation
  TH1D *hA = new TH1D("hA", "hA", nbins, 0, nbins);
  TH1D *hB = new TH1D("hB", "hB", nbins, 0, nbins);

  for (int i = 0; i < mean * nbins; ++i) {
    hA->Fill(gRandom->Uniform(0, nbins));
    hB->Fill(gRandom->Uniform(0, nbins));
  }

  // Methods Configuration
  struct Method { string name; string formula; void (*fcn)(int&, double*, double&, double*, int); };
  Method list[] =
	{
    {"default",  "#chi^{2}_{i} = (A_{i} - F #upoint B_{i})^{2} / (A_{i} #upoint (1+A_{i}/B_{i}))",             DefaultFCN},
    {"Pearson",  "#chi^{2}_{i} = (A_{i} - F #upoint B_{i})^{2} / (F #upoint B_{i} #upoint (1+A_{i}/B_{i}))",                 PearsonFCN},
    {"Yates",    "#chi^{2}_{i} = (|A_{i} - F #upoint B_{i}|-0.5)^{2} / (A_{i} #upoint (1+A_{i}/B_{i}))",       YatesFCN},
    {"YatesMod", "#chi^{2}_{i} = (|A_{i} - F #upoint B_{i}E|-0.5)^{2} / (F #upoint B_{i} #upoint (1+A_{i}/B_{i}))",           YatesPearsonFCN},
    {"corr",     "#chi^{2}_{i} = (A_{i} - F #upoint B_{i})^{2} / (A_{i} #upoint (1+A_{i}/B_{i}) + 0.5)",       NagyCsanadFCN},
    {"LL",       "#chi^{2}_{i} = (F #upoint B_{i} - A_{i} + A_{i} ln(A_{i}/(F #upoint B_{i})))", LikelihoodFCN}
  };

  for (auto &m : list)
	{
    ExecuteFit(hA, hB, m.fcn, m.name, m.formula, "poisson_vs_gauss_test_ratio_" + to_string(hitsM) + "M_" + m.name + ".png");
  }
}