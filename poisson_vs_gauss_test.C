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
TH1D* gHistoToFit = nullptr;

// --- Default chisquare with Neyman bias ---
void DefaultFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistoToFit->GetNbinsX(); i++) {
    double O_i = gHistoToFit->GetBinContent(i);
    if (O_i > 0) chi2 += pow(O_i - E, 2) / pow(gHistoToFit->GetBinError(i), 2); // Standard Neyman chi^2
  }
  f = chi2;
}

// --- Pearson chisquare ---
void PearsonFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistoToFit->GetNbinsX(); i++)
  {
    if (E > 0) chi2 += pow(gHistoToFit->GetBinContent(i) - E, 2) / E; // Expected variance
  }
  f = chi2;
}

// --- Yates chisquare ---
void YatesFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistoToFit->GetNbinsX(); i++)
  {
    double O_i = gHistoToFit->GetBinContent(i);
    if (O_i > 0) {
      double diff = abs(O_i - E) - 0.5;
      chi2 += (diff < 0 ? 0 : pow(diff, 2)) / O_i; // Correction with observed variance
    }
  }
  f = chi2;
}

// --- Yates & Pearson chisquare ---
void YatesPearsonFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistoToFit->GetNbinsX(); i++)
  {
    if (E > 0) {
    double O_i = gHistoToFit->GetBinContent(i);
      double diff = abs(O_i - E) - 0.5;
      chi2 += (diff < 0 ? 0 : pow(diff, 2)) / E; // Continuity correction
    }
  }
  f = chi2;
}

// --- Nagy, Csanad, et al. (Gaussian to Poisson best fit) chisquare ---
void NagyCsanadFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistoToFit->GetNbinsX(); i++)
  {
    double O_i = gHistoToFit->GetBinContent(i);
    chi2 += pow(O_i - E, 2) / (O_i + 0.5); // Modified variance for Gaussian to Poisson
  }
  f = chi2;
}

// --- Log Likelihood ---
void LikelihoodFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0];
  for (int i = 1; i <= gHistoToFit->GetNbinsX(); i++)
  {
    double O_i = gHistoToFit->GetBinContent(i);
    if (E > 0) chi2 += 2 * (E - O_i + (O_i > 0 ? O_i * log(O_i / E) : 0)); // Poisson LL
  }
  f = chi2;
}

// --- The "Neater" Runner Function ---
void ExecuteFit(TH1D* h, void (*fcn)(int&, double*, double&, double*, int), 
        string name, string formula, string saveName, double startVal)
{  
  gHistoToFit = h;
  TVirtualFitter *fitter = TVirtualFitter::Fitter(h);
  fitter->SetFCN(fcn);
  
  fitter->SetParameter(0, "Constant", startVal, 0.01, 0, 0);
  double args[2] = {5000, 0.01};
  fitter->ExecuteCommand("MIGRAD", args, 2);
  fitter->ExecuteCommand("HESSE", args, 2);

  double val = fitter->GetParameter(0);
  double err = fitter->GetParError(0);
  double chi2, edm, errdef;
  int nv, np;
  fitter->GetStats(chi2, edm, errdef, nv, np);
  int ndf = h->GetNbinsX() - 1;

  TCanvas *c = new TCanvas("c", "c", 800, 600);
  h->SetTitle(Form("%s: %s;Bin Index;Value", name.c_str(), formula.c_str()));
  h->SetMinimum(0.0*startVal);
  h->SetMaximum(2.0*startVal);
  if(startVal==1.0) h->SetMaximum(2.5*startVal);
  h->Draw("E");

  
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
}

// --- Main Test Function ---
void poisson_vs_gauss_test(int mean = 100)
{
  gStyle->SetOptStat(0);
  int nbins = 4e5;
    int hitsM = nbins*mean/1e6;

  // Data generation
  TH1D *h = new TH1D("h", "h", nbins, 0, nbins);

  for (int i = 0; i < mean * nbins; i++)
    h->Fill(gRandom->Uniform(0, nbins));
  h->Sumw2();

  // Methods Configuration
  struct Method { string name; string formula; void (*fcn)(int&, double*, double&, double*, int); };
  Method list[] =
  {
    {"default",  "#chi^{2}_{i} = (O_{i}-E)^{2} / O_{i}",             DefaultFCN},
    {"Pearson",  "#chi^{2}_{i} = (O_{i}-E)^{2} / E",                 PearsonFCN},
    {"Yates",    "#chi^{2}_{i} = (|O_{i}-E|-0.5)^{2} / O_{i}",       YatesFCN},
    {"YatesMod", "#chi^{2}_{i} = (|O_{i}-E|-0.5)^{2} / E",           YatesPearsonFCN},
    {"corr",     "#chi^{2}_{i} = (O_{i}-E)^{2} / (O_{i}+0.5)",       NagyCsanadFCN},
    {"LL",       "#chi^{2}_{i} = 2 (E - O_{i} + O_{i} ln(O_{i}/E))", LikelihoodFCN}
  };

  for (auto &m : list)
  {
    ExecuteFit(h, m.fcn, m.name, m.formula, "poisson_vs_gauss_test_" + to_string(hitsM) + "M_" + m.name + ".png", (double)mean);
  }
}#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include <iostream>

using namespace std;

TH1D* gHistoToFit = nullptr;

void Likelihood(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0]; // E: expected value
  for (int i = 1; i <= gHistoToFit->GetNbinsX(); i++)
	{
    double O_i = gHistoToFit->GetBinContent(i); // O_i i-th observed value 
    chi2 += 2*(E - O_i + O_i*log(O_i/E));
  }
  f = chi2;
}

void PearsonFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0]; // E: expected value
  for (int i = 1; i <= gHistoToFit->GetNbinsX(); i++)
	{
    double O_i = gHistoToFit->GetBinContent(i); // O_i i-th observed value 
    double sigma2 = E; // Pearson variance (expected E instead of observed O_i) to avoid Neyman bias
    
    if (sigma2 > 0) 
		{
      double diff = abs(O_i - E);
      chi2 += (diff * diff) / sigma2;
    }
  }
  f = chi2;
}

void YatesFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0]; // E: expected value
  for (int i = 1; i <= gHistoToFit->GetNbinsX(); i++)
	{
    double O_i = gHistoToFit->GetBinContent(i); // O_i i-th observed value 
    double sigma2 = E;
    
    if (sigma2 > 0) 
		{
      double diff = abs(O_i - E) - 0.5;
      if (diff < 0) diff = 0; // Yates' logic
      chi2 += (diff * diff) / sigma2;
    }
  }
  f = chi2;
}

void YatesModFCN(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chi2 = 0;
  double E = par[0]; // E: expected value
  for (int i = 1; i <= gHistoToFit->GetNbinsX(); i++)
	{
    double O_i = gHistoToFit->GetBinContent(i); // O_i i-th observed value 
    double sigma2 = O_i; // Pearson variance (expected E instead of observed O_i) to avoid Neyman bias
    
    if (sigma2 > 0) 
		{
      double diff = abs(O_i - E) - 0.5;
      if (diff < 0) diff = 0; // Yates' logic
      chi2 += (diff * diff) / sigma2;
    }
  }
  f = chi2;
}

void poisson_vs_gauss_test(int mean = 100) 
{
	const long int Nbins = 4e5;
	const long int Nhits = mean*Nbins;
	
  TH1D *h = new TH1D("h", "h", Nbins, 0, Nbins);

  for (int i = 0; i < Nhits; ++i)
    h->Fill(gRandom->Uniform(0, Nbins));

/******** DEFAULT METHOD ********/

  TFitResultPtr r = h->Fit("pol0", "SE");
  //TFitResultPtr r = h->Fit("pol0", "SL"); // Log likelihood
  //TFitResultPtr r = h->Fit("pol0", "SP"); // Pearson chi^2

  double parValue = r->Parameter(0);
	double parError = r->ParError(0);
	cerr << "err: " << r->ParError(0) << endl;
  double chi2 = r->Chi2();
  int ndf = r->Ndf();
  double prob = TMath::Prob(chi2,ndf); //r->Prob();

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1", "Fit Canvas", 1000, 700);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.5);
  h->SetMinimum(0.25*mean);
  h->SetMaximum(1.75*mean);
	h->SetTitle(Form("Random Uniform Distr., %dM hits, #chi_{i}^{2}=(O_{i}-E)^{2}/O_{i};Bin Index;Counts",(int)(Nhits/1e6)));
  h->Draw("E");

  TLatex l;
  l.SetNDC();
  l.SetTextSize(0.03);
  l.DrawLatex(0.15, 0.87, Form("Parameter: %.2f #pm %.2f", parValue, parError));
  l.DrawLatex(0.15, 0.83, Form("#chi^{2} / ndf: %.2f / %d", chi2, ndf));
  l.DrawLatex(0.15, 0.79, Form("Prob (CL): %.2e", prob));

  c1->Print(Form("poisson_vs_gauss_test_%dM.png",(int)(Nhits/1e6)));
	
/******** CORRECTING THE DIFFERENCE (YATES' CORRECTION) ********/

  gHistoToFit=h;
  TVirtualFitter *fitter = TVirtualFitter::Fitter(gHistoToFit);
  fitter->SetFCN(YatesFCN);
  fitter->SetParameter(0, "Constant", 100, 0.1, 0, 0);
  double fitarglist[2];
  fitarglist[0] = 5000; // number of function calls
  fitarglist[1] = 0.01; // tolerance
  fitter->ExecuteCommand("MIGRAD", fitarglist, 2);
  fitter->ExecuteCommand("HESSE", fitarglist, 2);
	
	h->SetTitle(Form("Random Uniform Distr., %dM hits, #chi_{i}^{2}=(|O_{i}-E|-0.5)^{2}/E;Bin Index;Counts",(int)(Nhits/1e6)));
  h->Draw("E");

  parValue = fitter->GetParameter(0);
  parError  = fitter->GetParError(0);
  double edm, errdef;
  int nvpar, nparx;
  fitter->GetStats(chi2, edm, errdef, nvpar, nparx);
  ndf = gHistoToFit->GetNbinsX() - 1;
  prob = TMath::Prob(chi2,ndf);
	
  l.DrawLatex(0.15, 0.87, Form("Parameter: %.2f #pm %.2f", parValue, parError));
  l.DrawLatex(0.15, 0.83, Form("#chi^{2} / ndf: %.2f / %d", chi2, ndf));
  l.DrawLatex(0.15, 0.79, Form("Prob (CL): %.2e", prob));

  c1->Print(Form("poisson_vs_gauss_test_%dM_Yates.png",(int)(Nhits/1e6)));
	
/******** CORRECTING SIGMA (PEARSON) ********/

  //fitter = TVirtualFitter::Fitter(gHistoToFit);
  fitter->SetFCN(PearsonFCN);
  fitter->SetParameter(0, "Constant", 100, 0.1, 0, 0);
  fitter->ExecuteCommand("MIGRAD", fitarglist, 2);
  fitter->ExecuteCommand("HESSE", fitarglist, 2);
	
	h->SetTitle(Form("Random Uniform Distr., %dM hits, #chi_{i}^{2}=(O_{i}-E)^{2}/E;Bin Index;Counts",(int)(Nhits/1e6)));
  h->Draw("E");

  parValue = fitter->GetParameter(0);
  parError  = fitter->GetParError(0);
  fitter->GetStats(chi2, edm, errdef, nvpar, nparx);
  ndf = gHistoToFit->GetNbinsX() - 1;
  prob = TMath::Prob(chi2,ndf);
	
  l.DrawLatex(0.15, 0.87, Form("Parameter: %.2f #pm %.2f", parValue, parError));
  l.DrawLatex(0.15, 0.83, Form("#chi^{2} / ndf: %.2f / %d", chi2, ndf));
  l.DrawLatex(0.15, 0.79, Form("Prob (CL): %.2e", prob));

  c1->Print(Form("poisson_vs_gauss_test_%dM_Pearson.png",(int)(Nhits/1e6)));
	
/******** YATES & PEARSON ********/

  gHistoToFit=h;
  fitter = TVirtualFitter::Fitter(gHistoToFit);
  fitter->SetFCN(YatesModFCN);
  fitter->SetParameter(0, "Constant", 100, 0.1, 0, 0);
  fitter->ExecuteCommand("MIGRAD", fitarglist, 2);
  fitter->ExecuteCommand("HESSE", fitarglist, 2);
	
	h->SetTitle(Form("Random Uniform Distr., %dM hits, #chi_{i}^{2}=(|O_{i}-E|-0.5)^{2}/O_{i};Bin Index;Counts",(int)(Nhits/1e6)));
  h->Draw("E");

  parValue = fitter->GetParameter(0);
  parError  = fitter->GetParError(0);
  fitter->GetStats(chi2, edm, errdef, nvpar, nparx);
  ndf = gHistoToFit->GetNbinsX() - 1;
  prob = TMath::Prob(chi2,ndf);
	
  l.DrawLatex(0.15, 0.87, Form("Parameter: %.2f #pm %.2f", parValue, parError));
  l.DrawLatex(0.15, 0.83, Form("#chi^{2} / ndf: %.2f / %d", chi2, ndf));
  l.DrawLatex(0.15, 0.79, Form("Prob (CL): %.2e", prob));

  c1->Print(Form("poisson_vs_gauss_test_%dM_YatesMod.png",(int)(Nhits/1e6)));
	
/******** CORRECTING SIGMA ********/
	
  for (int i = 1; i <= h->GetNbinsX(); i++) 
    h->SetBinError(i,sqrt(h->GetBinContent(i)+0.5));

  r = h->Fit("pol0", "SE");
  parValue = r->Parameter(0);
	parError = r->ParError(0);
  chi2 = r->Chi2();
  ndf = r->Ndf();
  prob = TMath::Prob(chi2,ndf); //r->Prob();
	
	h->SetTitle(Form("Random Uniform Distr., %dM hits, #chi_{i}^{2}=(O_{i}-E)^{2}/(O_{i}+0.5);Bin Index;Counts",(int)(Nhits/1e6)));
  h->Draw("E");

  l.DrawLatex(0.15, 0.87, Form("Parameter: %.2f #pm %.2f", parValue, parError));
  l.DrawLatex(0.15, 0.83, Form("#chi^{2} / ndf: %.2f / %d", chi2, ndf));
  l.DrawLatex(0.15, 0.79, Form("Prob (CL): %.2e", prob));

  c1->Print(Form("poisson_vs_gauss_test_%dM_corr.png",(int)(Nhits/1e6)));
	
/******** LIKELIHOOD ********/

  gHistoToFit=h;
  fitter = TVirtualFitter::Fitter(gHistoToFit);
  fitter->SetFCN(Likelihood);
  fitter->SetParameter(0, "Constant", 100, 0.1, 0, 0);
  fitter->ExecuteCommand("MIGRAD", fitarglist, 2);
  fitter->ExecuteCommand("HESSE", fitarglist, 2);
	
	h->SetTitle(Form("Random Uniform Distr., %dM hits, #chi_{i}^{2}=E - O_{i} + O_{i} log(O_{i}/E);Bin Index;Counts",(int)(Nhits/1e6)));
  h->Draw("E");

  parValue = fitter->GetParameter(0);
  parError  = fitter->GetParError(0);
  fitter->GetStats(chi2, edm, errdef, nvpar, nparx);
  ndf = gHistoToFit->GetNbinsX() - 1;
  prob = TMath::Prob(chi2,ndf);
	
  l.DrawLatex(0.15, 0.87, Form("Parameter: %.2f #pm %.2f", parValue, parError));
  l.DrawLatex(0.15, 0.83, Form("#chi^{2} / ndf: %.2f / %d", chi2, ndf));
  l.DrawLatex(0.15, 0.79, Form("Prob (CL): %.2e", prob));

  c1->Print(Form("poisson_vs_gauss_test_%dM_LL.png",(int)(Nhits/1e6)));
}
