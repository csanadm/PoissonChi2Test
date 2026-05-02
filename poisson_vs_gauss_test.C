#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include <iostream>

using namespace std;

TH1D* gHistoToFit = nullptr;

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

void poisson_vs_gauss_test() 
{
	const long int Nbins = 4e5;
	const int mean = 400;
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
}