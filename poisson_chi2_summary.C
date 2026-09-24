// poisson_chi2_summary.C
//
// Summary plot of the chi^2 / likelihood estimator tests, made from the "TABLE:" lines
// printed (to stderr) by poisson_vs_gauss_test.C and poisson_vs_gauss_test_ratio.C.
//
// Collect the table lines, e.g.:
//   for m in 100 200 400; do
//     root -l -b -q "poisson_vs_gauss_test.C($m)"       2>&1 | grep TABLE >> hist_table.log
//     root -l -b -q "poisson_vs_gauss_test_ratio.C($m)" 2>&1 | grep TABLE >> ratio_table.log
//   done
// then make the plot:
//   root -l -b -q 'poisson_chi2_summary.C("hist_table.log","ratio_table.log")'
//
// Layout (2x2 + legend column):
//   top row:    bias of the fitted value vs. lambda (single histogram | ratio)
//   bottom row: chi^2/NDF vs. lambda                 (single histogram | ratio)
// Markers: fit results. Dashed lines: analytic expectations. Grey bands: statistical spread
// (fit uncertainty in the top row, sqrt(2/NDF) in the bottom row).
// If the same method and lambda appear more than once in a file, the last line is used.

#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMath.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// ---------------- configuration ----------------
const double gNbins  = 4e5;   // number of bins used in the tests
const int    gFontPx = 30;    // font size in pixels (canvas is 1600 px wide);
                              // ~9 pt when the figure is printed 17 cm wide.
                              // Increase if the figure is printed smaller.
const double gXmin = 70, gXmax = 560;  // lambda range of the plots

struct Result { double val, err, chi2, chi2ndf; };
typedef map<double, Result> Series;        // lambda -> result (sorted, last line wins)
typedef map<string, Series> Table;         // method name -> series

struct MethodStyle { string label; int color; int marker; double offset; };

// Method names as printed by the test macros -> legend label, color, marker, x offset.
// The small multiplicative x offsets keep overlapping markers apart on the log axis.
const vector<pair<string, MethodStyle>> gMethods = {
  {"Classical",     {"Neyman",                       kRed+1,     20, 0.94}},
  {"Corrected",     {"Shifted variance", kOrange+7,  21, 0.96}},
  {"Yates",         {"Yates",                        kMagenta+1, 23, 0.98}},
  {"Pearson",       {"Pearson",                      kAzure+2,   24, 1.00}},
  {"Yates-Pearson", {"Yates-Pearson",                kCyan+2,    26, 1.02}},
  {"LogLikelihood", {"Log-likelihood",               kGreen+2,   33, 1.04}},
  {"CFLikelihood",  {"CF likelihood",                kGray+3,    34, 1.06}},
};

// ---------------- reading the TABLE lines ----------------
static string Trim(const string& s)
{
  size_t b = s.find_first_not_of(" \t\r\n");
  if (b == string::npos) return "";
  size_t e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

// Histogram lines: | name | 40M | 100 | val±err | chi2 | chi2/ndf | prob |
// Ratio lines:     | name | 40M | val±err | chi2 | chi2/ndf | prob |
Table ReadTable(const char* fname, bool isRatio)
{
  Table table;
  ifstream in(fname);
  if (!in) { cerr << "ERROR: cannot open " << fname << endl; return table; }
  string line;
  while (getline(in, line)) {
    size_t pos = line.find("TABLE:");
    if (pos == string::npos) continue;
    vector<string> f;
    stringstream ss(line.substr(pos + 6));
    string tok;
    while (getline(ss, tok, '|')) { tok = Trim(tok); if (!tok.empty()) f.push_back(tok); }
    const size_t nExpected = isRatio ? 6 : 7;
    if (f.size() != nExpected) { cerr << "WARNING: skipping line: " << line << endl; continue; }

    const int iVal = isRatio ? 2 : 3;
    const double hitsM = atof(f[1].c_str());                 // "40M" -> 40
    const double lambda = isRatio ? hitsM * 1e6 / gNbins : atof(f[2].c_str());
    const string& ve = f[iVal];
    const size_t pm = ve.find("\xC2\xB1");                    // UTF-8 "±"
    Result r;
    r.val     = atof(ve.substr(0, pm).c_str());
    r.err     = (pm == string::npos) ? 0 : atof(ve.substr(pm + 2).c_str());
    r.chi2    = atof(f[iVal + 1].c_str());
    r.chi2ndf = atof(f[iVal + 2].c_str());
    table[f[0]][lambda] = r;
  }
  cout << "Read " << table.size() << " methods from " << fname << endl;
  return table;
}

// ---------------- drawing helpers ----------------
static TH1F* DrawFrame(double ymin, double ymax, const char* ytitle)
{
  gPad->SetLogx();
  gPad->SetTicks(1, 1);
  TH1F* fr = gPad->DrawFrame(gXmin, ymin, gXmax, ymax, Form(";Occupancy #lambda;%s", ytitle));
  for (TAxis* ax : {fr->GetXaxis(), fr->GetYaxis()}) {
    ax->SetLabelFont(43); ax->SetLabelSize(gFontPx);
    ax->SetTitleFont(43); ax->SetTitleSize(gFontPx);
  }
  fr->GetXaxis()->SetMoreLogLabels();
  fr->GetXaxis()->SetNoExponent();
  fr->GetXaxis()->SetTitleOffset(1.25);
  fr->GetYaxis()->SetTitleOffset(1.9);
  fr->GetYaxis()->SetNdivisions(506);
  return fr;
}

// Grey band between center-halfwidth(x) and center+halfwidth(x)
static void DrawBand(TF1* halfwidth, double center)
{
  const int n = 100;
  TGraph* g = new TGraph(2 * n);
  for (int i = 0; i < n; i++) {
    double x = gXmin * pow(gXmax / gXmin, i / (n - 1.0));
    g->SetPoint(i, x, center + halfwidth->Eval(x));
    g->SetPoint(2 * n - 1 - i, x, center - halfwidth->Eval(x));
  }
  g->SetFillColorAlpha(kGray, 0.6);
  g->SetLineWidth(0);
  g->Draw("F");
}

static void DrawReference(double y)
{
  TLine* l = new TLine(gXmin, y, gXmax, y);
  l->SetLineStyle(3); l->SetLineColor(kGray + 2);
  l->Draw();
}

static void DrawCurve(const char* name, const char* formula, int color)
{
  TF1* f = new TF1(name, formula, gXmin, gXmax);
  f->SetNpx(200);
  f->SetLineColor(color); f->SetLineStyle(7); f->SetLineWidth(2);
  f->Draw("same");
}

// quantity: 0 = hist bias (val - lambda), 1 = ratio bias in % (100*(val-1)), 2 = chi2/NDF
static void DrawPoints(const Table& t, int quantity)
{
  for (auto& m : gMethods) {
    auto it = t.find(m.first);
    if (it == t.end()) continue;
    TGraph* g = new TGraph();
    for (auto& p : it->second) {
      const double lam = p.first;
      const Result& r = p.second;
      double y = (quantity == 0) ? r.val - lam : (quantity == 1) ? 100.0 * (r.val - 1.0) : r.chi2ndf;
      g->SetPoint(g->GetN(), lam * m.second.offset, y);
    }
    g->SetMarkerStyle(m.second.marker);
    g->SetMarkerColor(m.second.color);
    g->SetLineColor(m.second.color);
    g->SetMarkerSize(2.2);
    g->Draw("P");
  }
}

// Panel title in the top-left (top = true) or bottom-right (top = false) corner
static void PanelTitle(const char* text, bool top = true)
{
  TLatex l; l.SetNDC(); l.SetTextFont(43); l.SetTextSize(gFontPx);
  if (top) {
    l.SetTextAlign(13);
    l.DrawLatex(gPad->GetLeftMargin() + 0.03, 1.0 - gPad->GetTopMargin() - 0.04, text);
  } else {
    l.SetTextAlign(31);
    l.DrawLatex(1.0 - gPad->GetRightMargin() - 0.03, gPad->GetBottomMargin() + 0.04, text);
  }
}

static TPad* MakePad(const char* name, double x1, double y1, double x2, double y2)
{
  TPad* p = new TPad(name, name, x1, y1, x2, y2);
  p->SetLeftMargin(0.21); p->SetRightMargin(0.03);
  p->SetTopMargin(0.04);  p->SetBottomMargin(0.17);
  p->Draw();
  return p;
}

// ---------------- main ----------------
void poisson_chi2_summary(const char* histFile = "hist_table.log",
                          const char* ratioFile = "ratio_table.log",
                          const char* outName = "poisson_chi2_summary")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  Table hist  = ReadTable(histFile,  false);
  Table ratio = ReadTable(ratioFile, true);

  TCanvas* c = new TCanvas("csum", "csum", 1600, 1100);
  const double xl = 0.78;   // start of the legend column
  TPad* pHB = MakePad("pHB", 0.0,    0.5, xl / 2, 1.0);
  TPad* pRB = MakePad("pRB", xl / 2, 0.5, xl,     1.0);
  TPad* pHC = MakePad("pHC", 0.0,    0.0, xl / 2, 0.5);
  TPad* pRC = MakePad("pRC", xl / 2, 0.0, xl,     0.5);
  const double sigmaChi2 = sqrt(2.0 / (gNbins - 1));

  // (a) single histogram: absolute bias in counts
  pHB->cd();
  DrawFrame(-1.3, 0.9, "Bias of fitted value, #hat{E} #minus #lambda");
  DrawBand(new TF1("bHB", Form("sqrt(x/%g)", gNbins), gXmin, gXmax), 0.0);
  DrawReference(0.0);
  DrawCurve("aHN", "-1-1/x",        kRed + 1);      // Neyman, Eq. (shift_bias) with c=0
  DrawCurve("aHC", "-1-0.5/x",      kOrange + 7);   // shifted variance, c=0.5
  DrawCurve("aHP", "sqrt(x*x+x)-x", kAzure + 2);    // Pearson
  DrawPoints(hist, 0);
  PanelTitle("(a) Single histogram");
  gPad->RedrawAxis();

  // (b) ratio: relative bias in percent
  pRB->cd();
  DrawFrame(-3.6, 0.8, "Bias of fitted ratio, #hat{F} #minus 1 (%)");
  DrawBand(new TF1("bRB", Form("100*sqrt(2/(%g*x))", gNbins), gXmin, gXmax), 0.0);
  DrawReference(0.0);
  DrawCurve("aRN", "-300/x", kRed + 1);             // Neyman ratio, F = 1 - 3/lambda
  DrawPoints(ratio, 1);
  PanelTitle("(b) Ratio of two histograms");
  gPad->RedrawAxis();

  // (c) single histogram: chi2/NDF
  pHC->cd();
  DrawFrame(0.905, 1.025, "#chi^{2}/NDF");
  DrawBand(new TF1("bHC", Form("%g", sigmaChi2), gXmin, gXmax), 1.0);
  DrawReference(1.0);
  DrawCurve("aHY", "1-sqrt(2/(pi*x))+1/(4*x)", kCyan + 2);  // Yates, sigma^2 = lambda
  DrawPoints(hist, 2);
  PanelTitle("(c) Single histogram", false);
  gPad->RedrawAxis();

  // (d) ratio: chi2/NDF
  pRC->cd();
  DrawFrame(0.905, 1.025, "#chi^{2}/NDF");
  DrawBand(new TF1("bRC", Form("%g", sigmaChi2), gXmin, gXmax), 1.0);
  DrawReference(1.0);
  DrawCurve("aRY", "1-sqrt(1/(pi*x))+1/(8*x)", kCyan + 2);  // Yates, sigma^2 = 2 lambda
  DrawPoints(ratio, 2);
  PanelTitle("(d) Ratio of two histograms", false);
  gPad->RedrawAxis();

  // legend column
  c->cd();
  TPad* pL = new TPad("pL", "pL", xl, 0.0, 1.0, 1.0);
  pL->Draw(); pL->cd();
  TLegend* leg = new TLegend(0.02, 0.25, 0.98, 0.85);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->SetTextFont(43); leg->SetTextSize(gFontPx);
  for (auto& m : gMethods) {
    TGraph* g = new TGraph();
    g->SetMarkerStyle(m.second.marker); g->SetMarkerColor(m.second.color); g->SetMarkerSize(2.2);
    leg->AddEntry(g, m.second.label.c_str(), "p");
  }
  TF1* fDummy = new TF1("fDummy", "0", 0, 1);
  fDummy->SetLineStyle(7); fDummy->SetLineColor(kBlack); fDummy->SetLineWidth(2);
  leg->AddEntry(fDummy, "Analytic expectation", "l");
  TGraph* gBand = new TGraph();
  gBand->SetFillColorAlpha(kGray, 0.6); gBand->SetLineWidth(0);
  leg->AddEntry(gBand, "Statistical spread", "f");
  leg->Draw();

  c->Print(Form("%s.pdf", outName));
  c->Print(Form("%s.png", outName));
}
