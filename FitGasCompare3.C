// ================================================================
// FitGasCompare_3x3.cxx
//
// 3x3 summary canvas with ratio panels (same style as standalone plots)
//
// Rows    : NPE, SPR, sigma(theta_C)
// Columns : sipmDef, sipm75, sipmUV
//
// Each cell: upper pad (main plot) + lower pad (ratio/delta panel)
//            exactly like the standalone 900x900 canvases
//
// Canvas size: 1800 x 2700 (portrait, 3 cols x 600, 3 rows x 900)
//
// Run with: root -l -q FitGasCompare_3x3.cxx
// ================================================================

#ifndef __CINT__
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TSystem.h>
#include <iostream>
#include <cmath>
#include <cstring>
#endif

using namespace std;

// ============================================================
// GLOBAL CONFIGURATION
// ============================================================

const int  nCases = 4;
const char *caseDir[nCases]   = {"dataDef", "dataOpt1", "dataOpt2",  "dataOpt3"};
const char *caseLabel[nCases] = {"Default", "Abs+Ray",  "Abs. Len.", "Ray. Len."};

const int  nSiPM = 3;
const char *sipmDir[nSiPM]   = { "sipmDef",      "sipm75",   "sipmUV"};
const char *sipmLabel[nSiPM] = { "SiPM Default", "SiPM 75",  "SiPM UV"};

const char *sipmLabel2[3][nSiPM] = { 
			"<NPE> (default)",	"<NPE> (75 #mu m)",	"<NPE> (UV ext.)",	
			"SPR (default)",	"SPR (75 #mu m)",	"SPR (UV ext.)",	
			"#sigma(#theta_{C}^{trk}) (default)",	"#sigma(#theta_{C}^{trk}) (75 #mu m)", "#sigma(#theta_{C}^{trk}) (UV ext.)"			
			};

const int  nRad = 2;
const char *radName[nRad] = {"Gas", "Aerogel"};

const int nEtaPlot = 4;
const int nEtaFile = 5;

double etaPlot[nEtaPlot] = {1.8, 2.4, 2.8, 3.2};

const double etaBandMin = 1.7;
const double etaBandMax = 3.3;

// ============================================================
// STYLE
// ============================================================

void SetPublicationStyle(){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetTitleSize(0.050, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleOffset(1.50, "X");
    gStyle->SetTitleOffset(1.55, "Y");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLineWidth(1);
}

// ============================================================
// FILE / HISTOGRAM HELPERS
// ============================================================

TString GetFileName(int sipm, int dataCase, int fileIndex){
    TString filename;
    filename.Form("%s/%s/ana_211_30_%d_0.root", caseDir[dataCase],  sipmDir[sipm], fileIndex);
    return filename;
}

TString GetHistPath(const char *rad,  const char *histBase){
    TString path;
    path.Form("pid%s/%s_%s", rad, histBase, rad);
    return path;
}

// ============================================================
// GAUSSIAN FIT
// ============================================================

bool FitGaussian(TH1 *h,
                 double &mean,    double &meanErr,    double &sigma,    double &sigmaErr,
                 double fitMin = -1.0,    double fitMax = -1.0){
    mean = 0.0;    meanErr = 0.0;    sigma = 0.0;    sigmaErr = 0.0;
    if (!h) return false;
    if (h->GetEntries() <= 0) return false;

    double xmin = fitMin;
    double xmax = fitMax;
    if (fitMin < 0.0 || fitMax < 0.0) {
        double hmean = h->GetMean();
        double hrms  = h->GetRMS();
        if (hrms <= 0.0) return false;
        xmin = hmean - 3.0 * hrms;
        xmax = hmean + 3.0 * hrms;
    }
    double axisMin = h->GetXaxis()->GetXmin();
    double axisMax = h->GetXaxis()->GetXmax();
    if (xmin < axisMin) xmin = axisMin;
    if (xmax > axisMax) xmax = axisMax;
    if (xmax <= xmin) return false;

    TString fname;
    fname.Form("gaus_%s", h->GetName());
    TF1 fit(fname.Data(), "gaus", xmin, xmax);
    fit.SetParameter(1, h->GetMean());
    fit.SetParameter(2, h->GetRMS());
    int status = h->Fit(&fit, "RQ0");
    if (status != 0) return false;

    mean     = fit.GetParameter(1);
    meanErr  = fit.GetParError(1);
    sigma    = std::fabs(fit.GetParameter(2));
    sigmaErr = fit.GetParError(2);
    return true;
}

bool GetGaussianMean(TH1 *h, double &mean, double &meanErr,
                     double fitMin = -1.0, double fitMax = -1.0){
    double sigma = 0.0, sigmaErr = 0.0;
    return FitGaussian(h, mean, meanErr, sigma, sigmaErr, fitMin, fitMax);
}

bool GetGaussianSigma(TH1 *h, double &sigma, double &sigmaErr,
                      double fitMin = -1.0,  double fitMax = -1.0){
    double mean = 0.0, meanErr = 0.0;
    return FitGaussian(h, mean, meanErr, sigma, sigmaErr, fitMin, fitMax);
}

// ============================================================
// OBSERVABLE READERS
// ============================================================

bool GetNPE(TFile *f, const char *rad, double &value, double &error){
    value = 0.0; error = 0.0;
    if (!f) return false;
    TString path = GetHistPath(rad, "npe_dist");
    TH1 *h = dynamic_cast<TH1 *>(f->Get(path));
    if (!h){ cerr << "ERROR: Histogram not found: " << path << endl; return false; }
    return GetGaussianMean(h, value, error, 0.0, 70.0);
}

bool GetThetaMean(TFile *f, const char *rad, double &value, double &error){
    value = 0.0; error = 0.0;
    if (!f) return false;
    TString path = GetHistPath(rad, "theta_dist");
    TH1 *h = dynamic_cast<TH1 *>(f->Get(path));
    if (!h) { cerr << "ERROR: Histogram not found: " << path << endl; return false; }
    double fitMin, fitMax;
    if (strcmp(rad, "Gas") == 0) { fitMin = 10.0; fitMax = 40.0; }
    else if (strcmp(rad, "Aerogel") == 0) { fitMin = 150.0; fitMax = 250.0; }
    else { fitMin = -1.0; fitMax = -1.0; }
    return GetGaussianMean(h, value, error, fitMin, fitMax);
}

bool GetThetaResidualSigma(TFile *f, const char *rad, double &value, double &error){
    value = 0.0; error = 0.0;
    if (!f) return false;
    TString path = GetHistPath(rad, "thetaResid_dist");
    TH1 *h = dynamic_cast<TH1 *>(f->Get(path));
    if (!h){ cerr << "ERROR: Histogram not found: " << path << endl; return false; }
    return GetGaussianSigma(h, value, error);
}

bool GetSinglePhotonResolution(TFile *f, const char *rad, double &value, double &error){
    value = 0.0; error = 0.0;
    if (!f) return false;
    TString path;
    path.Form("pid%s/photonTheta_vs_photonPhi_%s", rad, rad);
    TH2 *h2 = dynamic_cast<TH2 *>(f->Get(path));
    if (!h2){ cerr << "ERROR: 2D histogram not found: " << path << endl; return false; }
    TString projName;
    projName.Form("photonThetaProjection_%s", rad);
    TH1D *hProj = h2->ProjectionY(projName.Data());
    if (!hProj) return false;
    hProj->SetDirectory(nullptr);
    double fitMin, fitMax;
    if (strcmp(rad, "Gas") == 0) { fitMin = 10.0; fitMax = 40.0; }
    else if (strcmp(rad, "Aerogel") == 0) { fitMin = 150.0; fitMax = 250.0; }
    else { fitMin = -1.0; fitMax = -1.0; }
    bool ok = GetGaussianSigma( hProj, value, error, fitMin, fitMax);
    delete hProj;
    return ok;
}

// ============================================================
// GRAPH CREATION / STYLE
// ============================================================

TGraphErrors *MakeGraph(double values[nEtaPlot], double errors[nEtaPlot]){
    TGraphErrors *g = new TGraphErrors(nEtaPlot);
    for (int i = 0; i < nEtaPlot; i++) {
        g->SetPoint(i, etaPlot[i], values[i]);
        g->SetPointError(i, 0.0, errors[i]);
    }
    return g;
}

void StyleGraph(TGraphErrors *g, int caseIndex, bool ratio = false){
    if (!g) return;
    int color[4] = {kBlack, kBlue + 1, kRed + 1, kGreen + 3};
    int mark[4]  = {20,     25,        24,        23};
    g->SetLineColor(color[caseIndex]);
    g->SetMarkerColor(color[caseIndex]);
    g->SetMarkerStyle(mark[caseIndex]);
    g->SetMarkerSize(1.5);
    g->SetLineWidth(2);
    if (caseIndex == 0) g->SetLineStyle(1);
    else if (caseIndex == 1) g->SetLineStyle(1);
    else if (caseIndex == 2) g->SetLineStyle(2);
    else g->SetLineStyle(3);
}

void DrawText(double x, double y, const char *text, double size = 0.045){
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(size);
    latex.DrawLatex(x, y, text);
}

// ============================================================
// RATIO / DELTA PANEL
// ============================================================

TGraphErrors *CalculateRatio( TGraphErrors *numerator, TGraphErrors *denominator){
    if (!numerator || !denominator) return nullptr;
    TGraphErrors *ratio = new TGraphErrors(nEtaPlot);
    for (int i = 0; i < nEtaPlot; i++){
        double xn, yn, xd, yd;
        numerator->GetPoint(i, xn, yn);
        denominator->GetPoint(i, xd, yd);
        double enY = numerator->GetErrorY(i);
        double edY = denominator->GetErrorY(i);
        double r = yn - yd;
        double er = std::sqrt( enY * enY + edY * edY );
        ratio->SetPoint(i, xn, r);
        ratio->SetPointError(i, 0.0, er);
    }
    return ratio;
}

// ============================================================
// DRAW RATIO PANEL (standalone version, for 3x3 cell)
// ============================================================

void DrawRatioPanelInPad( TGraphErrors *graphs[nCases], const char *observableName,
                          double ymin, double ymax){
    // Pad already created and cd()'d by caller
    TH1F *frame = new TH1F(Form("ratioFrame_%s", observableName), "", 100, 1.7, 3.3);
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);
    frame->GetYaxis()->SetLimits(ymin, ymax);
    frame->GetYaxis()->SetNdivisions(5, 0, 3);
    frame->GetXaxis()->SetTitle("#eta");
    TString ratioTitle;
    ratioTitle.Form("#Delta_{(%s)}", observableName);
    frame->GetYaxis()->SetTitle(ratioTitle);
    frame->GetXaxis()->SetLabelSize(0.13);
    frame->GetXaxis()->SetTitleSize(0.14);
    frame->GetXaxis()->SetTitleOffset(0.9);
    frame->GetYaxis()->SetLabelSize(0.11);
    frame->GetYaxis()->SetLabelOffset(0.01);    
    frame->GetYaxis()->SetTitleSize(0.11);
    frame->GetYaxis()->SetTitleOffset(0.6);
    frame->Draw();

    TGraphErrors *ratio[nCases] = {nullptr};
    ratio[0] = CalculateRatio(graphs[0], graphs[0]);
    for (int c = 1; c < nCases; c++){
        ratio[c] = CalculateRatio(graphs[c], graphs[0]);
    }
    for (int c = 1; c < nCases; c++) {
        if (!ratio[c]) continue;
        StyleGraph(ratio[c], c, true);
        ratio[c]->Draw("P SAME");
    }
}

// ============================================================
// OBSERVABLE + RATIO CANVAS (original per-SiPM per-rad plots)
// ============================================================

void MakeObservableRatioCanvas( TGraphErrors *graphs[nCases], const char *observableName,
    const char *yTitle, const char *title, const char *outputPDF, const char *outputPNG,
    const char *outputROOT, double upperYmin, double upperYmax, double deltaYmin, double deltaYmax){

    TCanvas *c = new TCanvas("cObservableRatio", title, 900, 900);

    TPad *upper = new TPad("upperPad", "upperPad", 0.0, 0.32, 1.0, 1.0);
    upper->SetTopMargin(0.07);
    upper->SetBottomMargin(0.015);
    upper->SetLeftMargin(0.14);
    upper->SetRightMargin(0.04);
    upper->Draw();
    upper->cd();

    double ymin = 1.0e30;
    double ymax = -1.0e30;
    for (int ccase = 0; ccase < nCases; ccase++) {
        if (!graphs[ccase]) continue;
        for (int i = 0; i < nEtaPlot; i++) {
            double x, y;
            graphs[ccase]->GetPoint( i, x, y );
            double ey = graphs[ccase]->GetErrorY(i);
            if (y <= 0.0) continue;
            ymin = std::min( ymin, y - ey);
            ymax = std::max( ymax, y + ey);
        }
    }
    if (ymin == 1.0e30 || ymax == -1.0e30) { ymin = 0.0; ymax = 1.0; }
    double range = ymax - ymin;
    if (range <= 0.0) range = 1.0;
    ymin -= 0.13 * range;
    ymax += 0.15 * range;
    if (ymin < 0.0) ymin = 0.0;

    TH1F *frame = new TH1F( "observableFrame", "", 100, 1.7, 3.3 );
    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);
    frame->GetXaxis()->SetTitle("");
    frame->GetYaxis()->SetTitle(yTitle);
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTitleSize(0.0);
    frame->GetYaxis()->SetLabelSize(0.06);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->Draw();

    for (int ccase = 0; ccase < nCases; ccase++) {
        if (!graphs[ccase]) continue;
        StyleGraph( graphs[ccase], ccase, false );
        graphs[ccase]->Draw("P SAME");
    }

    TLegend *leg = new TLegend( 0.68, 0.68, 0.92, 0.92);
    leg->SetTextSize(0.048);
    for (int ccase = 0; ccase < nCases; ccase++) {
        if (!graphs[ccase]) continue;
        leg->AddEntry( graphs[ccase], caseLabel[ccase], "lp");
    }
    leg->Draw();
    DrawText( 0.18, 0.87, title, 0.050);

    c->cd();
    TPad *ratioPad = new TPad("ratioPad", "ratioPad", 0.0, 0.0, 1.0, 0.33);
    ratioPad->SetTopMargin(0.025);
    ratioPad->SetBottomMargin(0.28);
    ratioPad->SetLeftMargin(0.14);
    ratioPad->SetRightMargin(0.04);
    ratioPad->Draw();
    ratioPad->cd();
    DrawRatioPanelInPad( graphs, observableName, deltaYmin, deltaYmax );

    c->SaveAs(outputPDF);
    c->SaveAs(outputPNG);
    c->SaveAs(outputROOT);
    delete c;
}

// ============================================================
// 3x3 SUMMARY CANVAS WITH RATIO PANELS
// ============================================================
// Rows    : NPE, SPR, sigma(theta_C)
// Columns : sipmDef, sipm75, sipmUV
//
// Each cell has upper pad (main) + lower pad (ratio)
// Canvas: 1800 (W) x 2700 (H)  -- 3 cols x 600, 3 rows x 900
// ============================================================

void Make3x3SummaryCanvas(
    TGraphErrors* gNPE[nSiPM][nCases],
    TGraphErrors* gSPR[nSiPM][nCases],
    TGraphErrors* gSigma[nSiPM][nCases],
    const char* rad,
    const char* outPDF,
    const char* outPNG)
{
    // Canvas: 1800 x 2700 (portrait)
    // Each cell: 600 x 900 (same proportions as standalone 900x900 canvas)
    TCanvas* c3 = new TCanvas("c3x3", Form("3x3 Summary %s", rad), 1800, 1800);
    c3->SetFixedAspectRatio(kFALSE);

    const char* rowTitle[3] = {"<NPE>", "SPR [mrad]", "#sigma(#theta_{C}) [mrad]"};
    double deltaYmin[3] = {-0.49, -0.19, -0.12};
    double deltaYmax[3] = { 2.45,  0.11,  0.12};

    // Cell dimensions in NDC
    double cellW = 1.0 / 3.0;  // 0.3333
    double cellH = 1.0 / 3.0;  // 0.3333

    for(int row=0; row<3; row++){
        for(int col=0; col<3; col++){
            // Select graphs for this cell
            TGraphErrors* graphs[nCases];
            if(row==0) for(int cc=0; cc<nCases; cc++) graphs[cc] = gNPE[col][cc];
            else if(row==1) for(int cc=0; cc<nCases; cc++) graphs[cc] = gSPR[col][cc];
            else for(int cc=0; cc<nCases; cc++) graphs[cc] = gSigma[col][cc];

            // Determine y-range for main plot
            double ymin = 1.0e30, ymax = -1.0e30;
            for(int cc=0; cc<nCases; cc++){
                if(!graphs[cc]) continue;
                
                for(int i=0; i<nEtaPlot; i++){
                    double x,y; graphs[cc]->GetPoint(i,x,y);
                    double ey = graphs[cc]->GetErrorY(i);
                    if(y<=0) continue;
                    ymin = std::min(ymin, y-ey);
                    ymax = std::max(ymax, y+ey);
                }
            }
            if(ymin==1.0e30 || ymax==-1.0e30){ ymin=0; ymax=1; }
            double range = ymax-ymin;
            if(range<=0) range=1;
            ymin -= 0.12*range;
            ymax += 0.15*range;
            if(ymin<0) ymin=0;

            // Cell position in NDC
            double x1 = col * cellW;
            double x2 = (col + 1) * cellW;
            double y1 = (2 - row) * cellH;  // row 0 at top
            double y2 = (3 - row) * cellH;

            // ---- Main (upper) pad: 68% of cell height ----
            double mainY1 = y1 + 0.32 * cellH;
            TPad* padMain = new TPad(
                Form("padMain_%d_%d", row, col),
                Form("padMain_%d_%d", row, col),
                x1, mainY1, x2, y2);
            padMain->SetTopMargin(0.07);
            padMain->SetBottomMargin(0.02);
            padMain->SetLeftMargin(0.14);
            padMain->SetRightMargin(0.04);
            c3->cd();
            padMain->Draw();
            padMain->cd();

            TH1F* frame = new TH1F(
                Form("frame_%d_%d", row, col), "",
                100, 1.7, 3.3);
            frame->SetMinimum(ymin);
            frame->SetMaximum(ymax);
            frame->GetXaxis()->SetTitle("");
            frame->GetYaxis()->SetTitle(rowTitle[row]);
            frame->GetXaxis()->SetLabelSize(0.0);
            frame->GetXaxis()->SetTitleSize(0.0);
            frame->GetYaxis()->SetLabelSize(0.06);
            frame->GetYaxis()->SetTitleSize(0.06);
            if(row==2)frame->GetYaxis()->SetNdivisions(5, 0, 3);
            frame->GetYaxis()->SetTitleOffset(1.1);
            frame->Draw();

            for(int cc=0; cc<nCases; cc++){
                if(!graphs[cc]) continue;
                StyleGraph(graphs[cc], cc, false);
                graphs[cc]->Draw("P SAME");
            }

            // Legend only on first column
            if(col>=0){
                TLegend* leg = new TLegend(0.64, 0.65, 0.89, 0.92);
                leg->SetTextSize(0.045);
                for(int cc=0; cc<nCases; cc++){
                    if(!graphs[cc]) continue;
                    leg->AddEntry(graphs[cc], caseLabel[cc], "lp");
                }
                leg->Draw();
            }

            // Column label (top of cell)
            // TLatex lat; lat.SetNDC(); lat.SetTextSize(0.055);
            // lat.DrawLatex(0.18, 0.87, sipmLabel[col]);
            
            TLatex lat; lat.SetNDC(); lat.SetTextSize(0.06);
            lat.DrawLatex(0.18, 0.86, sipmLabel2[row][col]);
            
            
            // ---- Ratio (lower) pad: 32% of cell height ----
            TPad* padRatio = new TPad(
                Form("padRatio_%d_%d", row, col),
                Form("padRatio_%d_%d", row, col),
                x1, y1, x2, mainY1);
            padRatio->SetTopMargin(0.02);
            padRatio->SetBottomMargin(0.28);
            padRatio->SetLeftMargin(0.14);
            padRatio->SetRightMargin(0.04);
            c3->cd();
            padRatio->Draw();
            padRatio->cd();

            DrawRatioPanelInPad(graphs, rowTitle[row], deltaYmin[row], deltaYmax[row]);
        }
    }

    c3->SaveAs(outPDF);
    c3->SaveAs(outPNG);
    delete c3;
}

// ============================================================
// MAIN
// ============================================================

void FitGasCompare3(){

    SetPublicationStyle();

    // Storage for 3x3 summary (filled during normal loop)
    TGraphErrors* summaryNPE[nSiPM][nCases];
    TGraphErrors* summarySPR[nSiPM][nCases];
    TGraphErrors* summarySigma[nSiPM][nCases];
    for(int s=0; s<nSiPM; s++)
        for(int c=0; c<nCases; c++){
            summaryNPE[s][c] = nullptr;
            summarySPR[s][c] = nullptr;
            summarySigma[s][c] = nullptr;
        }

    for (int sipm = 0; sipm < nSiPM; sipm++) {
        for (int irad = 0; irad < nRad; irad++) {
            const char *rad = radName[irad];

            cout << endl;
            cout << "=================================================" << endl;
            cout << "SiPM     : " << sipmLabel[sipm] << endl;
            cout << "Radiator : " << rad             << endl;
            cout << "=================================================" << endl;

            double npe[nCases][nEtaFile] = {};
            double npeErr[nCases][nEtaFile] = {};
            double thetaMean[nCases][nEtaPlot] = {};
            double thetaMeanErr[nCases][nEtaPlot] = {};
            double thetaResidual[nCases][nEtaPlot] = {};
            double thetaResidualErr[nCases][nEtaPlot] = {};
            double photonResolution[nCases][nEtaPlot] = {};
            double photonResolutionErr[nCases][nEtaPlot] = {};

            for (int ccase = 0; ccase < nCases; ccase++) {
                cout << endl;
                cout << "Case: " << caseLabel[ccase] << endl;
                for (int ifile = 0; ifile < nEtaFile; ifile++) {
                    TString filename =  GetFileName(sipm, ccase, ifile);
                    cout << "  Opening: " << filename << endl;
                    TFile *f = TFile::Open( filename, "READ" );
                    if (!f || f->IsZombie()) {
                        cerr << "ERROR: Cannot open " << filename << endl;
                        if (f) delete f;
                        continue;
                    }

                    bool okNPE = GetNPE( f, rad, npe[ccase][ifile], npeErr[ccase][ifile]);
                    if (!okNPE) cerr << "WARNING: NPE Gaussian fit failed: " << filename << endl;

                    if (ifile < nEtaPlot) {
                        bool okTheta = GetThetaMean(f, rad, thetaMean[ccase][ifile], thetaMeanErr[ccase][ifile] );
                        if (!okTheta) cerr << "WARNING: Theta Gaussian fit failed: " << filename << endl;
                        bool okResidual = GetThetaResidualSigma(f, rad, thetaResidual[ccase][ifile], thetaResidualErr[ccase][ifile]);
                        if (!okResidual) cerr << "WARNING: Theta residual Gaussian fit failed: " << filename << endl;
                        bool okPhoton = GetSinglePhotonResolution(f, rad, photonResolution[ccase][ifile], photonResolutionErr[ccase][ifile]);
                        if (!okPhoton) cerr << "WARNING: Photon resolution Gaussian fit failed: " << filename << endl;
                    }
                    f->Close();
                    delete f;
                }
            }

            cout << endl;
            cout << "-------------------------------------------------" << endl;
            cout << "NPE values" << endl;
            cout << "-------------------------------------------------" << endl;
            for (int ccase = 0; ccase < nCases; ccase++) {
                cout << caseLabel[ccase] << " : ";
                for (int i = 0; i < nEtaPlot; i++) {
                    cout << npe[ccase][i] << " +/- " << npeErr[ccase][i];
                    if (i < nEtaPlot - 1) cout << "   ";
                }
                cout << endl;
            }

            TString outDir;
            outDir.Form( "out/%s/%s", sipmDir[sipm], rad);
            TString mkdirCommand;
            mkdirCommand.Form("mkdir -p %s", outDir.Data());
            gSystem->Exec(mkdirCommand);

            TGraphErrors *gNPE[nCases] = {nullptr};
            TGraphErrors *gTheta[nCases] = {nullptr};
            TGraphErrors *gResidual[nCases] = {nullptr};
            TGraphErrors *gPhoton[nCases] = {nullptr};

            for (int ccase = 0; ccase < nCases; ccase++) {
                gNPE[ccase]      = MakeGraph( npe[ccase], npeErr[ccase] );
                gTheta[ccase]    = MakeGraph( thetaMean[ccase], thetaMeanErr[ccase]);
                gResidual[ccase] = MakeGraph( thetaResidual[ccase], thetaResidualErr[ccase]);
                gPhoton[ccase]   = MakeGraph( photonResolution[ccase], photonResolutionErr[ccase]);
            }

            // Save pointers for 3x3 summary
            for(int ccase=0; ccase<nCases; ccase++){
                summaryNPE[sipm][ccase]    = gNPE[ccase];
                summarySPR[sipm][ccase]    = gPhoton[ccase];
                summarySigma[sipm][ccase]  = gResidual[ccase];
            }

            TString npePDF, npePNG, npeROOT;
            TString thetaPDF, thetaPNG, thetaROOT;
            TString residualPDF, residualPNG, residualROOT;
            TString photonPDF, photonPNG, photonROOT;

            npePDF.Form( "%s/NPE.pdf",  outDir.Data());
            npePNG.Form( "%s/NPE.png",  outDir.Data());
            npeROOT.Form("%s/NPE.root", outDir.Data());

            thetaPDF.Form( "%s/ThetaMean.pdf",  outDir.Data());
            thetaPNG.Form( "%s/ThetaMean.png",  outDir.Data());
            thetaROOT.Form("%s/ThetaMean.root", outDir.Data());

            residualPDF.Form( "%s/ThetaResidualSigma.pdf", outDir.Data());
            residualPNG.Form( "%s/ThetaResidualSigma.png", outDir.Data());
            residualROOT.Form("%s/ThetaResidualSigma.root", outDir.Data());

            photonPDF.Form(   "%s/SinglePhotonResolution.pdf",  outDir.Data());
            photonPNG.Form(   "%s/SinglePhotonResolution.png",  outDir.Data());
            photonROOT.Form(  "%s/SinglePhotonResolution.root", outDir.Data());

            TString npeTitle;
            npeTitle.Form( "%s (%s)", sipmLabel[sipm], rad);
            MakeObservableRatioCanvas( gNPE, "<NPE>", "<NPE>", npeTitle.Data(),
                npePDF.Data(), npePNG.Data(), npeROOT.Data(),
                0.0, 70.0, -2.1, 2.1 );

            TString thetaTitle;
            thetaTitle.Form( "%s (%s)", sipmLabel[sipm], rad );
            MakeObservableRatioCanvas(gTheta, "<#theta_{C}> [mrad]", "<#theta_{C}> [mrad]", thetaTitle.Data(),
                thetaPDF.Data(), thetaPNG.Data(), thetaROOT.Data(),
                0.0, 300.0, -0.2, 0.2);

            TString residualTitle;
            residualTitle.Form("%s (%s)", sipmLabel[sipm], rad);
            MakeObservableRatioCanvas( gResidual, "#sigma_{#theta,res} [mrad]", "#sigma_{#theta,res} [mrad]",
                residualTitle.Data(), residualPDF.Data(), residualPNG.Data(), residualROOT.Data(),
                0.0, 2.0, -0.1, 0.1 );

            TString photonTitle;
            photonTitle.Form( "%s (%s)",  sipmLabel[sipm], rad );
            MakeObservableRatioCanvas( gPhoton, "SPR [mrad]", "SPR [mrad]",
                photonTitle.Data(), photonPDF.Data(), photonPNG.Data(), photonROOT.Data(),
                0.0, 5.0, -0.5, 0.5);

            cout << endl;
            cout << "Output written to: " << outDir << endl;
        }
    }

    // ============================================================
    // Produce 3x3 summary canvases (one per radiator)
    // ============================================================
    for(int irad=0; irad<nRad; irad++){
        const char* rad = radName[irad];
        TString sumPDF, sumPNG;
        sumPDF.Form("out/Summary_3x3_%s.pdf", rad);
        sumPNG.Form("out/Summary_3x3_%s.png", rad);
        Make3x3SummaryCanvas(summaryNPE, summarySPR, summarySigma, rad,
                             sumPDF.Data(), sumPNG.Data());
        cout << "3x3 summary saved: " << sumPDF << "  " << sumPNG << endl;
    }

    cout << endl;
    cout << "=================================================" << endl;
    cout << "FitGasCompare completed." << endl;
    cout << "=================================================" << endl;
}

