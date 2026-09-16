// ================================================================
// FitGasCompare.C
//
// Data configurations:
//   dataDef  = Default
//   dataOpt1 = Abs+Ray
//   dataOpt2 = Abs
//   dataOpt3 = Ray
//
// SiPM configurations:
//   sipmDef
//   sipm75
//   sipmUV
//
// Radiators:
//   Gas
//   Aerogel
//
// FILES:
//   ana_211_30_0_0.root  -> eta = 1.8
//   ana_211_30_1_0.root  -> eta = 2.4
//   ana_211_30_2_0.root  -> eta = 2.8
//   ana_211_30_3_0.root  -> eta = 3.2
//   ana_211_30_4_0.root  -> eta = 1.7 to 3.3
//
// IMPORTANT:
//   Index 4 is NOT plotted as an eta point.
//   It is ONLY for the horizontal band (if used).
//
// Gaussian fitting is performed for ALL four observables:
//
//   1. NPE:
//        npe_dist_<Rad>             -> Gaussian MEAN
//
//   2. Theta:
//        theta_dist_<Rad>           -> Gaussian MEAN
//
//   3. Theta residual:
//        thetaResid_dist_<Rad>     -> Gaussian SIGMA
//
//   4. Single photon:
//        photonTheta_vs_photonPhi_<Rad>
//        -> Y projection
//        -> Gaussian SIGMA
//
// Main NPE canvas:
//   Upper pad : NPE vs eta
//   Lower pad :
//       Abs+Ray / Default
//       Abs     / Default
//       Ray     / Default
//       ratio = 1 dotted line
//
// Output:
//   out/<SiPM>/<Radiator>/
//
// ================================================================


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

#include <iostream>
#include <cmath>
#include <cstring>

using namespace std;


// ============================================================
// GLOBAL CONFIGURATION
// ============================================================

const int  nCases = 4;
const char *caseDir[nCases]   = {"dataDef", "dataOpt1", "dataOpt2",  "dataOpt3"};
const char *caseLabel[nCases] = {"Default", "Abs+Ray",  "Abs. Len.", "Ray. Len."};

// ------------------------------------------------------------
// SiPM configurations
// ------------------------------------------------------------

const int  nSiPM = 3;
const char *sipmDir[nSiPM]   = { "sipmDef",      "sipm75",   "sipmUV"};
const char *sipmLabel[nSiPM] = { "SiPM Default", "SiPM 75",  "SiPM UV"};

// ------------------------------------------------------------
// Radiators
// ------------------------------------------------------------

const int  nRad = 2;
const char *radName[nRad] = {"Gas", "Aerogel"};

// ============================================================
// Y-RANGES FOR PLOTS
// ============================================================

// Upper panel Y-ranges
// Set useAuto = false if you want to explicitly specify the range.

struct YRange {
    double ymin;
    double ymax;
    bool useAuto;
};

// NPE
YRange npeRange = {0.0, 70.0, false};

// Theta mean
YRange thetaRange = {0.0, 300.0, false};

// Theta residual sigma
YRange residualRange = {0.0, 2.0, false};

// Single photon resolution
YRange photonRange = {0.0, 5.0, false};


// ============================================================
// BOTTOM Δ-PANEL Y-RANGES
// ============================================================

YRange npeDeltaRange      = {-2.1, 2.1, false};
YRange thetaDeltaRange    = {-0.2, 0.2, false};
YRange residualDeltaRange = {-0.1, 0.1, false};
YRange photonDeltaRange   = {-0.5, 0.5, false};


// ------------------------------------------------------------
// eta points
//
// File 0 -> eta = 1.8
// File 1 -> eta = 2.4
// File 2 -> eta = 2.8
// File 3 -> eta = 3.2
//
// File 4 -> full eta range 1.7--3.3
//            ONLY for NPE constant band (if used)
// ------------------------------------------------------------

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

    gStyle->SetLineWidth(2);
}


// ============================================================
// FILE NAME
// ============================================================

TString GetFileName(int sipm, int dataCase, int fileIndex){

    TString filename;
    filename.Form("%s/%s/ana_211_30_%d_0.root", caseDir[dataCase],  sipmDir[sipm], fileIndex);
    return filename;
}


// ============================================================
// HISTOGRAM PATH
// ============================================================

TString GetHistPath(const char *rad,  const char *histBase){

    TString path;
    path.Form("pid%s/%s_%s", rad, histBase, rad);
    return path;
}


// ============================================================
// GAUSSIAN FIT
//
// If fitMin and fitMax are not supplied, use
// mean +/- 2 RMS as the initial fit range.
// The range is also restricted to the histogram axis.
// ============================================================

bool FitGaussian(TH1 *h,
                 double &mean,		double &meanErr,	double &sigma,	double &sigmaErr,
                 double fitMin = -1.0,	double fitMax = -1.0){
                 
    mean    = 0.0;	meanErr = 0.0;	sigma   = 0.0;	sigmaErr = 0.0;

    if (!h)return false;
    if (h->GetEntries() <= 0)return false;

    // --------------------------------------------------------
    // Determine fit range
    // --------------------------------------------------------

    double xmin = fitMin;
    double xmax = fitMax;

    if (fitMin < 0.0 || fitMax < 0.0) {

        double hmean = h->GetMean();
        double hrms  = h->GetRMS();

        if (hrms <= 0.0) return false;

        xmin = hmean - 3.0 * hrms;
        xmax = hmean + 3.0 * hrms;
    }


    // Restrict to histogram axis range
    double axisMin = h->GetXaxis()->GetXmin();
    double axisMax = h->GetXaxis()->GetXmax();

    if (xmin < axisMin) xmin = axisMin;
    if (xmax > axisMax) xmax = axisMax;
    if (xmax <= xmin) return false;


    // --------------------------------------------------------
    // Gaussian
    // --------------------------------------------------------

    TString fname;
    fname.Form("gaus_%s", h->GetName());

    TF1 fit(fname.Data(), "gaus", xmin, xmax);

    fit.SetParameter(1, h->GetMean());
    fit.SetParameter(2, h->GetRMS());

    int status = h->Fit(&fit, "RQ0");

    if (status != 0) return false;


    // --------------------------------------------------------
    // Extract Gaussian parameters
    // --------------------------------------------------------

    mean     = fit.GetParameter(1);
    meanErr  = fit.GetParError(1);

    sigma    = std::fabs(fit.GetParameter(2));
    sigmaErr = fit.GetParError(2);

    return true;
}


// ============================================================
// GAUSSIAN MEAN
// ============================================================

bool GetGaussianMean(TH1 *h, double &mean, double &meanErr,
                     double fitMin = -1.0, double fitMax = -1.0){
                     
    double sigma = 0.0;
    double sigmaErr = 0.0;

    return FitGaussian(h, mean, meanErr, sigma, sigmaErr, fitMin, fitMax);
}


// ============================================================
// GAUSSIAN SIGMA
// ============================================================

bool GetGaussianSigma(TH1 *h, double &sigma, double &sigmaErr,
                      double fitMin = -1.0,  double fitMax = -1.0){
                      
    double mean = 0.0;
    double meanErr = 0.0;

    return FitGaussian(h, mean, meanErr, sigma, sigmaErr, fitMin, fitMax);
}


// ============================================================
// GET NPE
//
// Gaussian mean of:
// pid<Radiator>/npe_dist_<Radiator>
//
// NPE fit range: 0--70
// ============================================================

bool GetNPE(TFile *f, const char *rad, double &value, double &error){

    value = 0.0;
    error = 0.0;

    if (!f) return false;
    TString path = GetHistPath(rad, "npe_dist");

    TH1 *h = dynamic_cast<TH1 *>(f->Get(path));

    if (!h){
        cerr << "ERROR: Histogram not found: " << path << endl;
        return false;
    }

    return GetGaussianMean(h, value, error, 0.0, 70.0);
}


// ============================================================
// GET THETA MEAN
//
// Gaussian mean of:
// pid<Radiator>/theta_dist_<Radiator>
//
// Dynamic fit range = mean +/- 2 RMS
// ============================================================

bool GetThetaMean(TFile *f, const char *rad, double &value, double &error){

    value = 0.0;
    error = 0.0;

    if (!f) return false;

    TString path = GetHistPath(rad, "theta_dist");

    TH1 *h = dynamic_cast<TH1 *>(f->Get(path));

    if (!h) {
        cerr << "ERROR: Histogram not found: " << path << endl;
        return false;
    }

    double fitMin;
    double fitMax;

    if (strcmp(rad, "Gas") == 0) {
        fitMin = 10.0;
        fitMax = 40.0;
    }
    else if (strcmp(rad, "Aerogel") == 0) {
        fitMin = 150.0;
        fitMax = 250.0;
    }
    else {
        fitMin = -1.0;
        fitMax = -1.0;
    }

    return GetGaussianMean(h, value, error, fitMin, fitMax);
}



// ============================================================
// GET THETA RESIDUAL SIGMA
//
// Gaussian sigma of:
// pid<Radiator>/thetaResid_dist_<Radiator>
// ============================================================

bool GetThetaResidualSigma(TFile *f, const char *rad, double &value, double &error){

    value = 0.0;
    error = 0.0;

    if (!f) return false;
    TString path = GetHistPath(rad, "thetaResid_dist");

    TH1 *h = dynamic_cast<TH1 *>(f->Get(path));

    if (!h){
        cerr << "ERROR: Histogram not found: " << path << endl;
        return false;
    }

    return GetGaussianSigma(h, value, error);
}


// ============================================================
// GET SINGLE PHOTON RESOLUTION
//
// 2D histogram:
// pid<Radiator>/photonTheta_vs_photonPhi_<Radiator>
//
// ProjectionY() -> Gaussian fit -> sigma
// ============================================================

bool GetSinglePhotonResolution(TFile *f, const char *rad, double &value, double &error){

    value = 0.0;
    error = 0.0;

    if (!f) return false;

    TString path;
    path.Form("pid%s/photonTheta_vs_photonPhi_%s", rad, rad);


    TH2 *h2 = dynamic_cast<TH2 *>(f->Get(path));
    if (!h2){
        cerr << "ERROR: 2D histogram not found: " << path << endl;
        return false;
    }


    TString projName;
    projName.Form("photonThetaProjection_%s", rad);

    TH1D *hProj = h2->ProjectionY(projName.Data());

    if (!hProj) return false;
    hProj->SetDirectory(nullptr);


    double fitMin;
    double fitMax;

    if (strcmp(rad, "Gas") == 0) { fitMin = 10.0; fitMax = 40.0; }
    else if (strcmp(rad, "Aerogel") == 0) { fitMin = 150.0; fitMax = 250.0; }
    else { fitMin = -1.0; fitMax = -1.0; }

    bool ok = GetGaussianSigma( hProj, value, error, fitMin, fitMax);

    delete hProj;
    return ok;
}


// ============================================================
// GRAPH CREATION
// ============================================================

TGraphErrors *MakeGraph(double values[nEtaPlot], double errors[nEtaPlot]){

    TGraphErrors *g = new TGraphErrors(nEtaPlot);

    for (int i = 0; i < nEtaPlot; i++) {
        g->SetPoint(i, etaPlot[i], values[i]);
        g->SetPointError(i, 0.0, errors[i]);
    }

    return g;
}


// ============================================================
// GRAPH STYLE
// ============================================================

// "ratio corresponds to delta now"

void StyleGraph(TGraphErrors *g, int caseIndex, bool ratio = false){

    if (!g) return;

    // ROOT standard colors
    int color[4] = {kBlack, kBlue + 1, kRed + 1, kGreen + 3};
    int mark[4]  = {20, 	25, 	24, 	 23};

    g->SetLineColor(color[caseIndex]);
    g->SetMarkerColor(color[caseIndex]);
    g->SetMarkerStyle(mark[caseIndex]);
    //g->SetMarkerStyle(ratio ? 20 + caseIndex : 20 + caseIndex);
    g->SetMarkerSize(1.5);
    g->SetLineWidth(2);


    if (caseIndex == 0) g->SetLineStyle(1);
    else if (caseIndex == 1) g->SetLineStyle(1);
    else if (caseIndex == 2) g->SetLineStyle(2);
    else g->SetLineStyle(3);
}


// ============================================================
// DRAW LABEL
// ============================================================

void DrawText(double x, double y, const char *text, double size = 0.045){
    TLatex latex;

    latex.SetNDC();
    latex.SetTextSize(size);

    latex.DrawLatex(x, y, text);
}


// ============================================================
// CALCULATE RATIO
//
// ratio = numerator / denominator
//
// error propagation:
// dr/r = sqrt(
//    (dn/n)^2 +
//    (dd/d)^2
// )
// ============================================================

TGraphErrors *CalculateRatio( TGraphErrors *numerator, TGraphErrors *denominator){

    if (!numerator || !denominator) return nullptr;

    TGraphErrors *ratio = new TGraphErrors(nEtaPlot);

    for (int i = 0; i < nEtaPlot; i++){
    
        double xn, yn;
        double enX, enY;

        double xd, yd;
        double edX, edY;

        numerator->GetPoint(i, xn, yn);
        denominator->GetPoint(i, xd, yd);

        enX = numerator->GetErrorX(i);
        enY = numerator->GetErrorY(i);

        edX = denominator->GetErrorX(i);
        edY = denominator->GetErrorY(i);

        double r = 0.0;
        double er = 0.0;

        r = yn - yd;
        er = std::sqrt( enY * enY + edY * edY );
       
        ratio->SetPoint(i, xn, r);
        ratio->SetPointError(i, 0.0, er);
    }

    return ratio;
}


// ============================================================
// DRAW RATIO/DELTA PANEL
// ============================================================

void DrawRatioPanel( TGraphErrors *graphs[nCases], const char *observableName, 
                    double ymin, double ymax){

    // --------------------------------------------------------
    // Ratio/Delta pad
    // --------------------------------------------------------

    TPad *pad = new TPad("ratioPad", "ratioPad", 0.0, 0.0, 1.0, 0.33);

    pad->SetTopMargin(0.02);
    pad->SetBottomMargin(0.28);
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.04);

    pad->Draw();
    pad->cd();


    // --------------------------------------------------------
    // Ratio/Delta frame
    // --------------------------------------------------------

    // double ymin = -2.1;
    // double ymax = 2.1;

    TH1F *frame = new TH1F("ratioFrame", "", 100, 1.7, 3.3);

    frame->SetMinimum(ymin);
    frame->SetMaximum(ymax);

    frame->GetYaxis()->SetLimits(ymin, ymax);
    frame->GetYaxis()->SetNdivisions(5, 2, 5);
    
    frame->GetXaxis()->SetTitle("#eta");

    TString ratioTitle;
    ratioTitle.Form("#Delta_{(%s)}", observableName);
    frame->GetYaxis()->SetTitle(ratioTitle);

    frame->GetXaxis()->SetLabelSize(0.12);
    frame->GetXaxis()->SetTitleSize(0.13);
    frame->GetXaxis()->SetTitleOffset(0.9);

    frame->GetYaxis()->SetLabelSize(0.11);
    frame->GetYaxis()->SetTitleSize(0.12);
    frame->GetYaxis()->SetTitleOffset(0.6);
    
    frame->Draw();

    // --------------------------------------------------------
    // Ratio/Delta curves
    // --------------------------------------------------------

    TGraphErrors *ratio[nCases] = {nullptr};

    // Default / Default
    ratio[0] = CalculateRatio(graphs[0], graphs[0]);

    // Other three cases / Default
    for (int c = 1; c < nCases; c++){
        ratio[c] = CalculateRatio(graphs[c], graphs[0]);
    }


    // --------------------------------------------------------
    // Draw ratio curves
    // --------------------------------------------------------

    for (int c = 1; c < nCases; c++) {
        if (!ratio[c]) continue;

        StyleGraph(ratio[c], c, true);

        ratio[c]->Draw("P SAME");
    }

    pad->Modified();
}


// ============================================================
// COMMON OBSERVABLE + RATIO CANVAS
//
// Upper pad:
//     observable vs eta
//
// Lower pad:
//     observable / Default
//
// Used for:
//     NPE
//     Theta Mean
//     Theta Residual Sigma
//     Single Photon Resolution
// ============================================================

void MakeObservableRatioCanvas( TGraphErrors *graphs[nCases], const char *observableName,
    const char *yTitle, const char *title, const char *outputPDF, const char *outputPNG, 
    const char *outputROOT, double upperYmin, double upperYmax, double deltaYmin, double deltaYmax){
    
    // --------------------------------------------------------
    // Canvas
    // --------------------------------------------------------

    TCanvas *c = new TCanvas("cObservableRatio", title, 900, 900);


    // --------------------------------------------------------
    // Upper pad
    // --------------------------------------------------------

    TPad *upper = new TPad("upperPad", "upperPad", 0.0, 0.32, 1.0, 1.0);

    upper->SetTopMargin(0.07);
    upper->SetBottomMargin(0.015);
    upper->SetLeftMargin(0.14);
    upper->SetRightMargin(0.04);

    upper->Draw();
    upper->cd();


    // --------------------------------------------------------
    // Get y range
    // --------------------------------------------------------

    double ymin = 1.0e30;
    double ymax = -1.0e30;


    for (int ccase = 0; ccase < nCases; ccase++) {
        if (!graphs[ccase]) continue;

        for (int i = 0; i < nEtaPlot; i++) {
            double x;
            double y;

            graphs[ccase]->GetPoint( i, x, y );
            double ey = graphs[ccase]->GetErrorY(i);

            if (y <= 0.0) continue;

            ymin = std::min( ymin, y - ey);
            ymax = std::max( ymax, y + ey);
        }
    }

    if (ymin == 1.0e30 || ymax == -1.0e30) {
    
        ymin = 0.0;
        ymax = 1.0;
    }


    // Add margins
    double range = ymax - ymin;

    if (range <= 0.0) range = 1.0;

    ymin -= 0.12 * range;
    ymax += 0.15 * range;

    if (ymin < 0.0) ymin = 0.0;
    
    
    // double ymin = upperYmin;
    // double ymax = upperYmax;

    // --------------------------------------------------------
    // Upper frame
    // --------------------------------------------------------

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


    // --------------------------------------------------------
    // Observable curves
    // --------------------------------------------------------

    for (int ccase = 0; ccase < nCases; ccase++) {
    
        if (!graphs[ccase]) continue;

        StyleGraph( graphs[ccase], ccase, false );
        graphs[ccase]->Draw("P SAME");
    }


    // --------------------------------------------------------
    // Legend
    // --------------------------------------------------------

    TLegend *leg = new TLegend( 0.68, 0.68, 0.92, 0.92);
    leg->SetTextSize(0.048);

    for (int ccase = 0; ccase < nCases; ccase++) {
        if (!graphs[ccase]) continue;
        leg->AddEntry( graphs[ccase], caseLabel[ccase], "lp");
    }

    leg->Draw();


    // --------------------------------------------------------
    // Labels
    // --------------------------------------------------------

    DrawText( 0.18, 0.87, title, 0.050);


    // --------------------------------------------------------
    // Ratio pad
    // --------------------------------------------------------

    c->cd();

    DrawRatioPanel( graphs, observableName, deltaYmin, deltaYmax );


    // --------------------------------------------------------
    // Save
    // --------------------------------------------------------

    c->SaveAs(outputPDF);
    c->SaveAs(outputPNG);
    c->SaveAs(outputROOT);
    delete c;
}


// ============================================================
// MAIN ANALYSIS
// ============================================================

void FitGasCompare(){

    SetPublicationStyle();

    // ========================================================
    // Loop over SiPMs
    // ========================================================

    for (int sipm = 0; sipm < nSiPM; sipm++) {
    
        // ----------------------------------------------------
        // Loop over radiators
        // ----------------------------------------------------

        for (int irad = 0; irad < nRad; irad++) {
        
            const char *rad = radName[irad];

            cout << endl;
            cout << "=================================================" << endl;
            cout << "SiPM     : " << sipmLabel[sipm] << endl;
            cout << "Radiator : " << rad             << endl;
            cout << "=================================================" << endl;

            // =================================================
            // Arrays
            //
            // [case][eta]
            // =================================================

            double npe[nCases][nEtaFile] = {};
            double npeErr[nCases][nEtaFile] = {};

            double thetaMean[nCases][nEtaPlot] = {};
            double thetaMeanErr[nCases][nEtaPlot] = {};

            double thetaResidual[nCases][nEtaPlot] = {};
            double thetaResidualErr[nCases][nEtaPlot] = {};

            double photonResolution[nCases][nEtaPlot] = {};
            double photonResolutionErr[nCases][nEtaPlot] = {};

            // =================================================
            // Read ALL cases
            //
            // Gaussian fitting is performed independently
            // for every case.
            // =================================================

            for (int ccase = 0; ccase < nCases; ccase++) {
            
                cout << endl;
                cout << "Case: " << caseLabel[ccase] << endl;

                // ------------------------------------------------
                // Loop over all five files
                // ------------------------------------------------

                for (int ifile = 0; ifile < nEtaFile; ifile++) {
                
                    TString filename =  GetFileName(sipm, ccase, ifile);
                    cout << "  Opening: " << filename << endl;

                    TFile *f = TFile::Open( filename, "READ" );

                    if (!f || f->IsZombie()) {
                    
                        cerr << "ERROR: Cannot open " << filename << endl;
                        if (f) delete f;
                        continue;
                    }


                    // =================================================
                    // NPE
                    // Fit for ALL four cases.
                    // =================================================

                    bool okNPE = GetNPE( f, rad, npe[ccase][ifile], npeErr[ccase][ifile]);
                    if (!okNPE) cerr << "WARNING: NPE Gaussian fit failed: " << filename << endl;

                    // =================================================
                    // Other Observables.
                    // =================================================

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


            // ========================================================
            // Print values
            // ========================================================

            cout << endl;
            cout << "-------------------------------------------------" << endl;
            cout << "NPE values" 					<< endl;
            cout << "-------------------------------------------------" << endl;

            for (int ccase = 0; ccase < nCases; ccase++) {

                cout << caseLabel[ccase] << " : ";
                for (int i = 0; i < nEtaPlot; i++) {
                    cout << npe[ccase][i] << " +/- " << npeErr[ccase][i];
                    if (i < nEtaPlot - 1) cout << "   ";
                }
                cout << endl;
            }


            // ========================================================
            // Create output directory
            // ========================================================

            TString outDir;
            outDir.Form( "out/%s/%s", sipmDir[sipm], rad);
            TString mkdirCommand;
            mkdirCommand.Form("mkdir -p %s", outDir.Data());
            gSystem->Exec(mkdirCommand);


            // ========================================================
            // Create graphs
            // ========================================================

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

            // ========================================================
            // OUTPUT FILE NAMES
            // ========================================================

            TString npePDF;
            TString npePNG;
            TString npeROOT;

            TString thetaPDF;
            TString thetaPNG;
            TString thetaROOT;
            
            TString residualPDF;
            TString residualPNG;
            TString residualROOT;            

            TString photonPDF;
            TString photonPNG;
            TString photonROOT;

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

            // ========================================================
            // NPE
            //
            // Upper: NPE vs eta
            // Lower: NPE - Default
            //
            // ========================================================

            TString npeTitle;
            npeTitle.Form( "%s (%s)", sipmLabel[sipm], rad);
            MakeObservableRatioCanvas( gNPE, "<NPE>", "<NPE>", npeTitle.Data(), 
            	npePDF.Data(), npePNG.Data(), npeROOT.Data(),
                npeRange.ymin, npeRange.ymax, npeDeltaRange.ymin, npeDeltaRange.ymax );

            // ========================================================
            // THETA MEAN
            //
            // Same layout:
            // Upper: Theta mean vs eta
            // Lower: Theta mean - Default
            // ========================================================

            TString thetaTitle;
            thetaTitle.Form( "%s (%s)", sipmLabel[sipm], rad );
            MakeObservableRatioCanvas(gTheta, "<#theta_{C}> [mrad]", "<#theta_{C}> [mrad]", thetaTitle.Data(),
                thetaPDF.Data(), thetaPNG.Data(), thetaROOT.Data(), 
                thetaRange.ymin, thetaRange.ymax, thetaDeltaRange.ymin, thetaDeltaRange.ymax);

            // ========================================================
            // THETA RESIDUAL SIGMA
            //
            // Upper: sigma vs eta
            // Lower: sigma - Default
            // ========================================================

            TString residualTitle;
            residualTitle.Form("%s (%s)", sipmLabel[sipm], rad);
            MakeObservableRatioCanvas( gResidual, "#sigma_{#theta,res} [mrad]", "#sigma_{#theta,res} [mrad]",
                residualTitle.Data(), residualPDF.Data(), residualPNG.Data(), residualROOT.Data(),
                residualRange.ymin, residualRange.ymax, residualDeltaRange.ymin, residualDeltaRange.ymax );

            // ========================================================
            // SINGLE PHOTON RESOLUTION
            //
            // Upper: resolution vs eta
            // Lower: resolution - Default
            // ========================================================

            TString photonTitle;
            photonTitle.Form( "%s (%s)",  sipmLabel[sipm], rad );
            MakeObservableRatioCanvas( gPhoton, "SPR [mrad]", "SPR [mrad]",
                photonTitle.Data(), photonPDF.Data(), photonPNG.Data(), photonROOT.Data(),
                photonRange.ymin, photonRange.ymax, photonDeltaRange.ymin, photonDeltaRange.ymax);

            // ========================================================
            // Cleanup
            // ========================================================

            delete gNPE[0];
            delete gNPE[1];
            delete gNPE[2];
            delete gNPE[3];

            delete gTheta[0];
            delete gTheta[1];
            delete gTheta[2];
            delete gTheta[3];

            delete gResidual[0];
            delete gResidual[1];
            delete gResidual[2];
            delete gResidual[3];

            delete gPhoton[0];
            delete gPhoton[1];
            delete gPhoton[2];
            delete gPhoton[3];

            cout << endl;
            cout << "Output written to: " << outDir << endl;
        }
    }


    cout << endl;
    cout << "=================================================" << endl;
    cout << "FitGasCompare completed." << endl;
    cout << "=================================================" << endl;
}

