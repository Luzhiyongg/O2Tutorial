//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

// #include "FlowContainer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TGraphAsymmErrors.h"
#include <cstring>
#include <vector>
#include <map>
#include <array>
#include "../MyRun3Analysis/include/FlowContainerCalculation.h"

using namespace std;

vector<string> cumulnatList = {"c22","c22_gap10","c32","c32_gap10"};

vector<Color_t> colorsList = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kOrange, kViolet, kPink};

vector<EMarkerStyle> markerStylesList = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullDiamond, kFullCross, kFullCrossX, kFullThreeTriangles,
 kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown, kOpenStar, kOpenDiamond, kOpenCross, kOpenCrossX, kOpenThreeTriangles};

vector<ELineStyle> lineStylesList = {kSolid, kDashed, kDotted};

double MySqrt(double x) {
    if (x <= 0)
        return -1;
    else
        return sqrt(x);
}

double ErrorSqrt(double x, double ex) {
    if (x <= 0)
        return 10.;
    else
        return abs(ex * 0.5 / sqrt(x));
}

TH1D* FromCumulantToVn(TH1D* cumulant) {
    TH1D* vn = new TH1D(Form("%s_vn", cumulant->GetName()), Form("%s_vn", cumulant->GetName()), cumulant->GetNbinsX(), cumulant->GetXaxis()->GetXmin(), cumulant->GetXaxis()->GetXmax());
    for (int iBin = 1; iBin <= cumulant->GetNbinsX(); iBin++) {
        double cumulantValue = cumulant->GetBinContent(iBin);
        double vnValue = MySqrt(cumulantValue);
        double vnError = ErrorSqrt(cumulantValue, cumulant->GetBinError(iBin));
        // Printf("BincEnter: %f, cumulantValue: %f, vnValue: %f, vnError: %f", cumulant->GetBinCenter(iBin), cumulantValue, vnValue, vnError);
        vn->SetBinContent(iBin, vnValue);
        vn->SetBinError(iBin, vnError);
    }
    return vn;
}

TH1D* RunbyRunCalculateNSC32(TH1D* hCorr22_Full, TH1D* hCorr22_Gap, TH1D* hCorr32_Full, TH1D* hCorr32_Gap, TH1D* hCorr3232){
    TH1D* target = new TH1D(Form("%s",hCorr3232->GetName()),Form("%s",hCorr3232->GetName()),hCorr3232->GetNbinsX(),hCorr3232->GetXaxis()->GetXmin(),hCorr3232->GetXaxis()->GetXmax());
    for(int i=1;i<=hCorr3232->GetNbinsX();i++){
        target->SetBinContent(i,
        (hCorr3232->GetBinContent(i)-hCorr22_Full->GetBinContent(i)*hCorr32_Full->GetBinContent(i))
        /(hCorr22_Gap->GetBinContent(i)*hCorr32_Gap->GetBinContent(i))
        );
        double err = Error_NSC(hCorr3232->GetBinContent(i),hCorr3232->GetBinError(i),
        hCorr22_Full->GetBinContent(i),hCorr22_Full->GetBinError(i),
        hCorr32_Full->GetBinContent(i),hCorr32_Full->GetBinError(i),
        hCorr22_Gap->GetBinContent(i),hCorr22_Gap->GetBinError(i),
        hCorr32_Gap->GetBinContent(i),hCorr32_Gap->GetBinError(i)
        );
        target->SetBinError(i,err);
    }
    return target;
}

void DrawRunByRunFlow(){
    TFile *file = TFile::Open("./AnalysisResults/AnalysisResults_LHC23zzo_pass4_QC1_sampling_359365.root");
    if (!file) {
        std::cout << "Cannot open file" << std::endl;
        return;
    }
    gDirectory->Cd("flow-runby-run");
    // get the list of runs
    TList *runsList = gDirectory->GetListOfKeys();
    TDirectory* runDir = gDirectory;
    int nRuns = runsList->GetEntries(); // number of runs
    // vector<int> IgnoreRuns = {544640, 544913, 545117, 545295, 545311, 544091, 544652, 544694}; 
    vector<int> IgnoreRuns = {544116,544122,544653,544742}; // LHCzzh
    vector<int> SelectedRuns = {545367, 545291, 545223, 544917}; // Good ITS Runs
    
    TFile* publish_Run2pass2 = new TFile("./PublicData/PbPb_Run2Pass2.root","READ");
    if (!publish_Run2pass2) {
        std::cout << "Cannot open file" << std::endl;
        return;
    }
    TGraphAsymmErrors* g_v2 = (TGraphAsymmErrors*)publish_Run2pass2->Get("v2{2}_Gap10_TPCPileUp");
    SetMarkerAndLine(g_v2,kBlack,kOpenSquare,kSolid,1.0);
    TGraphAsymmErrors* g_v32 = (TGraphAsymmErrors*)publish_Run2pass2->Get("v3{2}_Gap10_TPCPileUp");
    SetMarkerAndLine(g_v32,kBlack,kOpenSquare,kSolid,1.0);

    // Vn
    for(string histName : cumulnatList){
        // create a canvas
        TCanvas *canvas = new TCanvas(Form("canvas_%s", histName.c_str()), Form("canvas_%s", histName.c_str()), 800, 600);
        vector<TLegend*> legends;
        legends.push_back(new TLegend());
        int iColor = 0;
        int iMarkerStyle = 0;
        int iLineStyle = 0;
        int iLegend = 0;
        // loop over the runs
        for (int iRun = 0; iRun < nRuns; iRun++) {
            // get the run number
            TString runName = runsList->At(iRun)->GetName();
            int runNumber = runName.Atoi();
            if (!(runNumber > 0))
                continue;
            if (find(IgnoreRuns.begin(), IgnoreRuns.end(), runNumber) != IgnoreRuns.end()) {
                Printf("Skipping run %d", runNumber);
                continue;
            }
            if (!SelectedRuns.empty()) {
                if (find(SelectedRuns.begin(), SelectedRuns.end(), runNumber) == SelectedRuns.end()) {
                    Printf("Skipping run %d", runNumber);
                    continue;
                }
            } 
            // Printf("Processing run %d", runNumber);
            // printf("%d, ", runNumber);
            // get the histogram
            TH1D *cumulant = (TH1D*)runDir->Get(Form("%d/%s", runNumber, histName.c_str()));
            if (!cumulant) {
                Printf("Histogram %s not found for run %d", histName.c_str(), runNumber);
                continue;
            }
            
            TH1D* hist = mergeCentralityToTargetBin(FromCumulantToVn(cumulant),targetCentralityBins);

            // set the color, marker style, and line style
            hist->SetLineStyle(lineStylesList[iLineStyle]);
            hist->SetMarkerStyle(markerStylesList[iMarkerStyle]);
            hist->SetMarkerColor(colorsList[iColor]);
            hist->SetLineColor(colorsList[iColor]);
            legends[iLegend]->AddEntry(hist, Form("%d", runNumber), "p");
            // draw the histogram
            if(iRun == 0){
                hist->SetTitle(Form("%s", hist->GetName()));
                hist->Draw("HIST");
                hist->GetYaxis()->SetRangeUser(0, 0.05);
                hist->Draw("p");
            }
            else{
                hist->Draw("p same");
            }

            if (SelectedRuns.empty()) {
                iMarkerStyle++;
                if (iMarkerStyle == markerStylesList.size()) {
                    iMarkerStyle = 0;
                    iColor++;
                    if (iColor == colorsList.size()) iColor = 0;
                }
            }
            else {
                // for small number of runs, use different color first
                iColor++;
                if (iColor == colorsList.size()){
                    iColor = 0;
                    iMarkerStyle++;
                    if (iMarkerStyle == markerStylesList.size()) iMarkerStyle = 0;
                }
            }

            if (iRun % 36 ==0 && iRun!= 0){
                iLegend++;
                legends.push_back(new TLegend());
            }
        }
        if (histName.find("c22") != string::npos || histName.find("c22_gap10") != string::npos) {
            g_v2->Draw("p same");
            legends[0]->AddEntry(g_v2, "Pb-Pb Run2 Pass2", "p");
        }

        if (histName.find("c32") != string::npos || histName.find("c32_gap10") != string::npos) {
            g_v32->Draw("p same");
            legends[0]->AddEntry(g_v32, "Pb-Pb Run2 Pass2", "p");
        }
        // add the legends to the canvas
        for(TLegend* legend : legends){
            legend->Draw();
        }
    }

    // SC
    {
        TCanvas *canvas = new TCanvas(Form("canvas_NSC"), Form("canvas_NSC"), 800, 600);
        vector<TLegend*> legends;
        legends.push_back(new TLegend());
        int iColor = 0;
        int iMarkerStyle = 0;
        int iLineStyle = 0;
        int iLegend = 0;
        // loop over the runs
        for (int iRun = 0; iRun < nRuns; iRun++) {
            // get the run number
            TString runName = runsList->At(iRun)->GetName();
            int runNumber = runName.Atoi();
            if (!(runNumber > 0))
                continue;
            if (find(IgnoreRuns.begin(), IgnoreRuns.end(), runNumber) != IgnoreRuns.end()) {
                Printf("Skipping run %d", runNumber);
                continue;
            }
            if (!SelectedRuns.empty()) {
                if (find(SelectedRuns.begin(), SelectedRuns.end(), runNumber) == SelectedRuns.end()) {
                    Printf("Skipping run %d", runNumber);
                    continue;
                }
            } 

            TH1D *hCorr22_Full = (TH1D*)runDir->Get(Form("%d/c22", runNumber));
            TH1D *hCorr22_Gap = (TH1D*)runDir->Get(Form("%d/c22_gap10", runNumber));
            TH1D *hCorr32_Full = (TH1D*)runDir->Get(Form("%d/c32", runNumber));
            TH1D *hCorr32_Gap = (TH1D*)runDir->Get(Form("%d/c32_gap10", runNumber));
            TH1D *hCorr3232 = (TH1D*)runDir->Get(Form("%d/c3232", runNumber));
            if(!hCorr22_Full ||!hCorr22_Gap ||!hCorr32_Full ||!hCorr32_Gap ||!hCorr3232){
                Printf("Histogram not found for run %d", runNumber);
                return;
            }

            // hCorr22_Full = mergeCentralityToTargetBin(hCorr22_Full, targetCentralityBins);
            // hCorr22_Gap = mergeCentralityToTargetBin(hCorr22_Gap, targetCentralityBins);
            // hCorr32_Full = mergeCentralityToTargetBin(hCorr32_Full, targetCentralityBins);
            // hCorr32_Gap = mergeCentralityToTargetBin(hCorr32_Gap, targetCentralityBins);
            // hCorr3232 = mergeCentralityToTargetBin(hCorr3232, targetCentralityBins);

            TH1D* hist = RunbyRunCalculateNSC32(hCorr22_Full, hCorr22_Gap, hCorr32_Full, hCorr32_Gap, hCorr3232);
            hist = mergeCentralityToTargetBin(hist, targetCentralityBins);

            // set the color, marker style, and line style
            hist->SetLineStyle(lineStylesList[iLineStyle]);
            hist->SetMarkerStyle(markerStylesList[iMarkerStyle]);
            hist->SetMarkerColor(colorsList[iColor]);
            hist->SetLineColor(colorsList[iColor]);
            legends[iLegend]->AddEntry(hist, Form("%d", runNumber), "p");
            // draw the histogram
            if(iRun == 0){
                hist->SetTitle(Form("%s", hist->GetName()));
                hist->Draw("HIST");
                hist->GetYaxis()->SetRangeUser(0, 0.05);
                hist->Draw("p");
            }
            else{
                hist->Draw("p same");
            }

            if (SelectedRuns.empty()) {
                iMarkerStyle++;
                if (iMarkerStyle == markerStylesList.size()) {
                    iMarkerStyle = 0;
                    iColor++;
                    if (iColor == colorsList.size()) iColor = 0;
                }
            }
            else {
                // for small number of runs, use different color first
                iColor++;
                if (iColor == colorsList.size()){
                    iColor = 0;
                    iMarkerStyle++;
                    if (iMarkerStyle == markerStylesList.size()) iMarkerStyle = 0;
                }
            }
            

            if (iRun % 36 ==0 && iRun!= 0){
                iLegend++;
                legends.push_back(new TLegend());
            }
        }
        // add the legends to the canvas
        for(TLegend* legend : legends){
            legend->Draw();
        }
    }
    
}
