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

using namespace std;

vector<string> histNamesList = {"hPhi","hEta","hVtxZ","hMult","hCent","c22_gap10"};

vector<Color_t> colorsList = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kOrange, kViolet, kPink};

vector<EMarkerStyle> markerStylesList = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullDiamond, kFullCross, kFullCrossX, kFullThreeTriangles,
 kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown, kOpenStar, kOpenDiamond, kOpenCross, kOpenCrossX, kOpenThreeTriangles};

vector<ELineStyle> lineStylesList = {kSolid, kDashed, kDotted};

void DrawQA(){
    TFile *file = TFile::Open("./AnalysisResults/AnalysisResults_LHC23_PbPb_pass4_288089.root");
    if (!file) {
        std::cout << "Cannot open file" << std::endl;
        return;
    }
    gDirectory->Cd("flow-runby-run");
    // get the list of runs
    TList *runsList = gDirectory->GetListOfKeys();
    int nRuns = runsList->GetEntries(); // number of runs
    vector<int> IgnoreRuns = {544640, 544913, 545117, 545295, 545311, 544091, 544652, 544694}; 
    
    for(string histName : histNamesList){
        // create a canvas
        TCanvas *canvas = new TCanvas(Form("canvas_%s", histName.c_str()), Form("canvas_%s", histName.c_str()), 800, 800);
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
            if (find(IgnoreRuns.begin(), IgnoreRuns.end(), runNumber) != IgnoreRuns.end()) {
                Printf("Skipping run %d", runNumber);
                continue;
            }
            Printf("Processing run %d", runNumber);
            // get the histogram
            TH1D *hist = (TH1D*)gDirectory->Get(Form("%d/%s", runNumber, histName.c_str()));
            if (!hist) {
                Printf("Histogram %s not found for run %d", histName.c_str(), runNumber);
                continue;
            }
            // normalize the histogram
            if(histName != "c22_gap10" && hist->Integral() > 0) hist->Scale(1.0 / hist->Integral());

            // set the color, marker style, and line style
            hist->SetLineStyle(lineStylesList[iLineStyle]);
            hist->SetMarkerStyle(markerStylesList[iMarkerStyle]);
            hist->SetMarkerColor(colorsList[iColor]);
            hist->SetLineColor(colorsList[iColor]);
            legends[iLegend]->AddEntry(hist, Form("%d", runNumber), "p");
            // draw the histogram
            if(iRun == 0){
                hist->SetTitle(Form("%s", histName.c_str()));
                hist->Draw("HIST");
                hist->GetYaxis()->SetRangeUser(0, 0.05);
                hist->Draw("p");
            }
            else{
                hist->Draw("p same");
            }

            iMarkerStyle++;
            if (iMarkerStyle == markerStylesList.size()) {
                iMarkerStyle = 0;
                iColor++;
                if (iColor == colorsList.size()) iColor = 0;
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
