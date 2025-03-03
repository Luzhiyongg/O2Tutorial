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

std::vector<std::string> histNamesList = {"hPhiWeighted","hEta","hVtxZ","hMult","hCent","c22_gap10"};

std::vector<Color_t> colorsList = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kOrange, kViolet, kPink};

std::vector<EMarkerStyle> markerStylesList = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullDiamond, kFullCross, kFullCrossX, kFullThreeTriangles,
 kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown, kOpenStar, kOpenDiamond, kOpenCross, kOpenCrossX, kOpenThreeTriangles};

std::vector<ELineStyle> lineStylesList = {kSolid, kDashed, kDotted};

void DrawRunByRunSubwagon() {
    // Compare different subwagons for each run
    TFile *file = TFile::Open("./AnalysisResults/AnalysisResults_LHC23zzh_pass4_361386.root");
    if (!file) {
        std::cout << "Cannot open file" << std::endl;
        return;
    }
    std::vector<std::string> SubwagonNames = {"_HighOcc","_LowOcc"};
    std::vector<std::string> SubwagonExplanations = {"3k < Occupancy","Occupancy < 3k"};
    std::vector<TDirectory*> SubwagonDirs;
    for (auto& SubwagonName : SubwagonNames) {
        SubwagonDirs.push_back((TDirectory*)file->Get(Form("flow-runby-run%s", SubwagonName.c_str())));
        // SubwagonDirs.back()->ls();
    }
    TList *runsList = SubwagonDirs[0]->GetListOfKeys(); // list of runs
    int nRuns = runsList->GetEntries(); // number of runs
    std::vector<int> IgnoreRuns = {544653,544742,544091}; 

    for (int iRun = 0; iRun < nRuns; iRun++) {
        TString runName = runsList->At(iRun)->GetName();
        int runNumber = runName.Atoi();
        if (!(runNumber > 0))
            continue;
        if (find(IgnoreRuns.begin(), IgnoreRuns.end(), runNumber) != IgnoreRuns.end()) {
            Printf("Skipping run %d", runNumber);
            continue;
        }


        int nHists = histNamesList.size();
        int nCols = 3;
        int nRows = nHists / nCols + (nHists % nCols > 0);
        //set canvas size based on the number of histNamesList
        TCanvas *canvas = new TCanvas(Form("c%d", runNumber), Form("Run %d", runNumber), 1000., 250.*nRows);
        //Divide canvas based on the number of histNamesList
        canvas->Divide(nCols, nRows);

        // turn off statistics box
        gStyle->SetOptStat(0);


        for(int iHist = 0; iHist < nHists; iHist++) {
            std::string histName = histNamesList[iHist];
            canvas->cd(iHist+1);
            TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9);
            bool isFirst = false;
            int iColor = 0;
            int iMarkerStyle = 0;
            int iLineStyle = 0;
            for (int iSubwagon = 0; iSubwagon < SubwagonDirs.size(); iSubwagon++) {
                auto SubwagonDir = SubwagonDirs[iSubwagon];
                TH1D *hist = (TH1D*)SubwagonDir->Get(Form("%s/%s", runName.Data(), histName.c_str()));
                if (!hist) {
                    Printf("Histogram %s not found for run %d", histName.c_str(), runNumber);
                    continue;
                }
                hist->SetLineStyle(lineStylesList[iLineStyle]);
                hist->SetMarkerStyle(markerStylesList[iMarkerStyle]);
                hist->SetMarkerColor(colorsList[iColor]);
                hist->SetLineColor(colorsList[iColor]);
                legend->AddEntry(hist, SubwagonExplanations[iSubwagon].c_str(), "lp");
                if (!isFirst) {
                    hist->SetTitle(Form("%s", histName.c_str()));
                    hist->Draw("HIST");
                    // hist->GetYaxis()->SetRangeUser(0, 0.05);
                    hist->Draw("p");
                    isFirst = true;
                }
                else {
                    hist->Draw("p same");
                }
                iMarkerStyle++;
                iLineStyle++;
                iColor++;
            }
            legend->Draw();
        }
    }

}
