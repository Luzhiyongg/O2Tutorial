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
#include "TLatex.h"
#include <cstring>
#include <vector>
#include <map>
#include <array>
#include "../MyRun3Analysis/include/FlowContainerCalculation.h"
#include "./include/RunByRunCommon.h"

std::vector<std::string> histNamesList = {"hPhiWeighted","hEta","hVtxZ","hMult","hCent","c22_gap10"};

void DrawRunByRunSubwagon() {
    // Compare different subwagons for each run
    TFile *file = TFile::Open("./AnalysisResults/AnalysisResults_LHC23zzh_pass4_361386.root");
    //LHC23zzm_pass4_calo_361088
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
    // std::vector<int> IgnoreRuns = {544653,544742,544091}; 

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
        // ========================
        // 1. 主标题（画布顶部）
        // ========================
        // 创建标题专用Pad（占据画布顶部5%空间）
        TPad* titlePad = new TPad("titlePad", "", 0, 0.95, 1, 1); // (x1,y1,x2,y2) 单位归一化
        titlePad->Draw();
        titlePad->cd(); // 进入标题Pad

        // 在标题Pad中央绘制主标题
        TLatex* mainTitle = new TLatex(0.5, 0.5, Form("Run %d", runNumber));
        mainTitle->SetTextAlign(22);   // 中心对齐
        mainTitle->SetTextSize(0.6);   // 文字大小（相对Pad高度）
        mainTitle->SetTextColor(kBlack); // 可自定义颜色
        mainTitle->Draw();

        // ========================
        // 2. 子图区域（带分割）
        // ========================
        canvas->cd(); // 返回主画布
        TPad* mainContent = new TPad("mainContent", "", 0, 0, 1, 0.95); // 占据剩余95%空间
        mainContent->Draw();
        mainContent->Divide(nCols, nRows); // 分割为子图

        // turn off statistics box
        gStyle->SetOptStat(0);


        for(int iHist = 0; iHist < nHists; iHist++) {
            std::string histName = histNamesList[iHist];
            mainContent->cd(iHist+1);
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
