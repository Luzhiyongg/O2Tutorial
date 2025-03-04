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
#include "./include/RunByRunCommon.h"

using namespace std;

vector<string> histNamesList = {"hPhi","hEta","hVtxZ","hMult","hCent","c22_gap10","c3232"};

void DrawQA(){
    // TFile *file = TFile::Open("./AnalysisResults/AnalysisResults_LHC23zzo_pass4_QC1_sampling_359365.root");
    TFile *file = TFile::Open("./AnalysisResults/AnalysisResults_LHC23_PbPb_pass4_356235.root");
    if (!file) {
        std::cout << "Cannot open file" << std::endl;
        return;
    }
    gDirectory->Cd("flow-runby-run");
    // get the list of runs
    TList *runsList = gDirectory->GetListOfKeys();
    int nRuns = runsList->GetEntries(); // number of runs
    TDirectory* HadronicRate = (TDirectory*)gDirectory->Get("HadronicRate");
    if (!HadronicRate) {
        std::cout << "Cannot find HadronicRate directory" << std::endl;
        return;
    }
    gStyle->SetOptStat(0); // turn off statistics box
    
    for(string histName : histNamesList){
        // create a canvas
        TCanvas *canvas = new TCanvas(Form("canvas_%s", histName.c_str()), Form("canvas_%s", histName.c_str()), 800, 800);
        vector<TLegend*> legends;
        legends.push_back(new TLegend());
        bool isFirst = false;
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
            TH1D *hist = (TH1D*)gDirectory->Get(Form("%d/%s", runNumber, histName.c_str()));
            if (!hist) {
                Printf("Histogram %s not found for run %d", histName.c_str(), runNumber);
                continue;
            }
            // normalize the histogram
            if(histName != "c22_gap10" && hist->Integral() > 0) hist->Scale(1.0 / hist->Integral());

            double meanHadronicRate = GetHadronicRate(runNumber, HadronicRate);
            if (cutHadronicRate && (meanHadronicRate < minHadronicRate || meanHadronicRate > maxHadronicRate)) {
                Printf("Skipping run %d, hadronic rate %0.1f kHz", runNumber, meanHadronicRate);
                continue;
            }

            // set the color, marker style, and line style
            hist->SetLineStyle(lineStylesList[iLineStyle]);
            hist->SetMarkerStyle(markerStylesList[iMarkerStyle]);
            hist->SetMarkerColor(colorsList[iColor]);
            hist->SetLineColor(colorsList[iColor]);
            legends[iLegend]->AddEntry(hist, Form("%d, %0.1f kHz", runNumber, meanHadronicRate), "p");
            // draw the histogram
            if(!isFirst){
                hist->SetTitle(Form("%s", histName.c_str()));
                hist->Draw("HIST");
                hist->GetYaxis()->SetRangeUser(0, 0.1);
                hist->Draw("p");
                isFirst = true;
            }
            else{
                hist->Draw("p same");
            }

            if (SelectedRuns.empty() && !cutHadronicRate) {
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
