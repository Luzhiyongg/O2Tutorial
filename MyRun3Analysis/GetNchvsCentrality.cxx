#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;

void GetNchvsCentrality(string FileNameSuffix = "LHC23_PbPb_pass4_328395"){
    TFile* f = new TFile(Form("./AnalysisResults/AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    vector<string> SubwagonNames = {"_kIsGoodITSLayersAll"};
    for(int i=0;i<SubwagonNames.size();i++){
        TH2D* hTrackCorrection2d = (TH2D*)f->Get(Form("flow-task%s/hTrackCorrection2d",SubwagonNames[i].c_str()));
        TF1* f1 = new TF1("f1","[0]*x+[1]",0,4000);
        hTrackCorrection2d->Fit(f1,"R");
        // hTrackCorrection2d->Draw("colz");

        TH2D* hglobalTrackvsCentT0C = (TH2D*)f->Get(Form("flow-task%s/globalTracks_centT0C",SubwagonNames[i].c_str()));
        TH1D* hNchvsCentrality = hglobalTrackvsCentT0C->ProjectionX("hNchvsCentrality");
        // at every X, merge all the bins of the same X together
        Int_t nbinsx = hglobalTrackvsCentT0C->GetNbinsX();
        for(int j=1;j<nbinsx;j++) {
            double NchSum = 0;
            double NSum = 0;
            for(int k=1;k<=hglobalTrackvsCentT0C->GetNbinsY();k++) {
                double uncorrectedTracks = hglobalTrackvsCentT0C->GetYaxis()->GetBinCenter(k);
                double correctedTracks = f1->Eval(uncorrectedTracks);
                NchSum += hglobalTrackvsCentT0C->GetBinContent(j,k) * correctedTracks;
                NSum += hglobalTrackvsCentT0C->GetBinContent(j,k);
            }
            if (NSum > 0) {
                double NchMean = NchSum/NSum;
                hNchvsCentrality->SetBinContent(j,NchMean);
            }
            else {
                hNchvsCentrality->SetBinContent(j,0);
            }
        }
        hNchvsCentrality->SetTitle(Form("Run 3: N_{ch} vs Centrality"));
        hNchvsCentrality->Draw();
    }
}

void GetRun2NchvsCentrality() {
    TFile* f = new TFile("./nch_vs_cent.root","READ");
    TH2D* hCentvsNch = (TH2D*)f->Get("hCentvsNch");
    TH1D* hNchvsCentrality = hCentvsNch->ProjectionY("hNchvsCentrality");
    Int_t nbinsy = hCentvsNch->GetNbinsY();
    for(int i=1;i<=nbinsy;i++) {
        double NchSum = 0;
        double NSum = 0;
        for(int j=1;j<=hCentvsNch->GetNbinsX();j++) {
            NchSum += hCentvsNch->GetBinContent(j,i) * hCentvsNch->GetXaxis()->GetBinCenter(j);
            NSum += hCentvsNch->GetBinContent(j,i);
        }
        if (NSum > 0) {
            double NchMean = NchSum/NSum;
            hNchvsCentrality->SetBinContent(i,NchMean);
        }
        else {
            hNchvsCentrality->SetBinContent(i,0);
        }
    }
    hNchvsCentrality->SetTitle("Run 2: N_{ch} vs Centrality");
    hNchvsCentrality->Draw();
}
