#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;

TH1D* GetRun3NchvsCentrality(string FileNameSuffix = "LHC23_PbPb_pass4_341269"){
    TFile* f = new TFile(Form("./AnalysisResults/AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    TH2D* hTrackCorrection2d = (TH2D*)f->Get(Form("flow-task/hTrackCorrection2d"));
    TF1* f1 = new TF1("f1","[0]*x+[1]",0,4000);
    hTrackCorrection2d->Fit(f1,"R");
    // hTrackCorrection2d->Draw("colz");

    TH2D* hglobalTrackvsCentT0C = (TH2D*)f->Get(Form("flow-task/globalTracks_centT0C"));
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
    hNchvsCentrality->SetTitle(Form("Run 3: <N_{ch}> vs Centrality"));
    // hNchvsCentrality->Draw();
    return hNchvsCentrality;
}

void GetRun2NchvsCentrality_Emil() {
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
    hNchvsCentrality->SetTitle("Run 2: <N_{ch}> vs Centrality");
    hNchvsCentrality->Draw();
}

TH1D* GetRun2NchvsCentrality() {
    TFile *f = new TFile("./hMultCent_Run2.root","READ");
    TH2F* hMultCent = (TH2F*)f->Get("hMultCent");
    TH1D* hNchvsCentrality = hMultCent->ProjectionX("hNchvsCentrality");
    Int_t nbinsx = hMultCent->GetNbinsX();
    for (int i=1; i<=nbinsx; i++) {
        double NchSum = 0;
        double NSum = 0;
        for (int j=1; j<=hMultCent->GetNbinsY(); j++) {
            NchSum += hMultCent->GetBinContent(i,j) * hMultCent->GetYaxis()->GetBinCenter(j);
            NSum += hMultCent->GetBinContent(i,j);
        }
        if (NSum > 0) {
            double NchMean = NchSum/NSum;
            hNchvsCentrality->SetBinContent(i,NchMean);
        }
        else {
            hNchvsCentrality->SetBinContent(i,0);
        }
    }
    hNchvsCentrality->SetTitle("Run 2: <N_{ch}> vs Centrality");
    // hNchvsCentrality->Draw();

    return hNchvsCentrality;
}


void GetNchvsCentrality(){
    TH1D* Run3 = GetRun3NchvsCentrality();
    if (!Run3) {
        cout << "Run 3 data not found." << endl;
        return;
    }
    TH1D* Run2 = GetRun2NchvsCentrality();
    if (!Run2) {
        cout << "Run 2 data not found." << endl;
        return;
    }
    Run3->Divide(Run2);
    Run3->SetTitle("Run 3/Run 2: <N_{ch}> vs Centrality");
    Run3->Draw();

    TFile* f = new TFile("./Run3overRun2.root","RECREATE");
    Run3->Write();
    f->Close();
}
