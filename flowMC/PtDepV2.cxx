#include "TFile.h"
#include "TH1D.h"
#include <cmath>

void ProduceCorrectionFactor(double bmin, double bmax, int Centmin, int Centmax){
    TFile* f = new TFile("./AnalysisResults_LHC24k2_David_327409.root","READ");

    double pTAxis[] = {0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f}; 

    TH3D* hBVsPtVsPhiGenerated = (TH3D*)f->Get("flow-test/hBVsPtVsPhiGenerated");
    if (!hBVsPtVsPhiGenerated) {
        std::cerr << "Error: hBVsPtVsPhiGenerated not found in file." << std::endl;
        return;
    }
    // hBVsPtVsPhiGenerated->Draw("colz");
    hBVsPtVsPhiGenerated->GetXaxis()->SetRangeUser(bmin,bmax);
    TH2D* hPtVsPhiGenerated = (TH2D*)hBVsPtVsPhiGenerated->Project3D("yz");
    // hPtVsPhiGenerated->GetXaxis()->SetRangeUser(0,3);
    // hPtVsPhiGenerated->Draw("col");

    TH1D* hV2Generated = new TH1D("hV2Generated","hV2Generated",40,pTAxis);
    // At every B, calculate average phi;
    for (int i=1; i<=hPtVsPhiGenerated->GetNbinsX(); ++i) {
        double pt = hPtVsPhiGenerated->GetXaxis()->GetBinCenter(i);
        double phiSum = 0;
        int nEntries = 0;
        for (int j=1; j<=hPtVsPhiGenerated->GetNbinsY(); ++j) {
            double phi = hPtVsPhiGenerated->GetYaxis()->GetBinCenter(j);
            double num = hPtVsPhiGenerated->GetBinContent(i,j);
            if (num > 0) {
                phiSum += cos(2*phi)*num;
                nEntries += num;
            }
        }
        double v2 = phiSum/nEntries;
        hV2Generated->SetBinContent(hV2Generated->GetXaxis()->FindBin(pt),v2);
    }
    // hV2Generated->Draw();

    TH3D* hBVsPtVsPhiGlobal = (TH3D*)f->Get("flow-test/hBVsPtVsPhiGlobal");
    if (!hBVsPtVsPhiGlobal) {
        std::cerr << "Error: hBVsPtVsPhiGlobal not found in file." << std::endl;
        return;
    }
    // hBVsPtVsPhiGlobal->Draw("colz");
    hBVsPtVsPhiGlobal->GetXaxis()->SetRangeUser(bmin,bmax);
    TH2D* hPtVsPhiGlobal = (TH2D*)hBVsPtVsPhiGlobal->Project3D("yz");
    // hPtVsPhiGlobal->GetXaxis()->SetRangeUser(0,3);
    // hPtVsPhiGlobal->Draw("col");

    TH1D* hV2Global = new TH1D("hV2Global","hV2Global",40,pTAxis);
    // At every B, calculate average phi;
    for (int i=1; i<=hPtVsPhiGlobal->GetNbinsX(); ++i) {
        double pt = hPtVsPhiGlobal->GetXaxis()->GetBinCenter(i);
        double phiSum = 0;
        int nEntries = 0;
        for (int j=1; j<=hPtVsPhiGlobal->GetNbinsY(); ++j) {
            double phi = hPtVsPhiGlobal->GetYaxis()->GetBinCenter(j);
            double num = hPtVsPhiGlobal->GetBinContent(i,j);
            if (num > 0) {
                phiSum += cos(2*phi)*num;
                nEntries += num;
            }
        }
        double v2 = phiSum/nEntries;
        hV2Global->SetBinContent(hV2Global->GetXaxis()->FindBin(pt),v2);
    }
    // hV2Global->Draw();

    hV2Global->Divide(hV2Generated);
    // hV2Global->Draw();
    hV2Global->SetName("hRatio");
    TFile* fOut = new TFile(Form("./MCcorrection/CorrectionFactor_Cent%dTo%d_LHC24k2_David_327409.root",Centmin,Centmax),"RECREATE");
    hV2Global->Write();
    fOut->Close();
    f->Close();
}

void PtDepV2(){
    vector<int> Cent = {0,5,10,20,30,40,50,60,70,80,90,100};
    vector<double> IP = {0, 2.9, 4.3, 6.1, 7.5, 8.7, 9.9, 10.9, 12.1, 13.1, 14.1, 20.};
    if (Cent.size() != IP.size()) {
        std::cerr << "Error: Cent and IP vectors have different sizes." << std::endl;
        return;
    }
    for(int i=0; i<IP.size()-1; i++) {
        ProduceCorrectionFactor(IP[i],IP[i+1],Cent[i],Cent[i+1]);
    }
}
