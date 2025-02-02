#include "TFile.h"
#include "TH1D.h"
#include <cmath>

TH1D* GetCentVsIP(TFile* f, const char* IPdistname) {
    TH1D* fIPdist = (TH1D*)f->Get(IPdistname);
    if (!fIPdist) {
        std::cerr << "Error: " << IPdistname << " not found in file." << std::endl;
        return nullptr;
    }
    TH1D* fNormalised = reinterpret_cast<TH1D*>(fIPdist->Clone("fNormalised"));
    fNormalised->Sumw2(kTRUE);
    fNormalised->Scale(1./fNormalised->Integral());
    TH1D* fCentIP = reinterpret_cast<TH1D*>(fIPdist->Clone("fCentrality"));
    for(auto i(1);i<=fIPdist->GetNbinsX();++i){
        double centrality = fNormalised->Integral(1,i);
        fCentIP->SetBinContent(i,centrality);
    }
    return fCentIP;
}

double GetCentrality(const TH1D* fCentIP, double b) {
    double cent = fCentIP->GetBinContent(fCentIP->GetXaxis()->FindBin(b));
    return 100*cent;
}

void CentDepV2(){
    TFile* f = new TFile("./AnalysisResults_LHC24k2_David_327409.root","READ");
    TH1D* fCentIP = GetCentVsIP(f,"flow-test/hImpactParameter");

    TH3D* hBVsPtVsPhiGenerated = (TH3D*)f->Get("flow-test/hBVsPtVsPhiGenerated");
    if (!hBVsPtVsPhiGenerated) {
        std::cerr << "Error: hBVsPtVsPhiGenerated not found in file." << std::endl;
        return;
    }
    // hBVsPtVsPhiGenerated->Draw("colz");
    double ptmin = 0.2;
    double ptmax = 3;
    hBVsPtVsPhiGenerated->GetZaxis()->SetRangeUser(ptmin,ptmax);
    TH2D* hBVsPhiGenerated = (TH2D*)hBVsPtVsPhiGenerated->Project3D("xy");
    // hBVsPhiGenerated->GetXaxis()->SetRangeUser(0,3);
    // hBVsPhiGenerated->Draw("col");

    TH1D* hV2Generated = new TH1D("hV2Generated","hV2Generated",10,0.,100.);
    // At every B, calculate average phi;
    for (int i=1; i<=hBVsPhiGenerated->GetNbinsY(); ++i) {
        double b = hBVsPhiGenerated->GetYaxis()->GetBinCenter(i);
        double phiSum = 0;
        int nEntries = 0;
        for (int j=1; j<=hBVsPhiGenerated->GetNbinsX(); ++j) {
            double phi = hBVsPhiGenerated->GetXaxis()->GetBinCenter(j);
            double num = hBVsPhiGenerated->GetBinContent(j,i);
            if (num > 0) {
                phiSum += cos(2*phi)*num;
                nEntries += num;
            }
        }
        double v2 = phiSum/nEntries;
        double cent = GetCentrality(fCentIP,b);
        hV2Generated->SetBinContent(hV2Generated->GetXaxis()->FindBin(cent),v2);
    }
    // hV2Generated->Draw();

    
    TH3D* hBVsPtVsPhiGlobal = (TH3D*)f->Get("flow-test/hBVsPtVsPhiGlobal");
    if (!hBVsPtVsPhiGlobal) {
        std::cerr << "Error: hBVsPtVsPhiGlobal not found in file." << std::endl;
        return;
    }
    // hBVsPtVsPhiGlobal->Draw("colz");
    hBVsPtVsPhiGlobal->GetZaxis()->SetRangeUser(ptmin,ptmax);
    TH2D* hBVsPhiGlobal = (TH2D*)hBVsPtVsPhiGlobal->Project3D("xy");
    // hBVsPhiGlobal->GetXaxis()->SetRangeUser(0,3);
    // hBVsPhiGlobal->Draw("col");

    TH1D* hV2Global = new TH1D("hV2Global","hV2Global",10,0.,100.);
    // At every B, calculate average phi;
    for (int i=1; i<=hBVsPhiGlobal->GetNbinsY(); ++i) {
        double b = hBVsPhiGlobal->GetYaxis()->GetBinCenter(i);
        double phiSum = 0;
        int nEntries = 0;
        for (int j=1; j<=hBVsPhiGlobal->GetNbinsX(); ++j) {
            double phi = hBVsPhiGlobal->GetXaxis()->GetBinCenter(j);
            double num = hBVsPhiGlobal->GetBinContent(j,i);
            if (num > 0) {
                phiSum += cos(2*phi)*num;
                nEntries += num;
            }
        }
        double v2 = phiSum/nEntries;
        double cent = GetCentrality(fCentIP,b);
        hV2Global->SetBinContent(hV2Global->GetXaxis()->FindBin(cent),v2);
    }
    // hV2Global->SetLineColor(kRed);
    // hV2Global->Draw("SAMES");

    hV2Global->Divide(hV2Generated);
    hV2Global->Draw();

}
