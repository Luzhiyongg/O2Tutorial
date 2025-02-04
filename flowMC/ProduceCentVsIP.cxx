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
        fCentIP->SetBinContent(i,100*centrality);
    }
    return fCentIP;
}

void ProduceCentVsIP(){
    TFile* f = new TFile("./AnalysisResults_LHC24k2_David_327409.root","READ");
    TH1D* fCentIP = GetCentVsIP(f,"flow-test/hImpactParameter");

    if(!fCentIP){
        Printf("can not get flow-test/hImpactParameter");
        return;
    }


    TH1D* hCentVsIP = (TH1D*) fCentIP->Clone();
    hCentVsIP->SetName("ccdb_object");
    hCentVsIP->SetTitle("Centrality vs impact parameter");
    hCentVsIP->GetYaxis()->SetTitle("Centrality (%)");
    hCentVsIP->GetXaxis()->SetTitle("Impact parameter (fm)");

    // Eff->Draw();

    TFile* fout = new TFile("./CentVsIP/CentVsIP_LHC24k2_David_327409.root","RECREATE");
    hCentVsIP->Write();
    fout->Close();
}
