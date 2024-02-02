//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"

// #include "FlowContainer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TProfile.h"


void ProduceEfficiency(){
    TFile* f = new TFile("./AnalysisResults_MC_LHC23k6d.root","READ");
    TH1D* Rec = (TH1D*)f->Get("flow-pt-efficiency/hPtMCRec");
    TH1D* Truth = (TH1D*)f->Get("flow-pt-efficiency/hPtMCGen");
    if(!Rec){
        Printf("can not get flow-pt-efficiency/hPtMCRec");
        return;
    }
    if(!Truth){
        Printf("can not get flow-pt-efficiency/hPtMCGen");
        return;
    }

    TH1D* Eff = (TH1D*) Rec->Clone();
    Eff->Divide(Truth);
    Eff->SetName("ccdb_object");
    Eff->SetTitle("Efficiency");

    // Eff->Draw();

    TFile* fout = new TFile("./Efficiency/Eff_LHC23k6d.root","RECREATE");
    Eff->Write();
    fout->Close();
    
}
