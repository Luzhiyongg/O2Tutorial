//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"

// #include "FlowContainer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "GFWWeights.h"


void ProduceAcceptanceWeight(){
    TFile* f = new TFile("./AnalysisResults_LHC23zzh_pass4_260874.root","READ");
    GFWWeights* W = (GFWWeights*)f->Get("flow-task/weights");
    if(!W){
        Printf("can not getflow-pb-pb-task/weights");
        return;
    }


    GFWWeights* Weights = (GFWWeights*) W->Clone();
    Weights->SetName("ccdb_object");
    Weights->SetTitle("NUA");

    // Eff->Draw();

    TFile* fout = new TFile("./Acceptance/NUA_LHC23zzh_pass4_260874.root","RECREATE");
    Weights->Write();
    fout->Close();
    
}
