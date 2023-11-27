//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"

#include "FlowContainer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TProfile.h"

void ProcessFlowContainer(){
    TFile* f = new TFile("AnalysisResults.root","READ");
    FlowContainer* fc = (FlowContainer*)f->Get("gfw-tutorial/FlowContainer");
    // fc->SetIDName("ChGap");
    // fc->SetPropagateErrors(kTRUE);
    // TH1D* hV22 = (TH1D*)fc->GetVN2VsMulti(2);
    // TH1D* hV22 = (TH1D*)fc->GetHistCorrXXVsMulti("22");
    // if(!hV22){
    //     Printf("Can't get hV22");
    //     return;
    // }
    // hV22->SetMarkerSize(1.0);
    // hV22->SetMarkerStyle(kFullCircle);
    // hV22->SetMarkerColor(kBlack);
    // hV22->SetLineColor(kBlack);
    // hV22->SetLineStyle(kSolid);
    // Printf("Error: %f",hV22->GetBinError(2));
    // hV22->Draw("E1");
    
    // TProfile2D* profiles = (TProfile2D*)fc->GetProfile();
    // profiles->Draw();
    // Printf("fIDNname: %s",fc->fIDName.Data());

    // TProfile* pV22 = GetCorrXXVsMulti(profiles,fc->fIDName,"22");
    // if(!pV22){
    //     Printf("Can't get pV22");
    //     return;
    // }
    // pV22->Draw();
    // GetHistCorrXXVsPt
    // GetVN2VsPt

    // TH1D* hV22pt = (TH1D*)fc->GetVN2VsPt(2,0,4.9);
    // if(!hV22pt){
    //     Printf("Can't get hV22");
    //     return;
    // }
    // hV22pt->Draw();

    fc->SetIDName("Ch04GapA");
    // fc->SetPropagateErrors(kTRUE);
    TH1D* hV22 = (TH1D*)fc->GetHistCorrXXVsMulti("422");
    if(!hV22){
        Printf("Can't get hV22");
        return;
    }
    hV22->SetMarkerSize(1.0);
    hV22->SetMarkerStyle(kFullCircle);
    hV22->SetMarkerColor(kBlack);
    hV22->SetLineColor(kBlack);
    hV22->SetLineStyle(kSolid);
    Printf("Error: %f",hV22->GetBinError(2));
    hV22->Draw("E1");
    
}