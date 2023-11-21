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
    fc->SetIDName("ChGap");
    // fc->SetPropagateErrors(kTRUE);
    TH1D* hV22 = (TH1D*)fc->GetVN2VsMulti(2);
    // TH1D* hV22 = (TH1D*)fc->GetHistCorrXXVsMulti("ChGap22");
    if(!hV22){
        Printf("Can't get hV22");
        return;
    }
    hV22->Draw();
    
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
    
}