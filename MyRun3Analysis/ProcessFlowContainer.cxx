//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"

#include "FlowContainer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"

void SetMarkerAndLine(TH1D* graph, Int_t Color=0, Int_t style=0, Int_t linestyle=0, Float_t size=1){
    if(!graph)return;
    if(Color){
        graph->SetMarkerColor(Color);
        graph->SetLineColor(Color);
    }
    if(style){
        graph->SetMarkerStyle(style);
    }
    if(linestyle){
        graph->SetLineStyle(kSolid);
    }
    if(size){
        graph->SetMarkerSize(1.0);
    }
}

void ProcessFlowContainer(){
    TFile* f = new TFile("AnalysisResults_Hyperloop.root","READ");
    FlowContainer* fc = (FlowContainer*)f->Get("flow-pb-pb-task/FlowContainer");
    if(!fc){
        Printf("can not get flow container");
        return;
    }

    TCanvas* canvas1 = new TCanvas("canvas1","canvas1",900,900);
    TH1D* Hist  = new TH1D("Hist","v_{n} in LHC22s_pass5",8,0,80);
    Hist->SetMinimum(0.);
    Hist->SetMaximum(0.15);
    Hist->SetXTitle("Centrality/%");
    Hist->SetYTitle("v_{n}");
    Hist->Draw();
    fc->SetIDName("Ch10Gap");
    fc->SetPropagateErrors(kTRUE);
    // TH1D* hV22 = (TH1D*)fc->GetVN2VsMulti(4);
    // // TH1D* hV22 = (TH1D*)fc->GetHistCorrXXVsMulti("22");
    // if(!hV22){
    //     Printf("Can't get hV22");
    //     return;
    // }
    // SetMarkerAndLine(hV22,kBlack,kFullCircle,kSolid,1.0);
    // hV22->Draw("ESames");
    TH1D* hVn[3] = {nullptr};
    for(int i=0;i<3;i++){
        hVn[i] = (TH1D*)fc->GetVN2VsMulti(i+2);
        if(!hVn[i]){
            Printf("Can't get v%d",i+2);
            return;
        }
        hVn[i]->Draw("ESames");
    }
    SetMarkerAndLine(hVn[0],kBlack,kFullCircle,kSolid,1.0);
    SetMarkerAndLine(hVn[1],kRed,kFullCircle,kSolid,1.0);
    SetMarkerAndLine(hVn[2],kBlue,kFullCircle,kSolid,1.0);
    TLegend* legend = new TLegend(0.2,0.7,0.5,0.9);
    for(int i=0;i<3;i++)legend->AddEntry(hVn[i],Form("v_{%d}{2} |#Delta#eta|>1",i+2));
    legend->Draw();
    
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

    TCanvas* canvas2 = new TCanvas("canvas2","canvas3",900,900);
    TH1D* Hist2  = new TH1D("Hist2","v_{n}(p_{T}) in LHC22s_pass5",10,0.2,3.0);
    Hist2->SetMinimum(0.);
    Hist2->SetMaximum(0.1);
    Hist2->SetXTitle("p_{T}");
    Hist2->SetYTitle("v_{n}");
    Hist2->Draw();
    fc->SetIDName("ChGap");
    TH1D* hV22pt = (TH1D*)fc->GetVN2VsPt(2,0,4.9);
    if(!hV22pt){
        Printf("Can't get hV22");
        return;
    }
    SetMarkerAndLine(hV22pt,kBlack,kFullCircle,kSolid,1.0);
    hV22pt->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.8,0.5,0.9);
    legend2->AddEntry(hV22pt,Form("v_{2}{2}(p_{T}) |#Delta#eta|>1 cent:0~5%%"));
    legend2->Draw();

    // fc->SetIDName("Ch04GapA");
    // // fc->SetPropagateErrors(kTRUE);
    // TH1D* hV22 = (TH1D*)fc->GetHistCorrXXVsMulti("422");
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
    
}
