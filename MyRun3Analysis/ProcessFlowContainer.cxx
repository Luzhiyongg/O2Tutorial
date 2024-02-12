//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"

// #include "FlowContainer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <cstring>

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

void Output_vn(string FileNameSuffix, FlowContainer* fc){
    TCanvas* canvas1 = new TCanvas("Canvass_vn","Canvas_vn",900,900);
    TH1D* Hist  = new TH1D(Form("v_{n} in %s",FileNameSuffix.c_str()),Form("v_{n} in %s",FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(0.);
    Hist->SetMaximum(0.15);
    Hist->SetXTitle("Centrality/%");
    Hist->SetYTitle("v_{n}");
    Hist->Draw();
    fc->SetIDName("Ch10Gap");
    fc->SetPropagateErrors(kTRUE);
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
}

void Output_ptDiffvn(string FileNameSuffix, FlowContainer* fc){
    TCanvas* canvas2 = new TCanvas("canvas_ptDiffvn","canvas_ptDiffvn",900,900);
    TH1D* Hist2  = new TH1D(Form("v_{n}(p_{T}) in %s",FileNameSuffix.c_str()),Form("v_{n}(p_{T}) in %s",FileNameSuffix.c_str()),20,0.2,10.0);
    Hist2->SetMinimum(0.);
    Hist2->SetMaximum(0.5);
    Hist2->SetXTitle("p_{T}");
    Hist2->SetYTitle("v_{n}");
    Hist2->Draw();
    fc->SetIDName("ChGap");
    TH1D* hV22pt = (TH1D*)fc->GetVN2VsPt(2,0,5.);
    if(!hV22pt){
        Printf("Can't get hV22");
        return;
    }
    SetMarkerAndLine(hV22pt,kBlack,kFullCircle,kSolid,1.0);
    hV22pt->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.8,0.5,0.9);
    legend2->AddEntry(hV22pt,Form("v_{2}{2}(p_{T}) |#Delta#eta|>1 cent:0~5%%"));
    legend2->Draw();
}

void Output_vNL(string FileNameSuffix, FlowContainer* fc){
    
    TCanvas* canvas1 = new TCanvas("Canvass_vNL","Canvas_vNL",900,900);
    TH1D* Hist  = new TH1D(Form("v_NL in %s",FileNameSuffix.c_str()),Form("v_NL in %s",FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(0.);
    Hist->SetMaximum(0.15);
    Hist->SetXTitle("Centrality/%");
    Hist->SetYTitle("v_{n}");
    Hist->Draw();

    // fc->SetPropagateErrors(kTRUE);
    fc->SetIDName("Ch10GapA");
    TH1D* hCorr422A = (TH1D*)fc->GetHistCorrXXVsMulti("422");
    if(!hCorr422A){
        Printf("Can't get hCorr422A");
        return;
    }
    TH1D* hCorr422 = (TH1D*)hCorr422A->Clone();
    hCorr422->SetName("hCorr42");

    fc->SetIDName("Ch10GapB");
    TH1D* hCorr422B = (TH1D*)fc->GetHistCorrXXVsMulti("422");
    if(!hCorr422B){
        Printf("Can't get hCorr422B");
        return;
    }
    fc->SetIDName("Ch10Gap");
    TH1D* hCorr24= (TH1D*)fc->GetHistCorrXXVsMulti("24");
    if(!hCorr24){
        Printf("Can't get hCorr24");
        return;
    }

    for(int i=1;i<=hCorr422->GetNbinsX();i++){
        // Printf("%f",hCorr422A->GetBinCenter(i));
        hCorr422->SetBinContent(i,(hCorr422A->GetBinContent(i)+hCorr422B->GetBinContent(i))/2.);
    }

    TH1D* hv422 = (TH1D*)hCorr422A->Clone();
    for(int i=1;i<=hCorr422->GetNbinsX();i++){
        if(hCorr24->GetBinContent(i)>0.){
            hv422->SetBinContent(i,hCorr422->GetBinContent(i)/sqrt(hCorr24->GetBinContent(i)));
        }
        else{
            hv422->SetBinContent(i,0.);
            Printf("Warning: in %f %%, Ch10Gap24 is negative", hCorr422A->GetBinCenter(i));
        }
    }

    SetMarkerAndLine(hv422,kBlack,kFullCircle,kSolid,1.0);
    hv422->Draw("E1");

}

void ProcessFlowContainer(string FileNameSuffix = "LHC23zzh_pass2_small"){
    TFile* f = new TFile(Form("./AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    FlowContainer* fc = (FlowContainer*)f->Get("flow-pb-pb-task/FlowContainer");
    if(!fc){
        Printf("can not get flow container");
        return;
    }

    Output_vn(FileNameSuffix, fc);
    Output_ptDiffvn(FileNameSuffix, fc);
    Output_vNL(FileNameSuffix, fc);
    return;

    TCanvas* canvas1 = new TCanvas("canvas1","canvas1",900,900);
    TH1D* Hist  = new TH1D("Hist","v_{n} in LHC23zzi_pass2_QC1_sampling",8,0,80);
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

    TCanvas* canvas2 = new TCanvas("canvas2","canvas2",900,900);
    TH1D* Hist2  = new TH1D("Hist2","v_{n}(p_{T}) in LHC22s_pass5",10,0.2,5.0);
    Hist2->SetMinimum(0.);
    Hist2->SetMaximum(0.3);
    Hist2->SetXTitle("p_{T}");
    Hist2->SetYTitle("v_{n}");
    Hist2->Draw();
    fc->SetIDName("ChGap");
    TH1D* hV22pt = (TH1D*)fc->GetVN2VsPt(2,0,5.);
    if(!hV22pt){
        Printf("Can't get hV22");
        return;
    }
    SetMarkerAndLine(hV22pt,kBlack,kFullCircle,kSolid,1.0);
    hV22pt->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.8,0.5,0.9);
    legend2->AddEntry(hV22pt,Form("v_{2}{2}(p_{T}) |#Delta#eta|>1 cent:0~5%%"));
    legend2->Draw();


    
}
