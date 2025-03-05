/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch) 
 * @Date: 2025-03-05 22:01:51 
 * @Last Modified by:   Zhiyong Lu 
 * @Last Modified time: 2025-03-05 22:01:51 
 */
//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

// #include "FlowContainer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TObjArray.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include <cstring>
#include <vector>
#include <map>
#include <array>
#include "include/ErrorPropagation.h"

using namespace std;

bool ComparewithPublish = true;
bool OutputPNG = false;

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

void SetMarkerAndLine(TGraphAsymmErrors* graph, Int_t Color=0, Int_t style=0, Int_t linestyle=0, Float_t size=1){
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

Color_t GetColor(Int_t index){
    switch(index){
        default: return kBlack;
        case 0: return kRed;
        case 1: return kBlue;
        case 2: return kGreen;    
        case 3: return kMagenta;
        case 4: return kCyan;
        case 5: return kYellow;
        case 6: return kOrange;
        case 7: return kSpring;
        case 8: return kTeal;
        case 9: return kAzure;
        case 10: return kViolet;
        case 11: return kPink;
        case 12: return kGray;
        break;
    }
}

Double_t x_v2[] = {0,5,10,20,30,40,50,60};
Double_t x_pTDiff[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5,5.5,6,7,8,9,10};
TH1D* GetRatio(Int_t NBin, TH1D* h1, TH1D* h2, Bool_t ispTDiff=false){
    if(!h1 ||!h2)return 0;
    TH1D* ratio = nullptr;
    if(!ispTDiff)ratio = new TH1D(Form("ratio_%s",h1->GetName()),Form("ratio_%s",h1->GetName()),NBin,x_v2);
    else ratio = new TH1D(Form("ratio_%s",h1->GetName()),Form("ratio_%s",h1->GetName()),NBin,x_pTDiff);
    for(int i=1;i<=NBin;i++){
        Double_t error1 = h1->GetBinError(i);
        Double_t error2 = h2->GetBinError(i);
        Double_t value1 = h1->GetBinContent(i);
        Double_t value2 = h2->GetBinContent(i);
        Double_t error = Error_Ratio(value1,error1,value2,error2,0);
        Double_t value = value1/value2;
        ratio->SetBinContent(i,value);
        ratio->SetBinError(i,error);
    }
    return ratio;
}
        
void ProcessFlowContainerDrawDiffRatio(){
    vector<string> FileNameSuffixs;
    vector<TFile*> resultsFiles_vn;
    vector<TFile*> resultsFiles_v422;
    vector<TFile*> resultsFiles_chi422;
    vector<TFile*> resultsFiles_rho422;
    vector<TFile*> resultsFiles_NSC;
    vector<TFile*> resultsFiles_pTDiffv2Cent0To5;
    vector<TFile*> resultsFiles_pTDiffv2Cent5To10;
    vector<TFile*> resultsFiles_pTDiffv2Cent10To20;
    vector<TFile*> resultsFiles_pTDiffv2Cent20To30;
    vector<TFile*> resultsFiles_pTDiffv2Cent30To40;
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306813");
    FileNameSuffixs.push_back("LHC23zzh_pass4_small_306992");
    FileNameSuffixs.push_back("LHC23zzh_pass4_small_306992_ITSMatch");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306755");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306755_NUECheck");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_279260");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_279260_midOcc");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_279260_highOcc");
    // FileNameSuffixs.push_back("LHC23zzh_PbPb_pass4_small_279536");

    for(auto& suffix : FileNameSuffixs){
        string fileName_vn = "./ProcessOutput/vn_" + suffix + ".root";
        string fileName_v422 = "./ProcessOutput/v422_" + suffix + ".root";
        string fileName_chi422 = "./ProcessOutput/chi422_" + suffix + ".root";
        string fileName_rho422 = "./ProcessOutput/rho422_" + suffix + ".root";
        string fileName_NSC = "./ProcessOutput/NSC_" + suffix + ".root";
        string fileName_pTDiffv2Cent0To5 = "./ProcessOutput/pTDiffv2Cent0To5_" + suffix + ".root";
        string fileName_pTDiffv2Cent5To10 = "./ProcessOutput/pTDiffv2Cent5To10_" + suffix + ".root";
        string fileName_pTDiffv2Cent10To20 = "./ProcessOutput/pTDiffv2Cent10To20_" + suffix + ".root";
        string fileName_pTDiffv2Cent20To30 = "./ProcessOutput/pTDiffv2Cent20To30_" + suffix + ".root";
        string fileName_pTDiffv2Cent30To40 = "./ProcessOutput/pTDiffv2Cent30To40_" + suffix + ".root";

        TFile* resultsFile = TFile::Open(fileName_vn.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_vn << endl;
            return;
        }
        resultsFiles_vn.push_back(resultsFile);
        resultsFile = TFile::Open(fileName_v422.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_v422 << endl;
            return;
        }
        resultsFiles_v422.push_back(resultsFile);
        resultsFile = TFile::Open(fileName_chi422.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_chi422 << endl;
            return;
        }
        resultsFiles_chi422.push_back(resultsFile);
        resultsFile = TFile::Open(fileName_rho422.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_rho422 << endl;
            return;
        }
        resultsFiles_rho422.push_back(resultsFile);
        resultsFile = TFile::Open(fileName_NSC.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_NSC << endl;
            return;
        }
        resultsFiles_NSC.push_back(resultsFile);
        resultsFile = TFile::Open(fileName_pTDiffv2Cent0To5.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_pTDiffv2Cent0To5 << endl;
            return;
        }
        resultsFiles_pTDiffv2Cent0To5.push_back(resultsFile);
        resultsFile = TFile::Open(fileName_pTDiffv2Cent5To10.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_pTDiffv2Cent5To10 << endl;
            return;
        }
        resultsFiles_pTDiffv2Cent5To10.push_back(resultsFile);
        resultsFile = TFile::Open(fileName_pTDiffv2Cent10To20.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_pTDiffv2Cent10To20 << endl;
            return;
        }
        resultsFiles_pTDiffv2Cent10To20.push_back(resultsFile);
        resultsFile = TFile::Open(fileName_pTDiffv2Cent20To30.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_pTDiffv2Cent20To30 << endl;
            return;
        }
        resultsFiles_pTDiffv2Cent20To30.push_back(resultsFile);
        resultsFile = TFile::Open(fileName_pTDiffv2Cent30To40.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_pTDiffv2Cent30To40 << endl;
            return;
        }
        resultsFiles_pTDiffv2Cent30To40.push_back(resultsFile);
    }

    // No histogram statistics box
    gStyle->SetOptStat(0); 

    // =================
    // vn{2}
    // =================
    int index = 0;
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 1200);
    c1->Divide(1,2);
    c1->cd(1);
    gPad->SetBottomMargin(0.05);
    TLegend* leg = new TLegend(0.2,0.7,0.8,0.9);
    TH1D* frame_vn = new TH1D("frame_vn", "frame_vn", 60,0,60);
    frame_vn->SetMaximum(0.15);
    frame_vn->SetMinimum(0.);
    frame_vn->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_v22 = (TH1D*)resultsFiles_vn[i]->Get("Corr_corr_22_hist");
        SetMarkerAndLine(h_v22,GetColor(index),kFullCircle,kSolid,1.0);
        h_v22->Draw("ESames");
        leg->AddEntry(h_v22,Form("v_{2}{2} (%s)",FileNameSuffixs[i].c_str()),"lp");
        index++;
    }
    leg->Draw();

    c1->cd(2);
    gPad->SetTopMargin(0.05);
    TH1D* frame_ratio = new TH1D("frame_ratio", "frame_ratio", 60,0,60);
    frame_ratio->SetMaximum(1.05);
    frame_ratio->SetMinimum(0.9);
    frame_ratio->SetYTitle("Ratio");
    frame_ratio->Draw("AXIS");
    TF1* One = new TF1("One","1",0,60);
    One->SetLineColor(kBlack);
    One->SetLineStyle(kDashed);
    One->Draw("sames");
    index=0;
    for(int i=1;i<FileNameSuffixs.size();i++){
        TH1D* ratio = GetRatio(7,(TH1D*)resultsFiles_vn[i]->Get("Corr_corr_22_hist"),(TH1D*)resultsFiles_vn[0]->Get("Corr_corr_22_hist"));
        SetMarkerAndLine(ratio,GetColor(i),kFullCircle,kSolid,1.0);
        ratio->Draw("ESAMES");
        index++;
    }

    // =================
    // v422
    // =================
    index = 0;
    TCanvas* c2 = new TCanvas("c2", "c2", 800, 1200);
    c2->Divide(1,2);
    c2->cd(1);
    gPad->SetBottomMargin(0.05);
    TLegend* leg2 = new TLegend(0.2,0.2,0.5,0.4);
    TH1D* frame_v422 = new TH1D("frame_v422", "frame_v422", 90,0,90);
    frame_v422->SetMaximum(0.015);
    frame_v422->SetMinimum(-0.01);
    frame_v422->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_v422 = (TH1D*)resultsFiles_v422[i]->Get("hCorr422_mean");
        SetMarkerAndLine(h_v422,GetColor(index),kFullCircle,kSolid,1.0);
        h_v422->Draw("ESAMES");
        leg2->AddEntry(h_v422,Form("v_{4,22} (%s)",FileNameSuffixs[i].c_str()),"lp");
        index+=1;
    }
    leg2->Draw();

    c2->cd(2);
    gPad->SetTopMargin(0.05);
    TH1D* frame_ratio_v422 = new TH1D("frame_ratio_v422", "frame_ratio_v422", 60,0,60);
    frame_ratio_v422->SetMaximum(1.1);
    frame_ratio_v422->SetMinimum(0.9);
    frame_ratio_v422->SetYTitle("Ratio");
    frame_ratio_v422->Draw("AXIS");
    One->Draw("sames");
    index=0;
    for(int i=1;i<FileNameSuffixs.size();i++){
        TH1D* ratio_v422 = GetRatio(7,(TH1D*)resultsFiles_v422[i]->Get("hCorr422_mean"),(TH1D*)resultsFiles_v422[0]->Get("hCorr422_mean"));
        SetMarkerAndLine(ratio_v422,GetColor(i),kFullCircle,kSolid,1.0);
        ratio_v422->Draw("ESAMES");
        index++;
    }

    // =================
    // chi422
    // =================
    index = 0;
    TCanvas* c4 = new TCanvas("c4", "c4", 800, 1200);
    c4->Divide(1,2);
    c4->cd(1);
    gPad->SetBottomMargin(0.05);
    TLegend* leg4 = new TLegend(0.2,0.2,0.5,0.4);
    TH1D* frame_chi422 = new TH1D("frame_chi422", "frame_chi422", 60,0,60);
    frame_chi422->SetMaximum(2.);
    frame_chi422->SetMinimum(-5.5);
    frame_chi422->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_chi422 = (TH1D*)resultsFiles_chi422[i]->Get("hCorr422_mean");
        SetMarkerAndLine(h_chi422,GetColor(index),kFullCircle,kSolid,1.0);
        h_chi422->Draw("ESAMES");
        leg4->AddEntry(h_chi422,Form("#chi_{4,22} (%s)",FileNameSuffixs[i].c_str()),"lp");
        index+=1;
    }
    leg4->Draw();

    c4->cd(2);
    gPad->SetTopMargin(0.05);
    TH1D* frame_ratio_chi422 = new TH1D("frame_ratio_chi422", "frame_ratio_chi422", 60,0,60);
    frame_ratio_chi422->SetMaximum(1.1);
    frame_ratio_chi422->SetMinimum(0.9);
    frame_ratio_chi422->SetYTitle("Ratio");
    frame_ratio_chi422->Draw("AXIS");
    One->Draw("sames");
    index=0;
    for(int i=1;i<FileNameSuffixs.size();i++){
        TH1D* ratio_chi422 = GetRatio(7,(TH1D*)resultsFiles_chi422[i]->Get("hCorr422_mean"),(TH1D*)resultsFiles_chi422[0]->Get("hCorr422_mean"));
        SetMarkerAndLine(ratio_chi422,GetColor(i),kFullCircle,kSolid,1.0);
        ratio_chi422->Draw("ESAMES");
        index++;
    }


    // =================
    // rho422
    // =================
    index = 0;
    TCanvas* c_rho422 = new TCanvas("c_rho422", "c_rho422", 800, 1200);
    c_rho422->Divide(1,2);
    c_rho422->cd(1);
    gPad->SetBottomMargin(0.05);
    TLegend* leg_rho422 = new TLegend(0.2,0.2,0.5,0.4);
    TH1D* frame_rho422 = new TH1D("frame_rho422", "frame_rho422", 60,0,60);
    frame_rho422->SetMaximum(1.);
    frame_rho422->SetMinimum(-0.5);
    frame_rho422->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_rho422 = (TH1D*)resultsFiles_rho422[i]->Get("hCorr422_mean");
        SetMarkerAndLine(h_rho422,GetColor(index),kFullCircle,kSolid,1.0);
        h_rho422->Draw("ESAMES");
        leg_rho422->AddEntry(h_rho422,Form("#rho_{4,22} (%s)",FileNameSuffixs[i].c_str()),"lp");
        index+=1;
    }
    leg_rho422->Draw();

    c_rho422->cd(2);
    gPad->SetTopMargin(0.05);
    TH1D* frame_ratio_rho422 = new TH1D("frame_ratio_rho422", "frame_ratio_rho422", 60,0,60);
    frame_ratio_rho422->SetMaximum(1.1);
    frame_ratio_rho422->SetMinimum(0.9);
    frame_ratio_rho422->SetYTitle("Ratio");
    frame_ratio_rho422->Draw("AXIS");
    One->Draw("sames");
    index=0;
    for(int i=1;i<FileNameSuffixs.size();i++){
        TH1D* ratio_rho422 = GetRatio(7,(TH1D*)resultsFiles_rho422[i]->Get("hCorr422_mean"),(TH1D*)resultsFiles_rho422[0]->Get("hCorr422_mean"));
        SetMarkerAndLine(ratio_rho422,GetColor(i),kFullCircle,kSolid,1.0);
        ratio_rho422->Draw("ESAMES");
        index++;
    }


    // =================
    // NSC(3,2)
    // =================
    index = 0;
    TCanvas* c3 = new TCanvas("c3", "c3", 800, 1200);
    c3->Divide(1,2);
    c3->cd(1);
    gPad->SetBottomMargin(0.05);
    TLegend* leg3 = new TLegend(0.2,0.7,0.5,0.9);
    TH1D* frame_NSC = new TH1D("frame_NSC", "frame_NSC", 60,0,60);
    frame_NSC->SetMaximum(1.);
    frame_NSC->SetMinimum(-1);
    frame_NSC->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_NSC = (TH1D*)resultsFiles_NSC[i]->Get("NSC32");
        SetMarkerAndLine(h_NSC,GetColor(index),kFullCircle,kSolid,1.0);
        h_NSC->Draw("ESAMES");
        leg3->AddEntry(h_NSC,Form("NSC(3,2) (%s)",FileNameSuffixs[i].c_str()),"lp");
        index+=1;
    }
    leg3->Draw();

    c3->cd(2);
    gPad->SetTopMargin(0.05);
    TH1D* frame_ratio_NSC = new TH1D("frame_ratio_NSC", "frame_ratio_NSC", 60,0,60);
    frame_ratio_NSC->SetMaximum(1.1);
    frame_ratio_NSC->SetMinimum(0.9);
    frame_ratio_NSC->SetYTitle("Ratio");
    frame_ratio_NSC->Draw("AXIS");
    One->Draw("sames");
    index=0;
    for(int i=1;i<FileNameSuffixs.size();i++){
        TH1D* ratio_NSC = GetRatio(7,(TH1D*)resultsFiles_NSC[i]->Get("NSC32"),(TH1D*)resultsFiles_NSC[0]->Get("NSC32"));
        SetMarkerAndLine(ratio_NSC,GetColor(i),kFullCircle,kSolid,1.0);
        ratio_NSC->Draw("ESAMES");
        index++;
    }


    // =================
    // pTDiffv2Cent0To5
    // =================
    index = 0;
    TCanvas* c6 = new TCanvas("c6", "c6", 800, 1200);
    c6->Divide(1,2);
    c6->cd(1);
    gPad->SetBottomMargin(0.05);
    TLegend* leg6 = new TLegend(0.2,0.7,0.8,0.9);
    TH1D* frame_pTDiffv2Cent0To5 = new TH1D("frame_pTDiffv2Cent0To5", "frame_pTDiffv2Cent0To5", 50,0,10.);
    frame_pTDiffv2Cent0To5->SetMaximum(0.15);
    frame_pTDiffv2Cent0To5->SetMinimum(0.);
    frame_pTDiffv2Cent0To5->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_pTDiffv2Cent0To5 = (TH1D*)resultsFiles_pTDiffv2Cent0To5[i]->Get("pTDiffv2");
        SetMarkerAndLine(h_pTDiffv2Cent0To5,GetColor(index),kFullCircle,kSolid,1.0);
        h_pTDiffv2Cent0To5->Draw("ESAMES");
        leg6->AddEntry(h_pTDiffv2Cent0To5,Form("v_{2}(p_{T}) Cent:0~5%% (%s)",FileNameSuffixs[i].c_str()),"lp");
        index+=1;
    }
    leg6->Draw();

    c6->cd(2);
    gPad->SetTopMargin(0.05);
    TH1D* frame_ratio_pt05 = new TH1D("frame_ratio_pt05", "frame_ratio_pt05",50,0.,10.);
    frame_ratio_pt05->SetMaximum(1.1);
    frame_ratio_pt05->SetMinimum(0.9);
    frame_ratio_pt05->SetYTitle("Ratio");
    frame_ratio_pt05->Draw("AXIS");
    One->Draw("sames");
    index=0;
    for(int i=1;i<FileNameSuffixs.size();i++){
        TH1D* ratio_pt05 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent0To5[i]->Get("pTDiffv2"),(TH1D*)resultsFiles_pTDiffv2Cent0To5[0]->Get("pTDiffv2"),true);
        SetMarkerAndLine(ratio_pt05,GetColor(i),kFullCircle,kSolid,1.0);
        ratio_pt05->Draw("ESAMES");
        index++;
    }

    // =================
    // pTDiffv2Cent5To10
    // =================
    index = 0;
    TCanvas* c5 = new TCanvas("c5", "c5", 800, 1200);
    c5->Divide(1,2);
    c5->cd(1);
    gPad->SetBottomMargin(0.05);
    TLegend* leg5 = new TLegend(0.3,0.15,0.9,0.4);
    TH1D* frame_pTDiffv2Cent5To10 = new TH1D("frame_pTDiffv2Cent5To10", "frame_pTDiffv2Cent5To10", 50,0,10.);
    frame_pTDiffv2Cent5To10->SetMaximum(0.15);
    frame_pTDiffv2Cent5To10->SetMinimum(0.);
    frame_pTDiffv2Cent5To10->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_pTDiffv2Cent5To10 = (TH1D*)resultsFiles_pTDiffv2Cent5To10[i]->Get("pTDiffv2");
        SetMarkerAndLine(h_pTDiffv2Cent5To10,GetColor(index),kFullCircle,kSolid,1.0);
        h_pTDiffv2Cent5To10->Draw("ESAMES");
        leg5->AddEntry(h_pTDiffv2Cent5To10,Form("v_{2}(p_{T}) Cent:5~10%% (%s)",FileNameSuffixs[i].c_str()),"lp");
        index+=1;
    }
    leg5->Draw();

    c5->cd(2);
    gPad->SetTopMargin(0.05);
    TH1D* frame_ratio_pt510 = new TH1D("frame_ratio_pt510", "frame_ratio_pt510",50,0.,10.);
    frame_ratio_pt510->SetMaximum(1.1);
    frame_ratio_pt510->SetMinimum(0.9);
    frame_ratio_pt510->SetYTitle("Ratio");
    frame_ratio_pt510->Draw("AXIS");
    One->Draw("sames");
    index=0;
    for(int i=1;i<FileNameSuffixs.size();i++){
        TH1D* ratio_pt510 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent5To10[i]->Get("pTDiffv2"),(TH1D*)resultsFiles_pTDiffv2Cent5To10[0]->Get("pTDiffv2"),true);
        SetMarkerAndLine(ratio_pt510,GetColor(i),kFullCircle,kSolid,1.0);
        ratio_pt510->Draw("ESAMES");
        index++;
    }


    // =================
    // pTDiffv2Cent10To20
    // =================
    index = 0;
    TCanvas* c10To20 = new TCanvas("c10To20", "c10To20", 800, 1200);
    c10To20->Divide(1,2);
    c10To20->cd(1);
    gPad->SetBottomMargin(0.05);
    TLegend* leg10To20 = new TLegend(0.3,0.15,0.9,0.4);
    TH1D* frame_pTDiffv2Cent10To20 = new TH1D("frame_pTDiffv2Cent10To20", "frame_pTDiffv2Cent10To20", 50,0,10.);
    frame_pTDiffv2Cent10To20->SetMaximum(0.2);
    frame_pTDiffv2Cent10To20->SetMinimum(0.);
    frame_pTDiffv2Cent10To20->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_pTDiffv2Cent10To20 = (TH1D*)resultsFiles_pTDiffv2Cent10To20[i]->Get("pTDiffv2");
        SetMarkerAndLine(h_pTDiffv2Cent10To20,GetColor(index),kFullCircle,kSolid,1.0);
        h_pTDiffv2Cent10To20->Draw("ESAMES");
        leg10To20->AddEntry(h_pTDiffv2Cent10To20,Form("v_{2}(p_{T}) Cent:10~20%% (%s)",FileNameSuffixs[i].c_str()),"lp");
        index+=1;
    }
    leg10To20->Draw();

    c10To20->cd(2);
    gPad->SetTopMargin(0.05);
    TH1D* frame_ratio_pt1020 = new TH1D("frame_ratio_pt1020", "frame_ratio_pt1020",50,0.,10.);
    frame_ratio_pt1020->SetMaximum(1.1);
    frame_ratio_pt1020->SetMinimum(0.9);
    frame_ratio_pt1020->SetYTitle("Ratio");
    frame_ratio_pt1020->Draw("AXIS");
    One->Draw("sames");
    index=0;
    for(int i=1;i<FileNameSuffixs.size();i++){
        TH1D* ratio_pt1020 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent10To20[i]->Get("pTDiffv2"),(TH1D*)resultsFiles_pTDiffv2Cent10To20[0]->Get("pTDiffv2"),true);
        SetMarkerAndLine(ratio_pt1020,GetColor(i),kFullCircle,kSolid,1.0);
        ratio_pt1020->Draw("ESAMES");
        index++;
    }

    // =================
    // pTDiffv2Cent20To30
    // =================
    index = 0;
    TCanvas* c20To30 = new TCanvas("c20To30", "c20To30", 800, 1200);
    c20To30->Divide(1,2);
    c20To30->cd(1);
    gPad->SetBottomMargin(0.05);
    TLegend* leg20To30 = new TLegend(0.3,0.15,0.9,0.4);
    TH1D* frame_pTDiffv2Cent20To30 = new TH1D("frame_pTDiffv2Cent20To30", "frame_pTDiffv2Cent20To30", 50,0,10.);
    frame_pTDiffv2Cent20To30->SetMaximum(0.3);
    frame_pTDiffv2Cent20To30->SetMinimum(0.);
    frame_pTDiffv2Cent20To30->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_pTDiffv2Cent20To30 = (TH1D*)resultsFiles_pTDiffv2Cent20To30[i]->Get("pTDiffv2");
        SetMarkerAndLine(h_pTDiffv2Cent20To30,GetColor(index),kFullCircle,kSolid,1.0);
        h_pTDiffv2Cent20To30->Draw("ESAMES");
        leg20To30->AddEntry(h_pTDiffv2Cent20To30,Form("v_{2}(p_{T}) Cent:20~30%% (%s)",FileNameSuffixs[i].c_str()),"lp");
        index+=1;
    }
    leg20To30->Draw();

    c20To30->cd(2);
    gPad->SetTopMargin(0.05);
    TH1D* frame_ratio_pt2030 = new TH1D("frame_ratio_pt2030", "frame_ratio_pt2030",50,0.,10.);
    frame_ratio_pt2030->SetMaximum(1.1);
    frame_ratio_pt2030->SetMinimum(0.9);
    frame_ratio_pt2030->SetYTitle("Ratio");
    frame_ratio_pt2030->Draw("AXIS");
    One->Draw("sames");
    index=0;
    for(int i=1;i<FileNameSuffixs.size();i++){
        TH1D* ratio_pt2030 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent20To30[i]->Get("pTDiffv2"),(TH1D*)resultsFiles_pTDiffv2Cent20To30[0]->Get("pTDiffv2"),true);
        SetMarkerAndLine(ratio_pt2030,GetColor(i),kFullCircle,kSolid,1.0);
        ratio_pt2030->Draw("ESAMES");
        index++;
    }



    // =================
    // pTDiffv2Cent30To40
    // =================
    index = 0;
    TCanvas* c7 = new TCanvas("c7", "c7", 800, 1200);
    c7->Divide(1,2);
    c7->cd(1);
    gPad->SetBottomMargin(0.05);
    TLegend* leg7 = new TLegend(0.3,0.15,0.9,0.4);
    TH1D* frame_pTDiffv2Cent30To40 = new TH1D("frame_pTDiffv2Cent30To40", "frame_pTDiffv2Cent30To40", 50,0,10.);
    frame_pTDiffv2Cent30To40->SetMaximum(0.3);
    frame_pTDiffv2Cent30To40->SetMinimum(0.);
    frame_pTDiffv2Cent30To40->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_pTDiffv2Cent30To40 = (TH1D*)resultsFiles_pTDiffv2Cent30To40[i]->Get("pTDiffv2");
        SetMarkerAndLine(h_pTDiffv2Cent30To40,GetColor(index),kFullCircle,kSolid,1.0);
        h_pTDiffv2Cent30To40->Draw("ESAMES");
        leg7->AddEntry(h_pTDiffv2Cent30To40,Form("v_{2}(p_{T}) Cent:30~40%% (%s)",FileNameSuffixs[i].c_str()),"lp");
        index+=1;
    }
    leg7->Draw();

    c7->cd(2);
    gPad->SetTopMargin(0.05);
    TH1D* frame_ratio_pT3040 = new TH1D("frame_ratio_pT3040", "frame_ratio_pT3040",50,0.,10.);
    frame_ratio_pT3040->SetMaximum(1.1);
    frame_ratio_pT3040->SetMinimum(0.9);
    frame_ratio_pT3040->SetYTitle("Ratio");
    frame_ratio_pT3040->Draw("AXIS");
    One->Draw("sames");
    index=0;
    for(int i=1;i<FileNameSuffixs.size();i++){
        TH1D* ratio_pT3040 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent30To40[i]->Get("pTDiffv2"),(TH1D*)resultsFiles_pTDiffv2Cent30To40[0]->Get("pTDiffv2"),true);
        SetMarkerAndLine(ratio_pT3040,GetColor(i),kFullCircle,kSolid,1.0);
        ratio_pT3040->Draw("ESAMES");
        index++;
    }



    return;
    

}
