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

Double_t x_v2[] = {0,5,10,20,30,40,50,60,70,80,90};
// vector<Double_t> x_pass2 = {0,5,10,20,30,40,50,60,70,80,90};
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

TH1D* GetAbsValue(Int_t NBin, TH1D* h1, TH1D* h2, Bool_t ispTDiff=false){
    if(!h1 ||!h2)return 0;
    TH1D* AbsValue = nullptr;
    if(!ispTDiff)AbsValue = new TH1D(Form("AbsValue_%s",h1->GetName()),Form("AbsValue_%s",h1->GetName()),NBin,x_v2);
    else AbsValue = new TH1D(Form("AbsValue_%s",h1->GetName()),Form("AbsValue_%s",h1->GetName()),NBin,x_pTDiff);
    for(int i=1;i<=NBin;i++){
        Double_t error1 = h1->GetBinError(i);
        Double_t error2 = h2->GetBinError(i);
        Double_t value1 = h1->GetBinContent(i);
        Double_t value2 = h2->GetBinContent(i);
        Double_t error = sqrt(error1*error1 + error2*error2);
        Double_t value = value1 - value2;
        AbsValue->SetBinContent(i,value);
        AbsValue->SetBinError(i,error);
    }
    return AbsValue;
}

enum kObservable{
    kVn,
    kV24,
    kV26,
    kV28,
    kV210,
    kV422,
    kChi422,
    kRho422,
    kNSC23,
    kNSC24,
    kNSC234,
    kNSC345,
    kpTDiffv2Cent0To5,
    kpTDiffv2Cent5To10,
    kpTDiffv2Cent10To20,
    kpTDiffv2Cent20To30,
    kpTDiffv2Cent30To40,
    kpTDiffv2Cent40To50,
    kpTDiffv2Cent50To60,
    kpTDiffv2Cent60To70,
    kNObservable
};

std::map<kObservable,bool> IfCheckObservable = {
    {kVn,true},
    {kV24,true},
    {kV26,true},
    {kV28,true},
    {kV210,false},
    {kV422,true},
    {kChi422,true},
    {kRho422,true},
    {kNSC23,true},
    {kNSC24,true},
    {kNSC234,true},
    {kNSC345,true},
    {kpTDiffv2Cent0To5,true},
    {kpTDiffv2Cent5To10,true},
    {kpTDiffv2Cent10To20,true},
    {kpTDiffv2Cent20To30,true},
    {kpTDiffv2Cent30To40,true},
    {kpTDiffv2Cent40To50,true},
    {kpTDiffv2Cent50To60,true},
    {kpTDiffv2Cent60To70,true}
};
        
void ProcessFlowContainerDrawDiffRatio2Run2(){
    vector<string> FileNameSuffixs;
    vector<string> legendNames;
    vector<TFile*> resultsFiles_vn;
    vector<TFile*> resultsFiles_v24;
    vector<TFile*> resultsFiles_v26;
    vector<TFile*> resultsFiles_v28;
    vector<TFile*> resultsFiles_v210;
    vector<TFile*> resultsFiles_v422;
    vector<TFile*> resultsFiles_chi422;
    vector<TFile*> resultsFiles_rho422;
    vector<TFile*> resultsFiles_NSC23;
    vector<TFile*> resultsFiles_NSC24;
    vector<TFile*> resultsFiles_NSC234;
    vector<TFile*> resultsFiles_NSC345;
    vector<TFile*> resultsFiles_pTDiffv2Cent0To5;
    vector<TFile*> resultsFiles_pTDiffv2Cent5To10;
    vector<TFile*> resultsFiles_pTDiffv2Cent10To20;
    vector<TFile*> resultsFiles_pTDiffv2Cent20To30;
    vector<TFile*> resultsFiles_pTDiffv2Cent30To40;
    vector<TFile*> resultsFiles_pTDiffv2Cent40To50;
    vector<TFile*> resultsFiles_pTDiffv2Cent50To60;
    vector<TFile*> resultsFiles_pTDiffv2Cent60To70;
    // FileNameSuffixs.push_back("LHC23zzh_pass4_305042");
    // legendNames.push_back("LHC23zzh_pass4_305042");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_285891");
    // legendNames.push_back("previous results: with TPC rejection");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327042");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_326471");
    // legendNames.push_back("default bad NUA");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_277299");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_277299_VtxZ6");
    // legendNames.push_back("VtxZ6");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_277299_VtxZ8");
    // legendNames.push_back("VtxZ8");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_308058");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_308058_pT5GeV");
    // legendNames.push_back("pT 0.2-5 GeV/c");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_308058_ITSMatch");
    // legendNames.push_back("this update: global + Run3ITSall7Layers");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306992");
    // legendNames.push_back("this update: default (global)");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306992_ITSMatch");
    // legendNames.push_back("this update: global + Run3ITSall7Layers");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306755");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306755_NUECheck");
    // legendNames.push_back("default change NUE");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_279260");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_279260_midOcc");
    // legendNames.push_back("mid-occupancy");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_279260_highOcc");
    // legendNames.push_back("high-occupancy");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_279536");
    // legendNames.push_back("remove kTVXinTRD");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_279536_Mult_Cor");
    // legendNames.push_back("remove kTVXinTRD w/o MultCor");
    // FileNameSuffixs.push_back("LHC23zzn_pass3_I_A11_218295");
    // legendNames.push_back("mixed interaction rate");
    // FileNameSuffixs.push_back("LHC23zzn_pass3_I_A11_small_219651");
    // legendNames.push_back("low interaction rate");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327042");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327042_kIsGoodITSLayersAll");
    // legendNames.push_back("default + kIsGoodITSLayersAll");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_306468_FineCentBin");
    // legendNames.push_back("full 2023 data");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_308058");
    // legendNames.push_back("2023 data");
    // FileNameSuffixs.push_back("LHC24ar_pass1_medium_318934");
    // legendNames.push_back("2024 data");
    FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815");
    legendNames.push_back("default");
    FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815_kIsGoodITSLayersAll");
    legendNames.push_back("default+kIsGoodITSLayersAll");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815_GoodZvtxFT0vsPV");
    // legendNames.push_back("default+GoodZvtxFT0vsPV");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815_kNoCollInTimeRangeStandard");
    // legendNames.push_back("default+kNoCollInTimeRangeStandard");


    if (FileNameSuffixs.empty() || legendNames.empty()) {
        Printf("Error: no input files specified.");
        return;
    }

    for(auto& suffix : FileNameSuffixs){
        string fileName_vn = "./ProcessOutput/vn_" + suffix + ".root";
        string fileName_v24 = "./ProcessOutput/v24_" + suffix + ".root";
        string fileName_v26 = "./ProcessOutput/v26_" + suffix + ".root";
        string fileName_v28 = "./ProcessOutput/v28_" + suffix + ".root";
        string fileName_v210 = "./ProcessOutput/v210_" + suffix + ".root";
        string fileName_v422 = "./ProcessOutput/v422_" + suffix + ".root";
        string fileName_chi422 = "./ProcessOutput/chi422_" + suffix + ".root";
        string fileName_rho422 = "./ProcessOutput/rho422_" + suffix + ".root";
        string fileName_SC23 = "./ProcessOutput/NSC23_" + suffix + ".root";
        string fileName_SC24 = "./ProcessOutput/NSC24_" + suffix + ".root";
        string fileName_SC234 = "./ProcessOutput/NSC234_" + suffix + ".root";
        string fileName_SC345 = "./ProcessOutput/NSC345_" + suffix + ".root";
        string fileName_pTDiffv2Cent0To5 = "./ProcessOutput/pTDiffv2Cent0To5_" + suffix + ".root";
        string fileName_pTDiffv2Cent5To10 = "./ProcessOutput/pTDiffv2Cent5To10_" + suffix + ".root";
        string fileName_pTDiffv2Cent10To20 = "./ProcessOutput/pTDiffv2Cent10To20_" + suffix + ".root";
        string fileName_pTDiffv2Cent20To30 = "./ProcessOutput/pTDiffv2Cent20To30_" + suffix + ".root";
        string fileName_pTDiffv2Cent30To40 = "./ProcessOutput/pTDiffv2Cent30To40_" + suffix + ".root";
        string fileName_pTDiffv2Cent40To50 = "./ProcessOutput/pTDiffv2Cent40To50_" + suffix + ".root";
        string fileName_pTDiffv2Cent50To60 = "./ProcessOutput/pTDiffv2Cent50To60_" + suffix + ".root";
        string fileName_pTDiffv2Cent60To70 = "./ProcessOutput/pTDiffv2Cent60To70_" + suffix + ".root";

        TFile* resultsFile = nullptr;
        if (IfCheckObservable[kVn]) {
            resultsFile = TFile::Open(fileName_vn.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_vn << endl;
                return;
            }
            resultsFiles_vn.push_back(resultsFile);
        }
        if (IfCheckObservable[kV24]) {
            resultsFile = TFile::Open(fileName_v24.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_v24 << endl;
                return;
            }
            resultsFiles_v24.push_back(resultsFile);
        }
        if (IfCheckObservable[kV26]) {
            resultsFile = TFile::Open(fileName_v26.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_v26 << endl;
                return;
            }
            resultsFiles_v26.push_back(resultsFile);
        }
        if (IfCheckObservable[kV28]) {
            resultsFile = TFile::Open(fileName_v28.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_v28 << endl;
                return;
            }
            resultsFiles_v28.push_back(resultsFile);
        }
        if (IfCheckObservable[kV210]) {
            resultsFile = TFile::Open(fileName_v210.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_v210 << endl;
                return;
            }
            resultsFiles_v210.push_back(resultsFile);
        }
        if (IfCheckObservable[kV422]) {
            resultsFile = TFile::Open(fileName_v422.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_v422 << endl;
                return;
            }
            resultsFiles_v422.push_back(resultsFile);
        }
        if (IfCheckObservable[kChi422]) {
            resultsFile = TFile::Open(fileName_chi422.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_chi422 << endl;
                return;
            }
            resultsFiles_chi422.push_back(resultsFile);
        }
        if (IfCheckObservable[kRho422]) {
            resultsFile = TFile::Open(fileName_rho422.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_rho422 << endl;
                return;
            }
            resultsFiles_rho422.push_back(resultsFile);
        }
        if (IfCheckObservable[kNSC23]) {
            resultsFile = TFile::Open(fileName_SC23.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_SC23 << endl;
                return;
            }
            resultsFiles_NSC23.push_back(resultsFile);
        }
        if (IfCheckObservable[kNSC24]) {
            resultsFile = TFile::Open(fileName_SC24.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_SC24 << endl;
                return;
            }
            resultsFiles_NSC24.push_back(resultsFile);
        }
        if (IfCheckObservable[kNSC234]) {
            resultsFile = TFile::Open(fileName_SC234.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_SC234 << endl;
                return;
            }
            resultsFiles_NSC234.push_back(resultsFile);
        }
        if (IfCheckObservable[kNSC345]) {
            resultsFile = TFile::Open(fileName_SC345.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_SC345 << endl;
                return;
            }
            resultsFiles_NSC345.push_back(resultsFile);
        }
        if (IfCheckObservable[kpTDiffv2Cent0To5]) {
            resultsFile = TFile::Open(fileName_pTDiffv2Cent0To5.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_pTDiffv2Cent0To5 << endl;
                return;
            }
            resultsFiles_pTDiffv2Cent0To5.push_back(resultsFile);
        }
        if (IfCheckObservable[kpTDiffv2Cent5To10]) {
            resultsFile = TFile::Open(fileName_pTDiffv2Cent5To10.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_pTDiffv2Cent5To10 << endl;
                return;
            }
            resultsFiles_pTDiffv2Cent5To10.push_back(resultsFile);
        }
        if (IfCheckObservable[kpTDiffv2Cent10To20]) {
            resultsFile = TFile::Open(fileName_pTDiffv2Cent10To20.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_pTDiffv2Cent10To20 << endl;
                return;
            }
            resultsFiles_pTDiffv2Cent10To20.push_back(resultsFile);
        }
        if (IfCheckObservable[kpTDiffv2Cent20To30]) {
            resultsFile = TFile::Open(fileName_pTDiffv2Cent20To30.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_pTDiffv2Cent20To30 << endl;
                return;
            }
            resultsFiles_pTDiffv2Cent20To30.push_back(resultsFile);
        }
        if (IfCheckObservable[kpTDiffv2Cent30To40]) {
            resultsFile = TFile::Open(fileName_pTDiffv2Cent30To40.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_pTDiffv2Cent30To40 << endl;
                return;
            }
            resultsFiles_pTDiffv2Cent30To40.push_back(resultsFile);
        }
        if (IfCheckObservable[kpTDiffv2Cent40To50]) {
            resultsFile = TFile::Open(fileName_pTDiffv2Cent40To50.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_pTDiffv2Cent40To50 << endl;
                return;
            }
            resultsFiles_pTDiffv2Cent40To50.push_back(resultsFile);
        }
        if (IfCheckObservable[kpTDiffv2Cent50To60]) {
            resultsFile = TFile::Open(fileName_pTDiffv2Cent50To60.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_pTDiffv2Cent50To60 << endl;
                return;
            }
            resultsFiles_pTDiffv2Cent50To60.push_back(resultsFile);
        }
        if (IfCheckObservable[kpTDiffv2Cent60To70]) {
            resultsFile = TFile::Open(fileName_pTDiffv2Cent60To70.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_pTDiffv2Cent60To70 << endl;
                return;
            }
            resultsFiles_pTDiffv2Cent60To70.push_back(resultsFile);
        }
    }

    // No histogram statistics box
    gStyle->SetOptStat(0); 

    TFile* publish = new TFile("./PublicData/HEPData-ins1778342-v1-root.root","READ");
    TFile* publish_1671792 = new TFile("./PublicData/HEPData-ins1671792-v1-root.root","READ");
    TFile* publish_Run2pass2 = new TFile("./PublicData/PbPb_Run2Pass2.root","READ");
    int index = 0;
    TF1* One = new TF1("One","1",0,100);
    // =================
    // v2{2}
    // =================
    if (IfCheckObservable[kVn]) {
        index = 0;
        TCanvas* c1 = new TCanvas("c1", "c1", 800, 1200);
        c1->Divide(1,2);
        c1->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg = new TLegend(0.2,0.7,0.8,0.9);
        TH1D* frame_vn = new TH1D("frame_vn", "frame_vn", 90,0,90);
        frame_vn->SetMaximum(0.15);
        frame_vn->SetMinimum(0.);
        frame_vn->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_v22 = (TH1D*)resultsFiles_vn[i]->Get("Corr_corr_22_hist");
            SetMarkerAndLine(h_v22,GetColor(index),kFullCircle,kSolid,1.0);
            h_v22->Draw("ESames");
            leg->AddEntry(h_v22,Form("v_{2}{2} (%s)",legendNames[i].c_str()),"lp");
            index++;
        }

        TGraphAsymmErrors* g_v2 = (TGraphAsymmErrors*)publish_Run2pass2->Get("v2{2}_Gap10_TPCPileUp");
        // TFile* publish_v2 = new TFile("./v2_Gap10_TPCPileUp.root","READ");
        // TGraphAsymmErrors* g_v2 = (TGraphAsymmErrors*)publish_v2->Get("v2{2}_Gap10_TPCPileUp");
        SetMarkerAndLine(g_v2,kBlack,kOpenSquare,kSolid,1.0);
        g_v2->Draw("PE");
        leg->AddEntry(g_v2,Form("v_{2}{2} arXiv:2409.04343"));
        leg->Draw();
        
        TH1D* pub_v22 = new TH1D("pub_v22","pub_v22",10,x_v2);
        for(int i=0;i<10;i++){
            pub_v22->SetBinContent(i+1,g_v2->GetPointY(i));
            pub_v22->SetBinError(i+1,g_v2->GetErrorY(i));
        }


        c1->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio = new TH1D("frame_ratio", "frame_ratio", 90,0,90);
        frame_ratio->SetMaximum(1.5);
        frame_ratio->SetMinimum(0.7);
        frame_ratio->SetYTitle("Ratio to Run 2");
        frame_ratio->Draw("AXIS");
        One->SetLineColor(kBlack);
        One->SetLineStyle(kDashed);
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio = GetRatio(10,(TH1D*)resultsFiles_vn[i]->Get("Corr_corr_22_hist"),pub_v22);
            SetMarkerAndLine(ratio,GetColor(index),kFullCircle,kSolid,1.0);
            ratio->Draw("ESAMES");
            index++;
        }
        c1->SaveAs(Form("./OutputPDF/v22_RatioToRun2.pdf"));

        // =================
        // v3{2}
        // =================
        index = 0;
        TCanvas* c_v32 = new TCanvas("c_v32", "c_v32", 800, 1200);
        c_v32->Divide(1,2);
        c_v32->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_v32 = new TLegend(0.2,0.7,0.8,0.9);
        TH1D* frame_v32 = new TH1D("frame_v32", "frame_v32", 90,0,90);
        frame_v32->SetMaximum(0.15);
        frame_v32->SetMinimum(0.);
        frame_v32->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_v32 = (TH1D*)resultsFiles_vn[i]->Get("Corr_corr_32_hist");
            SetMarkerAndLine(h_v32,GetColor(index),kFullCircle,kSolid,1.0);
            h_v32->Draw("ESAMES");
            leg_v32->AddEntry(h_v32,Form("v_{3}{2} (%s)",legendNames[i].c_str()),"lp");
            index++;
        }
        // TFile* publish = new TFile("./PublicData/HEPData-ins1778342-v1-root.root","READ");
        TGraphAsymmErrors* g_v32 = (TGraphAsymmErrors*)publish_Run2pass2->Get("v3{2}_Gap10_TPCPileUp");
        SetMarkerAndLine(g_v32,kBlack,kOpenSquare,kSolid,1.0);
        g_v32->Draw("PE");
        leg_v32->AddEntry(g_v32,Form("v_{3}{2} arXiv:2409.04343"));
        leg_v32->Draw();

        TH1D* pub_v32 = new TH1D("pub_v32","pub_v32",10,x_v2);
        for(int i=0;i<10;i++){
            pub_v32->SetBinContent(i+1,g_v32->GetPointY(i));
            pub_v32->SetBinError(i+1,g_v32->GetErrorY(i));
        }

        c_v32->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_v32 = new TH1D("frame_ratio_v32", "frame_ratio_v32", 90,0,90);
        frame_ratio_v32->SetMaximum(1.5);
        frame_ratio_v32->SetMinimum(0.7);
        frame_ratio_v32->SetYTitle("Ratio to Run 2");
        frame_ratio_v32->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_v32 = GetRatio(10,(TH1D*)resultsFiles_vn[i]->Get("Corr_corr_32_hist"),pub_v32);
            SetMarkerAndLine(ratio_v32,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_v32->Draw("ESAMES");
            index++;
        }
        c_v32->SaveAs(Form("./OutputPDF/v32_RatioToRun2.pdf"));

        // =================
        // v4{2}
        // =================
        index = 0;
        TCanvas* c_v42 = new TCanvas("c_v42", "c_v42", 800, 1200);
        c_v42->Divide(1,2);
        c_v42->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_v42 = new TLegend(0.2,0.7,0.8,0.9);
        TH1D* frame_v42 = new TH1D("frame_v42", "frame_v42", 90,0,90);
        frame_v42->SetMaximum(0.15);
        frame_v42->SetMinimum(0.);
        frame_v42->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_v42 = (TH1D*)resultsFiles_vn[i]->Get("Corr_corr_42_hist");
            SetMarkerAndLine(h_v42,GetColor(index),kFullCircle,kSolid,1.0);
            h_v42->Draw("ESAMES");
            leg_v42->AddEntry(h_v42,Form("v_{4}{2} (%s)",legendNames[i].c_str()),"lp");
            index++;
        }
        // TFile* publish = new TFile("./PublicData/HEPData-ins1778342-v1-root.root","READ");
        // TGraphAsymmErrors* g_v42 = (TGraphAsymmErrors*)publish->Get("v4/Graph1D_y1");
        TGraphAsymmErrors* g_v42 = (TGraphAsymmErrors*)publish_Run2pass2->Get("v4{2}_Gap10_TPCPileUp");
        SetMarkerAndLine(g_v42,kBlack,kOpenSquare,kSolid,1.0);
        g_v42->Draw("PE");
        leg_v42->AddEntry(g_v42,Form("v_{4}{2} arXiv:2409.04343"));
        leg_v42->Draw();

        TH1D* pub_v42 = new TH1D("pub_v42","pub_v42",10,x_v2);
        for(int i=0;i<10;i++){
            pub_v42->SetBinContent(i+1,g_v42->GetPointY(i));
            pub_v42->SetBinError(i+1,g_v42->GetErrorY(i));
        }

        c_v42->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_v42 = new TH1D("frame_ratio_v42", "frame_ratio_v42", 90,0,90);
        frame_ratio_v42->SetMaximum(1.5);
        frame_ratio_v42->SetMinimum(0.7);
        frame_ratio_v42->SetYTitle("Ratio to Run 2");
        frame_ratio_v42->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_v42 = GetRatio(10,(TH1D*)resultsFiles_vn[i]->Get("Corr_corr_42_hist"),pub_v42);
            SetMarkerAndLine(ratio_v42,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_v42->Draw("ESAMES");
            index++;
        }
        c_v42->SaveAs(Form("./OutputPDF/v42_RatioToRun2.pdf"));
    }

    TFile* publish_ins1666817 = new TFile("./PublicData/HEPData-ins1666817-v1-root.root","READ");
    // =================
    // v2{4},v2{6},v2{8}
    // =================
    if (IfCheckObservable[kV24]) {
        index = 0;
        TCanvas* c_v24 = new TCanvas("c_v24", "c_v24", 900, 1200);
        c_v24->Divide(1,2);
        c_v24->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_v24 = new TLegend(0.2,0.7,0.8,0.9);
        TH1D* frame_v24 = new TH1D("frame_v24", "frame_v24", 90,0,90);
        frame_v24->SetMaximum(0.15);
        frame_v24->SetMinimum(0.);
        frame_v24->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_v24 = (TH1D*)resultsFiles_v24[i]->Get("Corr_corr_24_hist");
            SetMarkerAndLine(h_v24,GetColor(index),kFullCircle,kSolid,1.0);
            h_v24->Draw("ESAMES");
            leg_v24->AddEntry(h_v24,Form("v_{2}{4} (%s)",legendNames[i].c_str()),"lp");
            index++;
        }
        if (IfCheckObservable[kV26]) 
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_v26 = (TH1D*)resultsFiles_v26[i]->Get("corr_26_hist");
            SetMarkerAndLine(h_v26,GetColor(index),kFullSquare,kSolid,1.0);
            h_v26->Draw("ESAMES");
            leg_v24->AddEntry(h_v26,Form("v_{2}{6} (%s)",legendNames[i].c_str()),"lp");
            index++;
        }
        if (IfCheckObservable[kV28]) 
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_v28 = (TH1D*)resultsFiles_v28[i]->Get("corr_28_hist");
            SetMarkerAndLine(h_v28,GetColor(index),kFullStar,kSolid,1.0);
            h_v28->Draw("ESAMES");
            leg_v24->AddEntry(h_v28,Form("v_{2}{8} (%s)",legendNames[i].c_str()),"lp");
            index++;
        }
        if (IfCheckObservable[kV210]) 
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_v210 = (TH1D*)resultsFiles_v210[i]->Get("cN10_corr_210_hist");
            SetMarkerAndLine(h_v210,GetColor(index),kFullTriangleUp,kSolid,1.0);
            h_v210->Draw("ESAMES");
            leg_v24->AddEntry(h_v210,Form("v_{2}{10} (%s)",legendNames[i].c_str()),"lp");
            index++;
        }
        // TGraphAsymmErrors* g_v24 = (TGraphAsymmErrors*)publish_ins1666817->Get("Table 194/Graph1D_y1");
        TGraphAsymmErrors* g_v24 = (TGraphAsymmErrors*)publish_Run2pass2->Get("v24_TPCPileUp");
        SetMarkerAndLine(g_v24,kBlack,kOpenSquare,kSolid,1.0);
        g_v24->Draw("PE");
        // leg_v24->AddEntry(g_v24,Form("v_{2}{4} JHEP 07 (2018) 103"));
        leg_v24->AddEntry(g_v24,Form("v_{2}{4} arXiv:2409.04343"));
        leg_v24->Draw();

        TH1D* pub_v24 = new TH1D("pub_v24","pub_v24",7,x_v2);
        for(int i=0;i<7;i++){
            pub_v24->SetBinContent(i+1,g_v24->GetPointY(i));
            pub_v24->SetBinError(i+1,g_v24->GetErrorY(i));
        }

        c_v24->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_v24 = new TH1D("frame_ratio_v24", "frame_ratio_v24", 90,0,90);
        frame_ratio_v24->SetMaximum(1.05);
        frame_ratio_v24->SetMinimum(0.7);
        frame_ratio_v24->SetYTitle("Ratio to v2{4}");
        frame_ratio_v24->Draw("AXIS");
        One->Draw("sames");
        TLegend* leg_ratio_v24 = new TLegend(0.2,0.7,0.8,0.9);
        index=1;
        if (IfCheckObservable[kV26]) 
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_v26 = GetRatio(10,(TH1D*)resultsFiles_v26[i]->Get("corr_26_hist"),(TH1D*)resultsFiles_v24[i]->Get("Corr_corr_24_hist"));
            SetMarkerAndLine(ratio_v26,GetColor(index),kFullSquare,kSolid,1.0);
            ratio_v26->Draw("ESAMES");
            leg_ratio_v24->AddEntry(ratio_v26,Form("v_{2}{6}/v_{2}{4}"),"lp");
            index++;
        }
        if (IfCheckObservable[kV28]) 
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_v28 = GetRatio(10,(TH1D*)resultsFiles_v28[i]->Get("corr_28_hist"),(TH1D*)resultsFiles_v24[i]->Get("Corr_corr_24_hist")); 
            SetMarkerAndLine(ratio_v28,GetColor(index),kFullStar,kSolid,1.0);
            ratio_v28->Draw("ESAMES");
            leg_ratio_v24->AddEntry(ratio_v28,Form("v_{2}{8}/v_{2}{4}"),"lp");
            index++;
        }
        if (IfCheckObservable[kV210]) 
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_v210 = GetRatio(10,(TH1D*)resultsFiles_v210[i]->Get("cN10_corr_210_hist"),(TH1D*)resultsFiles_v24[i]->Get("Corr_corr_24_hist"));
            SetMarkerAndLine(ratio_v210,GetColor(index),kFullTriangleUp,kSolid,1.0);
            ratio_v210->Draw("ESAMES");
            leg_ratio_v24->AddEntry(ratio_v210,Form("v_{2}{10}/v_{2}{4}"),"lp");
            index++;
        }
        leg_ratio_v24->Draw();
        c_v24->SaveAs(Form("./OutputPDF/v2m_RatioToRun2.pdf"));
    }


    // =================
    // v422
    // =================
    if (IfCheckObservable[kV422]) {
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
            leg2->AddEntry(h_v422,Form("v_{4,22} (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        // TFile* publish = new TFile("./PublicData/HEPData-ins1778342-v1-root.root","READ");
        // TGraphAsymmErrors* g_v422 = (TGraphAsymmErrors*)publish->Get("v422/Graph1D_y1");
        TGraphAsymmErrors* g_v422 = (TGraphAsymmErrors*)publish_Run2pass2->Get("v422_Gap10_TPCPileUp");
        SetMarkerAndLine(g_v422,kBlack,kOpenSquare,kSolid,1.0);
        g_v422->Draw("PE");
        // leg2->AddEntry(g_v422,Form("v_{4,22} JHEP 05 (2020) 085, 2020"));
        leg2->AddEntry(g_v422,Form("v_{4,22} arXiv:2409.04343"));
        leg2->Draw();
        
        TH1D* pub_v422 = new TH1D("pub_v422","pub_v422",7,x_v2);
        for(int i=0;i<7;i++){
            pub_v422->SetBinContent(i+1,g_v422->GetPointY(i));
            pub_v422->SetBinError(i+1,g_v422->GetErrorY(i));
        }

        c2->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_v422 = new TH1D("frame_ratio_v422", "frame_ratio_v422", 90,0,90);
        frame_ratio_v422->SetMaximum(1.5);
        frame_ratio_v422->SetMinimum(0.5);
        frame_ratio_v422->SetYTitle("Ratio to Run 2");
        frame_ratio_v422->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_v422 = GetRatio(7,(TH1D*)resultsFiles_v422[i]->Get("hCorr422_mean"),pub_v422);
            SetMarkerAndLine(ratio_v422,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_v422->Draw("ESAMES");
            index++;
        }
        c2->SaveAs(Form("./OutputPDF/v422_RatioToRun2.pdf"));
    }

    // =================
    // chi422
    // =================
    if (IfCheckObservable[kChi422]) {
        index = 0;
        TCanvas* c4 = new TCanvas("c4", "c4", 800, 1200);
        c4->Divide(1,2);
        c4->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg4 = new TLegend(0.2,0.2,0.5,0.4);
        TH1D* frame_chi422 = new TH1D("frame_chi422", "frame_chi422", 90,0,90);
        frame_chi422->SetMaximum(2.);
        frame_chi422->SetMinimum(-5.5);
        frame_chi422->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_chi422 = (TH1D*)resultsFiles_chi422[i]->Get("hCorr422_mean");
            SetMarkerAndLine(h_chi422,GetColor(index),kFullCircle,kSolid,1.0);
            h_chi422->Draw("ESAMES");
            leg4->AddEntry(h_chi422,Form("#chi_{4,22} (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg4->Draw();

        // TGraphAsymmErrors* g_chi422 = (TGraphAsymmErrors*)publish->Get("chi422/Graph1D_y1");
        TGraphAsymmErrors* g_chi422 = (TGraphAsymmErrors*)publish_Run2pass2->Get("chi422_Gap10_TPCPileUp");
        SetMarkerAndLine(g_chi422,kBlack,kOpenSquare,kSolid,1.0);
        g_chi422->Draw("PE");
        // leg4->AddEntry(g_chi422,Form("#chi_{4,22} JHEP 05 (2020) 085, 2020"));
        leg4->AddEntry(g_chi422,Form("#chi_{4,22} arXiv:2409.04343"));
        leg4->Draw();

        TH1D* pub_chi422 = new TH1D("pub_chi422","pub_chi422",7,x_v2);
        for(int i=0;i<7;i++){
            pub_chi422->SetBinContent(i+1,g_chi422->GetPointY(i));
            pub_chi422->SetBinError(i+1,g_chi422->GetErrorY(i));
        }

        c4->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_chi422 = new TH1D("frame_ratio_chi422", "frame_ratio_chi422", 90,0,90);
        frame_ratio_chi422->SetMaximum(1.5);
        frame_ratio_chi422->SetMinimum(0.5);
        frame_ratio_chi422->SetYTitle("Ratio to Run 2");
        frame_ratio_chi422->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_chi422 = GetRatio(7,(TH1D*)resultsFiles_chi422[i]->Get("hCorr422_mean"),pub_chi422);
            SetMarkerAndLine(ratio_chi422,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_chi422->Draw("ESAMES");
            index++;
        }
        c4->SaveAs(Form("./OutputPDF/chi422_RatioToRun2.pdf"));
    }


    // =================
    // rho422
    // =================
    if (IfCheckObservable[kRho422]) {
        index = 0;
        TCanvas* c_rho422 = new TCanvas("c_rho422", "c_rho422", 800, 1200);
        c_rho422->Divide(1,2);
        c_rho422->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_rho422 = new TLegend(0.2,0.2,0.5,0.4);
        TH1D* frame_rho422 = new TH1D("frame_rho422", "frame_rho422", 90,0,90);
        frame_rho422->SetMaximum(1.);
        frame_rho422->SetMinimum(-0.5);
        frame_rho422->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_rho422 = (TH1D*)resultsFiles_rho422[i]->Get("hCorr422_mean");
            SetMarkerAndLine(h_rho422,GetColor(index),kFullCircle,kSolid,1.0);
            h_rho422->Draw("ESAMES");
            leg_rho422->AddEntry(h_rho422,Form("#rho_{4,22} (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg_rho422->Draw();

        // TGraphAsymmErrors* g_rho422 = (TGraphAsymmErrors*)publish->Get("rho422/Graph1D_y1");
        TGraphAsymmErrors* g_rho422 = (TGraphAsymmErrors*)publish_Run2pass2->Get("rho422_Gap10_TPCPileUp");
        SetMarkerAndLine(g_rho422,kBlack,kOpenSquare,kSolid,1.0);
        g_rho422->Draw("PE");
        // leg_rho422->AddEntry(g_rho422,Form("#rho_{4,22} JHEP 05 (2020) 085, 2020"));
        leg_rho422->AddEntry(g_rho422,Form("#rho_{4,22} arXiv:2409.04343"));
        leg_rho422->Draw();

        TH1D* pub_rho422 = new TH1D("pub_rho422","pub_rho422",7,x_v2);
        for(int i=0;i<7;i++){
            pub_rho422->SetBinContent(i+1,g_rho422->GetPointY(i));
            pub_rho422->SetBinError(i+1,g_rho422->GetErrorY(i));
        }

        c_rho422->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_rho422 = new TH1D("frame_ratio_rho422", "frame_ratio_rho422", 90,0,90);
        frame_ratio_rho422->SetMaximum(1.5);
        frame_ratio_rho422->SetMinimum(0.5);
        frame_ratio_rho422->SetYTitle("Ratio to Run 2");
        frame_ratio_rho422->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_rho422 = GetRatio(7,(TH1D*)resultsFiles_rho422[i]->Get("hCorr422_mean"),pub_rho422);
            SetMarkerAndLine(ratio_rho422,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_rho422->Draw("ESAMES");
            index++;
        }
        c_rho422->SaveAs(Form("./OutputPDF/rho422_RatioToRun2.pdf"));
    }


    TFile* publish_1848215 = new TFile("./PublicData/HEPData-ins1848215-v1-root.root","READ");
    // =================
    // NSC(2,3)
    // =================
    if (IfCheckObservable[kNSC23]) {
        index = 0;
        TCanvas* c3 = new TCanvas("c3", "c3", 800, 1200);
        c3->Divide(1,2);
        c3->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg3 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_NSC = new TH1D("frame_NSC", "frame_NSC", 90,0,90);
        frame_NSC->SetMaximum(2.);
        frame_NSC->SetMinimum(-2);
        frame_NSC->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_NSC = (TH1D*)resultsFiles_NSC23[i]->Get("NSC23");
            SetMarkerAndLine(h_NSC,GetColor(index),kFullCircle,kSolid,1.0);
            h_NSC->Draw("ESAMES");
            leg3->AddEntry(h_NSC,Form("NSC(2,3) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg3->Draw();

        TGraphAsymmErrors* g_NSC23 = (TGraphAsymmErrors*)publish_1848215->Get("Table 1/Graph1D_y1");
        SetMarkerAndLine(g_NSC23,kBlack,kOpenSquare,kSolid,1.0);
        g_NSC23->Draw("PE");
        leg3->AddEntry(g_NSC23,Form("NSC(2,3) PLB 818 (2021) 136354, 2021"));
        leg3->Draw();

        TH1D* pub_NSC23 = new TH1D("pub_NSC23","pub_NSC23",7,x_v2);
        for(int i=0;i<7;i++){
            pub_NSC23->SetBinContent(i+1,g_NSC23->GetPointY(i));
            pub_NSC23->SetBinError(i+1,g_NSC23->GetErrorY(i));
        }

        c3->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_NSC23 = new TH1D("frame_ratio_NSC23", "frame_ratio_NSC23", 90,0,90);
        frame_ratio_NSC23->SetMaximum(1.5);
        frame_ratio_NSC23->SetMinimum(0.5);
        frame_ratio_NSC23->SetYTitle("Ratio to Run 2");
        frame_ratio_NSC23->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_NSC23 = GetRatio(7,(TH1D*)resultsFiles_NSC23[i]->Get("NSC23"),pub_NSC23);
            SetMarkerAndLine(ratio_NSC23,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_NSC23->Draw("ESAMES");
            index++;
        }
        c3->SaveAs(Form("./OutputPDF/NSC23_RatioToRun2.pdf"));
    }

    // =================
    // NSC(2,4)
    // =================
    if (IfCheckObservable[kNSC24]) {
        index = 0;
        TCanvas* c_NSC24 = new TCanvas("c_NSC24", "c_NSC24", 800, 1200);
        c_NSC24->Divide(1,2);
        c_NSC24->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_NSC24 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_NSC24 = new TH1D("frame_NSC24", "frame_NSC24", 90,0,90);
        frame_NSC24->SetMaximum(2.);
        frame_NSC24->SetMinimum(-2);
        frame_NSC24->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_NSC24 = (TH1D*)resultsFiles_NSC24[i]->Get("NSC24");
            SetMarkerAndLine(h_NSC24,GetColor(index),kFullCircle,kSolid,1.0);
            h_NSC24->Draw("ESAMES");
            leg_NSC24->AddEntry(h_NSC24,Form("NSC(2,4) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg_NSC24->Draw();

        TGraphAsymmErrors* g_NSC24 = (TGraphAsymmErrors*)publish_1848215->Get("Table 2/Graph1D_y1");
        SetMarkerAndLine(g_NSC24,kBlack,kOpenSquare,kSolid,1.0);
        g_NSC24->Draw("PE");
        leg_NSC24->AddEntry(g_NSC24,Form("NSC(2,4) PLB 818 (2021) 136354, 2021"));
        leg_NSC24->Draw();

        TH1D* pub_NSC24 = new TH1D("pub_NSC24","pub_NSC24",7,x_v2);
        for(int i=0;i<7;i++){
            pub_NSC24->SetBinContent(i+1,g_NSC24->GetPointY(i));
            pub_NSC24->SetBinError(i+1,g_NSC24->GetErrorY(i));
        }

        c_NSC24->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_NSC24 = new TH1D("frame_ratio_NSC24", "frame_ratio_NSC24", 90,0,90);
        frame_ratio_NSC24->SetMaximum(1.5);
        frame_ratio_NSC24->SetMinimum(0.5);
        frame_ratio_NSC24->SetYTitle("Ratio to Run 2");
        frame_ratio_NSC24->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_NSC24 = GetRatio(7,(TH1D*)resultsFiles_NSC24[i]->Get("NSC24"),pub_NSC24);
            SetMarkerAndLine(ratio_NSC24,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_NSC24->Draw("ESAMES");
            index++;
        }
        c_NSC24->SaveAs(Form("./OutputPDF/NSC24_RatioToRun2.pdf"));
    }

    // =================
    // NSC(2,3,4)
    // =================
    if (IfCheckObservable[kNSC234]) {
        index = 0;
        TCanvas* c_NSC234 = new TCanvas("c_NSC234", "c_NSC234", 800, 1200);
        c_NSC234->Divide(1,2);
        c_NSC234->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_NSC234 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_NSC234 = new TH1D("frame_NSC234", "frame_NSC234", 90,0,90);
        frame_NSC234->SetMaximum(2.);
        frame_NSC234->SetMinimum(-2);
        frame_NSC234->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_NSC234 = (TH1D*)resultsFiles_NSC234[i]->Get("NSC234");
            SetMarkerAndLine(h_NSC234,GetColor(index),kFullCircle,kSolid,1.0);
            h_NSC234->Draw("ESAMES");
            leg_NSC234->AddEntry(h_NSC234,Form("NSC(2,3,4) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg_NSC234->Draw();

        TGraphAsymmErrors* g_NSC234 = (TGraphAsymmErrors*)publish_1848215->Get("Table 4/Graph1D_y1");
        SetMarkerAndLine(g_NSC234,kBlack,kOpenSquare,kSolid,1.0);
        g_NSC234->Draw("PE");
        leg_NSC234->AddEntry(g_NSC234,Form("NSC(2,3,4) PLB 818 (2021) 136354, 2021"));
        leg_NSC234->Draw();

        TH1D* pub_NSC234 = new TH1D("pub_NSC234","pub_NSC234",7,x_v2);
        for(int i=0;i<7;i++){
            pub_NSC234->SetBinContent(i+1,g_NSC234->GetPointY(i));
            pub_NSC234->SetBinError(i+1,g_NSC234->GetErrorY(i));
        }

        c_NSC234->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_NSC234 = new TH1D("frame_ratio_NSC234", "frame_ratio_NSC234", 90,0,90);
        frame_ratio_NSC234->SetMaximum(1.5);
        frame_ratio_NSC234->SetMinimum(0.5);
        frame_ratio_NSC234->SetYTitle("Ratio to Run 2");
        frame_ratio_NSC234->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_NSC234 = GetRatio(7,(TH1D*)resultsFiles_NSC234[i]->Get("NSC234"),pub_NSC234);
            SetMarkerAndLine(ratio_NSC234,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_NSC234->Draw("ESAMES");
            index++;
        }
        c_NSC234->SaveAs(Form("./OutputPDF/NSC234_RatioToRun2.pdf"));
    }

    // =================
    // NSC(3,4,5)
    // =================
    if (IfCheckObservable[kNSC345]) {
        index = 0;
        TCanvas* c_NSC345 = new TCanvas("c_NSC345", "c_NSC345", 800, 1200);
        c_NSC345->Divide(1,2);
        c_NSC345->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_NSC345 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_NSC345 = new TH1D("frame_NSC345", "frame_NSC345", 90,0,90);
        frame_NSC345->SetMaximum(2.);
        frame_NSC345->SetMinimum(-2);
        frame_NSC345->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_NSC345 = (TH1D*)resultsFiles_NSC345[i]->Get("NSC345");
            SetMarkerAndLine(h_NSC345,GetColor(index),kFullCircle,kSolid,1.0);
            h_NSC345->Draw("ESAMES");
            leg_NSC345->AddEntry(h_NSC345,Form("NSC(3,4,5) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg_NSC345->Draw();

        // TGraphAsymmErrors* g_NSC345 = (TGraphAsymmErrors*)publish_1848215->Get("Table 6/Graph1D_y1");
        // SetMarkerAndLine(g_NSC345,kBlack,kOpenSquare,kSolid,1.0);
        // g_NSC345->Draw("PE");
        // leg_NSC345->AddEntry(g_NSC345,Form("NSC(3,4,5) PLB 818 (2021) 136354, 2021"));
        // leg_NSC345->Draw();

        // TH1D* pub_NSC345 = new TH1D("pub_NSC345","pub_NSC345",7,x_v2);
        // for(int i=0;i<7;i++){
        //     pub_NSC345->SetBinContent(i+1,g_NSC345->GetPointY(i));
        //     pub_NSC345->SetBinError(i+1,g_NSC345->GetErrorY(i));
        // }

        // c_NSC345->cd(2);
        // gPad->SetTopMargin(0.05);
        // TH1D* frame_ratio_NSC345 = new TH1D("frame_ratio_NSC345", "frame_ratio_NSC345", 90,0,90);
        // frame_ratio_NSC345->SetMaximum(1.5);
        // frame_ratio_NSC345->SetMinimum(0.5);
        // frame_ratio_NSC345->SetYTitle("Ratio to Run 2");
        // frame_ratio_NSC345->Draw("AXIS");
        // One->Draw("sames");
        // index=0;
        // for(int i=0;i<FileNameSuffixs.size();i++){
        //     TH1D* ratio_NSC345 = GetRatio(7,(TH1D*)resultsFiles_NSC345[i]->Get("NSC345"),pub_NSC345);
        //     SetMarkerAndLine(ratio_NSC345,GetColor(index),kFullCircle,kSolid,1.0);
        //     ratio_NSC345->Draw("ESAMES");
        //     index++;
        // }
        c_NSC345->SaveAs(Form("./OutputPDF/NSC345_RatioToRun2.pdf"));
    }



    // =================
    // pTDiffv2Cent0To5
    // =================
    if (IfCheckObservable[kpTDiffv2Cent0To5]) {
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
            leg6->AddEntry(h_pTDiffv2Cent0To5,Form("v_{2}{2,|#Delta#eta|>1.0}(p_{T}) Cent:0~5%% (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        // TFile* publish_ins1666817 = new TFile("./PublicData/HEPData-ins1419244-v2-root.root","READ");
        // TFile* publish_ins1666817 = new TFile("./PublicData/HEPData-ins1666817-v1-root.root","READ");

        TGraphAsymmErrors* g_pt05 = nullptr;
        // g_pt05=(TGraphAsymmErrors*)publish_ins1666817->Get("Table 7/Graph1D_y1");
        g_pt05=(TGraphAsymmErrors*)publish_ins1666817->Get("Table 15/Graph1D_y1");
        if(!g_pt05)return;
        SetMarkerAndLine(g_pt05,kBlack,kOpenSquare,kSolid,1.0);
        g_pt05->Draw("PE");
        leg6->AddEntry(g_pt05,Form("v_{2}{2,|#Delta#eta|>1.}(p_{T}) Cent:0~5%% JHEP 07 (2018) 103"));
        leg6->Draw();

        TH1D* pub_pt05 = new TH1D("pub_pt05","pub_pt05",28,x_pTDiff);
        for(int i=0;i<28;i++){
            pub_pt05->SetBinContent(i+1,g_pt05->GetPointY(i));
            pub_pt05->SetBinError(i+1,g_pt05->GetErrorY(i));
        }

        c6->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_pt05 = new TH1D("frame_ratio_pt05", "frame_ratio_pt05",50,0.,10.);
        frame_ratio_pt05->SetMaximum(1.5);
        frame_ratio_pt05->SetMinimum(0.5);
        frame_ratio_pt05->SetYTitle("Ratio to Run 2");
        frame_ratio_pt05->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_pt05 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent0To5[i]->Get("pTDiffv2"),pub_pt05,true);
            SetMarkerAndLine(ratio_pt05,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_pt05->Draw("ESAMES");
            index++;
        }
        c6->SaveAs(Form("./OutputPDF/pTDiffv2Cent0To5_RatioToRun2.pdf"));
    }

    // =================
    // pTDiffv2Cent5To10
    // =================
    if (IfCheckObservable[kpTDiffv2Cent5To10]) {
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
            leg5->AddEntry(h_pTDiffv2Cent5To10,Form("v_{2}{2,|#Delta#eta|>1.0}(p_{T}) Cent:5~10%% (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        TGraphAsymmErrors* g_pt510=(TGraphAsymmErrors*)publish_ins1666817->Get("Table 23/Graph1D_y1");
        SetMarkerAndLine(g_pt510,kBlack,kOpenSquare,kSolid,1.0);
        g_pt510->Draw("PE");
        leg5->AddEntry(g_pt510,Form("v_{2}{2,|#Delta#eta|>1.}(p_{T}) Cent:5~10%% JHEP 07 (2018) 103"));
        leg5->Draw();

        TH1D* pub_pt510 = new TH1D("pub_pt510","pub_pt510",28,x_pTDiff);
        for(int i=0;i<28;i++){
            pub_pt510->SetBinContent(i+1,g_pt510->GetPointY(i));
            pub_pt510->SetBinError(i+1,g_pt510->GetErrorY(i));
        }

        c5->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_pt510 = new TH1D("frame_ratio_pt510", "frame_ratio_pt510",50,0.,10.);
        frame_ratio_pt510->SetMaximum(2);
        frame_ratio_pt510->SetMinimum(0.);
        frame_ratio_pt510->SetYTitle("Ratio to Run 2");
        frame_ratio_pt510->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_pt510 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent5To10[i]->Get("pTDiffv2"),pub_pt510,true);
            SetMarkerAndLine(ratio_pt510,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_pt510->Draw("ESAMES");
            index++;
        }
        c5->SaveAs(Form("./OutputPDF/pTDiffv2Cent5To10_RatioToRun2.pdf"));
    }

    // =================
    // pTDiffv2Cent10To20
    // =================
    if (IfCheckObservable[kpTDiffv2Cent10To20]) {
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
            leg10To20->AddEntry(h_pTDiffv2Cent10To20,Form("v_{2}{2,|#Delta#eta|>1.0}(p_{T}) Cent:10~20%% (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        TGraphAsymmErrors* g_pt1020=(TGraphAsymmErrors*)publish_ins1666817->Get("Table 33/Graph1D_y1");
        SetMarkerAndLine(g_pt1020,kBlack,kOpenSquare,kSolid,1.0);
        g_pt1020->Draw("PE");
        leg10To20->AddEntry(g_pt1020,Form("v_{2}{2,|#Delta#eta|>1.}(p_{T}) Cent:10~20%% JHEP 07 (2018) 103"));
        leg10To20->Draw();

        TH1D* pub_pt1020 = new TH1D("pub_pt1020","pub_pt1020",28,x_pTDiff);
        for(int i=0;i<28;i++){
            pub_pt1020->SetBinContent(i+1,g_pt1020->GetPointY(i));
            pub_pt1020->SetBinError(i+1,g_pt1020->GetErrorY(i));
        }

        c10To20->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_pt1020 = new TH1D("frame_ratio_pt1020", "frame_ratio_pt1020",50,0.,10.);
        frame_ratio_pt1020->SetMaximum(2);
        frame_ratio_pt1020->SetMinimum(0.);
        frame_ratio_pt1020->SetYTitle("Ratio to Run 2");
        frame_ratio_pt1020->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_pt1020 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent10To20[i]->Get("pTDiffv2"),pub_pt1020,true);
            SetMarkerAndLine(ratio_pt1020,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_pt1020->Draw("ESAMES");
            index++;
        }
        c10To20->SaveAs(Form("./OutputPDF/pTDiffv2Cent10To20_RatioToRun2.pdf"));
    }


    // =================
    // pTDiffv2Cent20To30
    // =================
    if (IfCheckObservable[kpTDiffv2Cent20To30]) {
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
            leg20To30->AddEntry(h_pTDiffv2Cent20To30,Form("v_{2}{2,|#Delta#eta|>1.0}(p_{T}) Cent:20~30%% (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        TGraphAsymmErrors* g_pt2030=(TGraphAsymmErrors*)publish_ins1666817->Get("Table 43/Graph1D_y1");
        SetMarkerAndLine(g_pt2030,kBlack,kOpenSquare,kSolid,1.0);
        g_pt2030->Draw("PE");
        leg20To30->AddEntry(g_pt2030,Form("v_{2}{2,|#Delta#eta|>1.}(p_{T}) Cent:20~30%% JHEP 07 (2018) 103"));
        leg20To30->Draw();

        TH1D* pub_pt2030 = new TH1D("pub_pt2030","pub_pt2030",28,x_pTDiff);
        for(int i=0;i<28;i++){
            pub_pt2030->SetBinContent(i+1,g_pt2030->GetPointY(i));
            pub_pt2030->SetBinError(i+1,g_pt2030->GetErrorY(i));
        }

        c20To30->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_pt2030 = new TH1D("frame_ratio_pt2030", "frame_ratio_pt2030",50,0.,10.);
        frame_ratio_pt2030->SetMaximum(2);
        frame_ratio_pt2030->SetMinimum(0.);
        frame_ratio_pt2030->SetYTitle("Ratio to Run 2");
        frame_ratio_pt2030->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_pt2030 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent20To30[i]->Get("pTDiffv2"),pub_pt2030,true);
            SetMarkerAndLine(ratio_pt2030,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_pt2030->Draw("ESAMES");
            index++;
        }
        c20To30->SaveAs(Form("./OutputPDF/pTDiffv2Cent20To30_RatioToRun2.pdf"));
    }


    // =================
    // pTDiffv2Cent30To40
    // =================
    if (IfCheckObservable[kpTDiffv2Cent30To40]) {
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
            leg7->AddEntry(h_pTDiffv2Cent30To40,Form("v_{2}{2,|#Delta#eta|>1.0}(p_{T}) Cent:30~40%% (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        TGraphAsymmErrors* g_pT3040=(TGraphAsymmErrors*)publish_ins1666817->Get("Table 53/Graph1D_y1");
        SetMarkerAndLine(g_pT3040,kBlack,kOpenSquare,kSolid,1.0);
        g_pT3040->Draw("PE");
        leg7->AddEntry(g_pT3040,Form("v_{2}{2,|#Delta#eta|>1.}(p_{T}) Cent:30~40%% JHEP 07 (2018) 103"));
        leg7->Draw();

        TH1D* pub_pt3040 = new TH1D("pub_pt3040","pub_pt3040",28,x_pTDiff);
        for(int i=0;i<28;i++){
            pub_pt3040->SetBinContent(i+1,g_pT3040->GetPointY(i));
            pub_pt3040->SetBinError(i+1,g_pT3040->GetErrorY(i));
        }

        c7->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_pT3040 = new TH1D("frame_ratio_pT3040", "frame_ratio_pT3040",50,0.,10.);
        frame_ratio_pT3040->SetMaximum(2);
        frame_ratio_pT3040->SetMinimum(0.);
        frame_ratio_pT3040->SetYTitle("Ratio to Run 2");
        frame_ratio_pT3040->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_pT3040 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent30To40[i]->Get("pTDiffv2"),pub_pt3040,true);
            SetMarkerAndLine(ratio_pT3040,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_pT3040->Draw("ESAMES");
            index++;
        }
        c7->SaveAs(Form("./OutputPDF/pTDiffv2Cent30To40_RatioToRun2.pdf"));
    }


    // =================
    // pTDiffv2Cent40To50
    // =================
    if (IfCheckObservable[kpTDiffv2Cent40To50]) {
        index = 0;
        TCanvas* c_40To50 = new TCanvas("c_40To50", "c_40To50", 800, 1200);
        c_40To50->Divide(1,2);
        c_40To50->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_40To50 = new TLegend(0.3,0.15,0.9,0.4);
        TH1D* frame_pTDiffv2Cent40To50 = new TH1D("frame_pTDiffv2Cent40To50", "frame_pTDiffv2Cent40To50", 50,0,10.);
        frame_pTDiffv2Cent40To50->SetMaximum(0.3);
        frame_pTDiffv2Cent40To50->SetMinimum(0.);
        frame_pTDiffv2Cent40To50->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_pTDiffv2Cent40To50 = (TH1D*)resultsFiles_pTDiffv2Cent40To50[i]->Get("pTDiffv2");
            SetMarkerAndLine(h_pTDiffv2Cent40To50,GetColor(index),kFullCircle,kSolid,1.0);
            h_pTDiffv2Cent40To50->Draw("ESAMES");
            leg_40To50->AddEntry(h_pTDiffv2Cent40To50,Form("v_{2}{2,|#Delta#eta|>1.0}(p_{T}) Cent:40~50%% (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        TGraphAsymmErrors* g_pT4050=(TGraphAsymmErrors*)publish_ins1666817->Get("Table 63/Graph1D_y1");
        SetMarkerAndLine(g_pT4050,kBlack,kOpenSquare,kSolid,1.0);
        g_pT4050->Draw("PE");
        leg_40To50->AddEntry(g_pT4050,Form("v_{2}{2,|#Delta#eta|>1.}(p_{T}) Cent:40~50%% JHEP 07 (2018) 103"));
        leg_40To50->Draw();

        TH1D* pub_pt4050 = new TH1D("pub_pt4050","pub_pt4050",28,x_pTDiff);
        for(int i=0;i<28;i++){
            pub_pt4050->SetBinContent(i+1,g_pT4050->GetPointY(i));
            pub_pt4050->SetBinError(i+1,g_pT4050->GetErrorY(i));
        }

        c_40To50->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_pT4050 = new TH1D("frame_ratio_pT4050", "frame_ratio_pT4050",50,0.,10.);
        frame_ratio_pT4050->SetMaximum(2);
        frame_ratio_pT4050->SetMinimum(0.);
        frame_ratio_pT4050->SetYTitle("Ratio to Run 2");
        frame_ratio_pT4050->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_pT4050 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent40To50[i]->Get("pTDiffv2"),pub_pt4050,true);
            SetMarkerAndLine(ratio_pT4050,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_pT4050->Draw("ESAMES");
            index++;
        }
        c_40To50->SaveAs(Form("./OutputPDF/pTDiffv2Cent40To50_RatioToRun2.pdf"));
    }


    // =================
    // pTDiffv2Cent50To60
    // =================
    if (IfCheckObservable[kpTDiffv2Cent50To60]) {
        index = 0;
        TCanvas* c_50To60 = new TCanvas("c_50To60", "c_50To60", 800, 1200);
        c_50To60->Divide(1,2);
        c_50To60->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_50To60 = new TLegend(0.3,0.15,0.9,0.4);
        TH1D* frame_pTDiffv2Cent50To60 = new TH1D("frame_pTDiffv2Cent50To60", "frame_pTDiffv2Cent50To60", 50,0,10.);
        frame_pTDiffv2Cent50To60->SetMaximum(0.3);
        frame_pTDiffv2Cent50To60->SetMinimum(0.);
        frame_pTDiffv2Cent50To60->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_pTDiffv2Cent50To60 = (TH1D*)resultsFiles_pTDiffv2Cent50To60[i]->Get("pTDiffv2");
            SetMarkerAndLine(h_pTDiffv2Cent50To60,GetColor(index),kFullCircle,kSolid,1.0);
            h_pTDiffv2Cent50To60->Draw("ESAMES");
            leg_50To60->AddEntry(h_pTDiffv2Cent50To60,Form("v_{2}{2,|#Delta#eta|>1.0}(p_{T}) Cent:50~60%% (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        TGraphAsymmErrors* g_pT5060=(TGraphAsymmErrors*)publish_ins1666817->Get("Table 73/Graph1D_y1");
        SetMarkerAndLine(g_pT5060,kBlack,kOpenSquare,kSolid,1.0);
        g_pT5060->Draw("PE");
        leg_50To60->AddEntry(g_pT5060,Form("v_{2}{2,|#Delta#eta|>1.}(p_{T}) Cent:50~60%% JHEP 07 (2018) 103"));
        leg_50To60->Draw();

        TH1D* pub_pt5060 = new TH1D("pub_pt5060","pub_pt5060",28,x_pTDiff);
        for(int i=0;i<28;i++){
            pub_pt5060->SetBinContent(i+1,g_pT5060->GetPointY(i));
            pub_pt5060->SetBinError(i+1,g_pT5060->GetErrorY(i));
        }

        c_50To60->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_pT5060 = new TH1D("frame_ratio_pT5060", "frame_ratio_pT5060",50,0.,10.);
        frame_ratio_pT5060->SetMaximum(2);
        frame_ratio_pT5060->SetMinimum(0.);
        frame_ratio_pT5060->SetYTitle("Ratio to Run 2");
        frame_ratio_pT5060->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_pT5060 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent50To60[i]->Get("pTDiffv2"),pub_pt5060,true);
            SetMarkerAndLine(ratio_pT5060,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_pT5060->Draw("ESAMES");
            index++;
        }
        c_50To60->SaveAs(Form("./OutputPDF/pTDiffv2Cent50To60_RatioToRun2.pdf"));
    }

    // =================
    // pTDiffv2Cent60To70
    // =================
    if (IfCheckObservable[kpTDiffv2Cent60To70]) {
        index = 0;
        TCanvas* c_60To70 = new TCanvas("c_60To70", "c_60To70", 800, 1200);
        c_60To70->Divide(1,2);
        c_60To70->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_60To70 = new TLegend(0.3,0.15,0.9,0.4);
        TH1D* frame_pTDiffv2Cent60To70 = new TH1D("frame_pTDiffv2Cent60To70", "frame_pTDiffv2Cent60To70", 50,0,10.);
        frame_pTDiffv2Cent60To70->SetMaximum(0.3);
        frame_pTDiffv2Cent60To70->SetMinimum(0.);
        frame_pTDiffv2Cent60To70->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_pTDiffv2Cent60To70 = (TH1D*)resultsFiles_pTDiffv2Cent60To70[i]->Get("pTDiffv2");
            SetMarkerAndLine(h_pTDiffv2Cent60To70,GetColor(index),kFullCircle,kSolid,1.0);
            h_pTDiffv2Cent60To70->Draw("ESAMES");
            leg_60To70->AddEntry(h_pTDiffv2Cent60To70,Form("v_{2}{2,|#Delta#eta|>1.0}(p_{T}) Cent:60~70%% (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        TGraphAsymmErrors* g_pT6070=(TGraphAsymmErrors*)publish_ins1666817->Get("Table 83/Graph1D_y1");
        SetMarkerAndLine(g_pT6070,kBlack,kOpenSquare,kSolid,1.0);
        g_pT6070->Draw("PE");
        leg_60To70->AddEntry(g_pT6070,Form("v_{2}{2,|#Delta#eta|>1.}(p_{T}) Cent:60~70%% JHEP 07 (2018) 103"));
        leg_60To70->Draw();

        TH1D* pub_pt6070 = new TH1D("pub_pt6070","pub_pt6070",28,x_pTDiff);
        for(int i=0;i<28;i++){
            pub_pt6070->SetBinContent(i+1,g_pT6070->GetPointY(i));
            pub_pt6070->SetBinError(i+1,g_pT6070->GetErrorY(i));
        }

        c_60To70->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_pT6070 = new TH1D("frame_ratio_pT6070", "frame_ratio_pT6070",50,0.,10.);
        frame_ratio_pT6070->SetMaximum(2);
        frame_ratio_pT6070->SetMinimum(0.);
        frame_ratio_pT6070->SetYTitle("Ratio to Run 2");
        frame_ratio_pT6070->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_pT6070 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2Cent60To70[i]->Get("pTDiffv2"),pub_pt6070,true);
            SetMarkerAndLine(ratio_pT6070,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_pT6070->Draw("ESAMES");
            index++;
        }
        c_60To70->SaveAs(Form("./OutputPDF/pTDiffv2Cent60To70_RatioToRun2.pdf"));
    }


    return;
    

}
