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
    kV422,
    kChi422,
    kRho422,
    kNSC,
    kSC23,
    kSC24,
    kSC234,
    kSC235,
    kSC246,
    kSC345,
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
    {kV24,false},
    {kV26,false},
    {kV28,false},
    {kV422,false},
    {kChi422,false},
    {kRho422,false},
    {kNSC,false},
    {kSC23,false},
    {kSC24,false},
    {kSC234,false},
    {kSC235,false},
    {kSC246,false},
    {kSC345,false},
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
    vector<TFile*> resultsFiles_v422;
    vector<TFile*> resultsFiles_chi422;
    vector<TFile*> resultsFiles_rho422;
    vector<TFile*> resultsFiles_NSC;
    vector<TFile*> resultsFiles_SC23;
    vector<TFile*> resultsFiles_SC24;
    vector<TFile*> resultsFiles_SC234;
    vector<TFile*> resultsFiles_SC235;
    vector<TFile*> resultsFiles_SC246;
    vector<TFile*> resultsFiles_SC345;
    vector<TFile*> resultsFiles_pTDiffv2Cent0To5;
    vector<TFile*> resultsFiles_pTDiffv2Cent5To10;
    vector<TFile*> resultsFiles_pTDiffv2Cent10To20;
    vector<TFile*> resultsFiles_pTDiffv2Cent20To30;
    vector<TFile*> resultsFiles_pTDiffv2Cent30To40;
    vector<TFile*> resultsFiles_pTDiffv2Cent40To50;
    vector<TFile*> resultsFiles_pTDiffv2Cent50To60;
    vector<TFile*> resultsFiles_pTDiffv2Cent60To70;
    // FileNameSuffixs.push_back("LHC23zzh_pass4_305042");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_285891");
    // legendNames.push_back("previous results: with TPC rejection");
    FileNameSuffixs.push_back("LHC23zzh_pass4_small_308058");
    legendNames.push_back("2023 data");
    FileNameSuffixs.push_back("LHC24ar_pass1_medium_318934");
    legendNames.push_back("2024 data");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_308058_pT5GeV");
    // legendNames.push_back("this update: pT 0.2-5 GeV/c");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_308058_ITSMatch");
    // legendNames.push_back("this update: global + Run3ITSall7Layers");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306992");
    // legendNames.push_back("this update: default (global)");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306992_ITSMatch");
    // legendNames.push_back("this update: global + Run3ITSall7Layers");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306755");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_306755_NUECheck");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_279260_midOcc");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_279260_highOcc");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_279536");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_279536_Mult_Cor");

    for(auto& suffix : FileNameSuffixs){
        string fileName_vn = "./ProcessOutput/vn_" + suffix + ".root";
        string fileName_v24 = "./ProcessOutput/v24_" + suffix + ".root";
        string fileName_v26 = "./ProcessOutput/v26_" + suffix + ".root";
        string fileName_v28 = "./ProcessOutput/v28_" + suffix + ".root";
        string fileName_v422 = "./ProcessOutput/v422_" + suffix + ".root";
        string fileName_chi422 = "./ProcessOutput/chi422_" + suffix + ".root";
        string fileName_rho422 = "./ProcessOutput/rho422_" + suffix + ".root";
        string fileName_NSC = "./ProcessOutput/NSC_" + suffix + ".root";
        string fileName_SC23 = "./ProcessOutput/SC23_" + suffix + ".root";
        string fileName_SC24 = "./ProcessOutput/SC24_" + suffix + ".root";
        string fileName_SC234 = "./ProcessOutput/SC234_" + suffix + ".root";
        string fileName_SC235 = "./ProcessOutput/SC235_" + suffix + ".root";
        string fileName_SC246 = "./ProcessOutput/SC246_" + suffix + ".root";
        string fileName_SC345 = "./ProcessOutput/SC345_" + suffix + ".root";
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
        if (IfCheckObservable[kNSC]) {
            resultsFile = TFile::Open(fileName_NSC.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_NSC << endl;
                return;
            }
            resultsFiles_NSC.push_back(resultsFile);
        }
        if (IfCheckObservable[kSC23]) {
            resultsFile = TFile::Open(fileName_SC23.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_SC23 << endl;
                return;
            }
            resultsFiles_SC23.push_back(resultsFile);
        }
        if (IfCheckObservable[kSC24]) {
            resultsFile = TFile::Open(fileName_SC24.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_SC24 << endl;
                return;
            }
            resultsFiles_SC24.push_back(resultsFile);
        }
        if (IfCheckObservable[kSC234]) {
            resultsFile = TFile::Open(fileName_SC234.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_SC234 << endl;
                return;
            }
            resultsFiles_SC234.push_back(resultsFile);
        }
        if (IfCheckObservable[kSC235]) {
            resultsFile = TFile::Open(fileName_SC235.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_SC235 << endl;
                return;
            }
            resultsFiles_SC235.push_back(resultsFile);
        }
        if (IfCheckObservable[kSC246]) {
            resultsFile = TFile::Open(fileName_SC246.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_SC246 << endl;
                return;
            }
            resultsFiles_SC246.push_back(resultsFile);
        }
        if (IfCheckObservable[kSC345]) {
            resultsFile = TFile::Open(fileName_SC345.c_str(), "READ");
            if(!resultsFile || resultsFile->IsZombie()){
                cout << "Error: cannot open file " << fileName_SC345 << endl;
                return;
            }
            resultsFiles_SC345.push_back(resultsFile);
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
        frame_ratio->SetMaximum(1.05);
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
        frame_ratio_v32->SetMaximum(1.05);
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
        frame_ratio_v42->SetMaximum(1.05);
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
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_v26 = (TH1D*)resultsFiles_v26[i]->Get("corr_26_hist");
            SetMarkerAndLine(h_v26,GetColor(index),kFullSquare,kSolid,1.0);
            h_v26->Draw("ESAMES");
            leg_v24->AddEntry(h_v26,Form("v_{2}{6} (%s)",legendNames[i].c_str()),"lp");
            index++;
        }
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_v28 = (TH1D*)resultsFiles_v28[i]->Get("corr_28_hist");
            SetMarkerAndLine(h_v28,GetColor(index),kFullStar,kSolid,1.0);
            h_v28->Draw("ESAMES");
            leg_v24->AddEntry(h_v28,Form("v_{2}{8} (%s)",legendNames[i].c_str()),"lp");
            index++;
        }
        TGraphAsymmErrors* g_v24 = (TGraphAsymmErrors*)publish_ins1666817->Get("Table 194/Graph1D_y1");
        SetMarkerAndLine(g_v24,kBlack,kOpenSquare,kSolid,1.0);
        g_v24->Draw("PE");
        leg_v24->AddEntry(g_v24,Form("v_{2}{4} JHEP 07 (2018) 103"));
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
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_v26 = GetRatio(7,(TH1D*)resultsFiles_v26[i]->Get("corr_26_hist"),(TH1D*)resultsFiles_v24[i]->Get("Corr_corr_24_hist"));
            SetMarkerAndLine(ratio_v26,GetColor(index),kFullSquare,kSolid,1.0);
            ratio_v26->Draw("ESAMES");
            leg_ratio_v24->AddEntry(ratio_v26,Form("v_{2}{6}/v_{2}{4}"),"lp");
            index++;
        }
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_v28 = GetRatio(7,(TH1D*)resultsFiles_v28[i]->Get("corr_28_hist"),(TH1D*)resultsFiles_v24[i]->Get("Corr_corr_24_hist")); 
            SetMarkerAndLine(ratio_v28,GetColor(index),kFullStar,kSolid,1.0);
            ratio_v28->Draw("ESAMES");
            leg_ratio_v24->AddEntry(ratio_v28,Form("v_{2}{8}/v_{2}{4}"),"lp");
            index++;
        }
        leg_ratio_v24->Draw();
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
        TGraphAsymmErrors* g_v422 = (TGraphAsymmErrors*)publish->Get("v422/Graph1D_y1");
        SetMarkerAndLine(g_v422,kBlack,kOpenSquare,kSolid,1.0);
        g_v422->Draw("PE");
        leg2->AddEntry(g_v422,Form("v_{4,22} JHEP 05 (2020) 085, 2020"));
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

        TGraphAsymmErrors* g_chi422 = (TGraphAsymmErrors*)publish->Get("chi422/Graph1D_y1");
        SetMarkerAndLine(g_chi422,kBlack,kOpenSquare,kSolid,1.0);
        g_chi422->Draw("PE");
        leg4->AddEntry(g_chi422,Form("#chi_{4,22} JHEP 05 (2020) 085, 2020"));
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

        TGraphAsymmErrors* g_rho422 = (TGraphAsymmErrors*)publish->Get("rho422/Graph1D_y1");
        SetMarkerAndLine(g_rho422,kBlack,kOpenSquare,kSolid,1.0);
        g_rho422->Draw("PE");
        leg_rho422->AddEntry(g_rho422,Form("#rho_{4,22} JHEP 05 (2020) 085, 2020"));
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
    }


    TFile* publish_NSC32 = new TFile("./PublicData/HEPData-ins1848215-v1-root.root","READ");
    // =================
    // NSC(3,2)
    // =================
    if (IfCheckObservable[kNSC]) {
        index = 0;
        TCanvas* c3 = new TCanvas("c3", "c3", 800, 1200);
        c3->Divide(1,2);
        c3->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg3 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_NSC = new TH1D("frame_NSC", "frame_NSC", 90,0,90);
        frame_NSC->SetMaximum(1.);
        frame_NSC->SetMinimum(-1);
        frame_NSC->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_NSC = (TH1D*)resultsFiles_NSC[i]->Get("NSC32");
            SetMarkerAndLine(h_NSC,GetColor(index),kFullCircle,kSolid,1.0);
            h_NSC->Draw("ESAMES");
            leg3->AddEntry(h_NSC,Form("NSC(3,2) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg3->Draw();

        TGraphAsymmErrors* g_NSC32 = (TGraphAsymmErrors*)publish_NSC32->Get("Table 1/Graph1D_y1");
        SetMarkerAndLine(g_NSC32,kBlack,kOpenSquare,kSolid,1.0);
        g_NSC32->Draw("PE");
        leg3->AddEntry(g_NSC32,Form("NSC(3,2) PLB 818 (2021) 136354, 2021"));
        leg3->Draw();

        TH1D* pub_NSC32 = new TH1D("pub_NSC32","pub_NSC32",7,x_v2);
        for(int i=0;i<7;i++){
            pub_NSC32->SetBinContent(i+1,g_NSC32->GetPointY(i));
            pub_NSC32->SetBinError(i+1,g_NSC32->GetErrorY(i));
        }

        c3->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_ratio_NSC32 = new TH1D("frame_ratio_NSC32", "frame_ratio_NSC32", 90,0,90);
        frame_ratio_NSC32->SetMaximum(1.5);
        frame_ratio_NSC32->SetMinimum(0.5);
        frame_ratio_NSC32->SetYTitle("Ratio to Run 2");
        frame_ratio_NSC32->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* ratio_NSC32 = GetRatio(7,(TH1D*)resultsFiles_NSC[i]->Get("NSC32"),pub_NSC32);
            SetMarkerAndLine(ratio_NSC32,GetColor(index),kFullCircle,kSolid,1.0);
            ratio_NSC32->Draw("ESAMES");
            index++;
        }
    }

    TFile* publish_ins1452590 = new TFile("./PublicData/HEPData-ins1452590-v1-root.root","READ");
    // =================
    // SC(2,3)
    // =================
    if (IfCheckObservable[kSC23]) {
        index = 0;
        TCanvas* c_sc23 = new TCanvas("c_sc23","c_sc23", 800, 1200);
        c_sc23->Divide(1,2);
        c_sc23->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_sc23 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_sc23 = new TH1D("frame_sc23", "frame_sc23", 90,0,90);
        frame_sc23->SetMaximum(5e-6);
        frame_sc23->SetMinimum(-5e-6);
        frame_sc23->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_sc23 = (TH1D*)resultsFiles_SC23[i]->Get("SC23");
            SetMarkerAndLine(h_sc23,GetColor(index),kFullCircle,kSolid,1.0);
            h_sc23->Draw("ESAMES");
            leg_sc23->AddEntry(h_sc23,Form("SC(2,3) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg_sc23->Draw();

        TGraphAsymmErrors* g_sc23 = (TGraphAsymmErrors*)publish_ins1452590->Get("Table 1/Graph1D_y1");
        SetMarkerAndLine(g_sc23,kBlack,kOpenSquare,kSolid,1.0);
        g_sc23->Draw("PE");
        leg_sc23->AddEntry(g_sc23,Form("SC(2,3) PRL. 117 (2016) 182301"));
        leg_sc23->Draw();

        TH1D* pub_sc23 = new TH1D("pub_sc23","pub_sc23",7,x_v2);
        for(int i=0;i<7;i++){
            pub_sc23->SetBinContent(i+1,g_sc23->GetPointY(i));
            pub_sc23->SetBinError(i+1,g_sc23->GetErrorY(i));
        }

        c_sc23->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_AbsValue_sc23 = new TH1D("frame_AbsValue_sc23", "frame_AbsValue_sc23", 90,0,90);
        frame_AbsValue_sc23->SetMaximum(1e-6);
        frame_AbsValue_sc23->SetMinimum(-1e-6);
        frame_AbsValue_sc23->SetYTitle("Run3 - Run 2");
        frame_AbsValue_sc23->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* AbsValue_sc23 = GetAbsValue(7,(TH1D*)resultsFiles_SC23[i]->Get("SC23"),pub_sc23);
            SetMarkerAndLine(AbsValue_sc23,GetColor(index),kFullCircle,kSolid,1.0);
            AbsValue_sc23->Draw("ESAMES");
            index++;
        }
    }

    // =================
    // SC(2,4)
    // =================
    if (IfCheckObservable[kSC24]) {
        index = 0;
        TCanvas* c_sc24 = new TCanvas("c_sc24","c_sc24", 800, 1200);
        c_sc24->Divide(1,2);
        c_sc24->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_sc24 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_sc24 = new TH1D("frame_sc24", "frame_sc24", 90,0,90);
        frame_sc24->SetMaximum(5e-6);
        frame_sc24->SetMinimum(-5e-6);
        frame_sc24->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_sc24 = (TH1D*)resultsFiles_SC24[i]->Get("SC24");
            SetMarkerAndLine(h_sc24,GetColor(index),kFullCircle,kSolid,1.0);
            h_sc24->Draw("ESAMES");
            leg_sc24->AddEntry(h_sc24,Form("SC(2,4) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg_sc24->Draw();

        // TFile* publish_ins1452590 = new TFile("./PublicData/HEPData-ins1452590-v1-root.root","READ");
        TGraphAsymmErrors* g_sc24 = (TGraphAsymmErrors*)publish_ins1452590->Get("Table 1/Graph1D_y2");
        SetMarkerAndLine(g_sc24,kBlack,kOpenSquare,kSolid,1.0);
        g_sc24->Draw("PE");
        leg_sc24->AddEntry(g_sc24,Form("SC(2,4) PRL. 117 (2016) 182301"));
        leg_sc24->Draw();

        TH1D* pub_sc24 = new TH1D("pub_sc24","pub_sc24",7,x_v2);
        for(int i=0;i<7;i++){
            pub_sc24->SetBinContent(i+1,g_sc24->GetPointY(i));
            pub_sc24->SetBinError(i+1,g_sc24->GetErrorY(i));
        }

        c_sc24->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_AbsValue_sc24 = new TH1D("frame_AbsValue_sc24", "frame_AbsValue_sc24", 90,0,90);
        frame_AbsValue_sc24->SetMaximum(1e-6);
        frame_AbsValue_sc24->SetMinimum(-1e-6);
        frame_AbsValue_sc24->SetYTitle("Run3 - Run 2");
        frame_AbsValue_sc24->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* AbsValue_sc24 = GetAbsValue(7,(TH1D*)resultsFiles_SC24[i]->Get("SC24"),pub_sc24);
            SetMarkerAndLine(AbsValue_sc24,GetColor(index),kFullCircle,kSolid,1.0);
            AbsValue_sc24->Draw("ESAMES");
            index++;
        }
    }


    TFile* publish_ins1839720 = new TFile("./PublicData/HEPData-ins1839720-v1-root.root","READ");
    // =================
    // SC(2,3,4)
    // =================
    if (IfCheckObservable[kSC234]) {
        index = 0;
        TCanvas* c_sc234 = new TCanvas("c_sc234", "c_sc234", 800, 1200);
        c_sc234->Divide(1,2);
        c_sc234->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_sc234 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_sc234 = new TH1D("frame_sc234", "frame_sc234", 90,0,90);
        frame_sc234->SetMaximum(5e-9);
        frame_sc234->SetMinimum(-5e-9);
        frame_sc234->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_sc234 = (TH1D*)resultsFiles_SC234[i]->Get("SC234");
            SetMarkerAndLine(h_sc234,GetColor(index),kFullCircle,kSolid,1.0);
            h_sc234->Draw("ESAMES");
            leg_sc234->AddEntry(h_sc234,Form("SC(2,3,4) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg_sc234->Draw();

        TGraphAsymmErrors* g_sc234 = (TGraphAsymmErrors*)publish_ins1839720->Get("Table 1/Graph1D_y1");
        SetMarkerAndLine(g_sc234,kBlack,kOpenSquare,kSolid,1.0);
        g_sc234->Draw("PE");
        leg_sc234->AddEntry(g_sc234,Form("SC(2,3,4) PRL. 127 (2021) 092302"));
        leg_sc234->Draw();

        TH1D* pub_sc234 = new TH1D("pub_sc234","pub_sc234",7,x_v2);
        for(int i=0;i<7;i++){
            pub_sc234->SetBinContent(i+1,g_sc234->GetPointY(i));
            pub_sc234->SetBinError(i+1,g_sc234->GetErrorY(i));
        }

        c_sc234->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_AbsValue_sc234 = new TH1D("frame_AbsValue_sc234", "frame_AbsValue_sc234", 90,0,90);
        frame_AbsValue_sc234->SetMaximum(1e-9);
        frame_AbsValue_sc234->SetMinimum(-1e-9);
        frame_AbsValue_sc234->SetYTitle("Run3 - Run 2");
        frame_AbsValue_sc234->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* AbsValue_sc234 = GetAbsValue(7,(TH1D*)resultsFiles_SC234[i]->Get("SC234"),pub_sc234);
            SetMarkerAndLine(AbsValue_sc234,GetColor(index),kFullCircle,kSolid,1.0);
            AbsValue_sc234->Draw("ESAMES");
            index++;
        }
    }
    
    // =================
    // SC(2,3,5)
    // =================
    if (IfCheckObservable[kSC235]) {
        index = 0;
        TCanvas* c_sc235 = new TCanvas("c_sc235", "c_sc235", 800, 1200);
        c_sc235->Divide(1,2);
        c_sc235->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_sc235 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_sc235 = new TH1D("frame_sc235", "frame_sc235", 90,0,90);
        frame_sc235->SetMaximum(5e-9);
        frame_sc235->SetMinimum(-5e-9);
        frame_sc235->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_sc235 = (TH1D*)resultsFiles_SC235[i]->Get("SC235");
            SetMarkerAndLine(h_sc235,GetColor(index),kFullCircle,kSolid,1.0);
            h_sc235->Draw("ESAMES");
            leg_sc235->AddEntry(h_sc235,Form("SC(2,3,5) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg_sc235->Draw();

        TGraphAsymmErrors* g_sc235 = (TGraphAsymmErrors*)publish_ins1839720->Get("Table 2/Graph1D_y1");
        SetMarkerAndLine(g_sc235,kBlack,kOpenSquare,kSolid,1.0);
        g_sc235->Draw("PE");
        leg_sc235->AddEntry(g_sc235,Form("SC(2,3,5) PRL. 127 (2021) 092302"));
        leg_sc235->Draw();

        TH1D* pub_sc235 = new TH1D("pub_sc235","pub_sc235",7,x_v2);
        for(int i=0;i<7;i++){
            pub_sc235->SetBinContent(i+1,g_sc235->GetPointY(i));
            pub_sc235->SetBinError(i+1,g_sc235->GetErrorY(i));
        }

        c_sc235->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_AbsValue_sc235 = new TH1D("frame_AbsValue_sc235", "frame_AbsValue_sc235", 90,0,90);
        frame_AbsValue_sc235->SetMaximum(1e-9);
        frame_AbsValue_sc235->SetMinimum(-1e-9);
        frame_AbsValue_sc235->SetYTitle("Run3 - Run 2");
        frame_AbsValue_sc235->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* AbsValue_sc235 = GetAbsValue(7,(TH1D*)resultsFiles_SC235[i]->Get("SC235"),pub_sc235);
            SetMarkerAndLine(AbsValue_sc235,GetColor(index),kFullCircle,kSolid,1.0);
            AbsValue_sc235->Draw("ESAMES");
            index++;
        }
    }

    // =================
    // SC(2,4,6)
    // =================
    if (IfCheckObservable[kSC246]) {
        index = 0;
        TCanvas* c_sc246 = new TCanvas("c_sc246", "c_sc246", 800, 1200);
        c_sc246->Divide(1,2);
        c_sc246->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_sc246 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_sc246 = new TH1D("frame_sc246", "frame_sc246", 90,0,90);
        frame_sc246->SetMaximum(5e-9);
        frame_sc246->SetMinimum(-5e-9);
        frame_sc246->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_sc246 = (TH1D*)resultsFiles_SC246[i]->Get("SC246");
            SetMarkerAndLine(h_sc246,GetColor(index),kFullCircle,kSolid,1.0);
            h_sc246->Draw("ESAMES");
            leg_sc246->AddEntry(h_sc246,Form("SC(2,4,6) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg_sc246->Draw();

        TGraphAsymmErrors* g_sc246 = (TGraphAsymmErrors*)publish_ins1839720->Get("Table 3/Graph1D_y1");
        SetMarkerAndLine(g_sc246,kBlack,kOpenSquare,kSolid,1.0);
        g_sc246->Draw("PE");
        leg_sc246->AddEntry(g_sc246,Form("SC(2,4,6) PRL. 127 (2021) 092302"));
        leg_sc246->Draw();

        TH1D* pub_sc246 = new TH1D("pub_sc246","pub_sc246",7,x_v2);
        for(int i=0;i<7;i++){
            pub_sc246->SetBinContent(i+1,g_sc246->GetPointY(i));
            pub_sc246->SetBinError(i+1,g_sc246->GetErrorY(i));
        }

        c_sc246->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_AbsValue_sc246 = new TH1D("frame_AbsValue_sc246", "frame_AbsValue_sc246", 90,0,90);
        frame_AbsValue_sc246->SetMaximum(1e-9);
        frame_AbsValue_sc246->SetMinimum(-1e-9);
        frame_AbsValue_sc246->SetYTitle("Run3 - Run 2");
        frame_AbsValue_sc246->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* AbsValue_sc246 = GetAbsValue(7,(TH1D*)resultsFiles_SC246[i]->Get("SC246"),pub_sc246);
            SetMarkerAndLine(AbsValue_sc246,GetColor(index),kFullCircle,kSolid,1.0);
            AbsValue_sc246->Draw("ESAMES");
            index++;
        }
    }

    // =================
    // SC(3,4,5)
    // =================
    if (IfCheckObservable[kSC345]) {
        index = 0;
        TCanvas* c_sc345 = new TCanvas("c_sc345", "c_sc345", 800, 1200);
        c_sc345->Divide(1,2);
        c_sc345->cd(1);
        gPad->SetBottomMargin(0.05);
        TLegend* leg_sc345 = new TLegend(0.2,0.7,0.5,0.9);
        TH1D* frame_sc345 = new TH1D("frame_sc345", "frame_sc345", 90,0,90);
        frame_sc345->SetMaximum(5e-9);
        frame_sc345->SetMinimum(-5e-9);
        frame_sc345->Draw("AXIS");
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* h_sc345 = (TH1D*)resultsFiles_SC345[i]->Get("SC345");
            SetMarkerAndLine(h_sc345,GetColor(index),kFullCircle,kSolid,1.0);
            h_sc345->Draw("ESAMES");
            leg_sc345->AddEntry(h_sc345,Form("SC(3,4,5) (%s)",legendNames[i].c_str()),"lp");
            index+=1;
        }
        leg_sc345->Draw();

        TGraphAsymmErrors* g_sc345 = (TGraphAsymmErrors*)publish_ins1839720->Get("Table 4/Graph1D_y1");
        SetMarkerAndLine(g_sc345,kBlack,kOpenSquare,kSolid,1.0);
        g_sc345->Draw("PE");
        leg_sc345->AddEntry(g_sc345,Form("SC(3,4,5) PRL. 127 (2021) 092302"));
        leg_sc345->Draw();

        TH1D* pub_sc345 = new TH1D("pub_sc345","pub_sc345",7,x_v2);
        for(int i=0;i<7;i++){
            pub_sc345->SetBinContent(i+1,g_sc345->GetPointY(i));
            pub_sc345->SetBinError(i+1,g_sc345->GetErrorY(i));
        }

        c_sc345->cd(2);
        gPad->SetTopMargin(0.05);
        TH1D* frame_AbsValue_sc345 = new TH1D("frame_AbsValue_sc345", "frame_AbsValue_sc345", 90,0,90);
        frame_AbsValue_sc345->SetMaximum(1e-9);
        frame_AbsValue_sc345->SetMinimum(-1e-9);
        frame_AbsValue_sc345->SetYTitle("Run3 - Run 2");
        frame_AbsValue_sc345->Draw("AXIS");
        One->Draw("sames");
        index=0;
        for(int i=0;i<FileNameSuffixs.size();i++){
            TH1D* AbsValue_sc345 = GetAbsValue(7,(TH1D*)resultsFiles_SC345[i]->Get("SC345"),pub_sc345);
            SetMarkerAndLine(AbsValue_sc345,GetColor(index),kFullCircle,kSolid,1.0);
            AbsValue_sc345->Draw("ESAMES");
            index++;
        }
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
    }


    return;
    

}
