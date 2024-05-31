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
        
void ProcessFlowContainerDrawDiffRoot(){
    vector<string> FileNameSuffixs;
    vector<TFile*> resultsFiles_vn;
    vector<TFile*> resultsFiles_v422;
    vector<TFile*> resultsFiles_NSC;
    FileNameSuffixs.push_back("LHC23zzh_pass2_192487");
    FileNameSuffixs.push_back("LHC23zzn_pass3_I_A11_small_217816");

    for(auto& suffix : FileNameSuffixs){
        string fileName_vn = "./ProcessOutput/vn_" + suffix + ".root";
        string fileName_v422 = "./ProcessOutput/v422_" + suffix + ".root";
        string fileName_NSC = "./ProcessOutput/NSC_" + suffix + ".root";

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
        resultsFile = TFile::Open(fileName_NSC.c_str(), "READ");
        if(!resultsFile || resultsFile->IsZombie()){
            cout << "Error: cannot open file " << fileName_NSC << endl;
            return;
        }
        resultsFiles_NSC.push_back(resultsFile);
    }

    int index = 0;
    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 800);
    TLegend* leg = new TLegend(0.2,0.7,0.5,0.9);
    TH1D* frame_vn = new TH1D("frame_vn", "frame_vn", 90,0,90);
    frame_vn->SetMaximum(0.15);
    frame_vn->SetMinimum(0.);
    frame_vn->Draw("AXIS");
    for(int i=0;i<FileNameSuffixs.size();i++){
        TH1D* h_v22 = (TH1D*)resultsFiles_vn[i]->Get("Corr_corr_22_hist");
        SetMarkerAndLine(h_v22,GetColor(index),kFullCircle,kSolid,1.0);
        h_v22->Draw("ESames");
        leg->AddEntry(h_v22,Form("v_{2}{2} (%s)",FileNameSuffixs[i].c_str()),"lp");
        TH1D* h_v32 = (TH1D*)resultsFiles_vn[i]->Get("Corr_corr_32_hist");
        SetMarkerAndLine(h_v32,GetColor(index+1),kFullCircle,kSolid,1.0);
        h_v32->Draw("ESames");
        leg->AddEntry(h_v32,Form("v_{3}{2} (%s)",FileNameSuffixs[i].c_str()),"lp");
        // TH1D* h_v42 = (TH1D*)resultsFiles_vn[i]->Get("Corr_corr_42_hist");
        // SetMarkerAndLine(h_v42,GetColor(index+2),kFullCircle,kSolid,1.0);
        // leg->AddEntry(h_v42,Form("v_{4}{2} (%s)",FileNameSuffixs[i].c_str()),"lp");
        index+=3;
    }
    if(ComparewithPublish){
        TFile* publish = new TFile("./HEPData-ins1778342-v1-root.root","READ");
        TGraphAsymmErrors* g_v2 = (TGraphAsymmErrors*)publish->Get("v2/Graph1D_y1");
        SetMarkerAndLine(g_v2,kBlack,kOpenSquare,kSolid,1.0);
        TGraphAsymmErrors* g_v3 = (TGraphAsymmErrors*)publish->Get("v3/Graph1D_y1");
        SetMarkerAndLine(g_v3,kRed,kOpenSquare,kSolid,1.0);
        TGraphAsymmErrors* g_v4 = (TGraphAsymmErrors*)publish->Get("v4/Graph1D_y1");
        SetMarkerAndLine(g_v4,kBlue,kOpenSquare,kSolid,1.0);
        g_v2->Draw("PE");
        g_v3->Draw("PE");
        // g_v4->Draw("PE");
        leg->AddEntry(g_v2,Form("v_{2}{2} JHEP 05 (2020) 085, 2020"));
        leg->AddEntry(g_v3,Form("v_{3}{2} JHEP 05 (2020) 085, 2020"));
        // leg->AddEntry(g_v4,Form("v_{4}{2} JHEP 05 (2020) 085, 2020"));
    }
    leg->Draw();

    index = 0;
    TCanvas* c2 = new TCanvas("c2", "c2", 1200, 800);
    TLegend* leg2 = new TLegend(0.2,0.7,0.5,0.9);
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

     if(ComparewithPublish){
        TFile* publish = new TFile("./HEPData-ins1778342-v1-root.root","READ");
        TGraphAsymmErrors* g = nullptr;
        // if(observable==v422)
        g=(TGraphAsymmErrors*)publish->Get("v422/Graph1D_y1");
        // else if(observable==chi422)g=(TGraphAsymmErrors*)publish->Get("chi422/Graph1D_y1");
        // else if(observable==rho422)g=(TGraphAsymmErrors*)publish->Get("rho422/Graph1D_y1");
        if(!g)return;
        SetMarkerAndLine(g,kBlack,kOpenSquare,kSolid,1.0);
        g->Draw("PE");
        leg2->AddEntry(g,Form("v_{4,22} JHEP 05 (2020) 085, 2020"));
        // legend2->AddEntry(g,Form("%s JHEP 05 (2020) 085, 2020",ObservableName[observable]));
    }
    leg2->Draw();

    index = 0;
    TCanvas* c3 = new TCanvas("c3", "c3", 1200, 800);
    TLegend* leg3 = new TLegend(0.2,0.7,0.5,0.9);
    TH1D* frame_NSC = new TH1D("frame_NSC", "frame_NSC", 90,0,90);
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
    if(ComparewithPublish){
        TFile* publish = new TFile("./HEPData-ins1848215-v1-root.root","READ");
        TGraphAsymmErrors* g = nullptr;
        g=(TGraphAsymmErrors*)publish->Get("Table 1/Graph1D_y1");
        if(!g)return;
        SetMarkerAndLine(g,kBlack,kOpenSquare,kSolid,1.0);
        g->Draw("PE");
        leg3->AddEntry(g,Form("NSC(3,2) PLB 818 (2021) 136354, 2021"));
    }
    leg3->Draw();



    return;
    

}
