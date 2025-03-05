/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch) 
 * @Date: 2025-03-05 22:01:48 
 * @Last Modified by:   Zhiyong Lu 
 * @Last Modified time: 2025-03-05 22:01:48 
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
#include "TGraphErrors.h"
#include "TF1.h"
#include <cstring>
#include <vector>
#include <map>
#include <array>
#include "include/ErrorPropagation.h"
#include "include/ProcessDefine.h"

using namespace std;

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
        case 0: return kBlue;
        case 1: return kRed;
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

std::map<kObservable,bool> IfCheckObservable = {
    {kVn,false},
    {kV24,false},
    {kV26,false},
    {kV28,false},
    {kV210,false},
    {kV422,false},
    {kChi422,false},
    {kRho422,false},
    {kNSC23,false},
    {kNSC24,false},
    {kNSC234,false},
    {kNSC345,false},
    {kpTDiffv2,false},
    {kpTDiffv3,false},
    {kpTDiffv4,false},
    {kpTDiffv24ChFull,true},
    {kpTDiffv24Ch10Gap,false},
    {kpTDiffv26ChFull,false}
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
    vector<vector<TFile*>> resultsFiles_pTDiffv2(pTDiffCent.size());
    vector<vector<TFile*>> resultsFiles_pTDiffv3(pTDiffCent.size());
    vector<vector<TFile*>> resultsFiles_pTDiffv4(pTDiffCent.size());
    vector<vector<TFile*>> resultsFiles_pTDiffv24ChFull(pTDiffCent.size());
    vector<vector<TFile*>> resultsFiles_pTDiffv24Ch10Gap(pTDiffCent.size());
    vector<vector<TFile*>> resultsFiles_pTDiffv26ChFull(pTDiffCent.size());

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
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_299557");
    // legendNames.push_back("full 2023 data");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_308058");
    // legendNames.push_back("2023 data");
    // FileNameSuffixs.push_back("LHC24ar_pass1_medium_318934");
    // legendNames.push_back("2024 data");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815_kIsGoodITSLayersAll");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815_kIsGoodITSLayersAll");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815_kIsGoodITSLayersAllnonclosurecorrected");
    // legendNames.push_back("default+ non-closure corrected");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815_GoodZvtxFT0vsPV");
    // legendNames.push_back("default+GoodZvtxFT0vsPV");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815_kNoCollInTimeRangeStandard");
    // legendNames.push_back("default+kNoCollInTimeRangeStandard");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815_GoodZvtxFT0vsPV");
    // legendNames.push_back("default+GoodZvtxFT0vsPV");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_327815_kNoCollInTimeRangeStandard");
    // legendNames.push_back("default+kNoCollInTimeRangeStandard");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_328079");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_328079_SelSwitch");
    // legendNames.push_back("default_SelSwitch");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_328079_kIsGoodITSLayersAll");
    // legendNames.push_back("default+kIsGoodITSLayersAll");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_328079_ITS");
    // legendNames.push_back("default+kIsGoodITSLayersAll+ITSAll7Layers");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_332826_kIsGoodITSLayersAll");
    // legendNames.push_back("full 2023 data+kIsGoodITSLayersAll");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_334605_kIsGoodITSLayersAll");
    // legendNames.push_back("default+NoColl+ITSclu");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_340172_kColl");
    // legendNames.push_back("default+NoColl+ITSclu");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_340440");
    // legendNames.push_back("centrality estimator: FT0C");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_340440_FT0M");
    // legendNames.push_back("centrality estimator: FT0M");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_340440");
    // legendNames.push_back("default (small dataset)");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_341269");
    // legendNames.push_back("default (full 2023, NUA not fitted)");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_341269_kColl");
    // legendNames.push_back("add kColl cuts");
    // FileNameSuffixs.push_back("LHC23_PbPb_pass4_332826_kIsGoodITSLayersAll");
    // legendNames.push_back("Max occupancy 3k (full 2023, NUA not fitted)");
    FileNameSuffixs.push_back("LHC23_PbPb_pass4_344339");
    legendNames.push_back("default (full 2023)");
    // legendNames.push_back("default maxInactiveChips:{8,8,8,111,111,195,195}");
    // legendNames.push_back("Global track");
    FileNameSuffixs.push_back("LHC23_PbPb_pass4_344339nonclosurecorrected"); //only for pTDiffv22
    legendNames.push_back("with MC non-closure corrected");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_355881_AverageNUA");
    // legendNames.push_back("average NUA");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_355881_RunByRunNUA");
    // legendNames.push_back("run-by-run NUA");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_356028");
    // legendNames.push_back("default");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_356028_ITSclu0");
    // legendNames.push_back("ITSclu>0");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_356028_occupancy5k");
    // legendNames.push_back("occupancy<5k");
    // FileNameSuffixs.push_back("LHC24ar_pass2_QC4_360475");
    // legendNames.push_back("LHC24ar_pass2_QC4");

    // FileNameSuffixs.push_back("LHC23zzm_pass4_calo_359360_QualityITS");
    // legendNames.push_back("isQualityTrackITS");
    // FileNameSuffixs.push_back("LHC23zzm_pass4_calo_359360_ITS");
    // legendNames.push_back("isQualityTrackITS+Primary+Acceptance");

    // FileNameSuffixs.push_back("LHC23zzh_pass4_359613_ITSclu0");
    // legendNames.push_back("LHC23zzh_pass4");
    // legendNames.push_back("# of ITS clusters > 0");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_359613_NokIsGoodITSLayersAll");
    // legendNames.push_back("without kIsGoodITSLayersAll");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_359613_occupancy1w");
    // legendNames.push_back("occupancy < 10k");
    

    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_361919");
    // legendNames.push_back("maxInactiveChips:{8,8,8,80,80,150,150}");
    // FileNameSuffixs.push_back("LHC23zzh_pass4_small_361917");
    // legendNames.push_back("maxInactiveChips:{8,8,8,50,50,100,100}");




    if (FileNameSuffixs.empty() || legendNames.empty()) {
        Printf("Error: no input files specified.");
        return;
    }
    int pTCentIndex=0;
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
        vector<string> fileName_pTDiffv2(pTDiffCent.size()-1);
        vector<string> fileName_pTDiffv3(pTDiffCent.size()-1);
        vector<string> fileName_pTDiffv4(pTDiffCent.size()-1);
        vector<string> fileName_pTDiffv24ChFull(pTDiffCent.size()-1);
        vector<string> fileName_pTDiffv24Ch10Gap(pTDiffCent.size()-1);
        vector<string> fileName_pTDiffv26ChFull(pTDiffCent.size()-1);
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            string temp_s = Form("./ProcessOutput/pTDiffv2Cent%dTo%d_",(int)pTDiffCent[j],(int)pTDiffCent[j+1]) + suffix + ".root";
            fileName_pTDiffv2[j] = temp_s;
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            string temp_s = Form("./ProcessOutput/pTDiffv3Cent%dTo%d_",(int)pTDiffCent[j],(int)pTDiffCent[j+1]) + suffix + ".root";
            fileName_pTDiffv3[j] = temp_s;
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            string temp_s = Form("./ProcessOutput/pTDiffv4Cent%dTo%d_",(int)pTDiffCent[j],(int)pTDiffCent[j+1]) + suffix + ".root";
            fileName_pTDiffv4[j] = temp_s;
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            string temp_s = Form("./ProcessOutput/pTDiffv24ChFullCent%dTo%d_",(int)pTDiffCent[j],(int)pTDiffCent[j+1]) + suffix + ".root";
            fileName_pTDiffv24ChFull[j] = temp_s;
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            string temp_s = Form("./ProcessOutput/pTDiffv24Ch10GapCent%dTo%d_",(int)pTDiffCent[j],(int)pTDiffCent[j+1]) + suffix + ".root";
            fileName_pTDiffv24Ch10Gap[j] = temp_s;
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            string temp_s = Form("./ProcessOutput/pTDiffv26ChFullCent%dTo%d_",(int)pTDiffCent[j],(int)pTDiffCent[j+1]) + suffix + ".root";
            fileName_pTDiffv26ChFull[j] = temp_s;
        }

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
        if (IfCheckObservable[kpTDiffv2]) {
            for (uint j=0; j<pTDiffCent.size()-1; j++) {
                resultsFile = TFile::Open(fileName_pTDiffv2[j].c_str(), "READ");
                if(!resultsFile || resultsFile->IsZombie()){
                    cout << "Error: cannot open file " << fileName_pTDiffv2[j] << endl;
                    return;
                }
                resultsFiles_pTDiffv2[j].push_back(resultsFile);
            }
        }
        if (IfCheckObservable[kpTDiffv3]) {
            for (uint j=0; j<pTDiffCent.size()-1; j++) {
                resultsFile = TFile::Open(fileName_pTDiffv3[j].c_str(), "READ");
                if(!resultsFile || resultsFile->IsZombie()){
                    cout << "Error: cannot open file " << fileName_pTDiffv3[j] << endl;
                    return;
                }
                resultsFiles_pTDiffv3[j].push_back(resultsFile);
            }
        }
        if (IfCheckObservable[kpTDiffv4]) {
            for (uint j=0; j<pTDiffCent.size()-1; j++) {
                resultsFile = TFile::Open(fileName_pTDiffv4[j].c_str(), "READ");
                if(!resultsFile || resultsFile->IsZombie()){
                    cout << "Error: cannot open file " << fileName_pTDiffv4[j] << endl;
                    return;
                }
                resultsFiles_pTDiffv4[j].push_back(resultsFile);
            }
        }
        if (IfCheckObservable[kpTDiffv24ChFull]) {
            for (uint j=0; j<pTDiffCent.size()-1; j++) {
                resultsFile = TFile::Open(fileName_pTDiffv24ChFull[j].c_str(), "READ");
                if(!resultsFile || resultsFile->IsZombie()){
                    cout << "Error: cannot open file " << fileName_pTDiffv24ChFull[j] << endl;
                    return;
                }
                resultsFiles_pTDiffv24ChFull[j].push_back(resultsFile);
            }
        }
        if (IfCheckObservable[kpTDiffv24Ch10Gap]) {
            for (uint j=0; j<pTDiffCent.size()-1; j++) {
                resultsFile = TFile::Open(fileName_pTDiffv24Ch10Gap[j].c_str(), "READ");
                if(!resultsFile || resultsFile->IsZombie()){
                    cout << "Error: cannot open file " << fileName_pTDiffv24Ch10Gap[j] << endl;
                    return;
                }
                resultsFiles_pTDiffv24Ch10Gap[j].push_back(resultsFile);
            }
        }
        if (IfCheckObservable[kpTDiffv26ChFull]) {
            for (uint j=0; j<pTDiffCent.size()-1; j++) {
                resultsFile = TFile::Open(fileName_pTDiffv26ChFull[j].c_str(), "READ");
                if(!resultsFile || resultsFile->IsZombie()){
                    cout << "Error: cannot open file " << fileName_pTDiffv26ChFull[j] << endl;
                    return;
                }
                resultsFiles_pTDiffv26ChFull[j].push_back(resultsFile);
            }
        }
        pTCentIndex++;
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
        
        TH1D* pub_v422 = new TH1D("pub_v422","pub_v422",9,x_v2);
        for(int i=0;i<9;i++){
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
            TH1D* ratio_v422 = GetRatio(9,(TH1D*)resultsFiles_v422[i]->Get("hCorr422_mean"),pub_v422);
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

        TH1D* pub_chi422 = new TH1D("pub_chi422","pub_chi422",9,x_v2);
        for(int i=0;i<9;i++){
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
            TH1D* ratio_chi422 = GetRatio(9,(TH1D*)resultsFiles_chi422[i]->Get("hCorr422_mean"),pub_chi422);
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

        TH1D* pub_rho422 = new TH1D("pub_rho422","pub_rho422",9,x_v2);
        for(int i=0;i<9;i++){
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
            TH1D* ratio_rho422 = GetRatio(9,(TH1D*)resultsFiles_rho422[i]->Get("hCorr422_mean"),pub_rho422);
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
    // pTDiffv2
    // =================
    vector<int> PubTable_pTDiffv2 = {15, 23, 33, 43, 53, 63, 73, 83};
    if (IfCheckObservable[kpTDiffv2]) {
        for (uint j=0;j<pTDiffCent.size()-1;j++) {
            index = 0;
            TCanvas* cpTDiffv2 = new TCanvas(Form("cpTDiffv2_%d",j), Form("cpTDiffv2_%d",j), 800, 1200);
            cpTDiffv2->Divide(1,2);
            cpTDiffv2->cd(1);
            gPad->SetBottomMargin(0.05);
            TLegend* leg_pTDiffv2 = new TLegend(0.2,0.7,0.9,0.9);
            TH1D* frame_pTDiffv2 = new TH1D(Form("frame_pTDiffv2_%d",j), Form("frame_pTDiffv2_%d",j),200,-0.8,200);
            frame_pTDiffv2->SetMaximum(0.3);
            frame_pTDiffv2->SetMinimum(0.);
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_pTDiffv2->Draw("AXIS");
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* h_pTDiffv2 = (TH1D*)resultsFiles_pTDiffv2[j][i]->Get(Form("pTDiffv2"));
                SetMarkerAndLine(h_pTDiffv2,GetColor(index),kFullCircle,kSolid,1.0);
                h_pTDiffv2->Draw("ESAMES");
                leg_pTDiffv2->AddEntry(h_pTDiffv2,Form("v_{2}{2,|#Delta#eta|>1.0}(p_{T}) Cent:%d~%d%% (%s)",(int)pTDiffCent[j],(int)pTDiffCent[j+1],legendNames[i].c_str()),"lp");
                index+=1;
            }

            TGraphAsymmErrors* g_pTDiffv2 = nullptr;
            if (j<PubTable_pTDiffv2.size()) {
                g_pTDiffv2 = (TGraphAsymmErrors*)publish_ins1666817->Get(Form("Table %d/Graph1D_y1",PubTable_pTDiffv2[j]));
                if(!g_pTDiffv2) continue;
                SetMarkerAndLine(g_pTDiffv2,kBlack,kOpenSquare,kSolid,1.0);
                g_pTDiffv2->Draw("PE");
                leg_pTDiffv2->AddEntry(g_pTDiffv2,Form("JHEP 07 (2018) 103"));
            }
            leg_pTDiffv2->Draw();

            if(!g_pTDiffv2) continue;
            TH1D* pub_pTDiffv2 = new TH1D(Form("pub_pTDiffv2_%d",j),Form("pub_pTDiffv2_%d",j),28,x_pTDiff);
            for(int i=0;i<28;i++){
                pub_pTDiffv2->SetBinContent(i+1,g_pTDiffv2->GetPointY(i));
                pub_pTDiffv2->SetBinError(i+1,g_pTDiffv2->GetErrorY(i));
            }

            cpTDiffv2->cd(2);
            gPad->SetTopMargin(0.05);
            TH1D* frame_ratio_pTDiffv2 = new TH1D(Form("frame_ratio_pTDiffv2_%d",j), Form("frame_ratio_pTDiffv2_%d",j),200,-0.8,200);
            frame_ratio_pTDiffv2->SetMaximum(1.4);
            frame_ratio_pTDiffv2->SetMinimum(0.6);
            frame_ratio_pTDiffv2->SetYTitle("Ratio to Run 2");
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_ratio_pTDiffv2->Draw("AXIS");
            One->Draw("sames");
            index=0;
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* ratio_pTDiffv2 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv2[j][i]->Get(Form("pTDiffv2")),pub_pTDiffv2,true);
                SetMarkerAndLine(ratio_pTDiffv2,GetColor(index),kFullCircle,kSolid,1.0);
                ratio_pTDiffv2->Draw("ESAMES");
                index++;
            }
            cpTDiffv2->SaveAs(Form("./OutputPDF/pTDiffv2Cent%dTo%d_RatioToRun2.pdf",(int)pTDiffCent[j],(int)pTDiffCent[j+1]));
        }
    }


    // =================
    // pTDiffv3
    // =================
    vector<int> PubTable_pTDiffv3 = {17, 25, 35, 45, 55, 65, 75, 85};
    if (IfCheckObservable[kpTDiffv3]) {
        for (uint j=0;j<pTDiffCent.size()-1;j++) {
            index = 0;
            TCanvas* cpTDiffv3 = new TCanvas(Form("cpTDiffv3_%d",j), Form("cpTDiffv3_%d",j), 800, 1200);
            cpTDiffv3->Divide(1,2);
            cpTDiffv3->cd(1);
            gPad->SetBottomMargin(0.05);
            TLegend* leg_pTDiffv3 = new TLegend(0.2,0.7,0.9,0.9);
            TH1D* frame_pTDiffv3 = new TH1D(Form("frame_pTDiffv3_%d",j), Form("frame_pTDiffv3_%d",j), 200,-0.8,200);
            frame_pTDiffv3->SetMaximum(0.15);
            frame_pTDiffv3->SetMinimum(0.);
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_pTDiffv3->Draw("AXIS");
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* h_pTDiffv3 = (TH1D*)resultsFiles_pTDiffv3[j][i]->Get(Form("pTDiffv3"));
                SetMarkerAndLine(h_pTDiffv3,GetColor(index),kFullCircle,kSolid,1.0);
                h_pTDiffv3->Draw("ESAMES");
                leg_pTDiffv3->AddEntry(h_pTDiffv3,Form("v_{3}{2,|#Delta#eta|>1.0}(p_{T}) Cent:%d~%d%% (%s)",(int)pTDiffCent[j],(int)pTDiffCent[j+1],legendNames[i].c_str()),"lp");
                index+=1;
            }

            TGraphAsymmErrors* g_pTDiffv3 = nullptr;
            if (j<PubTable_pTDiffv3.size()) {
                g_pTDiffv3 = (TGraphAsymmErrors*)publish_ins1666817->Get(Form("Table %d/Graph1D_y1",PubTable_pTDiffv3[j]));
                if(!g_pTDiffv3) continue;
                SetMarkerAndLine(g_pTDiffv3,kBlack,kOpenSquare,kSolid,1.0);
                g_pTDiffv3->Draw("PE");
                leg_pTDiffv3->AddEntry(g_pTDiffv3,Form("JHEP 07 (2018) 103"));
            }
            leg_pTDiffv3->Draw();

            if(!g_pTDiffv3) continue;
            TH1D* pub_pTDiffv3 = new TH1D(Form("pub_pTDiffv3_%d",j),Form("pub_pTDiffv3_%d",j),20,x_pTDiff);
            for(int i=0;i<28;i++){
                pub_pTDiffv3->SetBinContent(i+1,g_pTDiffv3->GetPointY(i));
                pub_pTDiffv3->SetBinError(i+1,g_pTDiffv3->GetErrorY(i));
            }

            cpTDiffv3->cd(2);
            gPad->SetTopMargin(0.05);
            TH1D* frame_ratio_pTDiffv3 = new TH1D(Form("frame_ratio_pTDiffv3_%d",j), Form("frame_ratio_pTDiffv3_%d",j),200,-0.8,200);
            frame_ratio_pTDiffv3->SetMaximum(1.4);
            frame_ratio_pTDiffv3->SetMinimum(0.6);
            frame_ratio_pTDiffv3->SetYTitle("Ratio to Run 2");
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_ratio_pTDiffv3->Draw("AXIS");
            One->Draw("sames");
            index=0;
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* ratio_pTDiffv3 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv3[j][i]->Get(Form("pTDiffv3")),pub_pTDiffv3,true);
                SetMarkerAndLine(ratio_pTDiffv3,GetColor(index),kFullCircle,kSolid,1.0);
                ratio_pTDiffv3->Draw("ESAMES");
                index++;
            }
            cpTDiffv3->SaveAs(Form("./OutputPDF/pTDiffv3Cent%dTo%d_RatioToRun2.pdf",(int)pTDiffCent[j],(int)pTDiffCent[j+1]));
        }
    }

    // =================
    // pTDiffv4
    // =================
    vector<int> PubTable_pTDiffv4 = {19, 27, 37, 47, 57, 67, 77, 87};
    if (IfCheckObservable[kpTDiffv4]) {
        for (uint j=0;j<pTDiffCent.size()-1;j++) {
            index = 0;
            TCanvas* cpTDiffv4 = new TCanvas(Form("cpTDiffv4_%d",j), Form("cpTDiffv4_%d",j), 800, 1200);
            cpTDiffv4->Divide(1,2);
            cpTDiffv4->cd(1);
            gPad->SetBottomMargin(0.05);
            TLegend* leg_pTDiffv4 = new TLegend(0.2,0.7,0.9,0.9);
            TH1D* frame_pTDiffv4 = new TH1D(Form("frame_pTDiffv4_%d",j), Form("frame_pTDiffv4_%d",j), 200,-0.8,200);
            frame_pTDiffv4->SetMaximum(0.15);
            frame_pTDiffv4->SetMinimum(0.);
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_pTDiffv4->Draw("AXIS");
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* h_pTDiffv4 = (TH1D*)resultsFiles_pTDiffv4[j][i]->Get(Form("pTDiffv4"));
                SetMarkerAndLine(h_pTDiffv4,GetColor(index),kFullCircle,kSolid,1.0);
                h_pTDiffv4->Draw("ESAMES");
                leg_pTDiffv4->AddEntry(h_pTDiffv4,Form("v_{4}{2,|#Delta#eta|>1.0}(p_{T}) Cent:%d~%d%% (%s)",(int)pTDiffCent[j],(int)pTDiffCent[j+1],legendNames[i].c_str()),"lp");
                index+=1;
            }

            TGraphAsymmErrors* g_pTDiffv4 = nullptr;
            if (j<PubTable_pTDiffv4.size()) {
                g_pTDiffv4 = (TGraphAsymmErrors*)publish_ins1666817->Get(Form("Table %d/Graph1D_y1",PubTable_pTDiffv4[j]));
                if(!g_pTDiffv4) continue;
                SetMarkerAndLine(g_pTDiffv4,kBlack,kOpenSquare,kSolid,1.0);
                g_pTDiffv4->Draw("PE");
                leg_pTDiffv4->AddEntry(g_pTDiffv4,Form("JHEP 07 (2018) 103"));
            }
            leg_pTDiffv4->Draw();

            if(!g_pTDiffv4) continue;
            TH1D* pub_pTDiffv4 = new TH1D(Form("pub_pTDiffv4_%d",j),Form("pub_pTDiffv4_%d",j),20,x_pTDiff);
            for(int i=0;i<28;i++){
                pub_pTDiffv4->SetBinContent(i+1,g_pTDiffv4->GetPointY(i));
                pub_pTDiffv4->SetBinError(i+1,g_pTDiffv4->GetErrorY(i));
            }

            cpTDiffv4->cd(2);
            gPad->SetTopMargin(0.05);
            TH1D* frame_ratio_pTDiffv4 = new TH1D(Form("frame_ratio_pTDiffv4_%d",j), Form("frame_ratio_pTDiffv4_%d",j),200,-0.8,200);
            frame_ratio_pTDiffv4->SetMaximum(1.4);
            frame_ratio_pTDiffv4->SetMinimum(0.6);
            frame_ratio_pTDiffv4->SetYTitle("Ratio to Run 2");
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_ratio_pTDiffv4->Draw("AXIS");
            One->Draw("sames");
            index=0;
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* ratio_pTDiffv4 = GetRatio(28,(TH1D*)resultsFiles_pTDiffv4[j][i]->Get(Form("pTDiffv4")),pub_pTDiffv4,true);
                SetMarkerAndLine(ratio_pTDiffv4,GetColor(index),kFullCircle,kSolid,1.0);
                ratio_pTDiffv4->Draw("ESAMES");
                index++;
            }
            cpTDiffv4->SaveAs(Form("./OutputPDF/pTDiffv4Cent%dTo%d_RatioToRun2.pdf",(int)pTDiffCent[j],(int)pTDiffCent[j+1]));
        }
    }

    // =================
    // pTDiffv24ChFull
    // =================
    // TFile* publish_ins2093750 = new TFile("./PublicData/HEPData-ins2093750-v1-root.root","READ");
    // vector<int> PubTable_pTDiffv24 = {0, 0, 37, 46, 55, 64, 0, 0};
    vector<int> PubTable_pTDiffv24 = {0, 31, 41, 51, 61, 71, 80, 87};
    if (IfCheckObservable[kpTDiffv24ChFull]) {
        for (uint j=0;j<pTDiffCent.size()-1;j++) {
            index = 0;
            TCanvas* cpTDiffv24ChFull = new TCanvas(Form("cpTDiffv24ChFull_%d",j), Form("cpTDiffv24ChFull_%d",j), 800, 1200);
            cpTDiffv24ChFull->Divide(1,2);
            cpTDiffv24ChFull->cd(1);
            gPad->SetBottomMargin(0.05);
            TLegend* leg_pTDiffv24ChFull = new TLegend(0.2,0.7,0.9,0.9);
            TH1D* frame_pTDiffv24ChFull = new TH1D(Form("frame_pTDiffv24ChFull_%d",j), Form("frame_pTDiffv24ChFull_%d",j), 200,-0.8,200);
            frame_pTDiffv24ChFull->SetMaximum(0.3);
            frame_pTDiffv24ChFull->SetMinimum(0.);
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_pTDiffv24ChFull->Draw("AXIS");
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* h_pTDiffv24ChFull = (TH1D*)resultsFiles_pTDiffv24ChFull[j][i]->Get(Form("pTDiffv24ChFull"));
                SetMarkerAndLine(h_pTDiffv24ChFull,GetColor(index),kFullCircle,kSolid,1.0);
                h_pTDiffv24ChFull->Draw("ESAMES");
                leg_pTDiffv24ChFull->AddEntry(h_pTDiffv24ChFull,Form("v_{2}{4} (p_{T}) Cent:%d~%d%% (%s)",(int)pTDiffCent[j],(int)pTDiffCent[j+1],legendNames[i].c_str()),"lp");
                index+=1;
            }

            TGraphAsymmErrors* g_pTDiffv24ChFull = nullptr;
            if (j<PubTable_pTDiffv2.size()) {
                if (!PubTable_pTDiffv24[j]) {
                    Printf("No corresponding table for pTDiffv24ChFull, skip this plot");
                    leg_pTDiffv24ChFull->Draw();
                    continue;
                }
                g_pTDiffv24ChFull = (TGraphAsymmErrors*)publish_ins1666817->Get(Form("Table %d/Graph1D_y1",PubTable_pTDiffv24[j]));
                if(!g_pTDiffv24ChFull) continue;
                SetMarkerAndLine(g_pTDiffv24ChFull,kBlack,kOpenSquare,kSolid,1.0);
                g_pTDiffv24ChFull->Draw("PE");
                leg_pTDiffv24ChFull->AddEntry(g_pTDiffv24ChFull,Form("v2{4}(p_{T}) JHEP 07 (2018) 103"));
            }
            leg_pTDiffv24ChFull->Draw();

            if(!g_pTDiffv24ChFull) continue;
            TH1D* pub_pTDiffv24ChFull = new TH1D(Form("pub_pTDiffv24ChFull_%d",j),Form("pub_pTDiffv24ChFull_%d",j),20,x_pTDiff);
            for(int i=0;i<28;i++){
                pub_pTDiffv24ChFull->SetBinContent(i+1,g_pTDiffv24ChFull->GetPointY(i));
                pub_pTDiffv24ChFull->SetBinError(i+1,g_pTDiffv24ChFull->GetErrorY(i));
            }

            cpTDiffv24ChFull->cd(2);
            gPad->SetTopMargin(0.05);
            TH1D* frame_ratio_pTDiffv24ChFull = new TH1D(Form("frame_ratio_pTDiffv24ChFull_%d",j), Form("frame_ratio_pTDiffv24ChFull_%d",j),200,-0.8,200);
            frame_ratio_pTDiffv24ChFull->SetMaximum(1.4);
            frame_ratio_pTDiffv24ChFull->SetMinimum(0.6);
            frame_ratio_pTDiffv24ChFull->SetYTitle("Ratio to Run 2");
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_ratio_pTDiffv24ChFull->Draw("AXIS");
            One->Draw("sames");
            index=0;
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* ratio_pTDiffv24ChFull = GetRatio(28,(TH1D*)resultsFiles_pTDiffv24ChFull[j][i]->Get(Form("pTDiffv24ChFull")),pub_pTDiffv24ChFull,true);
                SetMarkerAndLine(ratio_pTDiffv24ChFull,GetColor(index),kFullCircle,kSolid,1.0);
                ratio_pTDiffv24ChFull->Draw("ESAMES");
                index++;
            }
            cpTDiffv24ChFull->SaveAs(Form("./OutputPDF/pTDiffv24ChFullCent%dTo%d_RatioToRun2.pdf",(int)pTDiffCent[j],(int)pTDiffCent[j+1]));
        }
    }

    // =================
    // pTDiffv24Ch10Gap
    // =================
    if (IfCheckObservable[kpTDiffv24Ch10Gap]) {
        for (uint j=0;j<pTDiffCent.size()-1;j++) {
            index = 0;
            TCanvas* cpTDiffv24Ch10Gap = new TCanvas(Form("cpTDiffv24Ch10Gap_%d",j), Form("cpTDiffv24Ch10Gap_%d",j), 800, 1200);
            cpTDiffv24Ch10Gap->Divide(1,2);
            cpTDiffv24Ch10Gap->cd(1);
            gPad->SetBottomMargin(0.05);
            TLegend* leg_pTDiffv24Ch10Gap = new TLegend(0.2,0.7,0.9,0.9);
            TH1D* frame_pTDiffv24Ch10Gap = new TH1D(Form("frame_pTDiffv24Ch10Gap_%d",j), Form("frame_pTDiffv24Ch10Gap_%d",j), 200,-0.8,200);
            frame_pTDiffv24Ch10Gap->SetMaximum(0.3);
            frame_pTDiffv24Ch10Gap->SetMinimum(0.);
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_pTDiffv24Ch10Gap->Draw("AXIS");
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* h_pTDiffv24Ch10Gap = (TH1D*)resultsFiles_pTDiffv24Ch10Gap[j][i]->Get(Form("pTDiffv24Ch10Gap"));
                SetMarkerAndLine(h_pTDiffv24Ch10Gap,GetColor(index),kFullCircle,kSolid,1.0);
                h_pTDiffv24Ch10Gap->Draw("ESAMES");
                leg_pTDiffv24Ch10Gap->AddEntry(h_pTDiffv24Ch10Gap,Form("v_{2}{4, |#Delta#eta|>1.0} (p_{T}) Cent:%d~%d%% (%s)",(int)pTDiffCent[j],(int)pTDiffCent[j+1],legendNames[i].c_str()),"lp");
                index+=1;
            }

            TGraphAsymmErrors* g_pTDiffv24Ch10Gap = nullptr;
            if (j<PubTable_pTDiffv2.size()) {
                g_pTDiffv24Ch10Gap = (TGraphAsymmErrors*)publish_ins1666817->Get(Form("Table %d/Graph1D_y1",PubTable_pTDiffv24[j]));
                if(!g_pTDiffv24Ch10Gap) continue;
                SetMarkerAndLine(g_pTDiffv24Ch10Gap,kBlack,kOpenSquare,kSolid,1.0);
                g_pTDiffv24Ch10Gap->Draw("PE");
                leg_pTDiffv24Ch10Gap->AddEntry(g_pTDiffv24Ch10Gap,Form("v2{4}(p_{T}) JHEP 07 (2018) 103"));
            }
            leg_pTDiffv24Ch10Gap->Draw();

            if(!g_pTDiffv24Ch10Gap) continue;
            TH1D* pub_pTDiffv24Ch10Gap = new TH1D(Form("pub_pTDiffv24Ch10Gap_%d",j),Form("pub_pTDiffv24Ch10Gap_%d",j),20,x_pTDiff);
            for(int i=0;i<28;i++){
                pub_pTDiffv24Ch10Gap->SetBinContent(i+1,g_pTDiffv24Ch10Gap->GetPointY(i));
                pub_pTDiffv24Ch10Gap->SetBinError(i+1,g_pTDiffv24Ch10Gap->GetErrorY(i));
            }

            cpTDiffv24Ch10Gap->cd(2);
            gPad->SetTopMargin(0.05);
            TH1D* frame_ratio_pTDiffv24Ch10Gap = new TH1D(Form("frame_ratio_pTDiffv24Ch10Gap_%d",j), Form("frame_ratio_pTDiffv24Ch10Gap_%d",j),200,-0.8,200);
            frame_ratio_pTDiffv24Ch10Gap->SetMaximum(1.4);
            frame_ratio_pTDiffv24Ch10Gap->SetMinimum(0.6);
            frame_ratio_pTDiffv24Ch10Gap->SetYTitle("Ratio to Run 2");
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_ratio_pTDiffv24Ch10Gap->Draw("AXIS");
            One->Draw("sames");
            index=0;
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* ratio_pTDiffv24Ch10Gap = GetRatio(28,(TH1D*)resultsFiles_pTDiffv24Ch10Gap[j][i]->Get(Form("pTDiffv24Ch10Gap")),pub_pTDiffv24Ch10Gap,true);
                SetMarkerAndLine(ratio_pTDiffv24Ch10Gap,GetColor(index),kFullCircle,kSolid,1.0);
                ratio_pTDiffv24Ch10Gap->Draw("ESAMES");
                index++;
            }
            cpTDiffv24Ch10Gap->SaveAs(Form("./OutputPDF/pTDiffv24Ch10GapCent%dTo%d_RatioToRun2.pdf",(int)pTDiffCent[j],(int)pTDiffCent[j+1]));
        }
    }


    // =================
    // pTDiffv26ChFull
    // =================
    if (IfCheckObservable[kpTDiffv26ChFull]) {
        for (uint j=0;j<pTDiffCent.size()-1;j++) {
            index = 0;
            TCanvas* cpTDiffv26ChFull = new TCanvas(Form("cpTDiffv26ChFull_%d",j), Form("cpTDiffv26ChFull_%d",j), 800, 1200);
            cpTDiffv26ChFull->Divide(1,2);
            cpTDiffv26ChFull->cd(1);
            gPad->SetBottomMargin(0.05);
            TLegend* leg_pTDiffv26ChFull = new TLegend(0.2,0.7,0.9,0.9);
            TH1D* frame_pTDiffv26ChFull = new TH1D(Form("frame_pTDiffv26ChFull_%d",j), Form("frame_pTDiffv26ChFull_%d",j), 200,-0.8,200);
            frame_pTDiffv26ChFull->SetMaximum(0.3);
            frame_pTDiffv26ChFull->SetMinimum(0.);
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_pTDiffv26ChFull->Draw("AXIS");
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* h_pTDiffv26ChFull = (TH1D*)resultsFiles_pTDiffv26ChFull[j][i]->Get(Form("pTDiffv26ChFull"));
                SetMarkerAndLine(h_pTDiffv26ChFull,GetColor(index),kFullCircle,kSolid,1.0);
                h_pTDiffv26ChFull->Draw("ESAMES");
                leg_pTDiffv26ChFull->AddEntry(h_pTDiffv26ChFull,Form("v_{2}{6} (p_{T}) Cent:%d~%d%% (%s)",(int)pTDiffCent[j],(int)pTDiffCent[j+1],legendNames[i].c_str()),"lp");
                index+=1;
            }

            if (IfCheckObservable[kpTDiffv24ChFull]) {
                TH1D* h_pTDiffv24ChFull = (TH1D*)resultsFiles_pTDiffv24ChFull[j][0]->Get(Form("pTDiffv24ChFull"));
                if (!h_pTDiffv24ChFull) {
                    Printf("No pTDiffv24ChFull found");
                    return;
                }
                SetMarkerAndLine(h_pTDiffv24ChFull,kBlack,kFullCircle,kSolid,1.0);
                h_pTDiffv24ChFull->Draw("sames");
                leg_pTDiffv26ChFull->AddEntry(h_pTDiffv24ChFull,Form("v_{2}{4} (p_{T}) Cent:%d~%d%% (%s)",(int)pTDiffCent[j],(int)pTDiffCent[j+1],legendNames[0].c_str()),"lp");
                leg_pTDiffv26ChFull->Draw();

                cpTDiffv26ChFull->cd(2);
                gPad->SetTopMargin(0.05);
                TH1D* frame_ratio_pTDiffv26ChFull = new TH1D(Form("frame_ratio_pTDiffv26ChFull_%d",j), Form("frame_ratio_pTDiffv26ChFull_%d",j),200,-0.8,200);
                frame_ratio_pTDiffv26ChFull->SetMaximum(1.4);
                frame_ratio_pTDiffv26ChFull->SetMinimum(0.6);
                frame_ratio_pTDiffv26ChFull->SetYTitle("Ratio to v_{2}{4}");
                if (kpTDiffLogx)
                    gPad->SetLogx();
                frame_ratio_pTDiffv26ChFull->Draw("AXIS");
                One->Draw("sames");
                index=0;
                for(int i=0;i<FileNameSuffixs.size();i++){
                    TH1D* ratio_pTDiffv26ChFull = GetRatio(28,(TH1D*)resultsFiles_pTDiffv26ChFull[j][i]->Get(Form("pTDiffv26ChFull")),(TH1D*)resultsFiles_pTDiffv24ChFull[j][0]->Get(Form("pTDiffv24ChFull")),true);
                    if (!ratio_pTDiffv26ChFull) {
                        Printf("No ratio_pTDiffv26ChFull found");
                        return;
                    }
                    SetMarkerAndLine(ratio_pTDiffv26ChFull,GetColor(index),kFullCircle,kSolid,1.0);
                    ratio_pTDiffv26ChFull->Draw("ESAMES");
                    index++;
                }

                continue;
            }

            TGraphAsymmErrors* g_pTDiffv26ChFull = nullptr;
            if (j<PubTable_pTDiffv2.size()) {
                if (!PubTable_pTDiffv24[j]) {
                    Printf("No corresponding table for pTDiffv24ChFull, skip this plot");
                    leg_pTDiffv26ChFull->Draw();
                    continue;
                }
                g_pTDiffv26ChFull = (TGraphAsymmErrors*)publish_ins1666817->Get(Form("Table %d/Graph1D_y1",PubTable_pTDiffv24[j]));
                if(!g_pTDiffv26ChFull) continue;
                SetMarkerAndLine(g_pTDiffv26ChFull,kBlack,kOpenSquare,kSolid,1.0);
                g_pTDiffv26ChFull->Draw("PE");
                leg_pTDiffv26ChFull->AddEntry(g_pTDiffv26ChFull,Form("v2{4}(p_{T}) JHEP 07 (2018) 103"));
            }
            leg_pTDiffv26ChFull->Draw();

            if(!g_pTDiffv26ChFull) continue;
            TH1D* pub_pTDiffv26ChFull = new TH1D(Form("pub_pTDiffv26ChFull_%d",j),Form("pub_pTDiffv26ChFull_%d",j),20,x_pTDiff);
            for(int i=0;i<28;i++){
                pub_pTDiffv26ChFull->SetBinContent(i+1,g_pTDiffv26ChFull->GetPointY(i));
                pub_pTDiffv26ChFull->SetBinError(i+1,g_pTDiffv26ChFull->GetErrorY(i));
            }

            cpTDiffv26ChFull->cd(2);
            gPad->SetTopMargin(0.05);
            TH1D* frame_ratio_pTDiffv26ChFull = new TH1D(Form("frame_ratio_pTDiffv26ChFull_%d",j), Form("frame_ratio_pTDiffv26ChFull_%d",j),200,-0.8,200);
            frame_ratio_pTDiffv26ChFull->SetMaximum(1.4);
            frame_ratio_pTDiffv26ChFull->SetMinimum(0.6);
            frame_ratio_pTDiffv26ChFull->SetYTitle("Ratio to Run 2");
            if (kpTDiffLogx)
                gPad->SetLogx();
            frame_ratio_pTDiffv26ChFull->Draw("AXIS");
            One->Draw("sames");
            index=0;
            for(int i=0;i<FileNameSuffixs.size();i++){
                TH1D* ratio_pTDiffv26ChFull = GetRatio(28,(TH1D*)resultsFiles_pTDiffv26ChFull[j][i]->Get(Form("pTDiffv26ChFull")),pub_pTDiffv26ChFull,true);
                SetMarkerAndLine(ratio_pTDiffv26ChFull,GetColor(index),kFullCircle,kSolid,1.0);
                ratio_pTDiffv26ChFull->Draw("ESAMES");
                index++;
            }
            cpTDiffv26ChFull->SaveAs(Form("./OutputPDF/pTDiffv26ChFullCent%dTo%d_RatioToRun2.pdf",(int)pTDiffCent[j],(int)pTDiffCent[j+1]));
        }
    }


    return;
    

}
