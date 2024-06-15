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

bool ComparewithPublish = true;
bool OutputPNG = false;
bool RebinpTDiff = true;
std::vector<Double_t> pTDiffOriginBinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
std::vector<Double_t> pTDiffTargetBinning = {0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.5, 1.7, 2.0, 2.4, 3.0, 3.5, 4.0, 5.0};

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

void CalculateBootstrapError(std::vector<std::vector<double>>& ValueArray, std::vector<std::vector<double>>& ValueErrorArray, std::vector<double>& ErrorArray){
    int Nsample = ValueArray.size();
    int Nbin = ValueArray[0].size();
    std::vector<int> Count;
    std::vector<double> Mean;
    std::vector<double> SumWeight;
    Count.resize(Nbin);
    Mean.resize(Nbin);
    SumWeight.resize(Nbin);
    for(int i=0;i<Nbin;i++){
        SumWeight[i]=0;
        Mean[i]=0;
        for(int j=0;j<Nsample;j++){
            Count[i]++;
            Mean[i] += ValueArray[j][i] * (1./(ValueErrorArray[j][i]*ValueErrorArray[j][i]));
            // Printf("ValueArray[%d][%d]=%f, ValueErrorArray[%d][%d]=%f",j,i,ValueArray[j][i],j,i,ValueErrorArray[j][i]);
            SumWeight[i] += (1./(ValueErrorArray[j][i]*ValueErrorArray[j][i]));
        }
        Mean[i] = Mean[i]/SumWeight[i];
    }
    for (Int_t i = 0; i < Nsample; i++)
    {
        for (Int_t j = 0; j < Nbin; j++)
        {
            ErrorArray[j] += (ValueArray[i][j] - Mean[j]) * (ValueArray[i][j] - Mean[j]);
        }
    }
    for (Int_t i = 0; i < Nbin; i++)
    {
        if (Count[i] > 1)
            ErrorArray[i] = TMath::Sqrt(ErrorArray[i] / (Count[i] - 1));
        // Printf("ErrorArray[%d]=%f",i,ErrorArray[i]);
    }
    
}

void ResizeValueArray(std::vector<std::vector<std::vector<double>>>& ValueArray,
    std::vector<std::vector<std::vector<double>>>& ValueErrorArray,
    std::vector<std::vector<double>>& ErrorArray,
    int Nobs=1, int NofSample=10, int Nbin=9){

    ValueArray.resize(Nobs);
    ValueErrorArray.resize(Nobs);
    ErrorArray.resize(Nobs);    
    for(int i=0;i<Nobs;i++){
        ErrorArray[i].resize(Nbin);
        ValueArray[i].resize(NofSample);
        ValueErrorArray[i].resize(NofSample);
        for(int j=0;j<NofSample;j++)ValueArray[i][j].resize(Nbin);
        for(int j=0;j<NofSample;j++)ValueErrorArray[i][j].resize(Nbin);
    }

}

void Output_vn(string FileNameSuffix, FlowContainer* fc){
    TCanvas* canvas1 = new TCanvas("Canvas_vn","Canvas_vn",900,900);
    TH1D* Hist  = new TH1D(Form("v_{n}{2} in %s",FileNameSuffix.c_str()),Form("v_{n}{2} in %s",FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(0.);
    Hist->SetMaximum(0.15);
    Hist->SetXTitle("Centrality (%)");
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
    }

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=3;//v22,v32,v42
    TObjArray* subsamples = fc->GetSubProfiles();
    int NofSample = subsamples->GetEntries();
    int Nbin = hVn[0]->GetNbinsX();
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    
    for(int sample=0;sample<NofSample;sample++){
        fc->OverrideMainWithSub(sample,false);
        for(int i=0;i<Nobs;i++){
            TH1D* temp = (TH1D*)fc->GetVN2VsMulti(i+2);
            if(!temp){
                Printf("Can't get v%d",i+2);
                return;
            }
            for(int j=0;j<temp->GetNbinsX();j++){
                ValueArray[i][sample][j] = temp->GetBinContent(j+1);
                ValueErrorArray[i][sample][j] = temp->GetBinError(j+1);
            }
        }
    }
    for(int i=0;i<Nobs;i++){
        CalculateBootstrapError(ValueArray[i],ValueErrorArray[i],ErrorArray[i]);
    }
    for(int i=0;i<Nbin;i++){
        hVn[0]->SetBinError(i+1, ErrorArray[0][i]);
        hVn[1]->SetBinError(i+1, ErrorArray[1][i]);
        hVn[2]->SetBinError(i+1, ErrorArray[2][i]);
    }

    SetMarkerAndLine(hVn[0],kBlack,kFullCircle,kSolid,1.0);
    SetMarkerAndLine(hVn[1],kRed,kFullCircle,kSolid,1.0);
    SetMarkerAndLine(hVn[2],kBlue,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    for(int i=0;i<Nobs;i++)hVn[i]->Draw("ESames");
    TLegend* legend = new TLegend(0.2,0.7,0.5,0.9);
    for(int i=0;i<Nobs;i++)legend->AddEntry(hVn[i],Form("v_{%d}{2} |#Delta#eta|>1",i+2));
    legend->Draw();

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
        g_v4->Draw("PE");
        legend->AddEntry(g_v2,Form("v_{2}{2} JHEP 05 (2020) 085, 2020"));
        legend->AddEntry(g_v3,Form("v_{3}{2} JHEP 05 (2020) 085, 2020"));
        legend->AddEntry(g_v4,Form("v_{4}{2} JHEP 05 (2020) 085, 2020"));
    }

    if(OutputPNG){
        canvas1->SaveAs("./ProcessOutput/vn.png");
    }
}

void Output_vn4(string FileNameSuffix, FlowContainer* fc){
    TCanvas* canvas1 = new TCanvas("Canvas_vn4","Canvas_vn4",900,900);
    TH1D* Hist  = new TH1D(Form("v_{n}{4} in %s",FileNameSuffix.c_str()),Form("v_{n}{4} in %s",FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(0.);
    Hist->SetMaximum(0.15);
    Hist->SetXTitle("Centrality (%)");
    Hist->SetYTitle("v_{n}{4}");
    Hist->Draw();
    fc->SetIDName("ChFull");
    fc->SetPropagateErrors(kTRUE);
    TH1D* hVn[3] = {nullptr};
    hVn[0] = (TH1D*)fc->GetVN4VsMulti(2);
    // TH1D* temp_22 = (TH1D*)fc->GetHistCorrXXVsMulti("22");
    // TH1D* temp_24 = (TH1D*)fc->GetHistCorrXXVsMulti("24");
    // hVn[0] = (TH1D*)temp_22->Clone();
    // hVn[0]->SetName("c24");
    // for(int i=1;i<=temp_22->GetNbinsX();i++){
    //     double corr22 = temp_22->GetBinContent(i);
    //     double corr24 = temp_24->GetBinContent(i);
    //     double c24 = 2*corr22*corr22 - corr24;
    //     Printf("corr24: %f, corr22: %f",corr24, corr22);
    //     if(c24>0){
    //         hVn[0]->SetBinContent(i,sqrt(sqrt(c24)));
    //         Printf("v24: %f",sqrt(sqrt(c24)));
    //     }
    //     else{
    //         hVn[0]->SetBinContent(i,0);
    //     }
    // }
    // delete temp_22;
    // delete temp_24;

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=1;//v24
    TObjArray* subsamples = fc->GetSubProfiles();
    int NofSample = subsamples->GetEntries();
    int Nbin = hVn[0]->GetNbinsX();
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    
    for(int sample=0;sample<NofSample;sample++){
        fc->OverrideMainWithSub(sample,false);
        for(int i=0;i<Nobs;i++){
            TH1D* temp = (TH1D*)fc->GetVN4VsMulti(i+2);
            if(!temp){
                Printf("Can't get v%d",i+2);
                return;
            }
            for(int j=0;j<temp->GetNbinsX();j++){
                ValueArray[i][sample][j] = temp->GetBinContent(j+1);
                ValueErrorArray[i][sample][j] = temp->GetBinError(j+1);
            }
        }
    }
    for(int i=0;i<Nobs;i++){
        CalculateBootstrapError(ValueArray[i],ValueErrorArray[i],ErrorArray[i]);
    }
    for(int i=0;i<Nbin;i++){
        hVn[0]->SetBinError(i+1, ErrorArray[0][i]);
    }

    SetMarkerAndLine(hVn[0],kBlack,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    for(int i=0;i<Nobs;i++)hVn[i]->Draw("ESames");
    TLegend* legend = new TLegend(0.2,0.8,0.5,0.9);
    for(int i=0;i<Nobs;i++)legend->AddEntry(hVn[i],Form("v_{%d}{4}",i+2));
    legend->Draw();

    if(ComparewithPublish){
        TFile* publish = new TFile("./HEPData-ins1666817-v1-root.root","READ");
        TGraphAsymmErrors* g_v2 = (TGraphAsymmErrors*)publish->Get("Table 2/Graph1D_y1");
        SetMarkerAndLine(g_v2,kRed,kOpenSquare,kSolid,1.0);
        g_v2->Draw("PE");
        legend->AddEntry(g_v2,Form("v_{2}{4} JHEP 07 (2018) 103"));
    }

    if(OutputPNG){
        canvas1->SaveAs("./ProcessOutput/v24.png");
    }
}

void Output_ptDiffvn(string FileNameSuffix, FlowContainer* fc, Double_t CentMin=0., Double_t CentMax=5.){
    TCanvas* canvas2 = nullptr;
    bool anotherCanvas = false;
    // if canvas_ptDiffvn exist, new a canvas with different name
    if(gROOT->FindObject("canvas_ptDiffvn"))anotherCanvas = true;
    if(!anotherCanvas)canvas2 = new TCanvas("canvas_ptDiffvn","canvas_ptDiffvn",900,900);
    else canvas2 = new TCanvas("canvas_ptDiffvn_2","canvas_ptDiffvn_2",900,900);
    TH1D* Hist2  = new TH1D(Form("v_{2}(p_{T}) in %s",FileNameSuffix.c_str()),Form("v_{n}(p_{T}) in %s",FileNameSuffix.c_str()),20,0.2,10.0);
    Hist2->SetMinimum(0.);
    Hist2->SetMaximum(0.3);
    Hist2->GetXaxis()->SetRangeUser(0.2,5.0);
    Hist2->SetXTitle("p_{T}");
    Hist2->SetYTitle("v_{n}");
    Hist2->Draw();
    fc->SetIDName("ChGap");
    TH1D* hV22pt = (TH1D*)fc->GetVN2VsPt(2,CentMin,CentMax);
    if(!hV22pt){
        Printf("Can't get hV22");
        return;
    }
    if(RebinpTDiff){
        TH1D* pTMerge = new TH1D("pTMerge","pTMerge",pTDiffTargetBinning.size(),pTDiffTargetBinning.data());
        for(int i=0;i<pTDiffTargetBinning.size()-1;i++){
            auto it = std::find(pTDiffOriginBinning.begin(), pTDiffOriginBinning.end(), pTDiffTargetBinning[i]);
            if (it == pTDiffOriginBinning.end()){
                Printf("Can't find bin %f in origin binning",pTDiffTargetBinning[i]);
                break;
            }
            Int_t StartBin = it - pTDiffOriginBinning.begin();
            auto it2 = std::find(pTDiffOriginBinning.begin(), pTDiffOriginBinning.end(), pTDiffTargetBinning[i+1]);
            Int_t EndBin = it2 - pTDiffOriginBinning.begin();
            Double_t value = 0;
            Double_t weight =0;
            for(Int_t bin=StartBin+1;bin<=EndBin;bin++){
                value += hV22pt->GetBinContent(bin)/(hV22pt->GetBinError(bin)*hV22pt->GetBinError(bin));
                weight += 1./(hV22pt->GetBinError(bin)*hV22pt->GetBinError(bin));
            }
            value /= weight;
            weight = 1./weight;
            pTMerge->SetBinContent(i+1,value);
            Printf("pTMerge bin %d: %f",i+1,value);
            pTMerge->SetBinError(i+1,sqrt(weight));
        }
        hV22pt = pTMerge;
    }
    SetMarkerAndLine(hV22pt,kBlack,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    hV22pt->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.8,0.5,0.9);
    legend2->AddEntry(hV22pt,Form("v_{2}{2}(p_{T}) |#Delta#eta|>1 cent:%d~%d%%",(int)CentMin,(int)CentMax));
    legend2->Draw();

    if(ComparewithPublish){
        Int_t table = 0;
        if(CentMin==0. && CentMax==5.)table = 7;
        else if(CentMin==30. && CentMax==40.)table = 8;
        if(table>0){
            TFile* publish = new TFile("./HEPData-ins1419244-v2-root.root","READ");
            TGraphAsymmErrors* g_v2 = (TGraphAsymmErrors*)publish->Get(Form("Table %d/Graph1D_y1",table));
            SetMarkerAndLine(g_v2,kRed,kOpenSquare,kSolid,1.0);
            g_v2->Draw("PE");
            legend2->AddEntry(g_v2,Form("v_{2}{2}(p_{T}) Phys.Rev.Lett. 116 (2016) 132302, 2016"));
        }
    }

    if(OutputPNG){
        if(!anotherCanvas)
        canvas2->SaveAs("./ProcessOutput/ptDiffvn.png");
        else
        canvas2->SaveAs("./ProcessOutput/ptDiffvn_2.png");
    }
}

void GetNonlinearHistogram(FlowContainer* fc, TH1D*& hCorr422, TH1D*& hCorr24, TH1D*& hCorr42, string GapName="Ch10Gap"){
    fc->SetIDName(Form("%sA",GapName.c_str()));
    TH1D* temp_422A = (TH1D*)fc->GetHistCorrXXVsMulti("422");
    if(!temp_422A){Printf("Can't get temp_422A");return;}
    TH1D* hCorr422A = (TH1D*)temp_422A->Clone();
    delete temp_422A;
    hCorr422A->SetName("hCorr422A");

    fc->SetIDName(Form("%sB",GapName.c_str()));
    TH1D* temp_422B = (TH1D*)fc->GetHistCorrXXVsMulti("422");
    if(!temp_422B){Printf("Can't get temp_422B");return;}
    TH1D* hCorr422B = (TH1D*)temp_422B->Clone();
    delete temp_422B;
    hCorr422B->SetName("hCorr422B");

    if(!hCorr422){
        hCorr422 = (TH1D*)hCorr422A->Clone();
        hCorr422->SetName("hCorr422_mean");
    }
    for(int i=1;i<=hCorr422->GetNbinsX();i++){
        // Printf("%f",hCorr422A->GetBinCenter(i));
        hCorr422->SetBinContent(i,(hCorr422A->GetBinContent(i)+hCorr422B->GetBinContent(i))/2.);
        hCorr422->SetBinError(i,(hCorr422A->GetBinError(i)+hCorr422B->GetBinError(i))/2.);
    }

    fc->SetIDName(Form("%s",GapName.c_str()));
    TH1D* temp_24= (TH1D*)fc->GetHistCorrXXVsMulti("24");
    if(!temp_24){Printf("Can't get temp_24");return;}
    hCorr24 = (TH1D*)temp_24->Clone();
    delete temp_24;

    fc->SetIDName(Form("%s",GapName.c_str()));
    TH1D* temp_42= (TH1D*)fc->GetHistCorrXXVsMulti("42");
    if(!temp_42){Printf("Can't get temp_42");return;}
    hCorr42 = (TH1D*)temp_42->Clone();
    delete temp_42;
}

enum ObservableEnum{
    v422,
    chi422,
    rho422
};

std::map<ObservableEnum, Char_t*> ObservableName = {
  {v422, "v_{4,22}"},
  {chi422, "#chi_{4,22}"},
  {rho422, "#rho_{4,22}"}
};

std::map<ObservableEnum, Char_t*> OutputObservableName = {
  {v422, "v422"},
  {chi422, "chi422"},
  {rho422, "rho422"}
};

std::map<ObservableEnum, std::array<double, 2>> UserRangeMap = {
  {v422, {-0.006,0.015}},
  {chi422, {-5.5,2.}},
  {rho422, {-0.5,1.}}
};

void SetNonlinearValue(TH1D*& target, TH1D*& hCorr422, TH1D*& hCorr24, TH1D*& hCorr42, ObservableEnum observable=v422){
    if(!target)target = (TH1D*)hCorr422->Clone();
    for(int i=1;i<=hCorr422->GetNbinsX();i++){
        if(observable==chi422){
            target->SetBinContent(i,hCorr422->GetBinContent(i)/hCorr24->GetBinContent(i));
            target->SetBinError(i,Error_Chi(hCorr422->GetBinContent(i), hCorr422->GetBinError(i), hCorr24->GetBinContent(i), hCorr24->GetBinError(i)));
        }
        else{
            if(observable==v422){
                if(hCorr24->GetBinContent(i)>0.){
                    target->SetBinContent(i,hCorr422->GetBinContent(i)/sqrt(hCorr24->GetBinContent(i)));
                    target->SetBinError(i,Error_vNL(hCorr422->GetBinContent(i), hCorr422->GetBinError(i), hCorr24->GetBinContent(i), hCorr24->GetBinError(i)));
                }
                else{
                    target->SetBinContent(i,0.);
                    target->SetBinError(i,0.);
                    Printf("Warning: in %f %%, Ch10Gap24 is negative", hCorr422->GetBinCenter(i));
                }
            }
            else if(observable==rho422){
                if(hCorr24->GetBinContent(i)>0.&&hCorr42->GetBinContent(i)>0.){
                    target->SetBinContent(i,hCorr422->GetBinContent(i)/(sqrt(hCorr24->GetBinContent(i))*sqrt(hCorr42->GetBinContent(i))));
                    target->SetBinError(i,Error_Rho(hCorr422->GetBinContent(i), hCorr422->GetBinError(i), hCorr24->GetBinContent(i), hCorr24->GetBinError(i),hCorr42->GetBinContent(i), hCorr42->GetBinError(i)));
                }
                else{
                    target->SetBinContent(i,0.);
                    target->SetBinError(i,0.);
                    Printf("Warning: in %f %%, Ch10Gap24 or Ch10Gap42 is negative", hCorr422->GetBinCenter(i));
                }
            }
        }
        
    }
    delete hCorr422;
    delete hCorr24;
    delete hCorr42;
}

void Output_Nonlinear(string FileNameSuffix, FlowContainer* fc, ObservableEnum observable=v422){
    
    TCanvas* canvas1 = new TCanvas(Form("Canvas_%s",ObservableName[observable]),Form("Canvas_%s",ObservableName[observable]),900,900);
    TH1D* Hist  = new TH1D(Form("%s in %s",ObservableName[observable],FileNameSuffix.c_str()),Form("%s in %s",ObservableName[observable],FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(UserRangeMap[observable][0]);
    Hist->SetMaximum(UserRangeMap[observable][1]);
    Hist->SetXTitle("Centrality (%)");
    Hist->SetYTitle(ObservableName[observable]);
    Hist->Draw();

    // fc->SetPropagateErrors(kTRUE);
    TH1D* hCorr422=nullptr;
    TH1D* hCorr24=nullptr;
    TH1D* hCorr42=nullptr;
    GetNonlinearHistogram(fc,hCorr422,hCorr24,hCorr42,"Ch10Gap");

    TH1D* hv422 = (TH1D*)hCorr422->Clone();
    SetNonlinearValue(hv422,hCorr422,hCorr24,hCorr42,observable);

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=1;//v422
    TObjArray* subsamples = fc->GetSubProfiles();
    int NofSample = subsamples->GetEntries();
    int Nbin = hv422->GetNbinsX();
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    
    for(int sample=0;sample<NofSample;sample++){
        fc->OverrideMainWithSub(sample,false);
        for(int i=0;i<Nobs;i++){
            TH1D* hCorr422=nullptr;
            TH1D* hCorr24=nullptr;
            TH1D* hCorr42=nullptr;
            GetNonlinearHistogram(fc,hCorr422,hCorr24,hCorr42,"Ch10Gap");
            TH1D* temp = (TH1D*)hCorr422->Clone();
            SetNonlinearValue(temp,hCorr422,hCorr24,hCorr42,observable);
            if(!temp){
                Printf("Can't get v422");
                return;
            }
            for(int j=0;j<temp->GetNbinsX();j++){
                ValueArray[i][sample][j] = temp->GetBinContent(j+1);
                ValueErrorArray[i][sample][j] = temp->GetBinError(j+1);
            }
        }
    }
    for(int i=0;i<Nobs;i++){
        CalculateBootstrapError(ValueArray[i],ValueErrorArray[i],ErrorArray[i]);
    }
    for(int i=0;i<Nbin;i++){
        hv422->SetBinError(i+1, ErrorArray[0][i]);
    }

    SetMarkerAndLine(hv422,kBlack,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    hv422->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.85,0.5,0.9);
    legend2->AddEntry(hv422,Form("%s |#Delta#eta|>1",ObservableName[observable]));
    legend2->Draw();

    if(ComparewithPublish){
        TFile* publish = new TFile("./HEPData-ins1778342-v1-root.root","READ");
        TGraphAsymmErrors* g = nullptr;
        if(observable==v422)g=(TGraphAsymmErrors*)publish->Get("v422/Graph1D_y1");
        else if(observable==chi422)g=(TGraphAsymmErrors*)publish->Get("chi422/Graph1D_y1");
        else if(observable==rho422)g=(TGraphAsymmErrors*)publish->Get("rho422/Graph1D_y1");
        if(!g)return;
        SetMarkerAndLine(g,kRed,kOpenSquare,kSolid,1.0);
        g->Draw("PE");
        legend2->AddEntry(g,Form("%s JHEP 05 (2020) 085, 2020",ObservableName[observable]));
    }

    if(OutputPNG){
        canvas1->SaveAs(Form("./ProcessOutput/%s.png",OutputObservableName[observable]));
    }

}

void CalculateNSC32(FlowContainer* fc, TH1D*& target){

    TH1D* hCorr3232=nullptr;
    TH1D* hCorr22_Full=nullptr;
    TH1D* hCorr32_Full=nullptr;
    TH1D* hCorr22_Gap=nullptr;
    TH1D* hCorr32_Gap=nullptr;
    fc->SetIDName(Form("ChFull"));
    TH1D* temp_3232= (TH1D*)fc->GetHistCorrXXVsMulti("3232");
    if(!temp_3232){Printf("Can't get temp_3232");return;}
    hCorr3232 = (TH1D*)temp_3232->Clone();
    delete temp_3232;

    TH1D* temp_22_Full= (TH1D*)fc->GetHistCorrXXVsMulti("22");
    if(!temp_22_Full){Printf("Can't get temp_22_Full");return;}
    hCorr22_Full = (TH1D*)temp_22_Full->Clone();
    hCorr22_Full->SetName("hCorr22_Full");
    delete temp_22_Full;

    TH1D* temp_32_Full= (TH1D*)fc->GetHistCorrXXVsMulti("32");
    if(!temp_32_Full){Printf("Can't get temp_32_Full");return;}
    hCorr32_Full = (TH1D*)temp_32_Full->Clone();
    hCorr32_Full->SetName("hCorr32_Full");
    delete temp_32_Full;

    fc->SetIDName(Form("Ch10Gap"));
    TH1D* temp_22_Gap= (TH1D*)fc->GetHistCorrXXVsMulti("22");
    if(!temp_22_Gap){Printf("Can't get temp_22_Gap");return;}
    hCorr22_Gap = (TH1D*)temp_22_Gap->Clone();
    delete temp_22_Gap;

    TH1D* temp_32_Gap= (TH1D*)fc->GetHistCorrXXVsMulti("32");
    if(!temp_32_Gap){Printf("Can't get temp_32_Gap");return;}
    hCorr32_Gap = (TH1D*)temp_32_Gap->Clone();
    delete temp_32_Gap;

    target = (TH1D*)hCorr3232->Clone();
    for(int i=1;i<=hCorr3232->GetNbinsX();i++){
        target->SetBinContent(i,
        (hCorr3232->GetBinContent(i)-hCorr22_Full->GetBinContent(i)*hCorr32_Full->GetBinContent(i))
        /(hCorr22_Gap->GetBinContent(i)*hCorr32_Gap->GetBinContent(i))
        );
        double err = Error_NSC(hCorr3232->GetBinContent(i),hCorr3232->GetBinError(i),
        hCorr22_Full->GetBinContent(i),hCorr22_Full->GetBinError(i),
        hCorr32_Full->GetBinContent(i),hCorr32_Full->GetBinError(i),
        hCorr22_Gap->GetBinContent(i),hCorr22_Gap->GetBinError(i),
        hCorr32_Gap->GetBinContent(i),hCorr32_Gap->GetBinError(i)
        );
        target->SetBinError(i,err);
    }

    delete hCorr3232;
    delete hCorr22_Full;
    delete hCorr32_Full;
    delete hCorr22_Gap;
    delete hCorr32_Gap;
    
}

void Output_NSC(string FileNameSuffix, FlowContainer* fc){
    
    TCanvas* canvas1 = new TCanvas(Form("Canvas_NSC"),Form("Canvas_NSC"),900,900);
    TH1D* Hist  = new TH1D(Form("NSC in %s",FileNameSuffix.c_str()),Form("NSC in %s",FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(-1.);
    Hist->SetMaximum(1.);
    Hist->SetXTitle("Centrality (%)");
    Hist->SetYTitle("NSC(3,2)");
    Hist->Draw();

    // fc->SetPropagateErrors(kTRUE);
    TH1D* NSC32 = nullptr;
    CalculateNSC32(fc, NSC32);
    NSC32->SetName("NSC32");
    

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=1;//NSC
    TObjArray* subsamples = fc->GetSubProfiles();
    int NofSample = subsamples->GetEntries();
    int Nbin = NSC32->GetNbinsX();
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    
    for(int sample=0;sample<NofSample;sample++){
        fc->OverrideMainWithSub(sample,false);
        for(int i=0;i<Nobs;i++){
            TH1D* temp = nullptr;
            CalculateNSC32(fc, temp);
            if(!temp){
                Printf("Can't get NSC32");
                return;
            }
            for(int j=0;j<temp->GetNbinsX();j++){
                ValueArray[i][sample][j] = temp->GetBinContent(j+1);
                ValueErrorArray[i][sample][j] = temp->GetBinError(j+1);
            }
            delete temp;
        }
    }
    for(int i=0;i<Nobs;i++){
        CalculateBootstrapError(ValueArray[i],ValueErrorArray[i],ErrorArray[i]);
    }
    for(int i=0;i<Nbin;i++){
        NSC32->SetBinError(i+1, ErrorArray[0][i]);
    }

    SetMarkerAndLine(NSC32,kBlack,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    NSC32->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.85,0.5,0.9);
    legend2->AddEntry(NSC32,Form("NSC(3,2)"));
    legend2->Draw();

    if(ComparewithPublish){
        TFile* publish = new TFile("./HEPData-ins1848215-v1-root.root","READ");
        TGraphAsymmErrors* g = nullptr;
        g=(TGraphAsymmErrors*)publish->Get("Table 1/Graph1D_y1");
        if(!g)return;
        SetMarkerAndLine(g,kRed,kOpenSquare,kSolid,1.0);
        g->Draw("PE");
        legend2->AddEntry(g,Form("NSC(3,2) PLB 818 (2021) 136354, 2021"));
    }

    if(OutputPNG){
        canvas1->SaveAs(Form("./ProcessOutput/NSC.png"));
    }
}

void CompareNonlinearCorr(string FileNameSuffix, FlowContainer* fc){
    TCanvas* canvas1 = new TCanvas(Form("Canvas_NSC"),Form("Canvas_NSC"),900,900);
    TH1D* Hist  = new TH1D(Form("<<4,-2,-2>> in %s",FileNameSuffix.c_str()),Form("<<4,-2,-2>> in %s",FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(-5e-5);
    Hist->SetMaximum(2e-4);
    Hist->SetXTitle("Centrality (%)");
    Hist->SetYTitle("Correlation");
    Hist->Draw();

    TH1D* hCorr422 = nullptr;
    fc->SetIDName(Form("Ch10GapA"));
    TH1D* temp_422A = (TH1D*)fc->GetHistCorrXXVsMulti("422");
    if(!temp_422A){Printf("Can't get temp_422A");return;}
    TH1D* hCorr422A = (TH1D*)temp_422A->Clone();
    delete temp_422A;
    hCorr422A->SetName("hCorr422A");

    fc->SetIDName(Form("Ch10GapB"));
    TH1D* temp_422B = (TH1D*)fc->GetHistCorrXXVsMulti("422");
    if(!temp_422B){Printf("Can't get temp_422B");return;}
    TH1D* hCorr422B = (TH1D*)temp_422B->Clone();
    delete temp_422B;
    hCorr422B->SetName("hCorr422B");

    if(!hCorr422){
        hCorr422 = (TH1D*)hCorr422A->Clone();
        hCorr422->SetName("hCorr422_mean");
    }
    for(int i=1;i<=hCorr422->GetNbinsX();i++){
        // Printf("%f",hCorr422A->GetBinCenter(i));
        hCorr422->SetBinContent(i,(hCorr422A->GetBinContent(i)+hCorr422B->GetBinContent(i))/2.);
        hCorr422->SetBinError(i,(hCorr422A->GetBinError(i)+hCorr422B->GetBinError(i))/2.);
    }

    SetMarkerAndLine(hCorr422,kBlack,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    hCorr422->Draw("ESames");

}

void CompareNSC32Corr(string FileNameSuffix, FlowContainer* fc){
    TCanvas* canvas1 = new TCanvas(Form("Canvas_NSC"),Form("Canvas_NSC"),900,900);
    TH1D* Hist  = new TH1D(Form("NSC in %s",FileNameSuffix.c_str()),Form("NSC in %s",FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(0.);
    Hist->SetMaximum(2e-5);
    Hist->SetXTitle("Centrality (%)");
    Hist->SetYTitle("Correlation");
    Hist->Draw();

    fc->SetIDName(Form("ChFull"));
    TH1D* temp_3232= (TH1D*)fc->GetHistCorrXXVsMulti("3232");
    if(!temp_3232){Printf("Can't get temp_3232");return;}
    TH1D* hCorr3232 = (TH1D*)temp_3232->Clone();
    delete temp_3232;

    TH1D* temp_22_Full= (TH1D*)fc->GetHistCorrXXVsMulti("22");
    if(!temp_22_Full){Printf("Can't get temp_22_Full");return;}
    TH1D* hCorr22_Full = (TH1D*)temp_22_Full->Clone();
    hCorr22_Full->SetName("hCorr22_Full");
    delete temp_22_Full;

    TH1D* temp_32_Full= (TH1D*)fc->GetHistCorrXXVsMulti("32");
    if(!temp_32_Full){Printf("Can't get temp_32_Full");return;}
    TH1D* hCorr32_Full = (TH1D*)temp_32_Full->Clone();
    hCorr32_Full->SetName("hCorr32_Full");
    delete temp_32_Full;

    TH1D* hCorr32and22 = (TH1D*)hCorr22_Full->Clone();
    for(int i=1;i<=hCorr32and22->GetNbinsX();i++){
        hCorr32and22->SetBinContent(i, hCorr22_Full->GetBinContent(i)*hCorr32_Full->GetBinContent(i));
        hCorr32and22->SetBinError(i, 0., sqrt(
            pow(hCorr22_Full->GetBinContent(i)*hCorr32_Full->GetBinError(i),2)
            +pow(hCorr32_Full->GetBinContent(i)*hCorr22_Full->GetBinError(i),2)
        ));
    }

    SetMarkerAndLine(hCorr3232,kBlack,kFullCircle,kSolid,1.0);
    SetMarkerAndLine(hCorr32and22,kRed,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    hCorr3232->Draw("ESames");
    hCorr32and22->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.85,0.5,0.9);
    legend2->AddEntry(hCorr3232,Form("<<3,2,-3,-2>>"));
    legend2->AddEntry(hCorr32and22,Form("<<3,-3>><<2,-2>>"));
    legend2->Draw();

}

void ProcessFlowContainer(string FileNameSuffix = "LHC23zzn_pass3_I_A11_small_219651"){
    TFile* f = new TFile(Form("./AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    FlowContainer* fc = (FlowContainer*)f->Get("flow-pb-pb-task/FlowContainer");
    if(!fc){
        Printf("can not get flow container");
        return;
    }
    TObjArray* Subsamples = fc->GetSubProfiles();

    Output_vn(FileNameSuffix, fc);
    Output_vn4(FileNameSuffix, fc);
    Output_ptDiffvn(FileNameSuffix, fc);
    Output_ptDiffvn(FileNameSuffix, fc,30,40);
    Output_Nonlinear(FileNameSuffix, fc, v422);
    Output_Nonlinear(FileNameSuffix, fc, chi422);
    Output_Nonlinear(FileNameSuffix, fc, rho422);
    Output_NSC(FileNameSuffix, fc);
    // CompareNSC32Corr(FileNameSuffix, fc);
    // CompareNonlinearCorr(FileNameSuffix, fc);
    return;
    
}
