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
#include <cstring>
#include <vector>
#include <map>
#include <array>
#include "include/ErrorPropagation.h"

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
    for(int i=0;i<3;i++)hVn[i]->Draw("ESames");
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

std::map<ObservableEnum, std::array<double, 2>> UserRangeMap = {
  {v422, {-0.005,0.015}},
  {chi422, {-2.,2.}},
  {rho422, {0.,1.}}
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
    Hist->SetXTitle("Centrality/%");
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
    hv422->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.85,0.5,0.9);
    legend2->AddEntry(hv422,Form("%s |#Delta#eta|>1",ObservableName[observable]));
    legend2->Draw();

}

void ProcessFlowContainer(string FileNameSuffix = "LHC23zzh_pass2"){
    TFile* f = new TFile(Form("./AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    FlowContainer* fc = (FlowContainer*)f->Get("flow-pb-pb-task/FlowContainer");
    if(!fc){
        Printf("can not get flow container");
        return;
    }
    TObjArray* Subsamples = fc->GetSubProfiles();

    Output_vn(FileNameSuffix, fc);
    // Output_ptDiffvn(FileNameSuffix, fc);
    Output_Nonlinear(FileNameSuffix, fc, v422);
    Output_Nonlinear(FileNameSuffix, fc, chi422);
    Output_Nonlinear(FileNameSuffix, fc, rho422);
    return;
    
}
