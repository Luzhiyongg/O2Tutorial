#ifndef FlowContainerCalculation_h
#define FlowContainerCalculation_h

#include "./BasicSetting.h"
#include "./GraphSetting.h"

TH1D* mergeCentralityToTargetBin(TH1D* h, const std::vector<float>& targetCentralityBins){
    // merge from fine binning to wide binning
    // target binning is defined in BasicSetting.h
    // target binning shoule not cross the origin binning, because this function detects the edge of bins.
    TH1D* h_target = new TH1D(Form("%s",h->GetName()),Form("%s",h->GetName()),targetCentralityBins.size()-1,targetCentralityBins.data());
    for(int i=1;i<targetCentralityBins.size();i++){    
        double content=0;
        double weightSum=0;
        int count=0;
        for(int j=1;j<=h->GetNbinsX();j++){
            double lowEdge = h->GetBinLowEdge(j);
            double highEdge = h->GetBinLowEdge(j+1);
            if(lowEdge>=targetCentralityBins[i-1] && highEdge<=targetCentralityBins[i]){
                double weight = 1. / (h->GetBinError(j) * h->GetBinError(j));
                content += h->GetBinContent(j) * weight;
                weightSum += weight;
                count++;
            }
        }
        if(count > 0 && !std::isinf(content) && content==content){
            content = content/weightSum;
            double error = sqrt(1./weightSum)/sqrt(count);
            h_target->SetBinContent(i,content);
            h_target->SetBinError(i,error);
            // Printf("Bin %f, content=%f, error=%f",(targetCentralityBins[i-1]+targetCentralityBins[i])/2,content,error);
        }
        else{
            h_target->SetBinContent(i,-1);
            h_target->SetBinError(i,10.);
            Printf("Bin %f, Set to content=-1, error=10",(targetCentralityBins[i-1]+targetCentralityBins[i])/2);
        }
    }
    return h_target;
}

void CalculateBootstrapError(const std::vector<std::vector<double>>& ValueArray, const std::vector<std::vector<double>>& ValueErrorArray, std::vector<double>& ErrorArray){
    int Nsample = ValueArray.size();
    int Nbin = ValueArray[0].size();
    std::vector<int> Count;
    std::vector<double> Mean;
    std::vector<double> SumWeight;
    std::vector<std::vector<bool>> MaskArray;
    Count.resize(Nbin);
    Mean.resize(Nbin);
    SumWeight.resize(Nbin);
    MaskArray.resize(Nsample);
    for(int i=0;i<Nsample;i++){
        MaskArray[i].resize(Nbin);
        for(int j=0;j<Nbin;j++){
            if (ValueArray[i][j]==-1 && ValueErrorArray[i][j]==10) {
                MaskArray[i][j] = false;
                continue;
            }
            if (isZero(ValueArray[i][j]) && isZero(ValueErrorArray[i][j])) {
                MaskArray[i][j] = false;
                continue;
            }
            MaskArray[i][j] = true;
        }
    }
    for(int i=0;i<Nbin;i++){
        SumWeight[i]=0;
        Mean[i]=0;
        for(int j=0;j<Nsample;j++){
            if (!MaskArray[j][i]) continue;
            Count[i]++;
            Mean[i] += ValueArray[j][i] * (1./(ValueErrorArray[j][i]*ValueErrorArray[j][i]));
            // Printf("ValueArray[%d][%d]=%f, ValueErrorArray[%d][%d]=%f",j,i,ValueArray[j][i],j,i,ValueErrorArray[j][i]);
            SumWeight[i] += (1./(ValueErrorArray[j][i]*ValueErrorArray[j][i]));
        }
        if(Count[i]>0)
            Mean[i] = Mean[i]/SumWeight[i];
    }
    for (Int_t i = 0; i < Nsample; i++)
    {
        for (Int_t j = 0; j < Nbin; j++)
        {
            if (!MaskArray[i][j]) continue;
            ErrorArray[j] += (ValueArray[i][j] - Mean[j]) * (ValueArray[i][j] - Mean[j]);
        }
    }
    for (Int_t i = 0; i < Nbin; i++)
    {
        if (Count[i] > 2){
            ErrorArray[i] = TMath::Sqrt(ErrorArray[i] / (Count[i] - 1));
        }
        else {
            ErrorArray[i] = 10;
            Printf("Bin %d has only %d samples, set error to 10",i,Count[i]);
        }
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

void Output_vn(string FileNameSuffix, FlowContainer* fc, string Subwagon=""){
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
        hVn[i] = mergeCentralityToTargetBin(hVn[i],targetCentralityBins);
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
            temp = mergeCentralityToTargetBin(temp,targetCentralityBins);
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

    if(OutputRoot){
        TFile* fout = new TFile(Form("./ProcessOutput/vn_%s%s.root",FileNameSuffix.c_str(),Subwagon.c_str()),"RECREATE");
        hVn[0]->Write();
        hVn[1]->Write();
        hVn[2]->Write();
        fout->Close();
        // canvas1->SaveAs("./ProcessOutput/vn.png");
    }
}

TH1D* GetVN10(FlowContainer* fc, Int_t n) {
    TH1D *corrN2, *corrN4, *corrN6, *corrN8, *corrN10;
    corrN2 = fc->GetHistCorrXXVsMulti(Form("%i2", n));
    corrN4 = fc->GetHistCorrXXVsMulti(Form("%i4", n));
    corrN6 = fc->GetHistCorrXXVsMulti(Form("%i6", n));
    corrN8 = fc->GetHistCorrXXVsMulti(Form("%i8", n));
    corrN10 = fc->GetHistCorrXXVsMulti(Form("%i10", n));
    TH1D* rethist = dynamic_cast<TH1D*>(corrN10->Clone(Form("cN10_%s", corrN10->GetName())));
    rethist->Reset();
    for (int i = 1; i <= corrN10->GetNbinsX(); i++) {
        double cor10v = corrN10->GetBinContent(i);
        double cor10e = corrN10->GetBinError(i);
        double cor8v = corrN8->GetBinContent(i);
        double cor8e = corrN8->GetBinError(i);
        double cor6v = corrN6->GetBinContent(i);
        double cor6e = corrN6->GetBinError(i);
        double cor4v = corrN4->GetBinContent(i);
        double cor4e = corrN4->GetBinError(i);
        double cor2v = corrN2->GetBinContent(i);
        double cor2e = corrN2->GetBinError(i);
        double cn10 = cor10v - 25.*cor8v*cor2v - 100.*cor6v*cor4v + 400.*cor6v*cor2v*cor2v
        + 900.*cor4v*cor4v*cor2v - 3600.*cor4v*cor2v*cor2v*cor2v
        + 2800.*cor2v*cor2v*cor2v*cor2v*cor2v;
        double cn10e = Error_CN10(cor10v, cor10e, cor8v, cor8e, cor6v, cor6e, cor4v, cor4e, cor2v, cor2e);
        double vn10 = 0.;
        double vn10e = 0.;
        if (cn10 > 0.) {
            vn10 = TMath::Power(cn10/456., 0.1);
            vn10e = TMath::Abs(TMath::Power(1./456., 0.1)*0.1*TMath::Power(cn10, -0.9)*cn10e);
        }
        rethist->SetBinContent(i, vn10);
        rethist->SetBinError(i, vn10e);
    }

    rethist->SetBinContent(0, 0);
    rethist->SetBinError(0, 0);
    return rethist;
}

TH1D* GetVnm(FlowContainer* fc, Int_t n, Int_t m_particle = 4){
    if (m_particle == 2)
        return (TH1D*)fc->GetVN2VsMulti(n);
    else if (m_particle == 4)
        return (TH1D*)fc->GetVN4VsMulti(n);
    else if (m_particle == 6)
        return (TH1D*)fc->GetVN6VsMulti(n);
    else if (m_particle == 8)
        return (TH1D*)fc->GetVN8VsMulti(n);
    else if (m_particle == 10){
        return GetVN10(fc, n);
    }
    
    return nullptr;
}

void Output_Vnm(string FileNameSuffix, FlowContainer* fc, Int_t n = 2, Int_t m_particle = 4, string IDName = "ChFull", string Subwagon=""){
    TCanvas* canvas1 = new TCanvas(Form("Canvas_Vnm_%d_%d", n, m_particle),Form("Canvas_Vnm_%d_%d", n, m_particle),900,900);
    TH1D* Hist  = new TH1D(Form("v_{%d}{%d} in %s", n, m_particle, FileNameSuffix.c_str()), Form("v_{%d}{%d} in %s", n, m_particle, FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(0.);
    Hist->SetMaximum(0.15);
    Hist->SetXTitle("Centrality (%)");
    Hist->SetYTitle(Form("v_{%d}{%d}", n, m_particle));
    Hist->Draw();
    fc->SetIDName(IDName.c_str());
    fc->SetPropagateErrors(kTRUE);
    TH1D* hVn[3] = {nullptr};
    hVn[0] = GetVnm(fc, n, m_particle);
    hVn[0] = mergeCentralityToTargetBin(hVn[0],targetCentralityBins);
    if (hVn[0] == nullptr){
        Printf("Can't get v_{%d}{%d}",n, m_particle);
        return;
    }

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
            TH1D* temp = GetVnm(fc, n, m_particle);
            temp = mergeCentralityToTargetBin(temp,targetCentralityBins);
            if(!temp){
                Printf("Can't get v_{%d}{%d}", n, m_particle);
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
    for(int i=0;i<Nobs;i++)legend->AddEntry(hVn[i],Form("v_{%d}{%d}", n, m_particle));
    legend->Draw();

    if(ComparewithPublish){
        TFile* publish = new TFile("./HEPData-ins1666817-v1-root.root","READ");
        TGraphAsymmErrors* g_v2 = (TGraphAsymmErrors*)publish->Get("Table 2/Graph1D_y1");
        SetMarkerAndLine(g_v2,kRed,kOpenSquare,kSolid,1.0);
        g_v2->Draw("PE");
        legend->AddEntry(g_v2,Form("v_{2}{4} JHEP 07 (2018) 103"));
    }

    if(OutputRoot){
        TFile* fout = new TFile(Form("./ProcessOutput/v%d%d_%s%s.root", n, m_particle, FileNameSuffix.c_str(),Subwagon.c_str()),"RECREATE");
        hVn[0]->Write();
        fout->Close();
    }
}

void Output_ptDiffvn(string FileNameSuffix, FlowContainer* fc, Int_t n=2, Double_t CentMin=0., Double_t CentMax=5., string Subwagon=""){
    TCanvas* canvas2 = nullptr;
    bool anotherCanvas = false;
    // if canvas_ptDiffvn exist, new a canvas with different name
    if(gROOT->FindObject("canvas_ptDiffvn"))anotherCanvas = true;
    if(!anotherCanvas)canvas2 = new TCanvas("canvas_ptDiffvn","canvas_ptDiffvn",900,900);
    else canvas2 = new TCanvas("canvas_ptDiffvn_2","canvas_ptDiffvn_2",900,900);
    TH1D* Hist2  = new TH1D(Form("v_{n}(p_{T}) in %s",FileNameSuffix.c_str()),Form("v_{n}(p_{T}) in %s",FileNameSuffix.c_str()),200,0,200.0);
    Hist2->SetMinimum(0.);
    Hist2->SetMaximum(0.3);
    Hist2->GetXaxis()->SetRangeUser(0.2,200.0);
    Hist2->SetXTitle("p_{T}");
    Hist2->SetYTitle("v_{n}");
    Hist2->Draw();
    fc->SetIDName("Ch10Gap");
    fc->SetPropagateErrors(kTRUE);
    TH1D* hV22pt = (TH1D*)fc->GetVN2VsPt(n,CentMin,CentMax);
    hV22pt->SetName(Form("pTDiffv%d",n));
    if(!hV22pt){
        Printf("Can't get hV22");
        return;
    }

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=1;//v22(pT)
    TObjArray* subsamples = fc->GetSubProfiles();
    int NofSample = subsamples->GetEntries();
    int Nbin = hV22pt->GetNbinsX();
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    
    for(int sample=0;sample<NofSample;sample++){
        fc->OverrideMainWithSub(sample,false);
        for(int i=0;i<Nobs;i++){
            TH1D* temp = (TH1D*)fc->GetVN2VsPt(n,CentMin,CentMax);
            temp->SetName(Form("pTDiffv2_%d",sample));
            if(!temp){
                Printf("Can't get pTDiffv2_%d",sample);
                return;
            }
            for(int j=0;j<temp->GetNbinsX();j++){
                ValueArray[i][sample][j] = temp->GetBinContent(j+1);
                ValueErrorArray[i][sample][j] = temp->GetBinError(j+1);
                // Printf("pTDiffv2_%d_%d = %f +- %f",sample,j,ValueArray[i][sample][j],ValueErrorArray[i][sample][j]);
            }
        }
    }
    for(int i=0;i<Nobs;i++){
        CalculateBootstrapError(ValueArray[i],ValueErrorArray[i],ErrorArray[i]);
    }
    for(int i=0;i<Nbin;i++){
        hV22pt->SetBinError(i+1, ErrorArray[0][i]);
    }

    SetMarkerAndLine(hV22pt,kBlack,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    hV22pt->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.8,0.5,0.9);
    legend2->AddEntry(hV22pt,Form("v_{%d}{2}(p_{T}) |#Delta#eta|>1 cent:%d~%d%%",n,(int)CentMin,(int)CentMax));
    legend2->Draw();

    if(OutputRoot){
        TFile* fout = new TFile(Form("./ProcessOutput/pTDiffv%dCent%dTo%d_%s%s.root",n,(int)CentMin,(int)CentMax,FileNameSuffix.c_str(),Subwagon.c_str()),"RECREATE");
        hV22pt->Write();
        fout->Close();
    }
}

void Output_ptDiffvn4(string FileNameSuffix, FlowContainer* fc, Double_t CentMin=0., Double_t CentMax=5., string Subwagon=""){
    TCanvas* canvas2 = nullptr;
    bool anotherCanvas = false;
    // if canvas_ptDiffvn exist, new a canvas with different name
    if(gROOT->FindObject("canvas_ptDiffvn4"))anotherCanvas = true;
    if(!anotherCanvas)canvas2 = new TCanvas("canvas_ptDiffvn4","canvas_ptDiffvn4",900,900);
    else canvas2 = new TCanvas("canvas_ptDiffvn4_2","canvas_ptDiffvn4_2",900,900);
    TH1D* Hist2  = new TH1D(Form("v_{n}{4}(p_{T}) in %s",FileNameSuffix.c_str()),Form("v_{n}{4}(p_{T}) in %s",FileNameSuffix.c_str()),100,0,100.0);
    Hist2->SetMinimum(-1);
    Hist2->SetMaximum(2.);
    Hist2->GetXaxis()->SetRangeUser(0.,100.0);
    Hist2->SetXTitle("p_{T}");
    Hist2->SetYTitle("v_{n}{4}");
    Hist2->Draw();
    fc->SetIDName("ChFull");
    fc->SetPropagateErrors(kTRUE);
    TH1D* hV24pt = (TH1D*)fc->GetVN4VsPt(2,CentMin,CentMax);
    hV24pt->SetName("pTDiffv24");
    if(!hV24pt){
        Printf("Can't get hV24");
        return;
    }

    SetMarkerAndLine(hV24pt,kBlack,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    hV24pt->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.8,0.5,0.9);
    legend2->AddEntry(hV24pt,Form("v_{2}{4}(p_{T}) |#Delta#eta|>1 cent:%d~%d%%",(int)CentMin,(int)CentMax));
    legend2->Draw();

    if(OutputRoot){
        TFile* fout = new TFile(Form("./ProcessOutput/pTDiffv24Cent%dTo%d_%s%s.root",(int)CentMin,(int)CentMax,FileNameSuffix.c_str(),Subwagon.c_str()),"RECREATE");
        hV24pt->Write();
        fout->Close();
    }
}

void GetNonlinearHistogram(FlowContainer* fc, TH1D*& hCorr422, TH1D*& hCorr24, TH1D*& hCorr42, string GapName="Ch10Gap"){
    fc->SetIDName(Form("%sA",GapName.c_str()));
    TH1D* temp_422A = (TH1D*)fc->GetHistCorrXXVsMulti("422");
    if(!temp_422A){Printf("Can't get temp_422A");return;}
    TH1D* hCorr422A = (TH1D*)temp_422A->Clone();
    delete temp_422A;
    hCorr422A = mergeCentralityToTargetBin(hCorr422A,targetCentralityBins);
    hCorr422A->SetName("hCorr422A");

    fc->SetIDName(Form("%sB",GapName.c_str()));
    TH1D* temp_422B = (TH1D*)fc->GetHistCorrXXVsMulti("422");
    if(!temp_422B){Printf("Can't get temp_422B");return;}
    TH1D* hCorr422B = (TH1D*)temp_422B->Clone();
    delete temp_422B;
    hCorr422B = mergeCentralityToTargetBin(hCorr422B,targetCentralityBins);
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
    hCorr24 = mergeCentralityToTargetBin(hCorr24,targetCentralityBins);

    fc->SetIDName(Form("%s",GapName.c_str()));
    TH1D* temp_42= (TH1D*)fc->GetHistCorrXXVsMulti("42");
    if(!temp_42){Printf("Can't get temp_42");return;}
    hCorr42 = (TH1D*)temp_42->Clone();
    delete temp_42;
    hCorr42 = mergeCentralityToTargetBin(hCorr42,targetCentralityBins);
}

void SetNonlinearValue(TH1D*& target, TH1D*& hCorr422, TH1D*& hCorr24, TH1D*& hCorr42, nonlinearObservableEnum observable=v422){
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
                    target->SetBinContent(i,-1.);
                    target->SetBinError(i,10.);
                    Printf("Warning: in %f %%, Ch10Gap24 is negative", hCorr422->GetBinCenter(i));
                }
            }
            else if(observable==rho422){
                if(hCorr24->GetBinContent(i)>0.&&hCorr42->GetBinContent(i)>0.){
                    target->SetBinContent(i,hCorr422->GetBinContent(i)/(sqrt(hCorr24->GetBinContent(i))*sqrt(hCorr42->GetBinContent(i))));
                    target->SetBinError(i,Error_Rho(hCorr422->GetBinContent(i), hCorr422->GetBinError(i), hCorr24->GetBinContent(i), hCorr24->GetBinError(i),hCorr42->GetBinContent(i), hCorr42->GetBinError(i)));
                }
                else{
                    target->SetBinContent(i,-1.);
                    target->SetBinError(i,10.);
                    Printf("Warning: in %f %%, Ch10Gap24 or Ch10Gap42 is negative", hCorr422->GetBinCenter(i));
                }
            }
        }
        
    }
    delete hCorr422;
    delete hCorr24;
    delete hCorr42;
}

void Output_Nonlinear(string FileNameSuffix, FlowContainer* fc, nonlinearObservableEnum observable=v422, string Subwagon=""){
    
    TCanvas* canvas1 = new TCanvas(Form("Canvas_%s",nonlinearObservableName[observable]),Form("Canvas_%s",nonlinearObservableName[observable]),900,900);
    TH1D* Hist  = new TH1D(Form("%s in %s",nonlinearObservableName[observable],FileNameSuffix.c_str()),Form("%s in %s",nonlinearObservableName[observable],FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(nonlinearUserRangeMap[observable][0]);
    Hist->SetMaximum(nonlinearUserRangeMap[observable][1]);
    Hist->SetXTitle("Centrality (%)");
    Hist->SetYTitle(nonlinearObservableName[observable]);
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
        if(ErrorArray[0][i]!=ErrorArray[0][i])
        hv422->SetBinError(i+1, 0);
        else
        hv422->SetBinError(i+1, ErrorArray[0][i]);
        Printf("%d: %f+-%f",i,hv422->GetBinContent(i+1),hv422->GetBinError(i+1));
    }

    SetMarkerAndLine(hv422,kBlack,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    hv422->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.85,0.5,0.9);
    legend2->AddEntry(hv422,Form("%s |#Delta#eta|>1",nonlinearObservableName[observable]));
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
        legend2->AddEntry(g,Form("%s JHEP 05 (2020) 085, 2020",nonlinearObservableName[observable]));
    }

    if(OutputRoot){
        TFile* fout = new TFile(Form("./ProcessOutput/%s_%s%s.root",nonlinearOutputObservableName[observable],FileNameSuffix.c_str(),Subwagon.c_str()),"RECREATE");
        hv422->Write();
        fout->Close();
        // canvas1->SaveAs(Form("./ProcessOutput/%s.png",nonlinearOutputObservableName[observable]));
    }

}

void Output_Nonlinear(string FileNameSuffix, FlowContainer* fc, FlowContainer* fc_com, nonlinearObservableEnum observable=v422){
    
    TCanvas* canvas1 = new TCanvas(Form("Canvas_%s",nonlinearObservableName[observable]),Form("Canvas_%s",nonlinearObservableName[observable]),900,900);
    TH1D* Hist  = new TH1D(Form("%s in %s",nonlinearObservableName[observable],FileNameSuffix.c_str()),Form("%s in %s",nonlinearObservableName[observable],FileNameSuffix.c_str()),8,0,80);
    Hist->SetMinimum(nonlinearUserRangeMap[observable][0]);
    Hist->SetMaximum(nonlinearUserRangeMap[observable][1]);
    Hist->SetXTitle("Centrality (%)");
    Hist->SetYTitle(nonlinearObservableName[observable]);
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
    legend2->AddEntry(hv422,Form("%s |#Delta#eta|>1",nonlinearObservableName[observable]));
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
        legend2->AddEntry(g,Form("%s JHEP 05 (2020) 085, 2020",nonlinearObservableName[observable]));
    }

    if(OutputRoot){
        TFile* fout = new TFile(Form("./ProcessOutput/%s_%s.root",nonlinearOutputObservableName[observable],FileNameSuffix.c_str()),"RECREATE");
        hv422->Write();
        fout->Close();
        // canvas1->SaveAs(Form("./ProcessOutput/%s.png",nonlinearOutputObservableName[observable]));
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

void Output_NSC(string FileNameSuffix, FlowContainer* fc, string Subwagon=""){
    
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
    NSC32 = mergeCentralityToTargetBin(NSC32, targetCentralityBins);
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
            temp = mergeCentralityToTargetBin(temp, targetCentralityBins);
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
    legend2->AddEntry(NSC32,Form("NSC(3,2) |#Delta#eta|>1"));
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

    if(OutputRoot){
        TFile* fout = new TFile(Form("./ProcessOutput/NSC_%s%s.root",FileNameSuffix.c_str(),Subwagon.c_str()),"RECREATE");
        NSC32->Write();
        fout->Close();
        // canvas1->SaveAs(Form("./ProcessOutput/NSC.png"));
    }
}

TH1D* GetSCklm(FlowContainer* fc, Int_t k, Int_t l, Int_t m){

    TH1D* hCorrklkl=nullptr;
    TH1D* hCorrk2_Full=nullptr;
    TH1D* hCorrl2_Full=nullptr;

    fc->SetIDName(Form("ChFull"));
    hCorrklkl = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d%d%d%d", k, l, k, l));
    if (!hCorrklkl) {Printf("Can't get hCorrklkl"); return nullptr;}
    hCorrk2_Full = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d2", k));
    if (!hCorrk2_Full) {Printf("Can't get hCorrk2_Full"); return nullptr;}
    hCorrl2_Full = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d2", l));
    if (!hCorrl2_Full) {Printf("Can't get hCorrl2_Full"); return nullptr;}
    if (m < 0) {
        TH1D* out = (TH1D*)hCorrklkl->Clone();
        for (int i = 1; i <= hCorrklkl->GetNbinsX(); i++) {
            out->SetBinContent(i, hCorrklkl->GetBinContent(i) - hCorrk2_Full->GetBinContent(i) * hCorrl2_Full->GetBinContent(i));
            double err = Error_SCnm(hCorrklkl->GetBinContent(i), hCorrklkl->GetBinError(i), hCorrk2_Full->GetBinContent(i), hCorrk2_Full->GetBinError(i), hCorrl2_Full->GetBinContent(i), hCorrl2_Full->GetBinError(i));
            out->SetBinError(i, err);
        }
        return out;
    }

    TH1D* hCorrklmklm=nullptr;
    TH1D* hCorrkmkm=nullptr;
    TH1D* hCorrlmlm=nullptr;
    TH1D* hCorrm2_Full=nullptr;
    hCorrklmklm = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d%d%d%d%d%d", k, l, m, k, l, m));
    if (!hCorrklmklm) {Printf("Can't get hCorrklmklm"); return nullptr;}
    hCorrkmkm = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d%d%d%d", k, m, k, m));
    if (!hCorrkmkm) {Printf("Can't get hCorrkmkm"); return nullptr;}
    hCorrlmlm = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d%d%d%d", l, m, l, m));
    if (!hCorrlmlm) {Printf("Can't get hCorrlmlm"); return nullptr;}
    hCorrm2_Full = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d2", m));
    if (!hCorrm2_Full) {Printf("Can't get hCorrm2_Full"); return nullptr;}

    TH1D* out = (TH1D*)hCorrklmklm->Clone();
    for (int i = 1; i <= hCorrklmklm->GetNbinsX(); i++) {
        out->SetBinContent(i, hCorrklmklm->GetBinContent(i) 
        - hCorrklkl->GetBinContent(i) * hCorrm2_Full->GetBinContent(i)
        - hCorrkmkm->GetBinContent(i) * hCorrl2_Full->GetBinContent(i)
        - hCorrlmlm->GetBinContent(i) * hCorrk2_Full->GetBinContent(i)
        + 2*hCorrk2_Full->GetBinContent(i)*hCorrl2_Full->GetBinContent(i)*hCorrm2_Full->GetBinContent(i)
        );
        double err = Error_SCklm(
            hCorrklmklm->GetBinContent(i), hCorrklmklm->GetBinError(i),
            hCorrklkl->GetBinContent(i), hCorrklkl->GetBinError(i),
            hCorrkmkm->GetBinContent(i), hCorrkmkm->GetBinError(i),
            hCorrlmlm->GetBinContent(i), hCorrlmlm->GetBinError(i),
            hCorrk2_Full->GetBinContent(i), hCorrk2_Full->GetBinError(i),
            hCorrl2_Full->GetBinContent(i), hCorrl2_Full->GetBinError(i),
            hCorrm2_Full->GetBinContent(i), hCorrm2_Full->GetBinError(i)
        );
        out->SetBinError(i, err);
    }
    return out;
    
}

void Output_SCklm(string FileNameSuffix, FlowContainer* fc, Int_t k, Int_t l, Int_t m=-1, string Subwagon=""){
    
    TCanvas* canvas1 = nullptr;
    if (m < 0) canvas1 = new TCanvas(Form("Canvas_SC%d%d",k,l),Form("Canvas_SC%d%d",k,l),900,900);
    else canvas1 = new TCanvas(Form("Canvas_SC%d%d%d",k,l,m),Form("Canvas_SC%d%d%d",k,l,m),900,900);
    TH1D* Hist  = nullptr;
    if (m < 0) Hist = new TH1D(Form("SC%d%d in %s",k,l,FileNameSuffix.c_str()),Form("SC%d%d in %s",k,l,FileNameSuffix.c_str()),8,0,80);
    else Hist = new TH1D(Form("SC%d%d%d in %s",k,l,m,FileNameSuffix.c_str()),Form("SC%d%d%d in %s",k,l,m,FileNameSuffix.c_str()),8,0,80);
    if (m<0){
        Hist->SetMinimum(-1e-5);
        Hist->SetMaximum(1e-5);
    }
    else{
        Hist->SetMinimum(-1e-8);
        Hist->SetMaximum(1e-8);
    }
    Hist->SetXTitle("Centrality (%)");
    if(m<0) Hist->SetYTitle(Form("SC(%d,%d)",k,l));
    else Hist->SetYTitle(Form("SC(%d,%d,%d)",k,l,m));
    Hist->Draw();

    fc->SetPropagateErrors(kTRUE);
    TH1D* SCklm = GetSCklm(fc, k, l, m);
    if (!SCklm) {Printf("Can't get SCklm"); return;}
    SCklm = mergeCentralityToTargetBin(SCklm, targetCentralityBins);
    if (m < 0) SCklm->SetName(Form("SC%d%d",k,l));
    else SCklm->SetName(Form("SC%d%d%d",k,l,m));
    

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=1;//NSC
    TObjArray* subsamples = fc->GetSubProfiles();
    int NofSample = subsamples->GetEntries();
    int Nbin = SCklm->GetNbinsX();
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    
    for(int sample=0;sample<NofSample;sample++){
        fc->OverrideMainWithSub(sample,false);
        for(int i=0;i<Nobs;i++){
            TH1D* temp = GetSCklm(fc, k, l, m);
            temp = mergeCentralityToTargetBin(temp, targetCentralityBins);
            temp->SetName(Form("SC%d%d%d_%d",k,l,m,sample));
            if(!temp){
                Printf("Can't get SC%d%d%d_%d",k,l,m,sample);
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
        SCklm->SetBinError(i+1, ErrorArray[0][i]);
    }

    SetMarkerAndLine(SCklm,kBlack,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    SCklm->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.85,0.5,0.9);
    if (m < 0) legend2->AddEntry(SCklm,Form("SC(%d,%d)",k,l));
    else legend2->AddEntry(SCklm,Form("SC(%d,%d,%d)",k,l,m));
    legend2->Draw();

    if(OutputRoot){
        TFile* fout = nullptr;
        if (m < 0) fout = new TFile(Form("./ProcessOutput/SC%d%d_%s%s.root",k,l,FileNameSuffix.c_str(),Subwagon.c_str()),"RECREATE");
        else fout = new TFile(Form("./ProcessOutput/SC%d%d%d_%s%s.root",k,l,m,FileNameSuffix.c_str(),Subwagon.c_str()),"RECREATE");
        SCklm->Write();
        fout->Close();
    }
}

TH1D* GetNSCklm(FlowContainer* fc, Int_t k, Int_t l, Int_t m){

    TH1D* hCorrklkl=nullptr;
    TH1D* hCorrk2_Full=nullptr;
    TH1D* hCorrl2_Full=nullptr;
    TH1D* hCorrk2_Gap=nullptr;
    TH1D* hCorrl2_Gap=nullptr;

    fc->SetIDName(Form("ChFull"));
    hCorrklkl = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d%d%d%d", k, l, k, l));
    if (!hCorrklkl) {Printf("Can't get hCorrklkl"); return nullptr;}

    TH1D* temp_k2_Full= (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d2", k));
    if(!temp_k2_Full){Printf("Can't get temp_k2_Full");return nullptr;}
    hCorrk2_Full = (TH1D*)temp_k2_Full->Clone();
    hCorrk2_Full->SetName("hCorrk2_Full");
    delete temp_k2_Full;

    TH1D* temp_l2_Full= (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d2", l));
    if(!temp_l2_Full){Printf("Can't get temp_l2_Full");return nullptr;}
    hCorrl2_Full = (TH1D*)temp_l2_Full->Clone();
    hCorrl2_Full->SetName("hCorrl2_Full");
    delete temp_l2_Full;

    fc->SetIDName(Form("Ch10Gap"));
    TH1D* temp_k2_Gap= (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d2", k));
    if(!temp_k2_Gap){Printf("Can't get temp_k2_Gap");return nullptr;}
    hCorrk2_Gap = (TH1D*)temp_k2_Gap->Clone();
    hCorrk2_Gap->SetName("hCorrk2_Gap");
    delete temp_k2_Gap;

    TH1D* temp_l2_Gap= (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d2", l));
    if(!temp_l2_Gap){Printf("Can't get temp_l2_Gap");return nullptr;}
    hCorrl2_Gap = (TH1D*)temp_l2_Gap->Clone();
    hCorrl2_Gap->SetName("hCorrl2_Gap");
    delete temp_l2_Gap;
    if (m < 0) {
        TH1D* out = (TH1D*)hCorrklkl->Clone();
        for (int i = 1; i <= hCorrklkl->GetNbinsX(); i++) {
            out->SetBinContent(i, 
            (hCorrklkl->GetBinContent(i) - hCorrk2_Full->GetBinContent(i) * hCorrl2_Full->GetBinContent(i))
            /(hCorrk2_Gap->GetBinContent(i)*hCorrl2_Gap->GetBinContent(i)));

            double err = Error_NSC(hCorrklkl->GetBinContent(i), hCorrklkl->GetBinError(i), 
            hCorrk2_Full->GetBinContent(i), hCorrk2_Full->GetBinError(i), 
            hCorrl2_Full->GetBinContent(i), hCorrl2_Full->GetBinError(i), 
            hCorrk2_Gap->GetBinContent(i), hCorrk2_Gap->GetBinError(i), 
            hCorrl2_Gap->GetBinContent(i), hCorrl2_Gap->GetBinError(i));

            out->SetBinError(i, err);
        }
        return out;
    }

    TH1D* hCorrklmklm=nullptr;
    TH1D* hCorrkmkm=nullptr;
    TH1D* hCorrlmlm=nullptr;
    TH1D* hCorrm2_Full=nullptr;
    TH1D* hCorrm2_Gap=nullptr;

    fc->SetIDName(Form("ChFull"));
    hCorrklmklm = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d%d%d%d%d%d", k, l, m, k, l, m));
    if (!hCorrklmklm) {Printf("Can't get hCorrklmklm"); return nullptr;}
    hCorrkmkm = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d%d%d%d", k, m, k, m));
    if (!hCorrkmkm) {Printf("Can't get hCorrkmkm"); return nullptr;}
    hCorrlmlm = (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d%d%d%d", l, m, l, m));
    if (!hCorrlmlm) {Printf("Can't get hCorrlmlm"); return nullptr;}
    TH1D* temp_m2_Full= (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d2", m));
    if(!temp_m2_Full){Printf("Can't get temp_m2_Full");return nullptr;}
    hCorrm2_Full = (TH1D*)temp_m2_Full->Clone();
    hCorrm2_Full->SetName("hCorrm2_Full");
    delete temp_m2_Full;

    fc->SetIDName(Form("Ch10Gap"));
    TH1D* temp_m2_Gap= (TH1D*)fc->GetHistCorrXXVsMulti(Form("%d2", m));
    if(!temp_m2_Gap){Printf("Can't get temp_m2_Gap");return nullptr;}
    hCorrm2_Gap = (TH1D*)temp_m2_Gap->Clone();
    hCorrm2_Gap->SetName("hCorrm2_Gap");
    delete temp_m2_Gap;

    TH1D* out = (TH1D*)hCorrklmklm->Clone();
    for (int i = 1; i <= hCorrklmklm->GetNbinsX(); i++) {
        double scklm = hCorrklmklm->GetBinContent(i) 
        - hCorrklkl->GetBinContent(i) * hCorrm2_Full->GetBinContent(i)
        - hCorrkmkm->GetBinContent(i) * hCorrl2_Full->GetBinContent(i)
        - hCorrlmlm->GetBinContent(i) * hCorrk2_Full->GetBinContent(i)
        + 2*hCorrk2_Full->GetBinContent(i)*hCorrl2_Full->GetBinContent(i)*hCorrm2_Full->GetBinContent(i);

        double bottom = hCorrk2_Gap->GetBinContent(i)*hCorrl2_Gap->GetBinContent(i)*hCorrm2_Gap->GetBinContent(i);
        out->SetBinContent(i, scklm/bottom);
        double err = Error_NSCklm(
            hCorrklmklm->GetBinContent(i), hCorrklmklm->GetBinError(i),
            hCorrklkl->GetBinContent(i), hCorrklkl->GetBinError(i),
            hCorrkmkm->GetBinContent(i), hCorrkmkm->GetBinError(i),
            hCorrlmlm->GetBinContent(i), hCorrlmlm->GetBinError(i),
            hCorrk2_Full->GetBinContent(i), hCorrk2_Full->GetBinError(i),
            hCorrl2_Full->GetBinContent(i), hCorrl2_Full->GetBinError(i),
            hCorrm2_Full->GetBinContent(i), hCorrm2_Full->GetBinError(i),
            hCorrk2_Gap->GetBinContent(i), hCorrk2_Gap->GetBinError(i),
            hCorrl2_Gap->GetBinContent(i), hCorrl2_Gap->GetBinError(i),
            hCorrm2_Gap->GetBinContent(i), hCorrm2_Gap->GetBinError(i)
        );
        out->SetBinError(i, err);
    }
    return out;
    
}


void Output_NSCklm(string FileNameSuffix, FlowContainer* fc, Int_t k, Int_t l, Int_t m=-1, string Subwagon=""){
    
    TCanvas* canvas1 = nullptr;
    if (m < 0) canvas1 = new TCanvas(Form("Canvas_NSC%d%d",k,l),Form("Canvas_NSC%d%d",k,l),900,900);
    else canvas1 = new TCanvas(Form("Canvas_NSC%d%d%d",k,l,m),Form("Canvas_NSC%d%d%d",k,l,m),900,900);
    TH1D* Hist  = nullptr;
    if (m < 0) Hist = new TH1D(Form("NSC%d%d in %s",k,l,FileNameSuffix.c_str()),Form("NSC%d%d in %s",k,l,FileNameSuffix.c_str()),8,0,80);
    else Hist = new TH1D(Form("NSC%d%d%d in %s",k,l,m,FileNameSuffix.c_str()),Form("NSC%d%d%d in %s",k,l,m,FileNameSuffix.c_str()),8,0,80);
    if (m<0){
        Hist->SetMinimum(-1e-5);
        Hist->SetMaximum(1e-5);
    }
    else{
        Hist->SetMinimum(-1e-8);
        Hist->SetMaximum(1e-8);
    }
    Hist->SetXTitle("Centrality (%)");
    if(m<0) Hist->SetYTitle(Form("NSC(%d,%d)",k,l));
    else Hist->SetYTitle(Form("NSC(%d,%d,%d)",k,l,m));
    Hist->Draw();

    fc->SetPropagateErrors(kTRUE);
    TH1D* NSCklm = GetNSCklm(fc, k, l, m);
    if (!NSCklm) {Printf("Can't get NSCklm"); return;}
    NSCklm = mergeCentralityToTargetBin(NSCklm, targetCentralityBins);
    if (m < 0) NSCklm->SetName(Form("NSC%d%d",k,l));
    else NSCklm->SetName(Form("NSC%d%d%d",k,l,m));
    

    std::vector<std::vector<std::vector<double>>> ValueArray;
    std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    std::vector<std::vector<double>> ErrorArray;
    int Nobs=1;//NSC
    TObjArray* subsamples = fc->GetSubProfiles();
    int NofSample = subsamples->GetEntries();
    int Nbin = NSCklm->GetNbinsX();
    ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    
    for(int sample=0;sample<NofSample;sample++){
        fc->OverrideMainWithSub(sample,false);
        for(int i=0;i<Nobs;i++){
            TH1D* temp = GetNSCklm(fc, k, l, m);
            temp = mergeCentralityToTargetBin(temp, targetCentralityBins);
            temp->SetName(Form("NSC%d%d%d_%d",k,l,m,sample));
            if(!temp){
                Printf("Can't get NSC%d%d%d_%d",k,l,m,sample);
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
        NSCklm->SetBinError(i+1, ErrorArray[0][i]);
    }

    SetMarkerAndLine(NSCklm,kBlack,kFullCircle,kSolid,1.0);
    gStyle->SetOptStat("");
    NSCklm->Draw("ESames");
    TLegend* legend2 = new TLegend(0.2,0.85,0.5,0.9);
    if (m < 0) legend2->AddEntry(NSCklm,Form("NSC(%d,%d)",k,l));
    else legend2->AddEntry(NSCklm,Form("NSC(%d,%d,%d)",k,l,m));
    legend2->Draw();

    if(OutputRoot){
        TFile* fout = nullptr;
        if (m < 0) fout = new TFile(Form("./ProcessOutput/NSC%d%d_%s%s.root",k,l,FileNameSuffix.c_str(),Subwagon.c_str()),"RECREATE");
        else fout = new TFile(Form("./ProcessOutput/NSC%d%d%d_%s%s.root",k,l,m,FileNameSuffix.c_str(),Subwagon.c_str()),"RECREATE");
        NSCklm->Write();
        fout->Close();
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

#endif /* FlowContainerCalculation_h */
