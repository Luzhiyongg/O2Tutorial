#include "TFile.h"
#include "TH1D.h"
#include <cmath>
#include <string>
#include "TCanvas.h"
#include "TLegend.h"
#include "../MyRun3Analysis/include/ErrorPropagation.h"
#include "../MyRun3Analysis/include/FlowContainerCalculation.h"

TH1D* GetMCptDiffvn(std::string FileNameSuffix, FlowContainer* fc, Int_t n=2, Double_t CentMin=0., Double_t CentMax=5., std::string Subwagon=""){

    fc->SetIDName("ChFull");
    fc->SetPropagateErrors(kTRUE);
    // TH1D* hV22pt = (TH1D*)fc->GetVN2VsMulti(n);
    TH1D* hV22pt = (TH1D*)fc->GetVN2VsPt(n,CentMin,CentMax);
    hV22pt->SetName(Form("pTDiffv%d",n));
    if(!hV22pt){
        Printf("Can't get hV22");
        return nullptr;
    }

    // hV22pt->Draw();

    // std::vector<std::vector<std::vector<double>>> ValueArray;
    // std::vector<std::vector<std::vector<double>>> ValueErrorArray;
    // std::vector<std::vector<double>> ErrorArray;
    // int Nobs=1;//v22(pT)
    // TObjArray* subsamples = fc->GetSubProfiles();
    // int NofSample = subsamples->GetEntries();
    // int Nbin = hV22pt->GetNbinsX();
    // ResizeValueArray(ValueArray,ValueErrorArray,ErrorArray,Nobs,NofSample,Nbin);
    
    // for(int sample=0;sample<NofSample;sample++){
    //     fc->OverrideMainWithSub(sample,false);
    //     for(int i=0;i<Nobs;i++){
    //         TH1D* temp = (TH1D*)fc->GetVN2VsPt(n,CentMin,CentMax);
    //         temp->SetName(Form("pTDiffv2_%d",sample));
    //         if(!temp){
    //             Printf("Can't get pTDiffv2_%d",sample);
    //             return nullptr;
    //         }
    //         for(int j=0;j<temp->GetNbinsX();j++){
    //             ValueArray[i][sample][j] = temp->GetBinContent(j+1);
    //             ValueErrorArray[i][sample][j] = temp->GetBinError(j+1);
    //             // Printf("pTDiffv2_%d_%d = %f +- %f",sample,j,ValueArray[i][sample][j],ValueErrorArray[i][sample][j]);
    //         }
    //     }
    // }
    // for(int i=0;i<Nobs;i++){
    //     CalculateBootstrapError(ValueArray[i],ValueErrorArray[i],ErrorArray[i]);
    // }
    // for(int i=0;i<Nbin;i++){
    //     hV22pt->SetBinError(i+1, ErrorArray[0][i]);
    // }

    return hV22pt;
}

void ProduceCorrectionFactor(std::string filenameSubfix, int CentMin, int CentMax, std::string SubwagonName){
    TFile* f = new TFile(Form("./AnalysisResults/AnalysisResults%s.root",filenameSubfix.c_str()),"READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: File not found." << std::endl;
        return;
    }

    FlowContainer* fcReco = (FlowContainer*)f->Get(Form("flow-pt-efficiency%s/FlowContainerReco",SubwagonName.c_str()));
    if (!fcReco) {
        std::cerr << "Error: FlowContainerReco not found." << std::endl;
        return;
    }
    TH1D* hV2Reco = GetMCptDiffvn(filenameSubfix, fcReco, 2, CentMin, CentMax, SubwagonName);

    hV2Reco->Draw();
    SetMarkerAndLine(hV2Reco,kBlack,kFullCircle,kSolid,1.0);

    FlowContainer* fcTrue = (FlowContainer*)f->Get(Form("flow-pt-efficiency%s/FlowContainerTrue",SubwagonName.c_str()));
    if (!fcTrue) {
        std::cerr << "Error: FlowContainerTrue not found." << std::endl;
        return;
    }
    TH1D* hV2True = GetMCptDiffvn(filenameSubfix, fcTrue, 2, CentMin, CentMax, SubwagonName);

    hV2True->Draw("SAME");
    SetMarkerAndLine(hV2True,kRed,kFullCircle,kSolid,1.0);

    TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(hV2Reco,"Reco","p");
    leg->AddEntry(hV2True,"True","p");
    leg->Draw();

    // hV2Reco->Divide(hV2True);
    // hV2Reco->Draw();


    // hV2Reco->SetName("hRatio");
    // TFile* fOut = new TFile(Form("./MCcorrection/CorrectionFactor_Cent%dTo%d%s.root",Centmin,Centmax,filenameSubfix.c_str()),"RECREATE");
    // hV2Reco->Write();
    // fOut->Close();
    // f->Close();
}

void CumulantPtDepV2(){
    vector<double> pTDiffCent={0,5,10,20,30,40,50,60,70,80};
    
    std::string filenameSubfix = "_LHC24k2_353179";
    std::string SubwagonName = "";

    ProduceCorrectionFactor(filenameSubfix, 10, 20, SubwagonName);
    // for (uint j=0; j<pTDiffCent.size()-1; j++) {
    //     fc = (FlowContainer*)f->Get(Form("flow-pt-efficiency%s/FlowContainer",SubwagonName.c_str()));
    //     ProduceCorrectionFactor(FileNameSuffix, fc, 2, pTDiffCent[j], pTDiffCent[j+1], SubwagonName);
    // }
}
