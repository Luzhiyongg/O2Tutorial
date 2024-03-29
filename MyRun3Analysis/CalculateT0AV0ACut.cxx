//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TList.h"
#include <string>
#include <vector>
using namespace std;


void CalculateT0AV0ACut(string FileNameSuffix = "LHC23zzh_pass2"){
    TFile* f = new TFile(Form("./AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    TH2D* V0AT0A = (TH2D*)f->Get("flow-pb-pb-task/multV0A_multT0A");
    // V0AT0A->Draw();
    Int_t N_v0a = V0AT0A->GetYaxis()->GetNbins();
    Int_t N_t0a = V0AT0A->GetXaxis()->GetNbins();
    Double_t Min_t0a = V0AT0A->GetXaxis()->GetXmin();
    Double_t Max_t0a = V0AT0A->GetXaxis()->GetXmax();
    vector<double> MeanDis(N_t0a);
    vector<double> SigmaDis(N_t0a);
    TH1D* V0AT0AMean = new TH1D("V0AT0AMean","V0AT0AMean", N_t0a, Min_t0a, Max_t0a);
    TH1D* V0AT0ASigma = new TH1D("V0AT0ASigma","V0AT0ASigma", N_t0a, Min_t0a, Max_t0a);
    // for every T0A, calculate the Gaussian mean and sigma
    for(Int_t i=1;i<=N_t0a;i++){

        Double_t Sum=0.;
        Double_t SquareSum=0.;
        Double_t SquareMean=0.;
        Double_t Counter=0.;
        for(Int_t j=1;j<=N_v0a;j++){
            if(V0AT0A->GetBinContent(i,j)>0){
                Printf("%f * %f",V0AT0A->GetYaxis()->GetBinCenter(j), V0AT0A->GetBinContent(i,j));
                Sum += V0AT0A->GetYaxis()->GetBinCenter(j) * V0AT0A->GetBinContent(i,j);
                SquareSum += V0AT0A->GetYaxis()->GetBinCenter(j) * V0AT0A->GetYaxis()->GetBinCenter(j) * V0AT0A->GetBinContent(i,j);
                Counter+=V0AT0A->GetBinContent(i,j);
            }
        }
        if(Counter>1.){
            MeanDis[i-1] = Sum / Counter;
            V0AT0AMean->SetBinContent(i,MeanDis[i-1]);
            SquareMean = SquareSum / Counter;
            SigmaDis[i-1] = sqrt( (Counter/(Counter-1.)) * (SquareMean - MeanDis[i-1]*MeanDis[i-1]) );
            V0AT0ASigma->SetBinContent(i,SigmaDis[i-1]);
        }
    }

    TF1* f1 = new TF1("f1","[0]*x+[1]",Min_t0a,Max_t0a);
    TF1* f2 = new TF1("f2","[0]*x+[1]",Min_t0a,Max_t0a);
    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    V0AT0AMean->Draw();
    V0AT0AMean->Fit(f1);
    // TCanvas* c2 = new TCanvas("c2","c2",600,600);
    // V0AT0ASigma->Draw();
    // V0AT0ASigma->Fit(f2);

    TH2D* afterCut = (TH2D*)V0AT0A->Clone();
    for(Int_t i=1;i<=N_t0a;i++){
        for(Int_t j=1;j<=N_v0a;j++){
            if(abs(afterCut->GetYaxis()->GetBinCenter(j)-f1->Eval(afterCut->GetXaxis()->GetBinCenter(i))) > 5*SigmaDis[i-1]){
                afterCut->SetBinContent(i,j,0);
            }
        }
    }
    TCanvas* c3 = new TCanvas("c3","c3",800,600);
    c3->SetLogz();
    gStyle->SetOptStat("");
    afterCut->Draw("colorZ");

    // save in root file
    // TFile* output = new TFile(Form("./T0AV0ACut/T0AV0ACut_%s.root",FileNameSuffix.c_str()),"RECREATE");
    // TDirectory* dir = output->mkdir(Form("%s",FileNameSuffix.c_str()));
    // dir->cd();
    // V0AT0AMean->Write();
    // V0AT0ASigma->Write();

    // specific the dataset when upload to ccdb (set the valid run number)
    TFile* output = new TFile(Form("./T0AV0ACut/T0AV0ACut_%s.root",FileNameSuffix.c_str()),"RECREATE");
    TList* outputL = new TList();
    outputL->Add(V0AT0AMean);
    outputL->Add(V0AT0ASigma);
    outputL->Write("ccdb_object",1);
    output->Close();

}
