//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
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
    TH1D* gMean = new TH1D("gMean","gMean", N_t0a, Min_t0a, Max_t0a);
    TH1D* gSigma = new TH1D("gSigma","gSigma", N_t0a, Min_t0a, Max_t0a);
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
            gMean->SetBinContent(i,MeanDis[i-1]);
            SquareMean = SquareSum / Counter;
            SigmaDis[i-1] = sqrt( (Counter/(Counter-1.)) * (SquareMean - MeanDis[i-1]*MeanDis[i-1]) );
            gSigma->SetBinContent(i,SigmaDis[i-1]);
        }
    }

    TF1* f1 = new TF1("f1","[0]*x+[1]",Min_t0a,Max_t0a);
    TF1* f2 = new TF1("f2","[0]*x+[1]",Min_t0a,Max_t0a);
    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    gMean->Draw();
    gMean->Fit(f1);
    TCanvas* c2 = new TCanvas("c2","c2",600,600);
    gSigma->Draw();
    gSigma->Fit(f2);
}
