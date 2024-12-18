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


void TestT0AV0ACut(string FileNameSuffix = "LHC23zzh_pass2_small"){
    TFile* f = new TFile(Form("./AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    TH2D* V0AT0A = (TH2D*)f->Get("flow-pb-pb-task/multV0A_multT0A");
    // V0AT0A->Draw();
    Int_t N_v0a = V0AT0A->GetYaxis()->GetNbins();
    Int_t N_t0a = V0AT0A->GetXaxis()->GetNbins();
    Double_t Min_t0a = V0AT0A->GetXaxis()->GetXmin();
    Double_t Max_t0a = V0AT0A->GetXaxis()->GetXmax();

    TF1* f1 = new TF1("f1","[0]+[1]*x",Min_t0a,Max_t0a);
    TF1* f2 = new TF1("f2","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",Min_t0a,Max_t0a);
    f1->SetParameters(-1601.0581, 9.417652e-01);
    f2->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);

    TH2D* afterCut = (TH2D*)V0AT0A->Clone();
    for(Int_t i=1;i<=N_t0a;i++){
        for(Int_t j=1;j<=N_v0a;j++){
            // if(abs(afterCut->GetYaxis()->GetBinCenter(j)-f1->Eval(afterCut->GetXaxis()->GetBinCenter(i))) > 5*SigmaDis[i-1]){
            //     afterCut->SetBinContent(i,j,0);
            // }
            if(abs(afterCut->GetYaxis()->GetBinCenter(j)-f1->Eval(afterCut->GetXaxis()->GetBinCenter(i))) > 5*f2->Eval(afterCut->GetXaxis()->GetBinCenter(i))){
                afterCut->SetBinContent(i,j,0);
            }
        }
    }
    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    c1->SetLogz();
    gStyle->SetOptStat("");
    V0AT0A->Draw("colorZ");

    TCanvas* c3 = new TCanvas("c3","c3",800,600);
    c3->SetLogz();
    gStyle->SetOptStat("");
    afterCut->Draw("colorZ");
}
