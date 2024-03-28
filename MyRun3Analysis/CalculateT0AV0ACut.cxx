//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "TH2D.h"
#include <string>
#include <vector>
using namespace std;


void CalculateT0AV0ACut(string FileNameSuffix = "LHC23zzh_pass2"){
    TFile* f = new TFile(Form("./AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    TH2D* V0AT0A = (TH2D*)f->Get("flow-pb-pb-task/multV0A_multT0A");
    V0AT0A->Draw();
    Int_t N_v0a = V0AT0A->GetYaxis()->GetNbins();
    Int_t N_t0a = V0AT0A->GetXaxis()->GetNbins();
    vector<double> MeanDis(N_t0a);
    vector<double> SigmaDis(N_t0a);
    // for every T0A, calculate the Gaussian mean and sigma
    for(Int_t i=1;i<=N_t0a;i++){

        Double_t Sum=0.;
        Double_t Counter=0.;
        for(Int_t j=1;j<=N_v0a;j++){
            Sum += N_v0a = V0AT0A->GetYaxis()->GetBinCenter(j) * V0AT0A->GetBinContent(i,j);
        }
    }
}
