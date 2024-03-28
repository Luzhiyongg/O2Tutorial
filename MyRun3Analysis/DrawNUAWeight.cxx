//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"
#pragma GCC diagnostic ignored "-Wwritable-strings"

#include "TFile.h"
#include "GFWWeights.h"
#include <string>
using namespace std;


void DrawNUAWeight(string FileNameSuffix = "LHC23zzh_pass2"){
    TFile* f = new TFile(Form("./AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    // FlowContainer* fc = (FlowContainer*)f->Get("flow-pb-pb-task/FlowContainer");
    GFWWeights* weight = (GFWWeights*)f->Get("flow-pb-pb-task/weights");
    TH3D* h3d = (TH3D*)weight->GetDataArray();
    // h3d->Draw();
    TH1D* phiDis = weight->GetdNdPhi();
    phiDis->Draw();
}
