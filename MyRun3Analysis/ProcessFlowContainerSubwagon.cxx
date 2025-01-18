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
#include "include/FlowContainerCalculation.h"

void ProcessFlowContainerSubwagon(string FileNameSuffix = "LHC23zzh_pass4_small_327815"){
    // Produce flow results in root file for each subwagon

    TFile* f = new TFile(Form("./AnalysisResults/AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    vector<string> SubwagonNames = {"","_kIsGoodITSLayersAll"};
    for(uint i=0; i < SubwagonNames.size(); i++){
        FlowContainer* fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        if(!fc){
            Printf("can not get flow container");
            return;
        }
        Output_vn(FileNameSuffix, fc, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_Vnm(FileNameSuffix, fc, 2, 4, "ChFull", SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_Vnm(FileNameSuffix, fc, 2, 6, "ChFull", SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_Vnm(FileNameSuffix, fc, 2, 8, "ChFull", SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_Vnm(FileNameSuffix, fc, 2, 10, "ChFull", SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_ptDiffvn(FileNameSuffix, fc, 0, 5, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_ptDiffvn(FileNameSuffix, fc, 30, 40, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_ptDiffvn(FileNameSuffix, fc, 5, 10, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_ptDiffvn(FileNameSuffix, fc, 10, 20, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_ptDiffvn(FileNameSuffix, fc, 20, 30, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_ptDiffvn(FileNameSuffix, fc, 30, 40, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_ptDiffvn(FileNameSuffix, fc, 40, 50, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_ptDiffvn(FileNameSuffix, fc, 50, 60, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_ptDiffvn(FileNameSuffix, fc, 60, 70, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_ptDiffvn(FileNameSuffix, fc, 70, 80, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_Nonlinear(FileNameSuffix, fc, v422, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_Nonlinear(FileNameSuffix, fc, chi422, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_Nonlinear(FileNameSuffix, fc, rho422, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_NSC(FileNameSuffix, fc, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_NSCklm(FileNameSuffix, fc, 2, 3, -1, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_NSCklm(FileNameSuffix, fc, 2, 4, -1, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_NSCklm(FileNameSuffix, fc, 2, 3, 4, SubwagonNames[i]);
        // fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        // Output_SCklm(FileNameSuffix, fc, 2, 3, 5, SubwagonNames[i]);
        // fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        // Output_SCklm(FileNameSuffix, fc, 2, 4, 6, SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_NSCklm(FileNameSuffix, fc, 3, 4, 5, SubwagonNames[i]);
    }

    return;
    
}
