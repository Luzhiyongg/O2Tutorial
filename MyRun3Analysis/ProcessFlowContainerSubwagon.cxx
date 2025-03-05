/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch) 
 * @Date: 2024-02-01
 * @Last Modified by: Zhiyong Lu
 * @Last Modified time: 2025-03-05 22:01:25
 */
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

void ProcessFlowContainerSubwagon(string FileNameSuffix = "LHC23_PbPb_pass4_344339"){
    // Produce flow results in root file for each subwagon

    TFile* f = new TFile(Form("./AnalysisResults/AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    vector<string> SubwagonNames = {""};
    vector<double> pTDiffCent={0,5,10,20,30,40,50,60,70,80};

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
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
            Output_ptDiffvn(FileNameSuffix, fc, 2, pTDiffCent[j], pTDiffCent[j+1], SubwagonNames[i]);
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
            Output_ptDiffvn(FileNameSuffix, fc, 3, pTDiffCent[j], pTDiffCent[j+1], SubwagonNames[i]);
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
            Output_ptDiffvn(FileNameSuffix, fc, 4, pTDiffCent[j], pTDiffCent[j+1], SubwagonNames[i]);
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
            Output_PtDiffVnm(FileNameSuffix, fc, "ChFull", 2, 4, pTDiffCent[j], pTDiffCent[j+1], SubwagonNames[i]);
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
            Output_PtDiffVnm(FileNameSuffix, fc, "Ch10Gap", 2, 4, pTDiffCent[j], pTDiffCent[j+1], SubwagonNames[i]);
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            fc = (FlowContainer*)f->Get(Form("flow-task%s/FlowContainer",SubwagonNames[i].c_str()));
            Output_PtDiffVnm(FileNameSuffix, fc, "ChFull", 2, 6, pTDiffCent[j], pTDiffCent[j+1], SubwagonNames[i]);
        }
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
