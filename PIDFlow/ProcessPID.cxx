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

void ProcessPID(std::string FileNameSuffix = "LHC23zzh_pass4_small_356376_PID"){
    // Produce flow results in root file for each subwagon

    TFile* f = new TFile(Form("./AnalysisResults/AnalysisResults_%s.root",FileNameSuffix.c_str()),"READ");
    std::vector<std::string> SubwagonNames = {""};
    std::vector<double> pTDiffCent={0,5,10,20,30,40,50,60,70,80};

    outputDir = "./PIDProcessOutput"; // change output dir

    for(uint i=0; i < SubwagonNames.size(); i++){
        FlowContainer* fc = (FlowContainer*)f->Get(Form("flow-pbpb-pikp%s/FlowContainer",SubwagonNames[i].c_str()));
        if(!fc){
            Printf("can not get flow container");
            return;
        }
        Output_VnmWithID(FileNameSuffix, fc, 2, 2, "Ch08Gap", SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-pbpb-pikp%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_VnmWithID(FileNameSuffix, fc, 2, 2, "Pi08Gap", SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-pbpb-pikp%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_VnmWithID(FileNameSuffix, fc, 2, 2, "Ka08Gap", SubwagonNames[i]);
        fc = (FlowContainer*)f->Get(Form("flow-pbpb-pikp%s/FlowContainer",SubwagonNames[i].c_str()));
        Output_VnmWithID(FileNameSuffix, fc, 2, 2, "Pr08Gap", SubwagonNames[i]);
        // fc = (FlowContainer*)f->Get(Form("flow-pbpb-pikp%s/FlowContainer",SubwagonNames[i].c_str()));
        // Output_Vnm(FileNameSuffix, fc, 2, 4, "ChFull", SubwagonNames[i]);
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            fc = (FlowContainer*)f->Get(Form("flow-pbpb-pikp%s/FlowContainer",SubwagonNames[i].c_str()));
            Output_PtDiffVnm(FileNameSuffix, fc, "Ch08Gap", 2, 2, pTDiffCent[j], pTDiffCent[j+1], SubwagonNames[i]);
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            fc = (FlowContainer*)f->Get(Form("flow-pbpb-pikp%s/FlowContainer",SubwagonNames[i].c_str()));
            Output_PtDiffVnm(FileNameSuffix, fc, "Pi08Gap", 2, 2, pTDiffCent[j], pTDiffCent[j+1], SubwagonNames[i]);
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            fc = (FlowContainer*)f->Get(Form("flow-pbpb-pikp%s/FlowContainer",SubwagonNames[i].c_str()));
            Output_PtDiffVnm(FileNameSuffix, fc, "Ka08Gap", 2, 2, pTDiffCent[j], pTDiffCent[j+1], SubwagonNames[i]);
        }
        for (uint j=0; j<pTDiffCent.size()-1; j++) {
            fc = (FlowContainer*)f->Get(Form("flow-pbpb-pikp%s/FlowContainer",SubwagonNames[i].c_str()));
            Output_PtDiffVnm(FileNameSuffix, fc, "Pr08Gap", 2, 2, pTDiffCent[j], pTDiffCent[j+1], SubwagonNames[i]);
        }
    }

    return;
    
}
