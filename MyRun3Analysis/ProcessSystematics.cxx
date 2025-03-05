/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch) 
 * @Date: 2025-03-05 22:05:29 
 * @Last Modified by:   Zhiyong Lu 
 * @Last Modified time: 2025-03-05 22:05:29 
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
#include "TF1.h"
#include <cstring>
#include <vector>
#include <map>
#include <array>
#include "include/BasicSetting.h"
#include "include/ErrorPropagation.h"
#include "include/ProcessDefine.h"

using namespace std;

std::map<int,bool> IfCheckObservable = {
    {kVn,true},
    {kV24,true},
    {kV26,true},
    {kV28,true},
    {kV210,true},
    {kV422,true},
    {kChi422,true},
    {kRho422,true},
    {kNSC23,true},
    {kNSC24,true},
    {kNSC234,true},
    {kNSC345,true},
    {kpTDiffv2,true},
    {kpTDiffv3,true},
    {kpTDiffv4,true},
    {kpTDiffv24ChFull,true},
    {kpTDiffv24Ch10Gap,false},
    {kpTDiffv26ChFull,true}
};

void CalculateSystematic(int kobs, string obsName, string DefaultFileNameSuffix, string SysFileNameSuffix, string SysDescription, int pTDiffCentMin=-1, int pTDiffCentMax=-1) {

    string DefaultFile = "./ProcessOutput/" + ObservableFilesMap[kobs];
    if (pTDiffCentMin >= 0 && pTDiffCentMax >= 0)
        DefaultFile += "Cent" + to_string(pTDiffCentMin) + "To" + to_string(pTDiffCentMax);
    DefaultFile += "_" + DefaultFileNameSuffix + ".root";
    TFile* fDefault = TFile::Open(DefaultFile.c_str(), "READ");
    if (!fDefault || fDefault->IsZombie()) {
        cout << "Error: cannot open file " << DefaultFile << endl;
        return;
    }
    TH1D* hDefault = (TH1D*)fDefault->Get(obsName.c_str());
    if (!hDefault) {
        cout << "Error: cannot find histogram " << obsName << " in file " << DefaultFile << endl;
        return;
    }
    // Printf("Default histogram: %s", hDefault->GetName());

    string SysFile = "./ProcessOutput/" + ObservableFilesMap[kobs];
    if (pTDiffCentMin >= 0 && pTDiffCentMax >= 0)
        SysFile += "Cent" + to_string(pTDiffCentMin) + "To" + to_string(pTDiffCentMax);
    SysFile += "_" + SysFileNameSuffix + ".root";
    TFile* fSys = TFile::Open(SysFile.c_str(), "READ");
    if (!fSys || fSys->IsZombie()) {
        cout << "Error: cannot open file " << SysFile << endl;
        return;
    }
    // Printf("Systematic file: %s", SysFile.c_str());
    TH1D* hSys = (TH1D*)fSys->Get(obsName.c_str());
    if (!hSys) {
        cout << "Error: cannot find histogram " << obsName << " in file " << SysFile << endl;
        return;
    }
    // Printf("Systematic histogram: %s", hSys->GetName());

    Int_t nbins = hDefault->GetNbinsX();
    TH1D* hUncertainty = (TH1D*)hDefault->Clone("hUncertainty");
    TH1D* hBarlow = (TH1D*)hDefault->Clone("hBarlow");
    for (int ibin=1; ibin<=nbins; ibin++) {
        double defaultBinContent = hDefault->GetBinContent(ibin);
        double sysBinContent = hSys->GetBinContent(ibin);
        if (kobs>=kpTDiffv2 && (hDefault->GetBinCenter(ibin)>MaxpT)) {
            hUncertainty->SetBinContent(ibin, -1);
            hUncertainty->SetBinError(ibin, 10);
            hBarlow->SetBinContent(ibin, -1);
            // Printf("skip pT at %f GeV/c", hDefault->GetBinCenter(ibin));
            continue;
        }
        if (defaultBinContent==-1 || sysBinContent==-1 || defaultBinContent!=defaultBinContent || sysBinContent!=sysBinContent || isZero(defaultBinContent) || isZero(sysBinContent)) {
            hUncertainty->SetBinContent(ibin, -1);
            hUncertainty->SetBinError(ibin, 10);
            hBarlow->SetBinContent(ibin, -1);
            // Printf("obs: %s, default: %s, sys: %s, bin: %d, Set failed point to -1", obsName.c_str(), DefaultFile.c_str(), SysFile.c_str(), ibin);
            continue;
        }
        double defaultBinError = hDefault->GetBinError(ibin);
        double sysBinError = hSys->GetBinError(ibin);

        double uncertainty = 100 * (defaultBinContent - sysBinContent) / defaultBinContent;
        double rho = 0;
        double error = Error_PercentSys(defaultBinContent, defaultBinError, sysBinContent, sysBinError, rho);
        hUncertainty->SetBinContent(ibin, uncertainty);
        hUncertainty->SetBinError(ibin, error);

        double barlow = std::fabs(defaultBinContent - sysBinContent) / sqrt(defaultBinError*defaultBinError + sysBinError*sysBinError - 2*rho*defaultBinError*sysBinError);
        hBarlow->SetBinContent(ibin, barlow);
    }

    // Save the systematics uncertainties and barlow check in root file
    string outputFileName = "./Systematics/Uncertainties/Uncert_" + obsName;
    if (pTDiffCentMin >= 0 && pTDiffCentMax >= 0)
        outputFileName += "Cent" + to_string(pTDiffCentMin) + "To" + to_string(pTDiffCentMax);
    outputFileName += "_" + DefaultFileNameSuffix + "_VS_" + SysFileNameSuffix + ".root";
    TFile* fOutput = TFile::Open(outputFileName.c_str(), "RECREATE");
    if (!fOutput || fOutput->IsZombie()) {
        cout << "Error: cannot open file " << outputFileName << endl;
        return;
    }
    hUncertainty->Write();
    hBarlow->Write();

    // Close all the files
    fDefault->Close();
    fSys->Close();
    fOutput->Close();

}


void BatchCalculateSystematics(string& DefaultFileNameSuffix, vector<vector<string>>& SysFileNameSuffixs, vector<vector<string>>& SysDescription) {
    // Calculate relative uncertainties for each observable
    // and barlow check

    for (int kobs=0; kobs<kNObservable; kobs++) {
        if (!IfCheckObservable[kobs]) continue; // skip the observables not to be checked
        for(auto& obsName : ObservableNamesMap[kobs])
        for(int kSys=0; kSys<SysFileNameSuffixs.size(); kSys++) 
        for(int kSubSys=0; kSubSys<SysFileNameSuffixs[kSys].size(); kSubSys++){
            if (kobs<kpTDiffv2) {
                CalculateSystematic(kobs, obsName, DefaultFileNameSuffix, SysFileNameSuffixs[kSys][kSubSys], SysDescription[kSys][kSubSys]);
            }
            else {
                for (uint j=0; j<pTDiffCent.size()-1; j++) {
                    CalculateSystematic(kobs, obsName, DefaultFileNameSuffix, SysFileNameSuffixs[kSys][kSubSys], SysDescription[kSys][kSubSys], (int)pTDiffCent[j], (int)pTDiffCent[j+1]);
                }
            }
        }
    }

}

void ProcessGroupSystematics(vector<double>& tabUncertAndDecision, vector<int>& barlowDetails, int kobs, string obsName, string DefaultFileNameSuffix, vector<string>& SysFileNameSuffixs, vector<string>& SysDescription, int pTDiffCentMin=-1, int pTDiffCentMax=-1) {
    
    vector<double> Uncert(SysFileNameSuffixs.size());
    vector<double> Barlow(SysFileNameSuffixs.size());
    vector<double> barlowNumSave(SysFileNameSuffixs.size());
    vector<double> binCountSave(SysFileNameSuffixs.size());
    for (int kSys=0; kSys<SysFileNameSuffixs.size(); kSys++) {
        string sysCalFileName = "./Systematics/Uncertainties/Uncert_" + obsName;
        if (pTDiffCentMin >= 0 && pTDiffCentMax >= 0)
            sysCalFileName += "Cent" + to_string(pTDiffCentMin) + "To" + to_string(pTDiffCentMax);
        sysCalFileName += "_" + DefaultFileNameSuffix + "_VS_" + SysFileNameSuffixs[kSys] + ".root";
        TFile* fSysCal = TFile::Open(sysCalFileName.c_str(), "READ");
        if (!fSysCal || fSysCal->IsZombie()) {
            cout << "Error: cannot open file " << sysCalFileName << endl;
            return;
        }
        TH1D* hUncertainty = (TH1D*)fSysCal->Get("hUncertainty");
        if (!hUncertainty) {
            cout << "Error: cannot find histogram hUncertainty in file " << sysCalFileName << endl;
            return;
        }
        TH1D* hBarlow = (TH1D*)fSysCal->Get("hBarlow");
        if (!hBarlow) {
            cout << "Error: cannot find histogram hBarlow in file " << sysCalFileName << endl;
            return;
        }

        double FitRange = hBarlow->GetXaxis()->GetXmax();
        if (kobs>=kpTDiffv2) FitRange = MaxpT;
        TF1 *f0 = new TF1(Form("f0_%s", obsName.c_str()), "[0]", 0, FitRange);
        hUncertainty->Fit(f0, "R");
        double uncertainty = std::fabs(f0->GetParameter(0));
        Uncert[kSys] = uncertainty;
        
        int barlowNum = 0;
        int binCount = 0;
        for(int iBarlow=1; iBarlow<=hBarlow->GetNbinsX(); iBarlow++){
            if (kobs>=kpTDiffv2 && (hBarlow->GetBinCenter(iBarlow)>MaxpT)) 
                continue;
            double binContent = hBarlow->GetBinContent(iBarlow);
            if (binContent > BarlowRequirement) {
                barlowNum++;
            }
            binCount++;
        }

        if (binCount > 0)
            Barlow[kSys] = (double)barlowNum / (double)binCount;

        barlowNumSave[kSys] = barlowNum;
        binCountSave[kSys] = binCount;
        
    }

    tabUncertAndDecision[0] = Uncert[0];
    tabUncertAndDecision[1] = Barlow[0];
    barlowDetails[0] = barlowNumSave[0];
    barlowDetails[1] = binCountSave[0];
    for (int kSys=1; kSys<SysFileNameSuffixs.size(); kSys++) {
        if (Barlow[kSys] > tabUncertAndDecision[1]) {
            tabUncertAndDecision[0] = Uncert[kSys];
            tabUncertAndDecision[1] = Barlow[kSys];
            barlowDetails[0] = barlowNumSave[kSys];
            barlowDetails[1] = binCountSave[kSys];
        }
    }

    // Make decision based on Barlow check
    if (tabUncertAndDecision[1]<0.5) {
        tabUncertAndDecision[1] = 0;
    }
    else {
        tabUncertAndDecision[1] = 1;
    }
}

void ProcessSingleObsSystematics(bool kSummaryTable, bool kRootFile, int kobs, int kobsName, string DefaultFileNameSuffix, vector<vector<string>>& SysFileNameSuffixs, vector<vector<string>>& SysDescription, int pTDiffCentMin=-1, int pTDiffCentMax=-1) {
    vector<vector<double>> tabUncertAndDecision;
    tabUncertAndDecision.resize(SysFileNameSuffixs.size());
    for(auto& row: tabUncertAndDecision) row.resize(2);

    vector<vector<int>> barlowDetails(SysFileNameSuffixs.size());
    for(auto& row: barlowDetails) row.resize(2);

    for(int kSys=0; kSys<SysFileNameSuffixs.size(); kSys++) {
        ProcessGroupSystematics(tabUncertAndDecision[kSys], barlowDetails[kSys], kobs, ObservableNamesMap[kobs][kobsName], DefaultFileNameSuffix, SysFileNameSuffixs[kSys], SysDescription[kSys], pTDiffCentMin, pTDiffCentMax);
    }

    double SysTotalError = 0.;
    for(int kSys=0; kSys<SysFileNameSuffixs.size(); kSys++){
        if (tabUncertAndDecision[kSys][1] > 0)
            SysTotalError += tabUncertAndDecision[kSys][0] * tabUncertAndDecision[kSys][0];
    }
    SysTotalError = sqrt(SysTotalError);

    if (kSummaryTable) {
        string tableName = "./Systematics/SummaryTable/SysTable_" + ObservableOutputNamesMap[kobs][kobsName];
        if (pTDiffCentMin >= 0 && pTDiffCentMax >= 0)
            tableName += "Cent" + to_string(pTDiffCentMin) + "To" + to_string(pTDiffCentMax);
        tableName += ".tex";
        ofstream ofs(tableName.c_str());
        ofs << "\\begin{table}[htbp]" << endl;
        ofs << "\\caption{" << "$" << ObservablePrintNamesMap[kobs][kobsName] << "$";
        if (pTDiffCentMin >= 0 && pTDiffCentMax >= 0)
            ofs << " (centrality: " << pTDiffCentMin << "-" << pTDiffCentMax << "\\%)";
        ofs << " Systematics}" << endl;
        ofs << "\\label{tab:Sys_" << ObservableOutputNamesMap[kobs][kobsName] << "}" << endl;
        ofs << "\\centering" << endl;
        ofs << "\\begin{tabular}{|c|c|c|c|}" << endl;
        ofs << "\\hline" << endl;
        ofs << "Systematic & Uncertainty(\\%) & Barlow Check & Decision \\\\" << endl;
        ofs << "\\hline" << endl;
        for(int kSys=0; kSys<SysFileNameSuffixs.size(); kSys++){
            ofs << SysDescription[kSys][0] << " & " << std::fixed << std::setprecision(2) << tabUncertAndDecision[kSys][0] << " & ";
            ofs << std::defaultfloat << barlowDetails[kSys][0] << "/" << barlowDetails[kSys][1] << "=" << (double)barlowDetails[kSys][0]/(double)barlowDetails[kSys][1] << " & " << tabUncertAndDecision[kSys][1] << " \\\\" << endl;
        }
        ofs << "\\hline" << endl;
        ofs << "Total Uncertainty & " << std::fixed << std::setprecision(2) << SysTotalError << " & " << "N/A" << " & " << "N/A" << " \\\\" << endl;
        ofs << "\\hline" << endl;
        ofs << "\\end{tabular}" << endl;
        ofs << "\\end{table}" << endl;

        ofs.close();
    }

    if (kRootFile) {
         string DefaultFile = "./ProcessOutput/" + ObservableFilesMap[kobs];
        if (pTDiffCentMin >= 0 && pTDiffCentMax >= 0)
            DefaultFile += "Cent" + to_string(pTDiffCentMin) + "To" + to_string(pTDiffCentMax);
        DefaultFile += "_" + DefaultFileNameSuffix + ".root";
        TFile* fDefault = TFile::Open(DefaultFile.c_str(), "READ");
        if (!fDefault || fDefault->IsZombie()) {
            cout << "Error: cannot open file " << DefaultFile << endl;
            return;
        }
        TH1D* hDefault = (TH1D*)fDefault->Get(ObservableNamesMap[kobs][kobsName].c_str());
        if (!hDefault) {
            cout << "Error: cannot find histogram " << ObservableNamesMap[kobs][kobsName] << " in file " << DefaultFile << endl;
            return;
        }

        TGraphErrors* sysBox = new TGraphErrors();
        sysBox->SetName(Form("Sys_%s", ObservableOutputNamesMap[kobs][kobsName].c_str()));
        int nbins = hDefault->GetNbinsX();
        for (int ibin=1; ibin<=nbins; ibin++) {
            double defaultBinContent = hDefault->GetBinContent(ibin);
            sysBox->SetPoint(ibin-1, hDefault->GetBinCenter(ibin), defaultBinContent);

            double xerror = 1.0;
            if (kobs>=kpTDiffv2 && ibin<nbins) {
                xerror = 0.5*(hDefault->GetBinCenter(ibin+1) - hDefault->GetBinCenter(ibin));
            }

            sysBox->SetPointError(ibin-1, xerror, defaultBinContent*SysTotalError*0.01);
            
        }
        sysBox->SetFillStyle(1001);
        sysBox->SetFillColorAlpha(kBlack,0.5);

        string outputFileName = "./Systematics/TotalSystematics/Sys_" + ObservableOutputNamesMap[kobs][kobsName];
        if (pTDiffCentMin >= 0 && pTDiffCentMax >= 0)
            outputFileName += "Cent" + to_string(pTDiffCentMin) + "To" + to_string(pTDiffCentMax);
        outputFileName += "_" + DefaultFileNameSuffix + ".root";
        TFile* fOutput = TFile::Open(outputFileName.c_str(), "RECREATE");
        if (!fOutput || fOutput->IsZombie()) {
            cout << "Error: cannot open file " << outputFileName << endl;
            return;
        }
        sysBox->Write();
        fOutput->Close();
    }
}

void ProcessTotalSystematics(string& DefaultFileNameSuffix, vector<vector<string>>& SysFileNameSuffixs, vector<vector<string>>& SysDescription, bool kSummaryTable, bool kRootFile) {
    // kSummaryTable: print summary table of systematics in Latex
    // kRootFile: save systematics uncertainties in root file

    for (int kobs=0; kobs<kNObservable; kobs++) {
        if (!IfCheckObservable[kobs]) continue; // skip the observables not to be checked
        for(int kobsName=0; kobsName<ObservableNamesMap[kobs].size(); kobsName++){
            if (kobs<kpTDiffv2) {
                ProcessSingleObsSystematics(kSummaryTable, kRootFile, kobs, kobsName, DefaultFileNameSuffix, SysFileNameSuffixs, SysDescription);
            }
            else {
                for (uint j=0; j<pTDiffCent.size()-1; j++) {
                    ProcessSingleObsSystematics(kSummaryTable, kRootFile, kobs, kobsName, DefaultFileNameSuffix, SysFileNameSuffixs, SysDescription, (int)pTDiffCent[j], (int)pTDiffCent[j+1]);
                }
            }
        }
        
    }

}

void ProcessSystematics() {
    string DefaultFileNameSuffix = "LHC23_PbPb_pass4_344339";
    vector<vector<string>> SysFileNameSuffixs;
    vector<vector<string>> SysDescription;

    // Event Systematics
    SysFileNameSuffixs.push_back({"LHC23zzh_pass4_small_340440_FT0M"});
    SysDescription.push_back({"Cent Estimator: FT0M"});
    SysFileNameSuffixs.push_back({"LHC23_PbPb_pass4_341269_kColl"});
    SysDescription.push_back({"kNoColl flags"});
    SysFileNameSuffixs.push_back({"LHC23zzh_pass4_small_346469_NoMultCorrelation"});
    SysDescription.push_back({"w/o Multiplicity correlation"});
    SysFileNameSuffixs.push_back({"LHC23zzh_pass4_359613_ITSclu0"});
    SysDescription.push_back({"\\# of ITS clusters $>$ 0"});
    SysFileNameSuffixs.push_back({"LHC23zzh_pass4_small_356028_occupancy5k","LHC23zzh_pass4_359613_occupancy1w"});
    SysDescription.push_back({"occupancy $<=$ 5k", "occupancy $<=$ 10k"});

    // Track Systematics
    SysFileNameSuffixs.push_back({"LHC23_PbPb_pass4_341269"});
    SysDescription.push_back({"w/o NUA correction"});
    SysFileNameSuffixs.push_back({"LHC23zzh_pass4_small_346469_Chi2"});
    SysDescription.push_back({"TPC Chi2 $<$ 4."});
    SysFileNameSuffixs.push_back({"LHC23zzh_pass4_small_346469_Vtxz8"});
    SysDescription.push_back({"Vertex Z $<$ 8 cm"});
    SysFileNameSuffixs.push_back({"LHC23zzh_pass4_small_346469_TPCclu90"});
    SysDescription.push_back({"\\# of TPC clusters $>$ 90"});
    SysFileNameSuffixs.push_back({"LHC23zzh_pass4_small_345386_DCAzPt_id24178"});
    SysDescription.push_back({"DCAz Pt dependence cut"});
    
    // BatchCalculateSystematics(DefaultFileNameSuffix, SysFileNameSuffixs, SysDescription);
    ProcessTotalSystematics(DefaultFileNameSuffix, SysFileNameSuffixs, SysDescription, true, true);
}
