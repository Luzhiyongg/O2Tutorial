#ifndef RUNBYRUNCOMMON_H
#define RUNBYRUNCOMMON_H

#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"

std::vector<Color_t> colorsList = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kOrange, kViolet, kPink};

std::vector<EMarkerStyle> markerStylesList = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullDiamond, kFullCross, kFullCrossX, kFullThreeTriangles,
 kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown, kOpenStar, kOpenDiamond, kOpenCross, kOpenCrossX, kOpenThreeTriangles};

std::vector<ELineStyle> lineStylesList = {kSolid, kDashed, kDotted};

std::vector<int> IgnoreRuns = {544640, 544913, 545117, 545295, 545311, 544091, 544652, 544694}; 
// std::vector<int> IgnoreRuns = {544640, 544913, 545117, 545295, 545311, 544091, 544652, 544694, 545296, 544964}; 
// std::vector<int> IgnoreRuns = {544116,544122,544653,544742}; // LHCzzh
// std::vector<int> IgnoreRuns = {-1}; 

std::vector<int> SelectedRuns = {}; 
// std::vector<int> SelectedRuns = {545367, 545291, 545223, 544917}; // Good ITS Runs
// std::vector<int> SelectedRuns = {544095,544098,544116,544121,544122,544123,544124}; //LHCzzh_pass4
// std::vector<int> SelectedRuns = {544614,544652,544672,544693,544694,544696,544739,544742,544754,544767}; //LHC23zzm_pass4_calo
// std::vector<int> SelectedRuns = {544868,544886,544887,544896,544913,544914,544917,544961,544964,544968,544992,545008,545009,545041,545042,545044,545047,545062,545063,545064,545066,545117,545171,545185,545210,545222,545223,545246,545249,5452625,45291,545294,545295,545296,545311,545312,545332,545367}; //LHC23zzo_apass4

bool cutHadronicRate = true;
double maxHadronicRate = 50; // in kHz
double minHadronicRate = 40; // in kHz

double GetHadronicRate(int runNumber, TDirectory* HadronicRate) {
    if (!(runNumber > 0)) {
        std::cout << "Invalid run number: " << runNumber << std::endl;
        return 0;
    }
    int nRunsHadronicRate = HadronicRate->GetListOfKeys()->GetEntries(); // number of runs in HadronicRate directory
    TH2D* h_HadronicRate = (TH2D*)HadronicRate->Get(Form("%d", runNumber));
    if (!h_HadronicRate) {
        std::cout << "Cannot find h_HadronicRate for run " << runNumber << std::endl;
        return 0;
    }
    double meanHadronicRate = 0;
    double count = 0;
    TH1D* h_HadronicRateY = h_HadronicRate->ProjectionY("h_HadronicRateY");
    for (int iBin = 1; iBin <= h_HadronicRateY->GetNbinsX(); iBin++) {
        double binContent = h_HadronicRateY->GetBinContent(iBin);
        if (binContent > 0) {
            meanHadronicRate += h_HadronicRateY->GetBinCenter(iBin) * binContent;
            count += binContent;
        }
    }
    if (count > 0) {
        meanHadronicRate /= count;
    }
    // Printf("Run %d: meanHadronicRate = %.2f", runNumber, meanHadronicRate);
    return meanHadronicRate;
}

#endif // RUNBYRUNCOMMON_H
