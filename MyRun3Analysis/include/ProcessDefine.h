#ifndef PROCESS_DEFINE_H
#define PROCESS_DEFINE_H

#include <cstring>
#include <vector>
#include <map>
#include <array>

std::vector<double> pTDiffCent={0,5,10,20,30,40,50,60,70};

double BarlowRequirement = 1.;
double MaxpT = 200.;

enum kObservable{
    kVn,
    kV24,
    kV26,
    kV28,
    kV210,
    kV422,
    kChi422,
    kRho422,
    kNSC23,
    kNSC24,
    kNSC234,
    kNSC345,
    kpTDiffv2,
    kpTDiffv3,
    kpTDiffv4,
    kpTDiffv24ChFull,
    kpTDiffv24Ch10Gap,
    kpTDiffv26ChFull,
    kNObservable
};

std::map<int, std::string> ObservableFilesMap = {
    {kVn, "vn"},
    {kV24, "v24"},
    {kV26, "v26"},
    {kV28, "v28"},
    {kV210, "v210"},
    {kV422, "v422"},
    {kChi422, "chi422"},
    {kRho422, "rho422"},
    {kNSC23, "NSC23"},
    {kNSC24, "NSC24"},
    {kNSC234, "NSC234"},
    {kNSC345, "NSC345"},
    {kpTDiffv2, "pTDiffv2"},
    {kpTDiffv3, "pTDiffv3"},
    {kpTDiffv4, "pTDiffv4"},
    {kpTDiffv24ChFull, "pTDiffv24ChFull"},
    {kpTDiffv24Ch10Gap, "pTDiffv24Ch10Gap"},
    {kpTDiffv26ChFull, "pTDiffv26ChFull"}
};

map<int, vector<string>> ObservableNamesMap = {
    {kVn, {"Corr_corr_22_hist","Corr_corr_32_hist","Corr_corr_42_hist"}},
    {kV24, {"Corr_corr_24_hist"}},
    {kV26, {"corr_26_hist"}},
    {kV28, {"corr_28_hist"}},
    {kV210, {"cN10_corr_210_hist"}},
    {kV422, {"hCorr422_mean"}},
    {kChi422, {"hCorr422_mean"}},
    {kRho422, {"hCorr422_mean"}},
    {kNSC23, {"NSC23"}},
    {kNSC24, {"NSC24"}},
    {kNSC234, {"NSC234"}},
    {kNSC345, {"NSC345"}},
    {kpTDiffv2, {"pTDiffv2"}},
    {kpTDiffv3, {"pTDiffv3"}},
    {kpTDiffv4, {"pTDiffv4"}},
    {kpTDiffv24ChFull, {"pTDiffv24ChFull"}},
    {kpTDiffv24Ch10Gap, {"pTDiffv24Ch10Gap"}},
    {kpTDiffv26ChFull, {"pTDiffv26ChFull"}}
};

map<int, vector<string>> ObservableOutputNamesMap = {
    {kVn, {"v22","v32","v42"}},
    {kV24, {"v24"}},
    {kV26, {"v26"}},
    {kV28, {"v28"}},
    {kV210, {"v210"}},
    {kV422, {"v422"}},
    {kChi422, {"chi422"}},
    {kRho422, {"rho422"}},
    {kNSC23, {"NSC23"}},
    {kNSC24, {"NSC24"}},
    {kNSC234, {"NSC234"}},
    {kNSC345, {"NSC345"}},
    {kpTDiffv2, {"pTDiffv2"}},
    {kpTDiffv3, {"pTDiffv3"}},
    {kpTDiffv4, {"pTDiffv4"}},
    {kpTDiffv24ChFull, {"pTDiffv24ChFull"}},
    {kpTDiffv24Ch10Gap, {"pTDiffv24Ch10Gap"}},
    {kpTDiffv26ChFull, {"pTDiffv26ChFull"}}
};

#endif // PROCESS_DEFINE_H
