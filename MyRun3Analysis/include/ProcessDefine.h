/*
 * @Author: Zhiyong Lu (zhiyong.lu@cern.ch) 
 * @Date: 2024-02-01
 * @Last Modified by:   Zhiyong Lu 
 * @Last Modified time: 2025-03-05 22:02:53 
 */
#ifndef PROCESS_DEFINE_H
#define PROCESS_DEFINE_H

#include <cstring>
#include <vector>
#include <map>
#include <array>

std::vector<double> pTDiffCent={0,5,10,20,30,40,50,60,70};

double BarlowRequirement = 1.;
double MaxpT = 200.;
bool kDrawFirstOneSystematics = true;
bool kpTDiffLogx = true;

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

map<int, vector<string>> ObservablePrintNamesMap = {
    {kVn, {"v_2\\{2\\}","v_3\\{2\\}","v_4\\{2\\}"}},
    {kV24, {"v_2\\{4\\}"}},
    {kV26, {"v_2\\{6\\}"}},
    {kV28, {"v_2\\{8\\}"}},
    {kV210, {"v_2\\{10\\}"}},
    {kV422, {"v_{4,22}"}},
    {kChi422, {"\\chi_{4,22}"}},
    {kRho422, {"\\rho_{4,22}"}},
    {kNSC23, {"NSC(2,3)"}},
    {kNSC24, {"NSC(2,4)"}},
    {kNSC234, {"NSC(2,3,4)"}},
    {kNSC345, {"NSC(3,4,5)"}},
    {kpTDiffv2, {"v_2\\{2\\}(p_{T})"}},
    {kpTDiffv3, {"v_3\\{2\\}(p_{T})"}},
    {kpTDiffv4, {"v_4\\{2\\}(p_{T})"}},
    {kpTDiffv24ChFull, {"v_2\\{4\\}(p_{T})"}},
    {kpTDiffv24Ch10Gap, {"v_2\\{4, |\\Delta\\eta|>1.0\\}(p_{T})"}},
    {kpTDiffv26ChFull, {"v_2\\{6\\}(p_{T})"}}
};

#endif // PROCESS_DEFINE_H
