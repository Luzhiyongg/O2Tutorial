// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "FlowContainer.h"
#include "TList.h"
#include <TProfile.h>
#include <TRandom3.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct GfwTutorial {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")


  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPtHist{"axisPtHist", {100, 0., 10.}, "pt axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // Define output
  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);

  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>>;

  void init(InitContext const&)
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    // Add some output objects to the histogram registry
    registry.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {axisEta}});
    registry.add("hPt", "", {HistType::kTH1D, {axisPtHist}});
    registry.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    registry.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    registry.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    registry.add("hMeanPt", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("hMeanPtWithinGap08", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c22", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c32", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c42", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c24", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c26", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c22_gap04", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c22_gap06", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c22_gap08", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c22_gap10", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c32_gap04", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c32_gap06", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c32_gap08", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c32_gap10", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c42_gap04", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c42_gap06", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c42_gap08", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c42_gap10", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c422", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c422_gapA04", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c422_gapB04", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c422_gapA10", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c422_gapB10", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c3232", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c4242", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c3232_gap04", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c4242_gap04", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c24_gap04", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c3232_gap10", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c4242_gap10", "", {HistType::kTProfile, {axisMultiplicity}});
    // registry.add("c24_gap10", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22_gap08_Weff", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("c22_gap08_trackMeanPt", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("PtVariance_partA_WithinGap08", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("PtVariance_partB_WithinGap08", "", {HistType::kTProfile, {axisMultiplicity}});

    registry.add("BootstrapContainer00/hMeanPtWithinGap08", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("BootstrapContainer00/c22_gap08_Weff", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("BootstrapContainer00/c22_gap08_trackMeanPt", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("BootstrapContainer00/PtVariance_partA_WithinGap08", "", {HistType::kTProfile, {axisMultiplicity}});
    registry.add("BootstrapContainer00/PtVariance_partB_WithinGap08", "", {HistType::kTProfile, {axisMultiplicity}});

    for(int i=1;i<10;i++){
      registry.addClone("BootstrapContainer00/", Form("BootstrapContainer0%d/",i));
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size()-1;
    double* PtBins= &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins,PtBins);

    //add in FlowContainer to Get boostrap sample automatically
    TObjArray* oba = new TObjArray();
    oba->Add(new TNamed("ChGap22", "ChGap22"));  
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      oba->Add(new TNamed(Form("ChGap22_pt_%i",i+1),"ChGap22_pTDiff"));
    oba->Add(new TNamed("ChFull22", "ChFull22"));
    oba->Add(new TNamed("ChFull32", "ChFull32"));
    oba->Add(new TNamed("ChFull42", "ChFull42"));
    oba->Add(new TNamed("ChFull24", "ChFull24"));
    oba->Add(new TNamed("ChFull26", "ChFull26"));
    oba->Add(new TNamed("Ch04Gap22", "Ch04Gap22"));
    oba->Add(new TNamed("Ch06Gap22", "Ch06Gap22"));
    oba->Add(new TNamed("Ch08Gap22", "Ch08Gap22"));
    oba->Add(new TNamed("Ch10Gap22", "Ch10Gap22"));
    oba->Add(new TNamed("Ch04Gap32", "Ch04Gap32"));
    oba->Add(new TNamed("Ch06Gap32", "Ch06Gap32"));
    oba->Add(new TNamed("Ch08Gap32", "Ch08Gap32"));
    oba->Add(new TNamed("Ch10Gap32", "Ch10Gap32"));
    oba->Add(new TNamed("Ch04Gap42", "Ch04Gap42"));
    oba->Add(new TNamed("Ch06Gap42", "Ch06Gap42"));
    oba->Add(new TNamed("Ch08Gap42", "Ch08Gap42"));
    oba->Add(new TNamed("Ch10Gap42", "Ch10Gap42"));
    oba->Add(new TNamed("ChFull422", "ChFull422"));
    oba->Add(new TNamed("Ch04GapA422", "Ch04GapA422"));
    oba->Add(new TNamed("Ch04GapB422", "Ch04GapB422"));
    oba->Add(new TNamed("Ch10GapA422", "Ch10GapA422"));
    oba->Add(new TNamed("Ch10GapB422", "Ch10GapB422"));
    oba->Add(new TNamed("ChFull3232", "ChFull3232"));
    oba->Add(new TNamed("ChFull4242", "ChFull4242"));
    oba->Add(new TNamed("Ch04Gap3232", "Ch04Gap3232"));
    oba->Add(new TNamed("Ch04Gap4242", "Ch04Gap4242"));
    oba->Add(new TNamed("Ch04Gap24", "Ch04Gap24"));
    oba->Add(new TNamed("Ch10Gap3232", "Ch10Gap3232"));
    oba->Add(new TNamed("Ch10Gap4242", "Ch10Gap4242"));
    oba->Add(new TNamed("Ch10Gap24", "Ch10Gap24"));
    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(oba, axisMultiplicity, cfgNbootstrap);
    delete oba;

    //eta region
    fGFW->AddRegion("full", -0.8, 0.8, 1, 1);
    fGFW->AddRegion("refN04", -0.8, -0.2, 1, 1);//gap4 negative region
    fGFW->AddRegion("refP04", 0.2, 0.8, 1, 1);//gap4 positve region
    fGFW->AddRegion("refN06", -0.8, -0.3, 1, 1);//gap6 negative region
    fGFW->AddRegion("refP06", 0.3, 0.8, 1, 1);//gap6 positve region
    fGFW->AddRegion("refN08", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP08", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("refN10", -0.8, -0.5, 1, 1);
    fGFW->AddRegion("refP10", 0.5, 0.8, 1, 1);
    fGFW->AddRegion("refP", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("refN", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("poiN", -0.8, -0.4, 1+fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("olN", -0.8, -0.4, 1, 4);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {3 -3}", "ChFull32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {4 -4}", "ChFull42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 2 -2 -2 -2}", "ChFull26", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {2} refP04 {-2}", "Ch04Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN06 {2} refP06 {-2}", "Ch06Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Ch08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {2} refP10 {-2}", "Ch10Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {3} refP04 {-3}", "Ch04Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN06 {3} refP06 {-3}", "Ch06Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {3} refP08 {-3}", "Ch08Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {3} refP10 {-3}", "Ch10Gap32", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {4} refP04 {-4}", "Ch04Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN06 {4} refP06 {-4}", "Ch06Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {4} refP08 {-4}", "Ch08Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {4} refP10 {-4}", "Ch10Gap42", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN {2} refP {-2}", "ChGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN refN | olN {2} refP {-2}", "ChGap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {4 -2 -2}", "ChFull422", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {-2 -2} refP04 {4}", "Ch04GapA422", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {4} refP04 {-2 -2}", "Ch04GapB422", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {-2 -2} refP10 {4}", "Ch10GapA422", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {4} refP10 {-2 -2}", "Ch10GapB422", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {3 2 -3 -2}", "ChFull3232", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {4 2 -4 -2}", "ChFull4242", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {3 2} refP04 {-3 -2}", "Ch04Gap3232", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {4 2} refP04 {-4 -2}", "Ch04Gap4242", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN04 {2 2} refP04 {-2 -2}", "Ch04Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {3 2} refP10 {-3 -2}", "Ch10Gap3232", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {4 2} refP10 {-4 -2}", "Ch10Gap4242", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN10 {2 2} refP10 {-2 -2}", "Ch10Gap24", kFALSE));
    fGFW->CreateRegions();
  }

  template <char... chars>
  void FillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        registry.fill(tarName, cent, val, dnx);
      return;
    }
    return;
  }

  template <char... chars, char... chars2>
  void FillpTvnProfile(const GFW::CorrConfig& corrconf, const double& sum_pt, const double& WeffEvent, const ConstStr<chars...>& vnWeff, const ConstStr<chars2...>& vnpT, const double& cent)
  {
    double meanPt = sum_pt/WeffEvent;
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1){
        registry.fill(vnWeff, cent, val, dnx*WeffEvent);
        registry.fill(vnpT, cent, val*meanPt, dnx*WeffEvent);
      }
      return;
    }
    return;
  }

   void FillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        fFC->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      return;
    }
    for (Int_t i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (TMath::Abs(val) < 1)
        fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
    }
    return;
  }

  template <int SampleIndex>
  void FillBootstrap(const double& cent, 
  const double& ptSum_Gap08, 
  const double& weffEvent_WithinGap08, 
  const double& sum_pt_wSquare_WithinGap08, 
  const double& sum_ptSquare_wSquare_WithinGap08, 
  const double& WeffEvent_diff_WithGap08)
  {
    static constexpr std::string_view subDir[] = {"BootstrapContainer00/","BootstrapContainer01/","BootstrapContainer02/","BootstrapContainer03/","BootstrapContainer04/","BootstrapContainer05/","BootstrapContainer06/","BootstrapContainer07/","BootstrapContainer08/","BootstrapContainer09/"};

    if(weffEvent_WithinGap08>1e-6)registry.fill(HIST(subDir[SampleIndex])+HIST("hMeanPtWithinGap08"), cent, ptSum_Gap08/weffEvent_WithinGap08, weffEvent_WithinGap08);
    if(weffEvent_WithinGap08>1e-6)FillpTvnProfile(corrconfigs.at(7), ptSum_Gap08, weffEvent_WithinGap08, HIST(subDir[SampleIndex])+HIST("c22_gap08_Weff"), HIST(subDir[SampleIndex])+HIST("c22_gap08_trackMeanPt"), cent);
    if(WeffEvent_diff_WithGap08>1e-6){
      registry.fill(HIST(subDir[SampleIndex])+HIST("PtVariance_partA_WithinGap08"), cent, 
        (ptSum_Gap08*ptSum_Gap08 - sum_ptSquare_wSquare_WithinGap08) / WeffEvent_diff_WithGap08, 
        WeffEvent_diff_WithGap08);
      registry.fill(HIST(subDir[SampleIndex])+HIST("PtVariance_partB_WithinGap08"), cent, 
        (weffEvent_WithinGap08*ptSum_Gap08 - sum_pt_wSquare_WithinGap08) / WeffEvent_diff_WithGap08, 
        WeffEvent_diff_WithGap08);
    }
    
  }

  void process(aodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks)
  {
    int Ntot = tracks.size();
    if (Ntot < 1)
      return;
    float l_Random = fRndm->Rndm();
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"), Ntot);
    registry.fill(HIST("hCent"),collision.centFT0C());
    fGFW->Clear();
    const auto cent = collision.centFT0C();
    float weff = 1, wacc = 1;
    double weffEvent=0, waccEvent=0;
    int TrackNum=0;
    double ptSum=0., ptSum_Gap08=0.;
    double weffEvent_WithinGap08=0., weffEventSquare_WithinGap08=0.;
    double sum_ptSquare_wSquare_WithinGap08=0., sum_pt_wSquare_WithinGap08=0.;

    for (auto& track : tracks) {
      double pt = track.pt();
      double eta = track.eta();
      bool WithinPtPOI = (cfgCutPtPOIMin<pt) && (pt<cfgCutPtPOIMax); //within POI pT range
      bool WithinPtRef  = (cfgCutPtMin<pt) && (pt<cfgCutPtMax);  //within RF pT range
      bool WithinEtaGap08 = (eta >= -0.4) && (eta <= 0.4);
      if(WithinPtRef) {
        registry.fill(HIST("hPhi"), track.phi());
        registry.fill(HIST("hEta"), track.eta());
        registry.fill(HIST("hPt"), pt);
        weffEvent+=weff;
        waccEvent+=wacc;
        ptSum+=weff*pt;
        TrackNum++;
        if(WithinEtaGap08){
          ptSum_Gap08 += weff*pt;
          sum_pt_wSquare_WithinGap08 += weff*weff*pt;
          sum_ptSquare_wSquare_WithinGap08 += weff*weff*pt*pt;
          weffEvent_WithinGap08 += weff;
          weffEventSquare_WithinGap08 += weff*weff;
        }
      }
      if(WithinPtRef) fGFW->Fill(track.eta(), fPtAxis->FindBin(pt)-1, track.phi(), wacc * weff, 1);
      if(WithinPtPOI) fGFW->Fill(track.eta(), fPtAxis->FindBin(pt)-1, track.phi(), wacc * weff, 2);
      if(WithinPtPOI && WithinPtRef) fGFW->Fill(track.eta(), fPtAxis->FindBin(pt)-1, track.phi(), wacc * weff, 4);
    }

    double WeffEvent_diff_WithGap08 = weffEvent_WithinGap08*weffEvent_WithinGap08 - weffEventSquare_WithinGap08;
    //Filling TProfile
    //MeanPt
    registry.fill(HIST("hMeanPt"), cent, ptSum/weffEvent, weffEvent);
    if(weffEvent_WithinGap08>1e-6)registry.fill(HIST("hMeanPtWithinGap08"), cent, ptSum_Gap08/weffEvent_WithinGap08, weffEvent_WithinGap08);
    //v22-Pt
    //c22_gap8 * pt_withGap8
    if(weffEvent_WithinGap08>1e-6)FillpTvnProfile(corrconfigs.at(7), ptSum_Gap08, weffEvent_WithinGap08, HIST("c22_gap08_Weff"), HIST("c22_gap08_trackMeanPt"), cent);
    //PtVariance
    if(WeffEvent_diff_WithGap08>1e-6){
      registry.fill(HIST("PtVariance_partA_WithinGap08"), cent, 
        (ptSum_Gap08*ptSum_Gap08 - sum_ptSquare_wSquare_WithinGap08) / WeffEvent_diff_WithGap08, 
        WeffEvent_diff_WithGap08);
      registry.fill(HIST("PtVariance_partB_WithinGap08"), cent, 
        (weffEvent_WithinGap08*ptSum_Gap08 - sum_pt_wSquare_WithinGap08) / WeffEvent_diff_WithGap08, 
        WeffEvent_diff_WithGap08);
    }

    //Filling Bootstrap samples
    int SampleIndex = (int)(cfgNbootstrap*l_Random);
    if(SampleIndex==0)FillBootstrap<0>(cent, ptSum_Gap08, weffEvent_WithinGap08, sum_pt_wSquare_WithinGap08, sum_ptSquare_wSquare_WithinGap08, WeffEvent_diff_WithGap08);
    else if(SampleIndex==1)FillBootstrap<1>(cent, ptSum_Gap08, weffEvent_WithinGap08, sum_pt_wSquare_WithinGap08, sum_ptSquare_wSquare_WithinGap08, WeffEvent_diff_WithGap08);
    else if(SampleIndex==2)FillBootstrap<2>(cent, ptSum_Gap08, weffEvent_WithinGap08, sum_pt_wSquare_WithinGap08, sum_ptSquare_wSquare_WithinGap08, WeffEvent_diff_WithGap08);
    else if(SampleIndex==3)FillBootstrap<3>(cent, ptSum_Gap08, weffEvent_WithinGap08, sum_pt_wSquare_WithinGap08, sum_ptSquare_wSquare_WithinGap08, WeffEvent_diff_WithGap08);
    else if(SampleIndex==4)FillBootstrap<4>(cent, ptSum_Gap08, weffEvent_WithinGap08, sum_pt_wSquare_WithinGap08, sum_ptSquare_wSquare_WithinGap08, WeffEvent_diff_WithGap08);
    else if(SampleIndex==5)FillBootstrap<5>(cent, ptSum_Gap08, weffEvent_WithinGap08, sum_pt_wSquare_WithinGap08, sum_ptSquare_wSquare_WithinGap08, WeffEvent_diff_WithGap08);
    else if(SampleIndex==6)FillBootstrap<6>(cent, ptSum_Gap08, weffEvent_WithinGap08, sum_pt_wSquare_WithinGap08, sum_ptSquare_wSquare_WithinGap08, WeffEvent_diff_WithGap08);
    else if(SampleIndex==7)FillBootstrap<7>(cent, ptSum_Gap08, weffEvent_WithinGap08, sum_pt_wSquare_WithinGap08, sum_ptSquare_wSquare_WithinGap08, WeffEvent_diff_WithGap08);
    else if(SampleIndex==8)FillBootstrap<8>(cent, ptSum_Gap08, weffEvent_WithinGap08, sum_pt_wSquare_WithinGap08, sum_ptSquare_wSquare_WithinGap08, WeffEvent_diff_WithGap08);
    else if(SampleIndex==9)FillBootstrap<9>(cent, ptSum_Gap08, weffEvent_WithinGap08, sum_pt_wSquare_WithinGap08, sum_ptSquare_wSquare_WithinGap08, WeffEvent_diff_WithGap08);
    
    
    // Filling c22 with ROOT TProfile
    // FillProfile(corrconfigs.at(0), HIST("c22"), cent);
    // FillProfile(corrconfigs.at(1), HIST("c32"), cent);
    // FillProfile(corrconfigs.at(2), HIST("c42"), cent);
    // FillProfile(corrconfigs.at(3), HIST("c24"), cent);
    // FillProfile(corrconfigs.at(4), HIST("c26"), cent);
    // FillProfile(corrconfigs.at(5), HIST("c22_gap04"), cent);
    // FillProfile(corrconfigs.at(6), HIST("c22_gap06"), cent);
    // FillProfile(corrconfigs.at(7), HIST("c22_gap08"), cent);
    // FillProfile(corrconfigs.at(8), HIST("c22_gap10"), cent);
    // FillProfile(corrconfigs.at(11), HIST("c422"), cent);
    // FillProfile(corrconfigs.at(12), HIST("c422_gapA04"), cent);
    // FillProfile(corrconfigs.at(13), HIST("c422_gapB04"), cent);
    // FillProfile(corrconfigs.at(kkk), HIST("c24_gap04"), cent);

    //Filling Flow Container
    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      FillFC(corrconfigs.at(l_ind), cent, l_Random);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<GfwTutorial>(cfgc)};
}
