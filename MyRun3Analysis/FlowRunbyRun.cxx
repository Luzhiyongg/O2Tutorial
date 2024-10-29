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
//
// code author: Zhiyong Lu (zhiyong.lu@cern.ch)
// jira: PWGCF-254
// Produce Run-by-Run QA plots and flow analysis for Run3

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include "GFWPowerArray.h"
#include "GFW.h"
#include "GFWCumulant.h"
#include "GFWWeights.h"
#include "TList.h"
#include <TProfile.h>
#include <TRandom3.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowRunbyRun {

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtRefMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for all tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 10.0f, "Maximal pT for all tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  Configurable<std::vector<uint>> cfgRunNumber{"cfgRunNumber", std::vector<uint>{544095,544098,544116,544121,544122,544123,544124}, "Preconfigured run numbers"};
  

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10}, "pt axis for histograms"};
  ConfigurableAxis axisIndependent{"axisIndependent", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "X axis for histograms"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) 
  && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) 
  && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // Define output
  HistogramRegistry registry{"registry"};

  // define global variables
  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;

  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults>>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;

  void init(InitContext const&)
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(nolaterthan.value);

    //Add some output objects to the histogram registry
    registry.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
    registry.add("hEta", "", {HistType::kTH1D, {axisEta}});
    registry.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    registry.add("hMult", "", {HistType::kTH1D, {{3000,0.5,3000.5}}});
    registry.add("hCent", "", {HistType::kTH1D, {{90,0,90}}});
    registry.add("c22", "", {HistType::kTProfile, {axisIndependent}});

    fGFW->AddRegion("full", -0.8, 0.8, 1, 1);
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    fGFW->CreateRegions();
  }

  template<char... chars>
  void FillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf,0,kTRUE).real();
    if(dnx==0) return;
    if(!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf,0,kFALSE).real()/dnx;
      if(TMath::Abs(val)<1)
        registry.fill(tarName,cent,val,dnx);
      return;
    };
    return;
  }

  void process(aodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, aodTracks const& tracks)
  {
    int Ntot = tracks.size();
    if (Ntot < 1)
      return;
    float vtxz = collision.posZ();
    registry.fill(HIST("hVtxZ"), vtxz);
    registry.fill(HIST("hMult"),Ntot);
    registry.fill(HIST("hCent"),collision.centFT0C());
    fGFW->Clear();
    const auto cent = collision.centFT0C();
    float weff = 1, wacc = 1;
    for (auto& track : tracks) {
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hEta"), track.eta());
      
      fGFW->Fill(track.eta(), 1, track.phi(), wacc * weff, 1);
    }

    //Filling c22 with ROOT TProfile
    FillProfile(corrconfigs.at(0), HIST("c22"), cent);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowRunbyRun>(cfgc)};
}
