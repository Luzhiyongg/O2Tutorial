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
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>, CERN

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DataModel/DerivedExampleTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h" // For Filter

struct DerivedBasicConsumer {
  /// Function to aid in calculating delta-phi
  /// \param phi1 first phi value
  /// \param phi2 second phi value
  Double_t ComputeDeltaPhi(Double_t phi1, Double_t phi2)
  {
    Double_t deltaPhi = phi1 - phi2;
    if (deltaPhi < -TMath::Pi() / 2.) {
      deltaPhi += 2. * TMath::Pi();
    }
    if (deltaPhi > 3 * TMath::Pi() / 2.) {
      deltaPhi -= 2. * TMath::Pi();
    }
    return deltaPhi;
  }

  SliceCache cache;

  //Filter
  Filter collZfilter = nabs(aod::collision::posZ)<10.0f;

  //Partition
  Partition<aod::DrTracks> associatedTracks = aod::exampleTrackSpace::pt < 6.0f && aod::exampleTrackSpace::pt > 4.0f;
  Partition<aod::DrTracks> triggerTracks = aod::exampleTrackSpace::pt > 6.0f;

  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // define axes you want to use
    const AxisSpec axisCounter{1, 0, +1, ""};
    const AxisSpec axisVertexZ{300, -15, +15, "z (cm)"};
    const AxisSpec axisPt{100, 0., +10., "p_{T} (GeV/c)"};
    const AxisSpec axisDeltaPhi{100, -0.5*TMath::Pi(), +1.5*TMath::Pi(), "#Delta#phi"};
    const AxisSpec axisDeltaEta{100, -1.0, +1.0, "#Delta#eta"};

    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("VertexZ", "VertexZ", kTH1F, {axisVertexZ});
    histos.add("ptAssoHistogram", "ptAssoHistogram", kTH1F, {axisPt});
    histos.add("ptTrigHistogram", "ptTrigHistogram", kTH1F, {axisPt});
    histos.add("ptHistogram", "ptHistogram", kTH1F, {axisPt});

  }

  void process(soa::Filtered<aod::DrCollisions>::iterator const& collision, aod::DrTracks const& tracks)
  {
    histos.fill(HIST("eventCounter"), 0.5);
    histos.fill(HIST("VertexZ"), collision.posZ());

    for (auto& track : tracks){
      histos.fill(HIST("ptHistogram"), track.pt());
    }

    auto assoTracksThisCollision = associatedTracks->sliceByCached(aod::exampleTrackSpace::drCollisionId, collision.globalIndex(), cache);
    auto trigTracksThisCollision = triggerTracks->sliceByCached(aod::exampleTrackSpace::drCollisionId, collision.globalIndex(), cache);
    for (auto& track : assoTracksThisCollision)
      histos.fill(HIST("ptAssoHistogram"), track.pt());
    for (auto& track : trigTracksThisCollision)
      histos.fill(HIST("ptTrigHistogram"), track.pt());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<DerivedBasicConsumer>(cfgc, TaskName{"derived-basic-consumer"})};
  return workflow;
}
