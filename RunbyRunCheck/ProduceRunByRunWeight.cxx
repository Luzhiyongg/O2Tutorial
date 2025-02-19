//put in the first lines to ignore the warning message
#pragma GCC diagnostic ignored "-Winconsistent-missing-override"

// #include "FlowContainer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "GFWWeights.h"
#include "TObjArray.h"
#include <vector>

void ProduceRunByRunWeight() {
    TFile* f = new TFile("./AnalysisResults/AnalysisResults_LHC23zzh_pass4_small_352370.root","READ");
    if (!f || f->IsZombie()) {
        std::cout << "Error: cannot open file." << std::endl;
        return;
    }
    gDirectory->Cd("flow-runby-run");
    // get the list of runs
    TList *runsList = gDirectory->GetListOfKeys();
    int nObj = runsList->GetEntries(); // number of runs
    std::vector<double> AxisPt = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5,5.5,6,7,8,9,10,15,20,25,30,50,70,100,200,300,400,500,600,700,800,900,1000};

    TObjArray* listOfWeight = new TObjArray();
    listOfWeight->SetOwner(kTRUE);

    for (int iObj = 0; iObj < nObj; iObj++){
      TString runName = runsList->At(iObj)->GetName();
      int runNumber = runName.Atoi();
      if (!(runNumber>0))
        continue;
      
      TH3D* hPhiEtaVtxz = (TH3D*)gDirectory->Get(Form("%d/hPhiEtaVtxz",runNumber));
      if (!hPhiEtaVtxz) {
        std::cout << "Error: cannot find " << runNumber << "/hPhiEtaVtxz." << std::endl;
        return;
      }
      // hPhiEtaVtxz->Draw("colz");

      GFWWeights* fWeights = new GFWWeights(Form("%d",runNumber));
      fWeights->setPtBins(AxisPt.size() - 1, AxisPt.data());
      fWeights->init(false, false);
      fWeights->setTH3D(hPhiEtaVtxz);
      fWeights->setDataFilled(true);

      listOfWeight->Add(fWeights);

    }
    
    TFile* fout = new TFile("./AcceptanceRunByRun/NUA_LHC23zzh_pass4_small_352370.root","RECREATE");
    listOfWeight->Write("ccdb_object", TObject::kSingleKey);
    fout->Close();
    f->Close();

}
