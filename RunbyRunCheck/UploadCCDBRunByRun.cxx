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
#include "TDirectory.h"
#include <vector>


void UploadCCDBRunByRun()
{
    std::string FileNameSuffix = "LHC23_PbPb_pass4_356235";
    TFile *file = TFile::Open(Form("./AnalysisResults/AnalysisResults_%s.root", FileNameSuffix.c_str()), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "Cannot open file" << std::endl;
        return;
    }
    file->Cd("flow-runby-run");
    // get the list of runs
    TDirectory* RunTopDir = gDirectory;
    TList *runsList = RunTopDir->GetListOfKeys();
    int nObj = runsList->GetEntries(); // number of runs
    vector<int> IgnoreRuns = {-1}; 
    std::vector<double> AxisPt = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.5,5,5.5,6,7,8,9,10,15,20,25,30,50,70,100,200,300,400,500,600,700,800,900,1000};
    
    // init ccdb api
    o2::ccdb::CcdbApi ccdb_api;
    ccdb_api.init("https://alice-ccdb.cern.ch");
    std::map<std::string, std::string> metadataRCT, headers;

    for (int iObj = 0; iObj < nObj; iObj++){
        TString runName = runsList->At(iObj)->GetName();
        int runNumber = runName.Atoi();
        if (!(runNumber > 0)) {
            Printf("\"%s\" is not a valid run number", runName.Data());
            continue;
        }
        if (find(IgnoreRuns.begin(), IgnoreRuns.end(), runNumber) != IgnoreRuns.end()) {
            Printf("Skipping run %d", runNumber);
            continue;
        }

        TH3D* hPhiEtaVtxz = (TH3D*)RunTopDir->Get(Form("%d/hPhiEtaVtxz",runNumber));
        if (!hPhiEtaVtxz) {
            std::cout << "Error: cannot find " << runNumber << "/hPhiEtaVtxz." << std::endl;
            return;
        }
        // Printf("Get %d/hPhiEtaVtxz", runNumber);

        GFWWeights* fWeights = new GFWWeights(Form("%d",runNumber));
        fWeights->setPtBins(AxisPt.size() - 1, AxisPt.data());
        fWeights->init(false, false);
        fWeights->setTH3D(hPhiEtaVtxz);
        fWeights->setDataFilled(true);
        Printf("Set GFWWeight for run %d", runNumber);

        // Get RCT run information
        headers = ccdb_api.retrieveHeaders(Form("RCT/Info/RunInformation/%i", runNumber), metadataRCT, -1);
        int64_t tsSOR = atol(headers["SOR"].c_str());
        int64_t tsEOR = atol(headers["EOR"].c_str());
        cout << "Run " << runNumber << " SOR " << tsSOR << " EOR " << tsEOR << endl;

        // cout << "Defining metadata for this run..." << endl;
        std::map<std::string, std::string> metadata; // can be empty
        metadata.insert(std::pair{"Description", Form("NUA correcrion for run %s", Form("%i", runNumber))});
        metadata.insert(std::pair{"Author", "Zhiyong Lu"});

        cout << "Attempting CCDB upload..." << endl;
        try
        {
            ccdb_api.storeAsTFileAny(fWeights, Form("Users/z/zhlu/AcceptanceRunByRun/NUA_%s", FileNameSuffix.c_str()), metadata, tsSOR, tsEOR);
        }
        catch (std::exception const &e)
        {
            LOG(fatal) << "Failed at CCDB submission!";
        }
        cout << "Finished with upload of run " << runNumber << " update! " << endl << endl;
    }
    cout << "Done!" << endl;
}
