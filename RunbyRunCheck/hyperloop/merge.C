#include "TString.h"
#include "TFileMerger.h"
#include "TFile.h"
#include <fstream>

void merge(TString outputfile="AnalysisResults_merge.root") {
  //注意resultList.txt最后一行留空
  std::ifstream fin("resultList.txt");
  TString output(outputfile);
  TFileMerger m;
  m.OutputFile(output);
  TString file_name;
  while (fin >> file_name) {
	TFile* f = new TFile(file_name.Data(),"READ");
	if (f->IsZombie()) continue;
	f->Close();
    m.AddFile(file_name);
  }
  m.Merge();
}
