#ifndef BASICSETTING_H
#define BASICSETTING_H

bool ComparewithPublish = false;
bool OutputRoot = true;
bool RebinpTDiff = false;
bool ApplynonclosureCorrection = false;
double epsilon = 1e-10;
std::vector<Double_t> pTDiffOriginBinning = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.5, 4, 5, 6, 8, 10};
std::vector<Double_t> pTDiffTargetBinning = {0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.5, 1.7, 2.0, 2.4, 3.0, 3.5, 4.0, 5.0};

vector<float> targetCentralityBins = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};

enum nonlinearObservableEnum{
    v422,
    chi422,
    rho422
};

std::map<nonlinearObservableEnum, Char_t*> nonlinearObservableName = {
  {v422, "v_{4,22}"},
  {chi422, "#chi_{4,22}"},
  {rho422, "#rho_{4,22}"}
};

std::map<nonlinearObservableEnum, Char_t*> nonlinearOutputObservableName = {
  {v422, "v422"},
  {chi422, "chi422"},
  {rho422, "rho422"}
};

std::map<nonlinearObservableEnum, std::array<double, 2>> nonlinearUserRangeMap = {
  {v422, {-0.006,0.015}},
  {chi422, {-5.5,2.}},
  {rho422, {-0.5,1.}}
};

bool isZero(double x) {
  return (x < epsilon && x > -epsilon);
}

#endif // BASICSETTING_H
