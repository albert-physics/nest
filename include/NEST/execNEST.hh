#ifndef __EXECNEST_H__
#define __EXECNEST_H__ 1

#include "NEST.hh"
#include "NESTObservableArray.hh"

using namespace std;
using namespace NEST;


vector<vector<double>> GetBand(vector<double> S1s, vector<double> S2s,
                               bool resol, int nFold);

void GetEnergyRes(vector<double> Es);

int execNEST(VDetector* detector, double numEvts, const string& type,
             double eMin, double eMax, double inField, string position,
             const string& posiMuon, double fPos, int seed, bool no_seed,
             double dayNumber);

NESTObservableArray runNESTvec(VDetector* detector,
                               INTERACTION_TYPE scatterType,
                               std::vector<double> eList,
                               std::vector<std::vector<double>> pos3dxyz,
                               double inField = -1.0, uint64_t seed = 0,   
                               std::vector<double> ERYieldsParam = default_ERYieldsParam,
                               std::vector<double> NRYieldsParam = default_NRYieldsParam,
                               std::vector<double> NRERWidthsParam = default_NRERWidthsParam,
                               S1CalculationMode s1mode = NEST::S1CalculationMode::Full,
                               S2CalculationMode s2mode = NEST::S2CalculationMode::Full,
                               bool calculate_times = false);
#endif
