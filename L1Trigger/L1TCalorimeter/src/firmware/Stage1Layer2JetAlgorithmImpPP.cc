///
/// \class l1t::Stage1Layer2JetAlgorithmImpHI
///
///
/// \author: R. Alex Barbieri MIT
///

// This example implements algorithm version 1 and 2.

#include "L1Trigger/L1TCalorimeter/interface/Stage1Layer2JetAlgorithmImp.h"
#include "L1Trigger/L1TCalorimeter/interface/JetFinderMethods.h"
#include "L1Trigger/L1TCalorimeter/interface/PUSubtractionMethods.h"
#include "L1Trigger/L1TCalorimeter/interface/JetCalibrationMethods.h"
#include "L1Trigger/L1TCalorimeter/interface/legacyGtHelper.h"
#include "L1Trigger/L1TCalorimeter/interface/HardwareSortingMethods.h"

using namespace std;
using namespace l1t;

Stage1Layer2JetAlgorithmImpPP::Stage1Layer2JetAlgorithmImpPP(CaloParamsHelper const* params) : params_(params) {}

void Stage1Layer2JetAlgorithmImpPP::processEvent(const std::vector<l1t::CaloRegion>& regions,
                                                 const std::vector<l1t::CaloEmCand>& EMCands,
                                                 std::vector<l1t::Jet>* jets,
                                                 std::vector<l1t::Jet>* preGtJets) {
  std::vector<l1t::CaloRegion> subRegions;
  std::vector<l1t::Jet> uncalibjets;
  std::vector<l1t::Jet> unSortedJets;
  std::vector<l1t::Jet> preGtEtaJets;

  double towerLsb = params_->towerLsbSum();
  int jetSeedThreshold = floor(params_->jetSeedThreshold() / towerLsb + 0.5);
  std::string jetCalibrationType = params_->jetCalibrationType();
  std::vector<double> jetCalibrationParams = params_->jetCalibrationParams();

  //Region Correction will return uncorrected subregions
  //if regionPUSType is set to None in the config
  RegionCorrection(regions, &subRegions, params_);

  slidingWindowJetFinder(jetSeedThreshold, &subRegions, &uncalibjets);

  //will return jets with no response corrections
  //if jetCalibrationType is set to None in the config
  JetCalibration(&uncalibjets, jetCalibrationParams, &unSortedJets, jetCalibrationType, towerLsb);

  SortJets(&unSortedJets, &preGtEtaJets);

  // takes input jets (using region scales/eta) and outputs jets using Gt scales/eta
  JetToGtEtaScales(params_, &preGtEtaJets, preGtJets);
  JetToGtPtScales(params_, preGtJets, jets);

  //the jets should be sorted, highest pT first.
  // do not truncate the tau list, GT converter handles that
  // auto comp = [&](l1t::Jet i, l1t::Jet j)-> bool {
  //   return (i.hwPt() < j.hwPt() );
  // };

  // std::sort(jets->begin(), jets->end(), comp);
  // std::reverse(jets->begin(), jets->end());
}
