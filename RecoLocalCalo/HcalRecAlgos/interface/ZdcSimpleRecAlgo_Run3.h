#ifndef ZDCSIMPLERECALGO_RUN3_H
#define ZDCSIMPLERECALGO_RUN3_H 1

#include "DataFormats/HcalDigi/interface/HBHEDataFrame.h"
#include "DataFormats/HcalDigi/interface/HFDataFrame.h"
#include "DataFormats/HcalDigi/interface/HODataFrame.h"
#include "DataFormats/HcalDigi/interface/ZDCDataFrame.h"
#include "DataFormats/HcalDigi/interface/HcalCalibDataFrame.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalRecHit/interface/ZDCRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalCalibRecHit.h"
#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrationWidths.h"
#include "CondFormats/HcalObjects/interface/HcalPedestal.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseContainmentCorrection.h"
#include <memory>


#include "DataFormats/HcalDigi/interface/QIE10DataFrame.h"

/** \class ZdcSimpleRecAlgo_Run3

   This class reconstructs RecHits from Digis for ZDC  by addition
   of selected time samples, pedestal subtraction, and gain application. The
   time of the hit is reconstructed using a weighted peak bin calculation.
    
   \author M. Nickel Based on ZDCSimpleRecAlgo
*/
class HcalTimeSlew;

class ZdcSimpleRecAlgo_Run3 {
public:
  /** Full featured constructor for ZDC */
  ZdcSimpleRecAlgo_Run3(bool correctForTimeslew,
                   bool correctForContainment,
                   float fixedPhaseNs,
                   int recoMethod);
  /** Simple constructor for PMT-based detectors */
  ZdcSimpleRecAlgo_Run3(int recoMethod);
  void initPulseCorr(int toadd, const HcalTimeSlew* hcalTimeSlew_delay);

// reco method currently used to match L1 Trigger LUT energy formula                        
  ZDCRecHit reconstruct(const QIE10DataFrame& digi,
                        const std::vector<unsigned int>& myNoiseTS,
                        const std::vector<unsigned int>& mySignalTS,
                        const HcalCoder& coder,
                        const HcalCalibrations& calibs,
                        const HcalPedestal& effPeds,
                        const bool& isRPD) const;

private:
  int recoMethod_;
  bool correctForTimeslew_;
  bool correctForPulse_;
  float phaseNS_;

  std::unique_ptr<HcalPulseContainmentCorrection> pulseCorr_;
};

#endif
