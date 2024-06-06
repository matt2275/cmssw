#include "RecoLocalCalo/HcalRecAlgos/interface/ZdcSimpleRecAlgo_Run3.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"
#include <algorithm>  // for "max"
#include <iostream>
#include <cmath>

#include "DataFormats/METReco/interface/HcalCaloFlagLabels.h"

constexpr double MaximumFractionalError = 0.0005;  // 0.05% error allowed from this source

ZdcSimpleRecAlgo_Run3::ZdcSimpleRecAlgo_Run3(
    bool correctForTimeslew, bool correctForPulse, float phaseNS, int recoMethod, int lowGainOffset, double lowGainFrac)
    : recoMethod_(recoMethod),
      correctForTimeslew_(correctForTimeslew),
      correctForPulse_(correctForPulse),
      phaseNS_(phaseNS),
      lowGainOffset_(lowGainOffset),
      lowGainFrac_(lowGainFrac) {}

ZdcSimpleRecAlgo_Run3::ZdcSimpleRecAlgo_Run3(int recoMethod) : recoMethod_(recoMethod), correctForTimeslew_(false) {}

void ZdcSimpleRecAlgo_Run3::initPulseCorr(int toadd, const HcalTimeSlew* hcalTimeSlew_delay) {
  if (correctForPulse_) {
    pulseCorr_ = std::make_unique<HcalPulseContainmentCorrection>(
        toadd, phaseNS_, false, MaximumFractionalError, hcalTimeSlew_delay);
  }
}
namespace ZdcHelper {
  inline float subPedestal(const float energy,
                      const float ped,
                      const float width) {
                         if( energy - ped > width)return (energy - ped);
                         else return(0);
                         
                      }
                      
                      
  inline float calcNoise(const float energy,
                      const float ped,
                      const float width) {
                         if( energy - ped > width)return (energy - ped);
                         else return(0);
                         
                      }
  template <class Digi>                     
  bool isSaturated(Digi& digi){
     unsigned int digi_size = digi.size();
     bool saturated= false;
     for(unsigned int i =0; i < digi_size; i++){
        if(digi[i].adc() >= 255){saturated=true; break;}
     }
     return saturated;
  }
  
int GetNormalLUTInt(float charge){
   return std::min(std::max(0, int(charge/(50))), 1023);
}
int GetOOTPULUTInt(float charge){
   return std::min(std::max(0, int(charge*97/(50*256))), 1023);
}
}
//static float timeshift_ns_zdc(float wpksamp);

namespace ZdcSimpleRecAlgoImpl {
  template <class Digi, class RecHit>
  inline RecHit reco1(const Digi& digi,
                      const HcalCoder& coder,
                      const HcalCalibrations& calibs,
                      const std::vector<unsigned int>& myNoiseTS,
                      const std::vector<unsigned int>& mySignalTS,
                      int lowGainOffset,
                      double lowGainFrac,
                      bool slewCorrect,
                      const HcalPulseContainmentCorrection* corr,
                      HcalTimeSlew::BiasSetting slewFlavor) {
    CaloSamples tool;
    coder.adc2fC(digi, tool);
    // Reads noiseTS and signalTS from database
    int ifirst = mySignalTS[0];
    //    int n = mySignalTS.size();
    double ampl = 0;
    int maxI = -1;
    double maxA = -1e10;
    double ta = 0;
    double fc_ampl = 0;
    double lowGEnergy = 0;
    double TempLGAmp = 0;
    double energySOIp1 = 0;
    double timeSOIp1 = 0;
    double chargeWeightedTime = 0;
    //  TS increment for regular energy to lowGainEnergy
    // Signal in higher TS (effective "low Gain") has a fraction of the whole signal
    // This constant for fC --> GeV is dervied from 2010 PbPb analysis of single neutrons
    // assumed similar fraction for EM and HAD sections
    // this variable converts from current assumed TestBeam values for fC--> GeV
    // to the lowGain TS region fraction value (based on 1N Had, assume EM same response)
    double Allnoise = 0;
    int noiseslices = 0;
    int CurrentTS = 0;
    double noise = 0;
    int digi_size = digi.size();
    // bool isSaturated = false;
    // regular energy (both use same noise)
    for (unsigned int iv = 0; iv < myNoiseTS.size(); ++iv) {
      CurrentTS = myNoiseTS[iv];
      int capid = digi[CurrentTS].capid();
      float ped = calibs.pedestal(capid);
      if (CurrentTS >= digi_size)
        continue;
      Allnoise += ZdcHelper::subPedestal(tool[CurrentTS],ped,0);
      // if(tool[CurrentTS]>ped) Allnoise += tool[CurrentTS];
      // Allnoise += tool[CurrentTS];
      noiseslices++;
    }
    if (noiseslices != 0) {
      noise = (Allnoise) / double(noiseslices);
    } else {
      noise = 0;
    }
    HcalZDCDetId cell = digi.id();
    if(cell.section()>2) noise =0;
    for (unsigned int ivs = 0; ivs < mySignalTS.size(); ++ivs) {
      CurrentTS = mySignalTS[ivs];
      if (CurrentTS >= digi_size)
        continue;
      int capid = digi[CurrentTS].capid();
      //       if(noise<0){
      //       // flag hit as having negative noise, and don't subtract anything, because
      //       // it will falsely increase the energy
      //          noisefactor=0.;
      //       }
      float ped = calibs.pedestal(capid);
      ta = ZdcHelper::subPedestal(tool[CurrentTS],ped,0) - noise;
      // ta = tool[CurrentTS] - noise;
      if(ta> 0) fc_ampl += ta-ped;
      ta *= calibs.respcorrgain(capid);  // fC --> GeV
      if(ta>0) ampl += ta;
      if (ta > maxA) {
        maxA = ta;
        maxI = CurrentTS;
      }
      // std::cout << "capid " << capid << " ta " << ta << " noise " << noise << " gains " << calibs.respcorrgain(capid) << std::endl;
    }
    // calculate low Gain Energy (in 2010 PbPb, signal TS 4,5,6, lowGain TS: 6,7,8)
    for (unsigned int iLGvs = 0; iLGvs < mySignalTS.size(); ++iLGvs) {
      CurrentTS = mySignalTS[iLGvs] + lowGainOffset;
      if (CurrentTS >= digi_size)
        continue;
      int capid = digi[CurrentTS].capid();
      TempLGAmp = tool[CurrentTS] - noise;
      TempLGAmp *= calibs.respcorrgain(capid);  // fC --> GeV
      TempLGAmp *= lowGainFrac;                 // TS (signalRegion) --> TS (lowGainRegion)
      lowGEnergy += TempLGAmp;
    }
    //    if(ta<0){
    //      // flag hits that have negative energy
    //    }

    double time = -9999;
    // Time based on regular energy (lowGainEnergy signal assumed to happen at same Time)
    ////Cannot calculate time value with max ADC sample at first or last position in window....
    if (maxI == 0 || maxI == (tool.size() - 1)) {
      LogDebug("HCAL Pulse") << "ZdcSimpleRecAlgo::reco2 :"
                             << " Invalid max amplitude position, "
                             << " max Amplitude: " << maxI << " first: " << ifirst << " last: " << (tool.size() - 1)
                             << std::endl;
    } else {
      int capid = digi[maxI - 1].capid();
      double Energy0 = ((tool[maxI - 1]) * calibs.respcorrgain(capid));
      // if any of the energies used in the weight are negative, make them 0 instead
      // these are actually QIE values, not energy
      if (Energy0 < 0) {
        Energy0 = 0.;
      }
      capid = digi[maxI].capid();
      double Energy1 = ((tool[maxI]) * calibs.respcorrgain(capid));
      if (Energy1 < 0) {
        Energy1 = 0.;
      }
      capid = digi[maxI + 1].capid();
      double Energy2 = ((tool[maxI + 1]) * calibs.respcorrgain(capid));
      if (Energy2 < 0) {
        Energy2 = 0.;
      }
      //
      double TSWeightEnergy = ((maxI - 1) * Energy0 + maxI * Energy1 + (maxI + 1) * Energy2);
      double EnergySum = Energy0 + Energy1 + Energy2;
      double AvgTSPos = 0.;
      if (EnergySum != 0)
        AvgTSPos = TSWeightEnergy / EnergySum;
      // If time is zero, set it to the "nonsensical" -99
      // Time should be between 75ns and 175ns (Timeslices 3-7)
      if (AvgTSPos == 0) {
        time = -99;
      } else {
        time = (AvgTSPos * 25.0);
      }
      // if (corr != nullptr) {
        // // Apply phase-based amplitude correction:
        // ampl *= corr->getCorrection(fc_ampl);
      // }
    }
    
   
    for (unsigned int ivs = 0; ivs < mySignalTS.size(); ++ivs) {
      CurrentTS = mySignalTS[ivs] + 1;
      if (CurrentTS >= digi_size)
        continue;
      int capid = digi[CurrentTS].capid();
      // ta = tool[CurrentTS] - noise;
      ta = tool[CurrentTS];
      ta *= calibs.respcorrgain(capid);  // fC --> GeV
      if(ta>0) energySOIp1 += ta;
    }
    
    if (maxI < (tool.size() - 2)){
      int capid = digi[maxI].capid();
      double Energy0 = ((tool[maxI]) * calibs.respcorrgain(capid));
      // if any of the energies used in the weight are negative, make them 0 instead
      // these are actually QIE values, not energy
      if (Energy0 < 0) {
        Energy0 = 0.;
      }
      capid = digi[maxI +1 ].capid();
      double Energy1 = ((tool[maxI+1]) * calibs.respcorrgain(capid));
      if (Energy1 < 0) {
        Energy1 = 0.;
      }
      capid = digi[maxI + 2].capid();
      double Energy2 = ((tool[maxI + 2]) * calibs.respcorrgain(capid));
      if (Energy2 < 0) {
        Energy2 = 0.;
      }
      //
      double TSWeightEnergy = ((maxI) * Energy0 + (maxI+1) * Energy1 + (maxI + 2) * Energy2);
      double EnergySum = Energy0 + Energy1 + Energy2;
      double AvgTSPos = 0.;
      if (EnergySum != 0)
        AvgTSPos = TSWeightEnergy / EnergySum;
      // If time is zero, set it to the "nonsensical" -99
      // Time should be between 75ns and 175ns (Timeslices 3-7)
      if (AvgTSPos == 0) {
        timeSOIp1 = -99;
      } else {
        timeSOIp1 = (AvgTSPos * 25.0);
      }
    } 
    
    double tmp_energy = 0;
    double tmp_TSWeight_Energy = 0;
    for (int ts = 0; ts < 6; ++ts) {
      if (CurrentTS >= digi_size)
        continue;
      int capid = digi[ts].capid();
      // ta = tool[CurrentTS] - noise;
      ta = tool[ts];
      ta *= calibs.respcorrgain(capid);  // fC --> GeV
      if(ta>0){ tmp_energy += ta; tmp_TSWeight_Energy += (ts+1)*ta;}
    }    

    chargeWeightedTime = (tmp_TSWeight_Energy/tmp_energy -1) * 25.0;
    // std::cout <<  "Energy " << ampl << " fc " << tool[3] << " digi " << digi[maxI].adc() <<  std::endl;
    auto rh = RecHit(digi.id(), ampl, time, lowGEnergy);
    if(maxI >= 0 && maxI <tool.size()) rh.setSaturated(digi[maxI].adc()>=255);
    rh.setenergySOIp1(energySOIp1);
    rh.setTimeSOIp1(timeSOIp1);
    if(ZdcHelper::isSaturated(digi)) rh.setFlagField(1, HcalCaloFlagLabels::ADCSaturationBit);
    float tmp_tdctime = 0;
    if(maxI >= 0 && maxI <tool.size()){
       if(digi[maxI].le_tdc()==62) tmp_tdctime = -10;
       else if(digi[maxI].le_tdc()==63) tmp_tdctime = -30;
       else tmp_tdctime = maxI* 25. + (digi[maxI].le_tdc() / 2);
    }
    // if(maxI >= 0 && maxI <tool.size())  std::cout <<  "le tdc " << digi[maxI].le_tdc() << std::endl;
    if(maxI >= 0 && maxI <tool.size())  rh.setTDCtime(tmp_tdctime);
    rh.setChargeWeightedTime(chargeWeightedTime);
    return rh;
  }
}  // namespace ZdcSi

//reco2 method should make the Energy Assosciated with the ZDC Trigger Energy as of June 2024
// only difference between this method and the trigger value LUT value is the signal and noise TS
// are bins into 1024 bins of 50 GeV bin widths. i.e. energyint = std::min(std::max(0, int(energy/(50))), 1023)
// this assumes the noise TS =2 and signalTS =2
namespace ZdcSimpleRecAlgoImpl {
  template <class Digi, class RecHit>
  inline RecHit reco2(const Digi& digi,
                      const HcalCoder& coder,
                      const HcalCalibrations& calibs,
                      const HcalPedestal& effPeds,
                      const std::vector<unsigned int>& myNoiseTS,
                      const std::vector<unsigned int>& mySignalTS,
                      int lowGainOffset,
                      double lowGainFrac,
                      bool slewCorrect,
                      const HcalPulseContainmentCorrection* corr,
                      HcalTimeSlew::BiasSetting slewFlavor) {
    CaloSamples tool;
    coder.adc2fC(digi, tool);
    // Reads noiseTS and signalTS from database
    int ifirst = mySignalTS[0];
    //    int n = mySignalTS.size();
    double ampl = 0;
    int maxI = -1;
    double maxA = -1e10;
    double ta = 0;
    double fc_ampl = 0;
    double lowGEnergy = 0;
    double TempLGAmp = 0;
    double energySOIp1 = 0;
    double timeSOIp1 = 0;
    double chargeWeightedTime = 0;
    
    double noiseFrac = 97.0/256.0;
    bool matchTriggerInteger = true;
    
    //  TS increment for regular energy to lowGainEnergy
    // Signal in higher TS (effective "low Gain") has a fraction of the whole signal
    // This constant for fC --> GeV is dervied from 2010 PbPb analysis of single neutrons
    // assumed similar fraction for EM and HAD sections
    // this variable converts from current assumed TestBeam values for fC--> GeV
    // to the lowGain TS region fraction value (based on 1N Had, assume EM same response)
    double Allnoise = 0;
    int noiseslices = 0;
    int CurrentTS = 0;
    double noise = 0;
    int digi_size = digi.size();
    // bool isSaturated = false;
    // regular energy (both use same noise)
    for (unsigned int iv = 0; iv < myNoiseTS.size(); ++iv) {
      CurrentTS = myNoiseTS[iv];
      int capid = digi[CurrentTS].capid();
      // float ped = calibs.pedestal(capid);
      float ped = effPeds.getValue(capid);
      float width = effPeds.getWidth(capid);
      if (CurrentTS >= digi_size)
        continue;
      Allnoise += ZdcHelper::subPedestal(tool[CurrentTS],ped,width)*noiseFrac;
      noiseslices++;
    }
    if (noiseslices != 0) {
      noise = (Allnoise) / double(noiseslices);
    } else {
      noise = 0;
    }
    HcalZDCDetId cell = digi.id();
    if(cell.section()>2) noise =0;
    for (unsigned int ivs = 0; ivs < mySignalTS.size(); ++ivs) {
      CurrentTS = mySignalTS[ivs];
      if (CurrentTS >= digi_size)
        continue;
      int capid = digi[CurrentTS].capid();

      // float ped = calibs.pedestal(capid);
      float ped = effPeds.getValue(capid);
      float width = effPeds.getWidth(capid);
      
    if(matchTriggerInteger){
      float gain = calibs.respcorrgain(capid);
      ta = std::min(std::max(0, int(ZdcHelper::subPedestal(tool[CurrentTS],ped,width)*gain/(50))), 1023) - std::min(std::max(0, int(noise*gain/(50))), 1023);
    }
    else{
      ta = ZdcHelper::subPedestal(tool[CurrentTS],ped,width) - noise;
      ta *= calibs.respcorrgain(capid);  // fC --> GeV
    }
      if(ta>0) ampl += ta;
      if (ta > maxA) {
        maxA = ta;
        maxI = CurrentTS;
      }
    }
    
    // LowGainEnergy not used currently
      lowGEnergy =-99;

    double time = -9999;
    // Time based on regular energy (lowGainEnergy signal assumed to happen at same Time)
    ////Cannot calculate time value with max ADC sample at first or last position in window....
    if (maxI == 0 || maxI == (tool.size() - 1)) {
      LogDebug("HCAL Pulse") << "ZdcSimpleRecAlgo::reco2 :"
                             << " Invalid max amplitude position, "
                             << " max Amplitude: " << maxI << " first: " << ifirst << " last: " << (tool.size() - 1)
                             << std::endl;
    } else {
      int capid = digi[maxI - 1].capid();
      double Energy0 = ((tool[maxI - 1]) * calibs.respcorrgain(capid));
      // if any of the energies used in the weight are negative, make them 0 instead
      // these are actually QIE values, not energy
      if (Energy0 < 0) {
        Energy0 = 0.;
      }
      capid = digi[maxI].capid();
      double Energy1 = ((tool[maxI]) * calibs.respcorrgain(capid));
      if (Energy1 < 0) {
        Energy1 = 0.;
      }
      capid = digi[maxI + 1].capid();
      double Energy2 = ((tool[maxI + 1]) * calibs.respcorrgain(capid));
      if (Energy2 < 0) {
        Energy2 = 0.;
      }
      //
      double TSWeightEnergy = ((maxI - 1) * Energy0 + maxI * Energy1 + (maxI + 1) * Energy2);
      double EnergySum = Energy0 + Energy1 + Energy2;
      double AvgTSPos = 0.;
      if (EnergySum != 0)
        AvgTSPos = TSWeightEnergy / EnergySum;
      // If time is zero, set it to the "nonsensical" -99
      // Time should be between 75ns and 175ns (Timeslices 3-7)
      if (AvgTSPos == 0) {
        time = -99;
      } else {
        time = (AvgTSPos * 25.0);
      }
      // if (corr != nullptr) {
        // // Apply phase-based amplitude correction:
        // ampl *= corr->getCorrection(fc_ampl);
      // }
    }
    
   // find Energy for SOI+1
    for (unsigned int ivs = 0; ivs < mySignalTS.size(); ++ivs) {
      CurrentTS = mySignalTS[ivs] + 1;
      if (CurrentTS >= digi_size)
        continue;
      int capid = digi[CurrentTS].capid();
      // ta = tool[CurrentTS] - noise;
      ta = tool[CurrentTS];
      ta *= calibs.respcorrgain(capid);  // fC --> GeV
      if(ta>0) energySOIp1 += ta;
    }
 
   // find Time for SOI+1 
    if (maxI < (tool.size() - 2)){
      int capid = digi[maxI].capid();
      double Energy0 = ((tool[maxI]) * calibs.respcorrgain(capid));
      // if any of the energies used in the weight are negative, make them 0 instead
      // these are actually QIE values, not energy
      if (Energy0 < 0) {
        Energy0 = 0.;
      }
      capid = digi[maxI +1 ].capid();
      double Energy1 = ((tool[maxI+1]) * calibs.respcorrgain(capid));
      if (Energy1 < 0) {
        Energy1 = 0.;
      }
      capid = digi[maxI + 2].capid();
      double Energy2 = ((tool[maxI + 2]) * calibs.respcorrgain(capid));
      if (Energy2 < 0) {
        Energy2 = 0.;
      }
      //
      double TSWeightEnergy = ((maxI) * Energy0 + (maxI+1) * Energy1 + (maxI + 2) * Energy2);
      double EnergySum = Energy0 + Energy1 + Energy2;
      double AvgTSPos = 0.;
      if (EnergySum != 0)
        AvgTSPos = TSWeightEnergy / EnergySum;
      // If time is zero, set it to the "nonsensical" -99
      // Time should be between 75ns and 175ns (Timeslices 3-7)
      if (AvgTSPos == 0) {
        timeSOIp1 = -99;
      } else {
        timeSOIp1 = (AvgTSPos * 25.0);
      }
    } 


    //finding the charge weighed time for first 6 TS
    double tmp_energy = 0;
    double tmp_TSWeight_Energy = 0;
    for (int ts = 0; ts < 6; ++ts) {
      if (CurrentTS >= digi_size)
        continue;
      int capid = digi[ts].capid();
      // ta = tool[CurrentTS] - noise;
      ta = tool[ts];
      ta *= calibs.respcorrgain(capid);  // fC --> GeV
      if(ta>0){ tmp_energy += ta; tmp_TSWeight_Energy += (ts+1)*ta;}
    }    

    chargeWeightedTime = (tmp_TSWeight_Energy/tmp_energy -1) * 25.0;
    
    auto rh = RecHit(digi.id(), ampl, time, lowGEnergy);
    if(maxI >= 0 && maxI <tool.size()) rh.setSaturated(digi[maxI].adc()>=255);
    rh.setenergySOIp1(energySOIp1);
    rh.setTimeSOIp1(timeSOIp1);
    if(ZdcHelper::isSaturated(digi)) rh.setFlagField(1, HcalCaloFlagLabels::ADCSaturationBit);
    float tmp_tdctime = 0;
    if(maxI >= 0 && maxI <tool.size()){
       if(digi[maxI].le_tdc()==62) tmp_tdctime = -10;
       else if(digi[maxI].le_tdc()==63) tmp_tdctime = -30;
       else tmp_tdctime = maxI* 25. + (digi[maxI].le_tdc() / 2);
    }
    if(maxI >= 0 && maxI <tool.size())  rh.setTDCtime(tmp_tdctime);
    rh.setChargeWeightedTime(chargeWeightedTime);
    return rh;
  }
}  



ZDCRecHit ZdcSimpleRecAlgo_Run3::reconstruct(const QIE10DataFrame& digi,
                                        const std::vector<unsigned int>& myNoiseTS,
                                        const std::vector<unsigned int>& mySignalTS,
                                        const HcalCoder& coder,
                                        const HcalCalibrations& calibs) const {
    return ZdcSimpleRecAlgoImpl::reco1<QIE10DataFrame, ZDCRecHit>(
        digi, coder, calibs, myNoiseTS, mySignalTS, lowGainOffset_, lowGainFrac_, false, nullptr, HcalTimeSlew::Fast);
  edm::LogError("ZDCSimpleRecAlgoImpl::reconstruct, recoMethod was not declared");
  throw cms::Exception("ZDCSimpleRecoAlgos::exiting process");
}

ZDCRecHit ZdcSimpleRecAlgo_Run3::reconstruct(const QIE10DataFrame& digi,
                                        const std::vector<unsigned int>& myNoiseTS,
                                        const std::vector<unsigned int>& mySignalTS,
                                        const HcalCoder& coder,
                                        const HcalCalibrations& calibs,
                                        const HcalPedestal& effPeds) const {
    return ZdcSimpleRecAlgoImpl::reco2<QIE10DataFrame, ZDCRecHit>(
        digi, coder, calibs, effPeds, myNoiseTS, mySignalTS, lowGainOffset_, lowGainFrac_, false, nullptr, HcalTimeSlew::Fast);
  edm::LogError("ZDCSimpleRecAlgoImpl::reconstruct, recoMethod was not declared");
  throw cms::Exception("ZDCSimpleRecoAlgos::exiting process");
}