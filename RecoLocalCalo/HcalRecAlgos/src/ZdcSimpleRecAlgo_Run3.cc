#include "RecoLocalCalo/HcalRecAlgos/interface/ZdcSimpleRecAlgo_Run3.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"
#include <algorithm>  // for "max"
#include <iostream>
#include <cmath>

#include "DataFormats/METReco/interface/HcalCaloFlagLabels.h"

constexpr double MaximumFractionalError = 0.0005;  // 0.05% error allowed from this source

ZdcSimpleRecAlgo_Run3::ZdcSimpleRecAlgo_Run3(
    bool correctForTimeslew, bool correctForPulse, float phaseNS, int recoMethod)
    : recoMethod_(recoMethod),
      correctForTimeslew_(correctForTimeslew),
      correctForPulse_(correctForPulse),
      phaseNS_(phaseNS) {}

ZdcSimpleRecAlgo_Run3::ZdcSimpleRecAlgo_Run3(int recoMethod) : recoMethod_(recoMethod), correctForTimeslew_(false) {}

void ZdcSimpleRecAlgo_Run3::initPulseCorr(int toadd, const HcalTimeSlew* hcalTimeSlew_delay) {
  if (correctForPulse_) {
    pulseCorr_ = std::make_unique<HcalPulseContainmentCorrection>(
        toadd, phaseNS_, false, MaximumFractionalError, hcalTimeSlew_delay);
  }
}
namespace ZdcHelper {
  inline double subPedestal(const float energy,
                      const float ped,
                      const float width) {
                         if( energy - ped > width)return (energy - ped);
                         else return(0);
                         
                      }
                      
   // assumes that if noise TS is from pileup, noise TS should be 
   // about 40% of previous energy and only 40% of noiseTS will be subtracted from signal
  inline double calcNoise(const float energy_0, 
                      const float energy_1) {
                         // if(energy_0<=0 || energy_1<=0)return (energy_1);
                         // // check if noise TS is less than half of previous TS
                         // else if( energy_0 >= 2*energy_1) return( .4*energy_1); 
                         // else return(.4*energy_1);
                         return(.4*energy_1);
                      }
  inline double calcNoiseRPD(const float energy_0, 
                      const float energy_1) {
                         // if(energy_0<=0 || energy_1<=0)return (energy_1);
                         // // check if noise TS is less than half of previous TS
                         // else if( energy_0 >= 2*energy_1) return( .4*energy_1); 
                         // else return(.4*energy_1);
                         return(.4*energy_1);
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
  
}

namespace ZdcSimpleRecAlgoImpl {
  template <class Digi, class RecHit>
  inline RecHit reco1(const Digi& digi,
                      const HcalCoder& coder,
                      const HcalCalibrations& calibs,
                      const HcalPedestal& effPeds,
                      const std::vector<unsigned int>& myNoiseTS,
                      const std::vector<unsigned int>& mySignalTS) {
    CaloSamples tool;
    coder.adc2fC(digi, tool);
    // Reads noiseTS and signalTS from database
    int ifirst = mySignalTS[0];
    //    int n = mySignalTS.size();
    double ampl = 0;
    int maxI = -1;
    double maxA = -1e10;
    double ta = 0;
    double lowGEnergy = 0;
    double energySOIp1 = 0;
    double ratioSOIp1 = -1.0;
    double chargeWeightedTime = 0;
    
    
    double Allnoise = 0;
    int noiseslices = 0;
    int CurrentTS = 0;
    double noise = 0;
    int digi_size = digi.size();
    
    // determining noise
    for (unsigned int iv = 0; iv < myNoiseTS.size(); ++iv) {
      CurrentTS = myNoiseTS[iv];
      int capid = digi[CurrentTS].capid();
      // float ped = calibs.pedestal(capid);
      float ped = effPeds.getValue(capid);
      float width = effPeds.getWidth(capid);
      if (CurrentTS >= digi_size)
        continue;
      float energy_1 = ZdcHelper::subPedestal(tool[CurrentTS],ped,width); 
      float energy_0 = 0;
      if(CurrentTS >0)energy_0 = ZdcHelper::subPedestal(tool[CurrentTS-1],ped,width);  
      Allnoise += ZdcHelper::calcNoise(energy_0,energy_1);
      noiseslices++;
    }
    if (noiseslices != 0) {
      noise = (Allnoise) / double(noiseslices);
    } else {
      noise = 0;
    }
    
    
    for (unsigned int ivs = 0; ivs < mySignalTS.size(); ++ivs) {
      CurrentTS = mySignalTS[ivs];
      if (CurrentTS >= digi_size)
        continue;
      float energy1 = -1;
      int capid = digi[CurrentTS].capid();
      // float ped = calibs.pedestal(capid);
      float ped = effPeds.getValue(capid);
      float width = effPeds.getWidth(capid);
      float energy0 = std::max(0.0,ZdcHelper::subPedestal(tool[CurrentTS],ped,width));
      if (CurrentTS  < digi_size-1) energy1 = std::max(0.0,ZdcHelper::subPedestal(tool[CurrentTS+1],ped,width));
      ratioSOIp1 = (energy0>0 && energy1 > 0) ? energy0/energy1 : -1.0;
      ta = energy0 - noise;
      if(ta>0) ampl += ta;
      
      if (ta > maxA) {
        maxA = ta;
        maxI = CurrentTS;
      }
    }
      
    // LowGainEnergy not used currently
      lowGEnergy =-99;

    double time = -9999;
    // Time based on regular energy 
    ////Cannot calculate time value with max ADC sample at first or last position in window....
    if (maxI == 0 || maxI == (tool.size() - 1)) {
      LogDebug("HCAL Pulse") << "ZdcSimpleRecAlgo::reco2 :"
                             << " Invalid max amplitude position, "
                             << " max Amplitude: " << maxI << " first: " << ifirst << " last: " << (tool.size() - 1)
                             << std::endl;
    } else {
      int capid = digi[maxI - 1].capid();
      double Energy0 = std::max(0.0,((tool[maxI - 1]) * calibs.respcorrgain(capid)));
      
      capid = digi[maxI].capid();
      double Energy1 = std::max(0.0,((tool[maxI]) * calibs.respcorrgain(capid)));
      capid = digi[maxI + 1].capid();
      double Energy2 = std::max(0.0,((tool[maxI + 1]) * calibs.respcorrgain(capid)));
      
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
    }
    
    // find energy for signal TS + 1
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
    
    
    
    double tmp_energy = 0;
    double tmp_TSWeight_Energy = 0;
    for (int ts = 0; ts < 6; ++ts) {
      if (CurrentTS >= digi_size)
        continue;
      int capid = digi[ts].capid();
      
      // max sure there are no negative values in time calculation
      ta = std::max(0.0,tool[ts]);
      ta *= calibs.respcorrgain(capid);  // fC --> GeV
      if(ta>0){ tmp_energy += ta; tmp_TSWeight_Energy += (ts+1)*ta;}
    }    

    chargeWeightedTime = (tmp_TSWeight_Energy/tmp_energy -1) * 25.0;
    auto rh = RecHit(digi.id(), ampl, time, lowGEnergy);
    if(maxI >= 0 && maxI <tool.size()) rh.setSaturated(digi[maxI].adc()>=255);
    rh.setEnergySOIp1(energySOIp1);
    
    
    if(ZdcHelper::isSaturated(digi)) rh.setFlagField(1, HcalCaloFlagLabels::ADCSaturationBit);
    float tmp_tdctime = 0;
    if(maxI >= 0 && maxI <tool.size()){
       // TDC error codes will be 60=-1, 61 = -2, 62 = -3, 63 = -4
       // if(digi[maxI].le_tdc()>=60) tmp_tdctime = -1*(digi[maxI].le_tdc() - 59)

       // should change tdc values for 63, 62 to more meaningful values usesful for quick analysis currently          
       if(digi[maxI].le_tdc()==62) tmp_tdctime = -10; 
       else if(digi[maxI].le_tdc()==63) tmp_tdctime = -30;
       else tmp_tdctime = maxI* 25. + (digi[maxI].le_tdc() / 2);
    }
    
    if(maxI >= 0 && maxI <tool.size())  rh.setTDCtime(tmp_tdctime);
    rh.setChargeWeightedTime(chargeWeightedTime);
    rh.setRatioSOIp1(ratioSOIp1);
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
                      const std::vector<unsigned int>& mySignalTS) {
    CaloSamples tool;
    coder.adc2fC(digi, tool);
    // Reads noiseTS and signalTS from database
    int ifirst = mySignalTS[0];
    //    int n = mySignalTS.size();
    double ampl = 0;
    int maxI = -1;
    double maxA = -1e10;
    double ta = 0;
    double lowGEnergy = 0;
    double energySOIp1 = 0;
    double ratioSOIp1 = -1;
    double chargeWeightedTime = 0;
    
    double noiseFrac = 97.0/256.0;
    bool matchTriggerInteger = false;
    

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
    for (unsigned int ivs = 0; ivs < mySignalTS.size(); ++ivs) {
      CurrentTS = mySignalTS[ivs];
      if(CurrentTS >= digi_size) continue;
      float energy1 = -1;
      int capid = digi[CurrentTS].capid();
      // float ped = calibs.pedestal(capid);
      float ped = effPeds.getValue(capid);
      float width = effPeds.getWidth(capid);
      float energy0 = std::max(0.0,ZdcHelper::subPedestal(tool[CurrentTS],ped,width));
      if (CurrentTS  < digi_size-1) energy1 = std::max(0.0,ZdcHelper::subPedestal(tool[CurrentTS+1],ped,width));
      ratioSOIp1 = (energy0>0 && energy1 > 0) ? energy0/energy1 : -1.0;
    if(matchTriggerInteger){
      float gain = calibs.respcorrgain(capid);
      ta = std::min(std::max(0, int(energy0*gain/(50))), 1023) - std::min(std::max(0, int(noise*gain/(50))), 1023);
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
    // Time based on regular energy 
    ////Cannot calculate time value with max ADC sample at first or last position in window....
    if (maxI == 0 || maxI == (tool.size() - 1)) {
      LogDebug("HCAL Pulse") << "ZdcSimpleRecAlgo::reco2 :"
                             << " Invalid max amplitude position, "
                             << " max Amplitude: " << maxI << " first: " << ifirst << " last: " << (tool.size() - 1)
                             << std::endl;
    } else {
      int capid = digi[maxI - 1].capid();
      double Energy0 = std::max(0.0,((tool[maxI - 1]) * calibs.respcorrgain(capid)));
      
      capid = digi[maxI].capid();
      double Energy1 = std::max(0.0,((tool[maxI]) * calibs.respcorrgain(capid)));
      capid = digi[maxI + 1].capid();
      double Energy2 = std::max(0.0,((tool[maxI + 1]) * calibs.respcorrgain(capid)));
      
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
    }
    // find energy for signal TS + 1
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
    
    
    
    double tmp_energy = 0;
    double tmp_TSWeight_Energy = 0;
    for (int ts = 0; ts < 6; ++ts) {
      if (CurrentTS >= digi_size)
        continue;
      int capid = digi[ts].capid();
      
      // max sure there are no negative values in time calculation
      ta = std::max(0.0,tool[ts]);
      ta *= calibs.respcorrgain(capid);  // fC --> GeV
      if(ta>0){ tmp_energy += ta; tmp_TSWeight_Energy += (ts+1)*ta;}
    }    

    chargeWeightedTime = (tmp_TSWeight_Energy/tmp_energy -1) * 25.0;
    auto rh = RecHit(digi.id(), ampl, time, lowGEnergy);
    if(maxI >= 0 && maxI <tool.size()) rh.setSaturated(digi[maxI].adc()>=255);
    rh.setEnergySOIp1(energySOIp1);
    
    
    if(ZdcHelper::isSaturated(digi)) rh.setFlagField(1, HcalCaloFlagLabels::ADCSaturationBit);
    float tmp_tdctime = 0;
    if(maxI >= 0 && maxI <tool.size()){
       // TDC error codes will be 60=-1, 61 = -2, 62 = -3, 63 = -4
       // if(digi[maxI].le_tdc()>=60) tmp_tdctime = -1*(digi[maxI].le_tdc() - 59)

       // should change tdc values for 63, 62 to more meaningful values usesful for quick analysis currently          
       if(digi[maxI].le_tdc()==62) tmp_tdctime = -10; 
       else if(digi[maxI].le_tdc()==63) tmp_tdctime = -30;
       else tmp_tdctime = maxI* 25. + (digi[maxI].le_tdc() / 2);
    }
    
    if(maxI >= 0 && maxI <tool.size())  rh.setTDCtime(tmp_tdctime);
    rh.setChargeWeightedTime(chargeWeightedTime);
    rh.setRatioSOIp1(ratioSOIp1);
    return rh;
  }
}  



// reco3 is used for RPD rechits
namespace ZdcSimpleRecAlgoImpl {
  template <class Digi, class RecHit>
  inline RecHit reco3(const Digi& digi,
                      const HcalCoder& coder,
                      const HcalCalibrations& calibs,
                      const HcalPedestal& effPeds,
                      const std::vector<unsigned int>& myNoiseTS,
                      const std::vector<unsigned int>& mySignalTS) {
    CaloSamples tool;
    coder.adc2fC(digi, tool);
    // Reads noiseTS and signalTS from database
    int ifirst = mySignalTS[0];
    //    int n = mySignalTS.size();
    double ampl = 0;
    int maxI = -1;
    double maxA = -1e10;
    double ta = 0;
    double lowGEnergy = 0;
    double energySOIp1 = 0;
    double ratioSOIp1 = -1;
    double chargeWeightedTime = 0;
    
    
    double Allnoise = 0;
    int noiseslices = 0;
    int CurrentTS = 0;
    double noise = 0;
    int digi_size = digi.size();
    
    
    // determining noise
    for (unsigned int iv = 0; iv < myNoiseTS.size(); ++iv) {
      CurrentTS = myNoiseTS[iv];
      int capid = digi[CurrentTS].capid();
      // float ped = calibs.pedestal(capid);
      float ped = effPeds.getValue(capid);
      float width = effPeds.getWidth(capid);
      if (CurrentTS >= digi_size)
        continue;
      float energy_1 = ZdcHelper::subPedestal(tool[CurrentTS],ped,width); 
      float energy_0 = 0;
      if(CurrentTS >0) energy_0 = ZdcHelper::subPedestal(tool[CurrentTS-1],ped,width);  
      Allnoise += ZdcHelper::calcNoiseRPD(energy_0,energy_1);
      noiseslices++;
    }
    if (noiseslices != 0) {
      noise = (Allnoise) / double(noiseslices);
    } else {
      noise = 0;
    }
    
    
    for (unsigned int ivs = 0; ivs < mySignalTS.size(); ++ivs) {
      CurrentTS = mySignalTS[ivs];
      if(CurrentTS >= digi_size) continue;
      float energy1 = -1;
      int capid = digi[CurrentTS].capid();
      // float ped = calibs.pedestal(capid);
      float ped = effPeds.getValue(capid);
      float width = effPeds.getWidth(capid);
      float energy0 = std::max(0.0,ZdcHelper::subPedestal(tool[CurrentTS],ped,width));
      if (CurrentTS  < digi_size-1) energy1 = std::max(0.0,ZdcHelper::subPedestal(tool[CurrentTS+1],ped,width));
      ratioSOIp1 = (energy0>0 && energy1 > 0) ? energy0/energy1 : -1.0;
      ta = energy0 - noise;
      ta *= calibs.respcorrgain(capid);  // fC --> GeV 
      if(ta>0) ampl += ta;
      
      if (ta > maxA) {
        maxA = ta;
        maxI = CurrentTS;
      }
    }
      
      
    // LowGainEnergy not used currently
      lowGEnergy =-99;

    double time = -9999;
    // Time based on regular energy 
    ////Cannot calculate time value with max ADC sample at first or last position in window....
    if (maxI == 0 || maxI == (tool.size() - 1)) {
      LogDebug("HCAL Pulse") << "ZdcSimpleRecAlgo::reco2 :"
                             << " Invalid max amplitude position, "
                             << " max Amplitude: " << maxI << " first: " << ifirst << " last: " << (tool.size() - 1)
                             << std::endl;
    } else {
      int capid = digi[maxI - 1].capid();
      double Energy0 = std::max(0.0,((tool[maxI - 1]) * calibs.respcorrgain(capid)));
      
      capid = digi[maxI].capid();
      double Energy1 = std::max(0.0,((tool[maxI]) * calibs.respcorrgain(capid)));
      capid = digi[maxI + 1].capid();
      double Energy2 = std::max(0.0,((tool[maxI + 1]) * calibs.respcorrgain(capid)));
      
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
    }
    
    // find energy for signal TS + 1
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
    
    
    
    double tmp_energy = 0;
    double tmp_TSWeight_Energy = 0;
    for (int ts = 0; ts < 6; ++ts) {
      if (CurrentTS >= digi_size)
        continue;
      int capid = digi[ts].capid();
      
      // max sure there are no negative values in time calculation
      ta = std::max(0.0,tool[ts]);
      ta *= calibs.respcorrgain(capid);  // fC --> GeV
      if(ta>0){ tmp_energy += ta; tmp_TSWeight_Energy += (ts+1)*ta;}
    }    

    chargeWeightedTime = (tmp_TSWeight_Energy/tmp_energy -1) * 25.0;
    auto rh = RecHit(digi.id(), ampl, time, lowGEnergy);
    if(maxI >= 0 && maxI <tool.size()) rh.setSaturated(digi[maxI].adc()>=255);
    rh.setEnergySOIp1(energySOIp1);
    
    
    if(ZdcHelper::isSaturated(digi)) rh.setFlagField(1, HcalCaloFlagLabels::ADCSaturationBit);
    float tmp_tdctime = 0;
    if(maxI >= 0 && maxI <tool.size()){
       // TDC error codes will be 60=-1, 61 = -2, 62 = -3, 63 = -4
       // if(digi[maxI].le_tdc()>=60) tmp_tdctime = -1*(digi[maxI].le_tdc() - 59)

       // should change tdc values for 63, 62 to more meaningful values usesful for quick analysis currently          
       if(digi[maxI].le_tdc()==62) tmp_tdctime = -10; 
       else if(digi[maxI].le_tdc()==63) tmp_tdctime = -30;
       else tmp_tdctime = maxI* 25. + (digi[maxI].le_tdc() / 2);
    }
    
    if(maxI >= 0 && maxI <tool.size())  rh.setTDCtime(tmp_tdctime);
    rh.setChargeWeightedTime(chargeWeightedTime);
    rh.setRatioSOIp1(ratioSOIp1);
    return rh;
  }
}  // namespace ZdcSi



ZDCRecHit ZdcSimpleRecAlgo_Run3::reconstruct(const QIE10DataFrame& digi,
                                        const std::vector<unsigned int>& myNoiseTS,
                                        const std::vector<unsigned int>& mySignalTS,
                                        const HcalCoder& coder,
                                        const HcalCalibrations& calibs,
                                        const HcalPedestal& effPeds,
                                        const bool& isRPD) const {
    if(isRPD) return ZdcSimpleRecAlgoImpl::reco3<QIE10DataFrame, ZDCRecHit>(
        digi, coder, calibs, effPeds, myNoiseTS, mySignalTS);
    else{
       if(recoMethod_==1) return ZdcSimpleRecAlgoImpl::reco1<QIE10DataFrame, ZDCRecHit>(
        digi, coder, calibs, effPeds, myNoiseTS, mySignalTS);
       if(recoMethod_==2) return ZdcSimpleRecAlgoImpl::reco2<QIE10DataFrame, ZDCRecHit>(
        digi, coder, calibs, effPeds, myNoiseTS, mySignalTS);
    }
   edm::LogError("ZDCSimpleRecAlgoImpl::reconstruct, recoMethod was not declared");
  throw cms::Exception("ZDCSimpleRecoAlgos::exiting process");
}