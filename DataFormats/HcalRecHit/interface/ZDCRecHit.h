#ifndef DATAFORMATS_HCALRECHIT_ZDCRECHIT_H
#define DATAFORMATS_HCALRECHIT_ZDCRECHIT_H 1

#include "DataFormats/HcalDetId/interface/HcalZDCDetId.h"
#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"

/** \class ZDCRecHit
 *  
 *\author J. Mans - Minnesota
 */
class ZDCRecHit : public CaloRecHit {
public:
  typedef HcalZDCDetId key_type;

  ZDCRecHit();
  ZDCRecHit(const HcalZDCDetId& id, float energy, float time, float lowGainEnergy);
  /// get the id
  HcalZDCDetId id() const { return HcalZDCDetId(detid()); }
  // follow EcalRecHit method of adding variable flagBits_ to CaloRecHit
  float lowGainEnergy() const { return lowGainEnergy_; };
  
  
  constexpr inline void setSaturated(const bool sat ) { saturated_ = sat;};
  constexpr inline bool saturated() const { return saturated_; };
  constexpr inline void setenergySOIp1(const float en) { energySOIp1_ = en;};
  constexpr inline float energySOIp1() const { return energySOIp1_;};
  constexpr inline void setTimeSOIp1(const float time) { timeSOIp1_ = time;};
  constexpr inline float timeSOIp1() const { return timeSOIp1_;};
  constexpr inline void setTDCtime(const float time) { TDCtime_ = time;};
  constexpr inline float TDCtime() const { return TDCtime_;};
  constexpr inline void setChargeWeightedTime(const float time) { chargeWeightedTime_ = time;};
  constexpr inline float chargeWeightedTime() const { return chargeWeightedTime_;}; 

 
private:
  float lowGainEnergy_;
  float saturated_;
  float energySOIp1_;
  float timeSOIp1_;
  float TDCtime_;
  float chargeWeightedTime_;
};

std::ostream& operator<<(std::ostream& s, const ZDCRecHit& hit);

#endif
