#include "DataFormats/HcalRecHit/interface/ZDCRecHit.h"

ZDCRecHit::ZDCRecHit() : CaloRecHit(), lowGainEnergy_()
    ,saturated_(false) , energySOIp1_(-99) , timeSOIp1_(-99), TDCtime_(-99), chargeWeightedTime_(-99)    {}

ZDCRecHit::ZDCRecHit(const HcalZDCDetId& id, float energy, float time, float lowGainEnergy)
    : CaloRecHit(id, energy, time), lowGainEnergy_(lowGainEnergy)
    ,saturated_(false) , energySOIp1_(-99) , timeSOIp1_(-99), TDCtime_(-99), chargeWeightedTime_(-99)   {}

std::ostream& operator<<(std::ostream& s, const ZDCRecHit& hit) {
  return s << hit.id() << ": " << hit.energy() << " GeV, " << hit.time() << " ns";
}
