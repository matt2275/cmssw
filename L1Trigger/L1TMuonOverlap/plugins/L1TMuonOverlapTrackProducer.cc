#include <iostream>
#include <vector>

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/DataRecord/interface/L1TMuonOverlapParamsRcd.h"
#include "CondFormats/L1TObjects/interface/L1TMuonOverlapParams.h"

#include "L1Trigger/L1TMuonOverlap/plugins/L1TMuonOverlapTrackProducer.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFProcessor.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlap/interface/OMTFConfiguration.h"
#include "L1Trigger/L1TMuonOverlap/interface/XMLConfigWriter.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"

L1TMuonOverlapTrackProducer::L1TMuonOverlapTrackProducer(const edm::ParameterSet& cfg)
    : m_Reconstruction(cfg, consumesCollector()) {
  produces<l1t::RegionalMuonCandBxCollection>("OMTF");

  inputTokenDTPh = consumes<L1MuDTChambPhContainer>(cfg.getParameter<edm::InputTag>("srcDTPh"));
  inputTokenDTTh = consumes<L1MuDTChambThContainer>(cfg.getParameter<edm::InputTag>("srcDTTh"));
  inputTokenCSC = consumes<CSCCorrelatedLCTDigiCollection>(cfg.getParameter<edm::InputTag>("srcCSC"));
  inputTokenRPC = consumes<RPCDigiCollection>(cfg.getParameter<edm::InputTag>("srcRPC"));
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
L1TMuonOverlapTrackProducer::~L1TMuonOverlapTrackProducer() {}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapTrackProducer::beginJob() { m_Reconstruction.beginJob(); }
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapTrackProducer::endJob() { m_Reconstruction.endJob(); }
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapTrackProducer::beginRun(edm::Run const& run, edm::EventSetup const& iSetup) {
  m_Reconstruction.beginRun(run, iSetup);
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
void L1TMuonOverlapTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& evSetup) {
  std::ostringstream str;

  std::unique_ptr<l1t::RegionalMuonCandBxCollection> candidates = m_Reconstruction.reconstruct(iEvent, evSetup);

  iEvent.put(std::move(candidates), "OMTF");
}
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TMuonOverlapTrackProducer);
