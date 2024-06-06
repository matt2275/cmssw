#include "ZdcHitReconstructorRunThree.h"
#include "DataFormats/Common/interface/EDCollection.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"

#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"

#include <iostream>

#include "CalibFormats/HcalObjects/interface/HcalCalibrationWidths.h"
#include "CondFormats/HcalObjects/interface/HcalPedestal.h"

/*  Zdc Hit reconstructor allows for CaloRecHits with status words */
ZdcHitReconstructorRunThree::ZdcHitReconstructorRunThree(edm::ParameterSet const& conf)

    : reco_(conf.getParameter<bool>("correctForTimeslew"),
            conf.getParameter<bool>("correctForPhaseContainment"),
            conf.getParameter<double>("correctionPhaseNS"),
            conf.getParameter<int>("recoMethod"),
            conf.getParameter<int>("lowGainOffset"),
            conf.getParameter<double>("lowGainFrac")),
      saturationFlagSetter_(nullptr),
      det_(DetId::Hcal),
      correctTiming_(conf.getParameter<bool>("correctTiming")),
      setNoiseFlags_(conf.getParameter<bool>("setNoiseFlags")),
      setHSCPFlags_(conf.getParameter<bool>("setHSCPFlags")),
      setSaturationFlags_(conf.getParameter<bool>("setSaturationFlags")),
      setTimingTrustFlags_(conf.getParameter<bool>("setTimingTrustFlags")),
      dropZSmarkedPassed_(conf.getParameter<bool>("dropZSmarkedPassed")),
      ignoreRPD_(conf.getParameter<bool>("ignoreRPD")),
      matchTrigger_(conf.getParameter<bool>("matchTrigger")),
      AuxTSvec_(conf.getParameter<std::vector<int> >("AuxTSvec")) {
         
  // commented out hcal and  castor digis since not used in RunThree ZDC    
  // tok_input_hcal = consumes<ZDCDigiCollection>(conf.getParameter<edm::InputTag>("digiLabelhcal"));
  // tok_input_castor = consumes<ZDCDigiCollection>(conf.getParameter<edm::InputTag>("digiLabelcastor"));
  tok_input_QIE10 = consumes<QIE10DigiCollection>(conf.getParameter<edm::InputTag>("digiLabelQIE10ZDC"));

  std::sort(AuxTSvec_.begin(), AuxTSvec_.end());  // sort vector in ascending TS order
  std::string subd = conf.getParameter<std::string>("Subdetector");

  if (setSaturationFlags_) {
    const edm::ParameterSet& pssat = conf.getParameter<edm::ParameterSet>("saturationParameters");
    saturationFlagSetter_ = new HcalADCSaturationFlag(pssat.getParameter<int>("maxADCvalue"));
  }
  if (!strcasecmp(subd.c_str(), "ZDC")) {
    det_ = DetId::Calo;
    subdet_ = HcalZDCDetId::SubdetectorId;
    produces<ZDCRecHitCollection>();
  } else if (!strcasecmp(subd.c_str(), "CALIB")) {
    subdet_ = HcalOther;
    subdetOther_ = HcalCalibration;
    produces<HcalCalibRecHitCollection>();
  } else {
    std::cout << "ZdcHitReconstructorRunThree is not associated with a specific subdetector!" << std::endl;
  }

  // ES tokens
  htopoToken_ = esConsumes<HcalTopology, HcalRecNumberingRecord, edm::Transition::BeginRun>();
  paramsToken_ = esConsumes<HcalLongRecoParams, HcalLongRecoParamsRcd, edm::Transition::BeginRun>();
  conditionsToken_ = esConsumes<HcalDbService, HcalDbRecord>();
  qualToken_ = esConsumes<HcalChannelQuality, HcalChannelQualityRcd>(edm::ESInputTag("", "withTopo"));
  sevToken_ = esConsumes<HcalSeverityLevelComputer, HcalSeverityLevelComputerRcd>();
}

ZdcHitReconstructorRunThree::~ZdcHitReconstructorRunThree() { delete saturationFlagSetter_; }

void ZdcHitReconstructorRunThree::beginRun(edm::Run const& r, edm::EventSetup const& es) {
  const HcalTopology& htopo = es.getData(htopoToken_);
  const HcalLongRecoParams& p = es.getData(paramsToken_);
  longRecoParams_ = std::make_unique<HcalLongRecoParams>(p);
  longRecoParams_->setTopo(&htopo);
}

void ZdcHitReconstructorRunThree::endRun(edm::Run const& r, edm::EventSetup const& es) {}

void ZdcHitReconstructorRunThree::produce(edm::Event& e, const edm::EventSetup& eventSetup) {
  // get conditions
  const HcalDbService* conditions = &eventSetup.getData(conditionsToken_);
  const HcalChannelQuality* myqual = &eventSetup.getData(qualToken_);
  const HcalSeverityLevelComputer* mySeverity = &eventSetup.getData(sevToken_);

  // define vectors to pass noiseTS and signalTS
  std::vector<unsigned int> mySignalTS;
  std::vector<unsigned int> myNoiseTS;

  if (det_ == DetId::Calo && subdet_ == HcalZDCDetId::SubdetectorId) {
    edm::Handle<QIE10DigiCollection> digi;
    e.getByToken(tok_input_QIE10, digi);
   


    // create empty output
    auto rec = std::make_unique<ZDCRecHitCollection>();
    rec->reserve(digi->size());
    
    
    // testing QEI10 conditions 
        for (auto it = digi->begin(); it != digi->end(); it++) {
      QIE10DataFrame QIE10_i = static_cast<QIE10DataFrame>(*it);     
      HcalZDCDetId cell = QIE10_i.id();
      if(cell.section() != 1 && cell.section()!= 2 && ignoreRPD_) continue;
      if(cell.section() == 1 && cell.channel()> 5) continue; // ignore extra EM channels 
      
      DetId detcell = (DetId)cell; // temporarly removed to avoid issue with RPD
      
      // check on cells to be ignored and dropped: (rof,20.Feb.09)
      const HcalChannelStatus* mydigistatus = myqual->getValues(detcell.rawId());
      if (mySeverity->dropChannel(mydigistatus->getValue()))
        continue;
      if (dropZSmarkedPassed_)
        if (QIE10_i.zsMarkAndPass())
          continue;
       
      const HcalCalibrations& calibrations = conditions->getHcalCalibrations(cell);
      const HcalQIECoder* channelCoder = conditions->getHcalCoder(cell);
      const HcalQIEShape* shape = conditions->getHcalShape(channelCoder);
      HcalCoderDb coder(*channelCoder, *shape);



      const HcalLongRecoParam* myParams = longRecoParams_->getValues(detcell);
      mySignalTS.clear();
      myNoiseTS.clear();
      mySignalTS = myParams->signalTS();
      myNoiseTS = myParams->noiseTS();
      
      if(matchTrigger_){
        const HcalPedestal* Peds = conditions->getPedestal(cell);
        const HcalPedestal* effPeds = conditions->getEffectivePedestal(cell);
        // for(int capid =0; capid<4; capid++) std::cout<<calibWidths.pedestal(capid)<<std::endl;
        rec->push_back(reco_.reconstruct(QIE10_i, myNoiseTS, mySignalTS, coder, calibrations,*effPeds)); 
      }
      
      else rec->push_back(reco_.reconstruct(QIE10_i, myNoiseTS, mySignalTS, coder, calibrations));
      
      
      //// saturationFlagSetter_ doesn't work with QIE10 so saturation is set in ZDCRecAlgo
      // (rec->back()).setFlags(0);
      // if (setSaturationFlags_)
        // saturationFlagSetter_->setSaturationFlag(rec->back(), QIE10_i);

     }
    // return result
    e.put(std::move(rec));
  }  // else if (det_==DetId::Calo...)

}  // void HcalHitReconstructor::produce(...)

//define this as a plug-in
DEFINE_FWK_MODULE(ZdcHitReconstructorRunThree);