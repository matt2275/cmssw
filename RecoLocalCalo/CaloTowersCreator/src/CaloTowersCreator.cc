/** \class CaloTowersCreator
  *  
  * Original author: J. Mans - Minnesota
  */

// Now we allow for the creation of towers from
// rejected hists as well: requested by the MET group
// for studies of the effect of noise clean up.
#include "CaloTowersCreationAlgo.h"
#include "EScales.h"

#include "CommonTools/Utils/interface/StringToEnumValue.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTowerTopology.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "CondFormats/DataRecord/interface/HcalPFCutsRcd.h"
#include "CondTools/Hcal/interface/HcalPFCutsHandler.h"
#include "CondFormats/EcalObjects/interface/EcalPFRecHitThresholds.h"
#include "CondFormats/DataRecord/interface/EcalPFRecHitThresholdsRcd.h"

class CaloTowersCreator : public edm::stream::EDProducer<> {
public:
  explicit CaloTowersCreator(const edm::ParameterSet& ps);
  ~CaloTowersCreator() override {}
  void produce(edm::Event& e, const edm::EventSetup& c) override;
  void beginRun(edm::Run const&, edm::EventSetup const&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  double EBEScale, EEEScale, HBEScale, HESEScale;
  double HEDEScale, HOEScale, HF1EScale, HF2EScale;

private:
  static const std::vector<double>& getGridValues();

  CaloTowersCreationAlgo algo_;
  edm::EDGetTokenT<HBHERecHitCollection> tok_hbhe_;
  edm::EDGetTokenT<HORecHitCollection> tok_ho_;
  edm::EDGetTokenT<HFRecHitCollection> tok_hf_;
  std::vector<edm::InputTag> ecalLabels_;
  std::vector<edm::EDGetTokenT<EcalRecHitCollection> > toks_ecal_;
  bool allowMissingInputs_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> tok_geom_;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> tok_topo_;
  edm::ESGetToken<CaloTowerTopology, HcalRecNumberingRecord> tok_cttopo_;
  edm::ESGetToken<CaloTowerConstituentsMap, CaloGeometryRecord> tok_ctmap_;
  edm::ESGetToken<EcalChannelStatus, EcalChannelStatusRcd> tok_ecalChStatus_;
  edm::ESGetToken<HcalChannelQuality, HcalChannelQualityRcd> tok_hcalChStatus_;
  edm::ESGetToken<HcalSeverityLevelComputer, HcalSeverityLevelComputerRcd> tok_hcalSevComputer_;
  edm::ESGetToken<EcalSeverityLevelAlgo, EcalSeverityLevelAlgoRcd> tok_ecalSevAlgo_;

  // more compact flags: all HCAL are combined

  unsigned int theHcalAcceptSeverityLevel_;
  std::vector<int> theEcalSeveritiesToBeExcluded_;

  // flag to use recovered hits
  bool theRecoveredHcalHitsAreUsed_;
  bool theRecoveredEcalHitsAreUsed_;

  // paramaters for creating towers from rejected hits

  bool useRejectedHitsOnly_;
  unsigned int theHcalAcceptSeverityLevelForRejectedHit_;
  //  for ECAL we have a list of problem flags
  std::vector<int> theEcalSeveritiesToBeUsedInBadTowers_;

  // Flags wheteher to use recovered hits for production of
  // "bad towers".
  // If the recoverd hits were already used for good towers,
  // these flags have no effect.
  // Note: These flags have no effect on the default tower reconstruction.
  bool useRejectedRecoveredHcalHits_;
  bool useRejectedRecoveredEcalHits_;

  edm::ESWatcher<HcalSeverityLevelComputerRcd> hcalSevLevelWatcher_;
  edm::ESWatcher<HcalChannelQualityRcd> hcalChStatusWatcher_;
  edm::ESWatcher<IdealGeometryRecord> caloTowerConstituentsWatcher_;
  edm::ESWatcher<EcalSeverityLevelAlgoRcd> ecalSevLevelWatcher_;
  EScales eScales_;

  edm::ESGetToken<HcalPFCuts, HcalPFCutsRcd> hcalCutsToken_;
  bool cutsFromDB;
  HcalPFCuts const* paramPF = nullptr;

  // Ecal noise thresholds
  edm::ESGetToken<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd> ecalPFRechitThresholdsToken_;
  bool ecalRecHitThresh_;
  EcalPFRecHitThresholds const* ecalThresholds = nullptr;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(CaloTowersCreator);

CaloTowersCreator::CaloTowersCreator(const edm::ParameterSet& conf)
    : algo_(conf.getParameter<double>("EBThreshold"),
            conf.getParameter<double>("EEThreshold"),

            conf.getParameter<bool>("UseEtEBTreshold"),
            conf.getParameter<bool>("UseEtEETreshold"),
            conf.getParameter<bool>("UseSymEBTreshold"),
            conf.getParameter<bool>("UseSymEETreshold"),

            conf.getParameter<double>("HcalThreshold"),
            conf.getParameter<double>("HBThreshold"),
            conf.getParameter<double>("HBThreshold1"),
            conf.getParameter<double>("HBThreshold2"),
            conf.getParameter<double>("HESThreshold"),
            conf.getParameter<double>("HESThreshold1"),
            conf.getParameter<double>("HEDThreshold"),
            conf.getParameter<double>("HEDThreshold1"),
            conf.getParameter<double>("HOThreshold0"),
            conf.getParameter<double>("HOThresholdPlus1"),
            conf.getParameter<double>("HOThresholdMinus1"),
            conf.getParameter<double>("HOThresholdPlus2"),
            conf.getParameter<double>("HOThresholdMinus2"),
            conf.getParameter<double>("HF1Threshold"),
            conf.getParameter<double>("HF2Threshold"),
            conf.getParameter<std::vector<double> >("EBGrid"),
            conf.getParameter<std::vector<double> >("EBWeights"),
            conf.getParameter<std::vector<double> >("EEGrid"),
            conf.getParameter<std::vector<double> >("EEWeights"),
            conf.getParameter<std::vector<double> >("HBGrid"),
            conf.getParameter<std::vector<double> >("HBWeights"),
            conf.getParameter<std::vector<double> >("HESGrid"),
            conf.getParameter<std::vector<double> >("HESWeights"),
            conf.getParameter<std::vector<double> >("HEDGrid"),
            conf.getParameter<std::vector<double> >("HEDWeights"),
            conf.getParameter<std::vector<double> >("HOGrid"),
            conf.getParameter<std::vector<double> >("HOWeights"),
            conf.getParameter<std::vector<double> >("HF1Grid"),
            conf.getParameter<std::vector<double> >("HF1Weights"),
            conf.getParameter<std::vector<double> >("HF2Grid"),
            conf.getParameter<std::vector<double> >("HF2Weights"),
            conf.getParameter<double>("EBWeight"),
            conf.getParameter<double>("EEWeight"),
            conf.getParameter<double>("HBWeight"),
            conf.getParameter<double>("HESWeight"),
            conf.getParameter<double>("HEDWeight"),
            conf.getParameter<double>("HOWeight"),
            conf.getParameter<double>("HF1Weight"),
            conf.getParameter<double>("HF2Weight"),
            conf.getParameter<double>("EcutTower"),
            conf.getParameter<double>("EBSumThreshold"),
            conf.getParameter<double>("EESumThreshold"),
            conf.getParameter<bool>("UseHO"),
            // (for momentum reconstruction algorithm)
            conf.getParameter<int>("MomConstrMethod"),
            conf.getParameter<double>("MomHBDepth"),
            conf.getParameter<double>("MomHEDepth"),
            conf.getParameter<double>("MomEBDepth"),
            conf.getParameter<double>("MomEEDepth"),
            conf.getParameter<int>("HcalPhase")),

      ecalLabels_(conf.getParameter<std::vector<edm::InputTag> >("ecalInputs")),
      allowMissingInputs_(conf.getParameter<bool>("AllowMissingInputs")),

      theHcalAcceptSeverityLevel_(conf.getParameter<unsigned int>("HcalAcceptSeverityLevel")),

      theRecoveredHcalHitsAreUsed_(conf.getParameter<bool>("UseHcalRecoveredHits")),
      theRecoveredEcalHitsAreUsed_(conf.getParameter<bool>("UseEcalRecoveredHits")),

      // paramaters controlling the use of rejected hits

      useRejectedHitsOnly_(conf.getParameter<bool>("UseRejectedHitsOnly")),

      theHcalAcceptSeverityLevelForRejectedHit_(
          conf.getParameter<unsigned int>("HcalAcceptSeverityLevelForRejectedHit")),

      useRejectedRecoveredHcalHits_(conf.getParameter<bool>("UseRejectedRecoveredHcalHits")),
      useRejectedRecoveredEcalHits_(conf.getParameter<bool>("UseRejectedRecoveredEcalHits")),
      cutsFromDB(conf.getParameter<bool>("usePFThresholdsFromDB")),
      ecalRecHitThresh_(conf.getParameter<bool>("EcalRecHitThresh")) {
  algo_.setMissingHcalRescaleFactorForEcal(conf.getParameter<double>("missingHcalRescaleFactorForEcal"));

  // register for data access
  tok_hbhe_ = consumes<HBHERecHitCollection>(conf.getParameter<edm::InputTag>("hbheInput"));
  tok_ho_ = consumes<HORecHitCollection>(conf.getParameter<edm::InputTag>("hoInput"));
  tok_hf_ = consumes<HFRecHitCollection>(conf.getParameter<edm::InputTag>("hfInput"));
  tok_geom_ = esConsumes<CaloGeometry, CaloGeometryRecord>();
  tok_topo_ = esConsumes<HcalTopology, HcalRecNumberingRecord>();
  tok_cttopo_ = esConsumes<CaloTowerTopology, HcalRecNumberingRecord>();
  tok_ctmap_ = esConsumes<CaloTowerConstituentsMap, CaloGeometryRecord>();
  tok_ecalChStatus_ = esConsumes<EcalChannelStatus, EcalChannelStatusRcd>();
  tok_hcalChStatus_ = esConsumes<HcalChannelQuality, HcalChannelQualityRcd>(edm::ESInputTag("", "withTopo"));
  tok_hcalSevComputer_ = esConsumes<HcalSeverityLevelComputer, HcalSeverityLevelComputerRcd>();
  tok_ecalSevAlgo_ = esConsumes<EcalSeverityLevelAlgo, EcalSeverityLevelAlgoRcd>();

  if (cutsFromDB) {
    hcalCutsToken_ = esConsumes<HcalPFCuts, HcalPFCutsRcd, edm::Transition::BeginRun>(edm::ESInputTag("", "withTopo"));
  }

  if (ecalRecHitThresh_) {
    ecalPFRechitThresholdsToken_ =
        esConsumes<EcalPFRecHitThresholds, EcalPFRecHitThresholdsRcd, edm::Transition::BeginRun>();
  }

  const unsigned nLabels = ecalLabels_.size();
  for (unsigned i = 0; i != nLabels; i++)
    toks_ecal_.push_back(consumes<EcalRecHitCollection>(ecalLabels_[i]));

  EBEScale = eScales_.EBScale;
  EEEScale = eScales_.EEScale;
  HBEScale = eScales_.HBScale;
  HESEScale = eScales_.HESScale;
  HEDEScale = eScales_.HEDScale;
  HOEScale = eScales_.HOScale;
  HF1EScale = eScales_.HF1Scale;
  HF2EScale = eScales_.HF2Scale;

  // get the Ecal severities to be excluded
  const std::vector<std::string> severitynames =
      conf.getParameter<std::vector<std::string> >("EcalRecHitSeveritiesToBeExcluded");

  theEcalSeveritiesToBeExcluded_ = StringToEnumValue<EcalSeverityLevel::SeverityLevel>(severitynames);

  // get the Ecal severities to be used for bad towers
  theEcalSeveritiesToBeUsedInBadTowers_ = StringToEnumValue<EcalSeverityLevel::SeverityLevel>(
      conf.getParameter<std::vector<std::string> >("EcalSeveritiesToBeUsedInBadTowers"));

  if (eScales_.instanceLabel.empty())
    produces<CaloTowerCollection>();
  else
    produces<CaloTowerCollection>(eScales_.instanceLabel);

#ifdef EDM_ML_DEBUG
  std::cout << "VI Producer " << (useRejectedHitsOnly_ ? "use rejectOnly " : " ")
            << (allowMissingInputs_ ? "allowMissing " : " ") << nLabels << ' ' << severitynames.size() << std::endl;
#endif
}

void CaloTowersCreator::beginRun(const edm::Run& run, const edm::EventSetup& es) {
  if (cutsFromDB) {
    paramPF = &es.getData(hcalCutsToken_);
  }

  if (ecalRecHitThresh_) {
    ecalThresholds = &es.getData(ecalPFRechitThresholdsToken_);
  }

  algo_.setThresFromDB(ecalThresholds, paramPF);
}

void CaloTowersCreator::produce(edm::Event& e, const edm::EventSetup& c) {
  // get the necessary event setup objects...
  edm::ESHandle<CaloGeometry> pG = c.getHandle(tok_geom_);
  edm::ESHandle<HcalTopology> htopo = c.getHandle(tok_topo_);
  edm::ESHandle<CaloTowerTopology> cttopo = c.getHandle(tok_cttopo_);
  edm::ESHandle<CaloTowerConstituentsMap> ctmap = c.getHandle(tok_ctmap_);

  // ECAL channel status map ****************************************
  edm::ESHandle<EcalChannelStatus> ecalChStatus = c.getHandle(tok_ecalChStatus_);
  const EcalChannelStatus* dbEcalChStatus = ecalChStatus.product();

  // HCAL channel status map ****************************************
  edm::ESHandle<HcalChannelQuality> hcalChStatus = c.getHandle(tok_hcalChStatus_);

  const HcalChannelQuality* dbHcalChStatus = hcalChStatus.product();

  // Assignment of severity levels **********************************
  edm::ESHandle<HcalSeverityLevelComputer> hcalSevLvlComputerHndl = c.getHandle(tok_hcalSevComputer_);
  const HcalSeverityLevelComputer* hcalSevLvlComputer = hcalSevLvlComputerHndl.product();

  edm::ESHandle<EcalSeverityLevelAlgo> ecalSevLvlAlgoHndl = c.getHandle(tok_ecalSevAlgo_);
  const EcalSeverityLevelAlgo* ecalSevLvlAlgo = ecalSevLvlAlgoHndl.product();

  algo_.setEBEScale(EBEScale);
  algo_.setEEEScale(EEEScale);
  algo_.setHBEScale(HBEScale);
  algo_.setHESEScale(HESEScale);
  algo_.setHEDEScale(HEDEScale);
  algo_.setHOEScale(HOEScale);
  algo_.setHF1EScale(HF1EScale);
  algo_.setHF2EScale(HF2EScale);
  algo_.setGeometry(cttopo.product(), ctmap.product(), htopo.product(), pG.product());

  // for treatment of problematic and anomalous cells

  algo_.setHcalChStatusFromDB(dbHcalChStatus);
  algo_.setEcalChStatusFromDB(dbEcalChStatus);

  algo_.setHcalAcceptSeverityLevel(theHcalAcceptSeverityLevel_);
  algo_.setEcalSeveritiesToBeExcluded(theEcalSeveritiesToBeExcluded_);

  algo_.setRecoveredHcalHitsAreUsed(theRecoveredHcalHitsAreUsed_);
  algo_.setRecoveredEcalHitsAreUsed(theRecoveredEcalHitsAreUsed_);

  algo_.setHcalSevLvlComputer(hcalSevLvlComputer);
  algo_.setEcalSevLvlAlgo(ecalSevLvlAlgo);

  algo_.setUseRejectedHitsOnly(useRejectedHitsOnly_);

  algo_.setHcalAcceptSeverityLevelForRejectedHit(theHcalAcceptSeverityLevelForRejectedHit_);
  algo_.SetEcalSeveritiesToBeUsedInBadTowers(theEcalSeveritiesToBeUsedInBadTowers_);

  algo_.setUseRejectedRecoveredHcalHits(useRejectedRecoveredHcalHits_);
  algo_.setUseRejectedRecoveredEcalHits(useRejectedRecoveredEcalHits_);

#ifdef EDM_ML_DEBUG
  std::cout << "VI Produce: " << (useRejectedHitsOnly_ ? "use rejectOnly " : " ")
            << (allowMissingInputs_ ? "allowMissing " : " ")
            << (theRecoveredEcalHitsAreUsed_ ? "use RecoveredEcal " : " ") << toks_ecal_.size() << ' '
            << theEcalSeveritiesToBeExcluded_.size() << ' ' << theEcalSeveritiesToBeUsedInBadTowers_.size()
            << std::endl;
#endif

  algo_.begin();  // clear the internal buffer

  // can't chain these in a big OR statement, or else it'll
  // get triggered for each of the first three events
  bool check1 = hcalSevLevelWatcher_.check(c);
  bool check2 = hcalChStatusWatcher_.check(c);
  bool check3 = caloTowerConstituentsWatcher_.check(c);
  if (check1 || check2 || check3) {
    algo_.makeHcalDropChMap();
  }

  // check ecal SevLev
  if (ecalSevLevelWatcher_.check(c))
    algo_.makeEcalBadChs();

  // ----------------------------------------------------------
  // For ecal error handling need to
  // have access to the EB and EE collections at the end of
  // tower reconstruction.

  edm::Handle<EcalRecHitCollection> ebHandle;
  edm::Handle<EcalRecHitCollection> eeHandle;

  for (std::vector<edm::EDGetTokenT<EcalRecHitCollection> >::const_iterator i = toks_ecal_.begin();
       i != toks_ecal_.end();
       i++) {
    edm::Handle<EcalRecHitCollection> ec_tmp;

    if (!e.getByToken(*i, ec_tmp))
      continue;
    if (ec_tmp->empty())
      continue;

    // check if this is EB or EE
    if ((ec_tmp->begin()->detid()).subdetId() == EcalBarrel) {
      ebHandle = ec_tmp;
    } else if ((ec_tmp->begin()->detid()).subdetId() == EcalEndcap) {
      eeHandle = ec_tmp;
    }
  }

  algo_.setEbHandle(ebHandle);
  algo_.setEeHandle(eeHandle);

  //-----------------------------------------------------------

  bool present;

  // Step A/C: Get Inputs and process (repeatedly)
  edm::Handle<HBHERecHitCollection> hbhe;
  present = e.getByToken(tok_hbhe_, hbhe);
  if (present || !allowMissingInputs_)
    algo_.process(*hbhe);

  edm::Handle<HORecHitCollection> ho;
  present = e.getByToken(tok_ho_, ho);
  if (present || !allowMissingInputs_)
    algo_.process(*ho);

  edm::Handle<HFRecHitCollection> hf;
  present = e.getByToken(tok_hf_, hf);
  if (present || !allowMissingInputs_)
    algo_.process(*hf);

  std::vector<edm::EDGetTokenT<EcalRecHitCollection> >::const_iterator i;
  for (i = toks_ecal_.begin(); i != toks_ecal_.end(); i++) {
    edm::Handle<EcalRecHitCollection> ec;
    present = e.getByToken(*i, ec);
    if (present || !allowMissingInputs_)
      algo_.process(*ec);
  }

  // Step B: Create empty output
  auto prod = std::make_unique<CaloTowerCollection>();

  // Step C: Process
  algo_.finish(*prod);

#ifdef EDM_ML_DEBUG
  int totc = 0;
  float totE = 0;
  reco::LeafCandidate::LorentzVector totP4;
  for (auto const& tw : (*prod)) {
    totc += tw.constituents().size();
    totE += tw.energy();
    totP4 += tw.p4();
    std::cout << "CaloTowerCreator: " << tw.id() << " with E " << tw.energy() << " and " << tw.constituents().size()
              << " constituents\n";
  }
  std::cout << "VI " << (*prod).size() << " " << totc << " " << totE << " " << totP4 << std::endl;
#endif

  // Step D: Put into the event
  if (eScales_.instanceLabel.empty())
    e.put(std::move(prod));
  else
    e.put(std::move(prod), eScales_.instanceLabel);
}

void CaloTowersCreator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<double>("EBSumThreshold", 0.2);
  desc.add<double>("HF2Weight", 1.0);
  desc.add<double>("EBWeight", 1.0);
  desc.add<double>("EESumThreshold", 0.45);
  desc.add<double>("HOThreshold0", 1.1);
  desc.add<double>("HOThresholdPlus1", 3.5);
  desc.add<double>("HOThresholdMinus1", 3.5);
  desc.add<double>("HOThresholdPlus2", 3.5);
  desc.add<double>("HOThresholdMinus2", 3.5);
  desc.add<double>("HBThreshold", 0.7);
  desc.add<double>("HBThreshold1", 0.7);
  desc.add<double>("HBThreshold2", 0.7);
  desc.add<double>("HF1Threshold", 0.5);
  desc.add<double>("HEDWeight", 1.0);
  desc.add<double>("EEWeight", 1.0);
  desc.add<double>("HESWeight", 1.0);
  desc.add<double>("HF1Weight", 1.0);
  desc.add<double>("HOWeight", 1.0);
  desc.add<double>("EBThreshold", 0.07);
  desc.add<double>("EEThreshold", 0.3);
  desc.add<double>("HcalThreshold", -1000.0);
  desc.add<double>("HF2Threshold", 0.85);
  desc.add<double>("HESThreshold", 0.8);
  desc.add<double>("HESThreshold1", 0.8);
  desc.add<double>("HEDThreshold", 0.8);
  desc.add<double>("HEDThreshold1", 0.8);
  desc.add<double>("EcutTower", -1000.0);
  desc.add<double>("HBWeight", 1.0);
  desc.add<double>("MomHBDepth", 0.2);
  desc.add<double>("MomHEDepth", 0.4);
  desc.add<double>("MomEBDepth", 0.3);
  desc.add<double>("MomEEDepth", 0.0);
  desc.add<bool>("UseHO", true);
  desc.add<bool>("UseEtEBTreshold", false);
  desc.add<bool>("UseSymEBTreshold", true);
  desc.add<bool>("UseEtEETreshold", false);
  desc.add<bool>("UseSymEETreshold", true);
  desc.add<bool>("UseHcalRecoveredHits", true);
  desc.add<bool>("UseEcalRecoveredHits", false);
  desc.add<bool>("UseRejectedHitsOnly", false);
  desc.add<bool>("UseRejectedRecoveredHcalHits", true);
  desc.add<bool>("UseRejectedRecoveredEcalHits", false);
  desc.add<double>("missingHcalRescaleFactorForEcal", 0.0);
  desc.add<bool>("AllowMissingInputs", false);
  desc.add<std::vector<double> >("HBGrid", {-1.0, 1.0, 10.0, 100.0, 1000.0});
  desc.add<std::vector<double> >("EEWeights", {1.0, 1.0, 1.0, 1.0, 1.0});
  desc.add<std::vector<double> >("HF2Weights", {1.0, 1.0, 1.0, 1.0, 1.0});
  desc.add<std::vector<double> >("HOWeights", {1.0, 1.0, 1.0, 1.0, 1.0});
  desc.add<std::vector<double> >("EEGrid", {-1.0, 1.0, 10.0, 100.0, 1000.0});
  desc.add<std::vector<double> >("HBWeights", {1.0, 1.0, 1.0, 1.0, 1.0});
  desc.add<std::vector<double> >("HF2Grid", {-1.0, 1.0, 10.0, 100.0, 1000.0});
  desc.add<std::vector<double> >("HEDWeights", {1.0, 1.0, 1.0, 1.0, 1.0});
  desc.add<std::vector<double> >("HF1Grid", {-1.0, 1.0, 10.0, 100.0, 1000.0});
  desc.add<std::vector<double> >("EBWeights", {1.0, 1.0, 1.0, 1.0, 1.0});
  desc.add<std::vector<double> >("HF1Weights", {1.0, 1.0, 1.0, 1.0, 1.0});
  desc.add<std::vector<double> >("HESGrid", {-1.0, 1.0, 10.0, 100.0, 1000.0});
  desc.add<std::vector<double> >("HESWeights", {1.0, 1.0, 1.0, 1.0, 1.0});
  desc.add<std::vector<double> >("HEDGrid", {-1.0, 1.0, 10.0, 100.0, 1000.0});
  desc.add<std::vector<double> >("HOGrid", {-1.0, 1.0, 10.0, 100.0, 1000.0});
  desc.add<std::vector<double> >("EBGrid", {-1.0, 1.0, 10.0, 100.0, 1000.0});
  desc.add<edm::InputTag>("hfInput", edm::InputTag("hfreco"));
  desc.add<edm::InputTag>("hbheInput", edm::InputTag("hbhereco"));
  desc.add<edm::InputTag>("hoInput", edm::InputTag("horeco"));
  desc.add<std::vector<edm::InputTag> >(
      "ecalInputs", {edm::InputTag("ecalRecHit", "EcalRecHitsEB"), edm::InputTag("ecalRecHit", "EcalRecHitsEE")});
  desc.add<int>("MomConstrMethod", 1);
  desc.add<unsigned int>("HcalAcceptSeverityLevel", 9);
  desc.add<std::vector<std::string> >("EcalRecHitSeveritiesToBeExcluded", {"kTime", "kWeird", "kBad"});
  desc.add<unsigned int>("HcalAcceptSeverityLevelForRejectedHit", 9999);
  desc.add<std::vector<std::string> >("EcalSeveritiesToBeUsedInBadTowers", {});
  desc.add<int>("HcalPhase", 0);
  desc.add<bool>("usePFThresholdsFromDB", true);
  desc.add<bool>("EcalRecHitThresh", false);
  descriptions.addDefault(desc);
}
