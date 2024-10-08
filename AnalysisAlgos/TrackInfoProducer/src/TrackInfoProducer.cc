//
// Package:    RecoTracker/TrackInfoProducer
// Class:      TrackInfoProducer
//
//
// Description: Produce TrackInfo from Trajectory
//
//
// Original Author:  Chiara Genta
//         Created:
//

// system include files
#include <memory>

// user include files
#include "AnalysisAlgos/TrackInfoProducer/interface/TrackInfoProducerAlgorithm.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfo.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfoTrackAssociation.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"

class TrackInfoProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackInfoProducer(const edm::ParameterSet& iConfig);

  ~TrackInfoProducer() override {}

  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  TrackInfoProducerAlgorithm theAlgo_;
  edm::EDGetTokenT<std::vector<Trajectory> > TrajectoryToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackCollectionToken_;
  edm::EDGetTokenT<TrajTrackAssociationCollection> assoMapToken_;
  std::string forwardPredictedStateTag_, backwardPredictedStateTag_, updatedStateTag_, combinedStateTag_;
};

TrackInfoProducer::TrackInfoProducer(const edm::ParameterSet& iConfig)
    : tkGeomToken_(esConsumes()),
      theAlgo_(iConfig),
      TrajectoryToken_(consumes<std::vector<Trajectory> >(iConfig.getParameter<edm::InputTag>("cosmicTracks"))),
      trackCollectionToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("cosmicTracks"))),
      assoMapToken_(consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("cosmicTracks"))) {
  produces<reco::TrackInfoCollection>();
  produces<reco::TrackInfoTrackAssociationCollection>();
}

void TrackInfoProducer::produce(edm::Event& theEvent, const edm::EventSetup& setup) {
  //
  // create empty output collections
  //
  std::unique_ptr<reco::TrackInfoCollection> outputColl(new reco::TrackInfoCollection);

  const TrackerGeometry* tracker = &setup.getData(tkGeomToken_);

  edm::Handle<std::vector<Trajectory> > TrajectoryCollection;
  edm::Handle<reco::TrackCollection> trackCollection;
  edm::Handle<TrajTrackAssociationCollection> assoMap;

  theEvent.getByToken(TrajectoryToken_, TrajectoryCollection);
  theEvent.getByToken(trackCollectionToken_, trackCollection);
  theEvent.getByToken(assoMapToken_, assoMap);

  //
  //run the algorithm
  //
  reco::TrackInfo output;

  edm::LogInfo("TrackInfoProducer") << "Loop on trajectories";
  std::map<reco::TrackRef, unsigned int> trackid;
  int i = 0;

  for (TrajTrackAssociationCollection::const_iterator it = assoMap->begin(); it != assoMap->end(); ++it) {
    const edm::Ref<std::vector<Trajectory> > traj = it->key;
    const reco::TrackRef track = it->val;
    trackid.insert(make_pair(track, i));
    i++;
    theAlgo_.run(traj, track, output, tracker);
    outputColl->push_back(output);
  }

  //put everything in the event
  edm::OrphanHandle<reco::TrackInfoCollection> rTrackInfo;

  //     if(forwardPredictedStateTag_!="") rTrackInfof = theEvent.put(std::move(outputFwdColl),forwardPredictedStateTag_ );
  //     if(backwardPredictedStateTag_!="") rTrackInfob =   theEvent.put(std::move(outputBwdColl),backwardPredictedStateTag_);
  //     if(updatedStateTag_!="") rTrackInfou =   theEvent.put(std::move(outputUpdatedColl),updatedStateTag_ );
  //     if(combinedStateTag_!="") rTrackInfoc =   theEvent.put(std::move(outputCombinedColl),combinedStateTag_ );
  rTrackInfo = theEvent.put(std::move(outputColl));
  std::unique_ptr<reco::TrackInfoTrackAssociationCollection> TIassociationColl(
      new reco::TrackInfoTrackAssociationCollection(assoMap->refProd().val, rTrackInfo));

  for (std::map<reco::TrackRef, unsigned int>::iterator ref_iter = trackid.begin(); ref_iter != trackid.end();
       ++ref_iter) {
    TIassociationColl->insert(ref_iter->first, edm::Ref<reco::TrackInfoCollection>(rTrackInfo, ref_iter->second));
  }

  theEvent.put(std::move(TIassociationColl));
}

DEFINE_FWK_MODULE(TrackInfoProducer);
