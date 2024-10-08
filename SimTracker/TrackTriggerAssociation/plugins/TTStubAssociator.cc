/*! \brief   
 *  \details Here, in the source file, the methods which do depend
 *           on the specific type <T> that can fit the template.
 *
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 19
 *
 */

#include "SimTracker/TrackTriggerAssociation/plugins/TTStubAssociator.h"

/// Implement the producer
template <>
void TTStubAssociator<Ref_Phase2TrackerDigi_>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  /// Exit if real data
  if (iEvent.isRealData())
    return;

  /// Exit if the vectors are uncorrectly dimensioned
  if (ttClusterTruthInputTags_.size() != ttStubsInputTags_.size()) {
    edm::LogError("TTStubAsso ") << "E R R O R! the InputTag vectors have different size!";
    return;
  }

  int ncont1 = 0;

  const TrackerGeometry* const theTrackerGeom = &iSetup.getData(theTrackerGeometryToken_);
  const TrackerTopology* const tTopo = &iSetup.getData(theTrackerTopologyToken_);

  /// Loop over the InputTags to handle multiple collections

  for (const auto& iTag : ttStubsTokens_) {
    /// Prepare output
    auto associationMapForOutput = std::make_unique<TTStubAssociationMap<Ref_Phase2TrackerDigi_>>();

    /// Get the Stubs already stored away
    edm::Handle<TTStubDetSetVec> ttStubHandle;
    iEvent.getByToken(iTag, ttStubHandle);

    /// Get the Cluster MC truth
    edm::Handle<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>> ttClusterAssociationMapHandle;
    iEvent.getByToken(ttClusterTruthTokens_.at(ncont1), ttClusterAssociationMapHandle);

    /// Prepare the necessary maps
    std::map<TTStubRef, TrackingParticlePtr> stubToTrackingParticleMap;
    std::map<TrackingParticlePtr, std::vector<TTStubRef>> trackingParticleToStubVectorMap;

    /// Loop over the input Stubs

    if (not ttStubHandle->empty()) {
      for (const auto& gd : theTrackerGeom->dets()) {
        DetId detid = gd->geographicalId();
        if (detid.subdetId() != StripSubdetector::TOB && detid.subdetId() != StripSubdetector::TID)
          continue;  // only run on OT

        if (!tTopo->isLower(detid))
          continue;  // loop on the stacks: choose the lower arbitrarily

        DetId stackDetid = tTopo->stack(detid);  // Stub module detid

        if (ttStubHandle->find(stackDetid) == ttStubHandle->end())
          continue;

        /// Get the DetSets of the Clusters
        edmNew::DetSet<TTStub<Ref_Phase2TrackerDigi_>> stubs = (*ttStubHandle)[stackDetid];

        for (auto contentIter = stubs.begin(); contentIter != stubs.end(); ++contentIter) {
          /// Make the reference to be put in the map
          TTStubRef tempStubRef = edmNew::makeRefTo(ttStubHandle, contentIter);

          /// Fill the inclusive map which is careless of the stub classification
          for (unsigned int ic = 0; ic < 2; ic++) {
            const std::vector<TrackingParticlePtr>& tempTPs =
                ttClusterAssociationMapHandle->findTrackingParticlePtrs(tempStubRef->clusterRef(ic));

            for (const TrackingParticlePtr& testTP : tempTPs) {
              if (testTP.isNull())
                continue;

              /// Prepare the maps wrt TrackingParticle
              if (trackingParticleToStubVectorMap.find(testTP) == trackingParticleToStubVectorMap.end()) {
                std::vector<TTStubRef> stubVector;
                trackingParticleToStubVectorMap.emplace(testTP, stubVector);
              }
              trackingParticleToStubVectorMap.find(testTP)->second.push_back(tempStubRef);  /// Fill the auxiliary map
            }
          }

          /// GENUINE for clusters means not combinatoric and
          /// not unknown: same MC truth content MUST be found
          /// in both clusters composing the stub
          if (ttClusterAssociationMapHandle->isUnknown(tempStubRef->clusterRef(0)) ||
              ttClusterAssociationMapHandle->isUnknown(tempStubRef->clusterRef(1))) {
            /// If at least one cluster is unknown, it means
            /// either unknown, either combinatoric
            /// Do nothing, and go to the next Stub
            continue;
          } else {
            /// Here both are clusters are genuine/combinatoric
            /// If both clusters have some known SimTrack content
            /// they must be compared to each other
            if (ttClusterAssociationMapHandle->isGenuine(tempStubRef->clusterRef(0)) &&
                ttClusterAssociationMapHandle->isGenuine(tempStubRef->clusterRef(1))) {
              /// If both clusters are genuine, they must be associated to the same TrackingParticle
              /// in order to return a genuine stub. Period. Note we can perform safely
              /// this comparison because, if both clusters are genuine, their TrackingParticle shall NEVER be NULL
              if (ttClusterAssociationMapHandle->findTrackingParticlePtr(tempStubRef->clusterRef(0)).get() ==
                  ttClusterAssociationMapHandle->findTrackingParticlePtr(tempStubRef->clusterRef(1)).get()) {
                /// Two genuine clusters with same SimTrack content mean genuine
                const TrackingParticlePtr& testTP =
                    ttClusterAssociationMapHandle->findTrackingParticlePtr(tempStubRef->clusterRef(0));

                /// Fill the map: by construction, this will always be the first time the
                /// stub is inserted into the map: no need for "find"
                stubToTrackingParticleMap.emplace(tempStubRef, testTP);

                /// At this point, go to the next Stub
                continue;
              } else {
                /// It means combinatoric
                continue;
              }
            }  /// End of two genuine clusters
            else {
              /// Here, at least one cluster is combinatoric
              TrackingParticle* prevTPAddress = nullptr;
              unsigned int whichTP = 0;

              const std::vector<TrackingParticlePtr>& trackingParticles0 =
                  ttClusterAssociationMapHandle->findTrackingParticlePtrs(tempStubRef->clusterRef(0));
              const std::vector<TrackingParticlePtr>& trackingParticles1 =
                  ttClusterAssociationMapHandle->findTrackingParticlePtrs(tempStubRef->clusterRef(1));

              bool escape = false;

              for (unsigned int i = 0; i < trackingParticles0.size() && !escape; i++) {
                /// Skip NULL pointers
                if (trackingParticles0.at(i).isNull())
                  continue;

                for (unsigned int k = 0; k < trackingParticles1.size() && !escape; k++) {
                  /// Skip NULL pointers
                  if (trackingParticles1.at(k).isNull())
                    continue;

                  if (trackingParticles0.at(i).get() == trackingParticles1.at(k).get()) {
                    /// Same SimTrack is present in both clusters
                    if (prevTPAddress == nullptr) {
                      prevTPAddress = const_cast<TrackingParticle*>(trackingParticles1.at(k).get());
                      whichTP = k;
                    }

                    /// If two different SimTracks are found in both clusters,
                    /// then the stub is for sure combinatoric
                    if (prevTPAddress != const_cast<TrackingParticle*>(trackingParticles1.at(k).get())) {
                      escape = true;
                      continue;
                    }
                  }
                }
              }  /// End of double loop over SimTracks of both clusters

              /// If two different SimTracks are found in both clusters,
              /// go to the next stub
              if (escape)
                continue;

              if (prevTPAddress == nullptr) {
                /// No SimTracks were found to be in both clusters
                continue;
              } else {
                /// Only one SimTrack was found to be present in both clusters
                /// even if one of the clusters (or both) are combinatoric:
                /// this means there is only one track that participates in
                /// both clusters, hence the stub is genuine
                TrackingParticlePtr testTP = trackingParticles1.at(whichTP);

                /// Fill the map: by construction, this will always be the first time the
                /// stub is inserted into the map: no need for "find"
                stubToTrackingParticleMap.emplace(tempStubRef, testTP);

                /// At this point, go to the next Stub
                continue;
              }  /// End of one single SimTrack in both clusters
            }  /// End of "at least one cluster is combinatoric"
          }  /// End of "both clusters are known, somehow..."
        }
      }  /// End of loop over Stubs
    }

    /// Clean the only map that needs cleaning
    /// Prepare the output map wrt TrackingParticle
    for (auto& p : trackingParticleToStubVectorMap) {
      /// Get the vector of references to TTStub
      std::vector<TTStubRef>& tempVector = p.second;

      /// Sort and remove duplicates
      std::sort(tempVector.begin(), tempVector.end());
      tempVector.erase(std::unique(tempVector.begin(), tempVector.end()), tempVector.end());
    }

    /// Also, create the pointer to the TTClusterAssociationMap
    edm::RefProd<TTClusterAssociationMap<Ref_Phase2TrackerDigi_>> theCluAssoMap(ttClusterAssociationMapHandle);

    /// Put the maps in the association object
    associationMapForOutput->setTTStubToTrackingParticleMap(stubToTrackingParticleMap);
    associationMapForOutput->setTrackingParticleToTTStubsMap(trackingParticleToStubVectorMap);
    associationMapForOutput->setTTClusterAssociationMap(theCluAssoMap);

    /// Put output in the event
    iEvent.put(std::move(associationMapForOutput), ttStubsInputTags_.at(ncont1).instance());

    ++ncont1;
  }  /// End of loop over InputTags
}
