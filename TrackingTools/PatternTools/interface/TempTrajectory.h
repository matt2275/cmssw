#ifndef CkfPattern_TempTrajectory_H
#define CkfPattern_TempTrajectory_H

#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "DataFormats/TrackCandidate/interface/TrajectoryStopReasons.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "FWCore/Utilities/interface/Visibility.h"

#include <vector>
#include <algorithm>
#include <limits>
#include "TrackingTools/PatternTools/interface/bqueue.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"

/** A class for detailed particle trajectory representation.
 *  It is used during trajectory building to "grow" a trajectory.
 *  The trajectory is represented as an ordered sequence of 
 *  TrajectoryMeasurement objects with a stack-like interface.
 *  The measurements are added to the Trajectory in the order of
 *  increasing precision: each new TrajectoryMeasurement is assumed to improve
 *  the precision of the last one, normally by adding a constraint from 
 *  a new RecHit.
 *     However the Trajectory class does not have the means to verify
 *  that measurements are added in the correct order, and thus cannot
 *  guarantee the order, which is the responsibility of the 
 *  TrajectoryBuilder. The Trajectory provides some security by
 *  allowing to add or remove measurements only on one of it's ends,
 *  with push(TM) and pop() methods. The last measurement in a Trajectory
 *  can thus be either the innermost (closest to the interaction point)
 *  or the outermost, depending on the way the Trajectory was built.
 *  The direction of building is represented as a PropagationDirection,
 *  which has two possible values: alongMomentum (outwards) and
 *  oppositeToMomentum (inwards), and is accessed with the direction()
 *  method.
 */

class TempTrajectory {
public:
  typedef cmsutils::bqueue<TrajectoryMeasurement> DataContainer;
  typedef TrackingRecHit::ConstRecHitContainer ConstRecHitContainer;
  typedef ConstRecHitContainer RecHitContainer;

  struct Payload {
    float theChiSquared = 0;
    float theDPhiCache = 0;
    float theCCCThreshold_ = std::numeric_limits<float>::max();

    uint8_t theNumberOfFoundHits = 0;
    uint8_t theNumberOfFoundPixelHits = 0;
    uint8_t theNumberOfLostHits = 0;
    uint8_t theNumberOfTrailingFoundHits = 0;
    uint8_t theNumberOfCCCBadHits_ = 0;

    // PropagationDirection
    int8_t theDirection = anyDirection;

    uint8_t theNHseed = 0;
    uint8_t theNLoops = 0;
    StopReason stopReason_ = StopReason::UNINITIALIZED;
  };

  /** Default constructor of an empty trajectory with undefined seed and 
   * undefined direction. This constructor is necessary in order to transiently
   * copy vector<Trajectory> in the edm::Event
   */

  TempTrajectory() {}

  /** Constructor of an empty trajectory with defined direction.
   *  No check is made in the push method that measurements are
   *  added in the correct direction.
   */
  TempTrajectory(PropagationDirection dir, unsigned char nhseed) : thePayload(std::make_unique<Payload>()) {
    thePayload->theDirection = dir;
    thePayload->theNHseed = nhseed;
  }

  TempTrajectory(TempTrajectory const& rh)
      : theData(rh.theData), thePayload(rh.thePayload ? std::make_unique<Payload>(*rh.thePayload) : nullptr) {}

  TempTrajectory& operator=(TempTrajectory const& rh) {
    theData = rh.theData;
    thePayload = rh.thePayload ? std::make_unique<Payload>(*rh.thePayload) : nullptr;
    return *this;
  }

  TempTrajectory(TempTrajectory&& rh) noexcept : theData(std::move(rh.theData)), thePayload(std::move(rh.thePayload)) {}

  TempTrajectory& operator=(TempTrajectory&& rh) noexcept {
    using std::swap;
    swap(theData, rh.theData);
    swap(thePayload, rh.thePayload);
    return *this;
  }

  void swap(TempTrajectory& rh) noexcept {
    using std::swap;
    swap(theData, rh.theData);
    swap(thePayload, rh.thePayload);
  }

  /// construct TempTrajectory from standard Trajectory
  explicit TempTrajectory(Trajectory&& traj);

  /// destruct a TempTrajectory
  ~TempTrajectory() = default;

  /** Add a new measurement to a Trajectory.
   *  The Chi2 of the trajectory is incremented by the value
   *  of tm.estimate() . 
   */
  void push(const TrajectoryMeasurement& tm) { push(tm, tm.estimate()); }

  void push(TrajectoryMeasurement&& tm) { push(std::forward<TrajectoryMeasurement>(tm), tm.estimate()); }

  template <typename... Args>
  void emplace(Args&&... args) {
    theData.emplace_back(std::forward<Args>(args)...);
    pushAux(theData.back().estimate());
  }

  /** Add a new sets of measurements to a Trajectory
   *  The sorting of hits in the other trajectory must match the one
   *  inside this trajectory (that is, both along or both opposite to momentum)
   *  (the input segment will be reset to an empty one)
   */
  void push(TempTrajectory const& segment);

  /** Add a new sets of measurements to a Trajectory
   *  Exactly like push(TempTrajectory), but it doesn't copy the data
   *  (the input segment will be reset to an empty one)
   */
  void join(TempTrajectory& segment);

  /** same as the one-argument push, but the trajectory Chi2 is incremented 
   *  by chi2Increment. Useful e.g. in trajectory smoothing.
   */
  void push(const TrajectoryMeasurement& tm, double chi2Increment) {
    theData.push_back(tm);
    pushAux(chi2Increment);
  }

  void push(TrajectoryMeasurement&& tm, double chi2Increment) {
    theData.push_back(std::move(tm));
    pushAux(chi2Increment);
  }

  template <typename... Args>
  void emplace(double chi2Increment, Args&&... args) {  // works only because the first Arg is never a double!
    theData.emplace_back(std::forward<Args>(args)...);
    pushAux(chi2Increment);
  }

  /** Remove the last measurement from the trajectory.
   */
  void pop();

  /** Access to the last measurement.
   *  It's the most precise one in a trajectory before smoothing.
   *  It's the outermost measurement if direction() == alongMomentum,
   *  the innermost one if direction() == oppositeToMomentum.
   */
  const TrajectoryMeasurement& lastMeasurement() const {
    check();
    return theData.back();
  }

  /** Access to the first measurement.
   *  It is the least precise one in a trajectory before smoothing.
   *  It is precise in a smoothed trajectory. 
   *  It's the innermost measurement if direction() == alongMomentum,
   *  the outermost one if direction() == oppositeToMomentum.
   */
  const TrajectoryMeasurement& firstMeasurement() const {
    check();
    return theData.front();
  }

  /** Return all measurements in a container.
   */
  const DataContainer& measurements() const { return theData; }

  /** Number of valid RecHits used to determine the trajectory.
   *  Can be less than the number of measurements in data() since
   *  detector layers crossed without using RecHits from them are also 
   *  stored as measurements.
   */
  int foundHits() const { return thePayload->theNumberOfFoundHits; }

  /** Number of valid pixel RecHits used to determine the trajectory.
   */
  int foundPixelHits() const { return thePayload->theNumberOfFoundPixelHits; }

  /** Number of detector layers crossed without valid RecHits.
   *  Used mainly as a criteria for abandoning a trajectory candidate
   *  during trajectory building.
   */
  int lostHits() const { return thePayload->theNumberOfLostHits; }

  /** Number of valid RecHits at the end of the trajectory after last lost hit.
   */
  int trailingFoundHits() const { return thePayload->theNumberOfTrailingFoundHits; }

  /** Number of hits that are not compatible with the CCC used during
   *  patter recognition. Used mainly as a criteria for abandoning a
   *  trajectory candidate during trajectory building.
   */
  int cccBadHits() const { return thePayload->theNumberOfCCCBadHits_; }

  //number of hits in seed
  unsigned int seedNHits() const { return thePayload->theNHseed; }

  /// True if trajectory has no measurements.
  bool empty() const { return theData.empty(); }

  /// Value of the raw Chi2 of the trajectory, not normalised to the N.D.F.
  float chiSquared() const { return thePayload->theChiSquared; }

  /** Direction of "growing" of the trajectory. 
   *  Possible values are alongMomentum (outwards) and 
   *  oppositeToMomentum (inwards).
   */
  PropagationDirection direction() const;

  /** Returns true if the Trajectory is valid.
   *  Trajectories are invalidated e.g. during ambiguity resolution.
   */
  bool isValid() const { return bool(thePayload); }

  /// Method to invalidate a trajectory. Useful during ambiguity resolution.
  void invalidate() { thePayload.reset(); }

  /** Definition of inactive Det from the Trajectory point of view.
   */
  static bool inactive(  //const Det& det
  ) {
    return false;
  }  //FIXME

  /// Redundant method, returns the layer of lastMeasurement() .
  const DetLayer* lastLayer() const {
    check();
    return theData.back().layer();
  }

  /// Convert to a standard Trajectory
  Trajectory toTrajectory() const;

  /// Pops out all the invalid hits on the tail
  void popInvalidTail();

  /// accessor to the delta phi angle betweem the directions of the two measurements on the last
  /// two layers crossed by the trajectory
  float dPhiCacheForLoopersReconstruction() const { return thePayload->theDPhiCache; }

  /// method to set the delta phi angle betweem the directions of the two measurements on the last
  /// two layers crossed by the trajectory
  void setDPhiCacheForLoopersReconstruction(float dphi) { thePayload->theDPhiCache = dphi; }

  bool isLooper() const { return (thePayload->theNLoops > 0); }
  signed char nLoops() const { return thePayload->theNLoops; }

  void setNLoops(int8_t value) { thePayload->theNLoops = value; }
  void incrementLoops() { thePayload->theNLoops++; }

  StopReason stopReason() const { return thePayload->stopReason_; }
  void setStopReason(StopReason s) { thePayload->stopReason_ = s; }

  int numberOfCCCBadHits(float ccc_threshold);

  static bool lost(const TrackingRecHit& hit) dso_internal;

private:
  /** Definition of what it means for a hit to be "lost".
   *  This definition is also used by the TrajectoryBuilder.
   */
  bool badForCCC(const TrajectoryMeasurement& tm) dso_internal;
  void updateBadForCCC(float ccc_threshold) dso_internal;

  void pushAux(double chi2Increment);

private:
  DataContainer theData;
  std::unique_ptr<Payload> thePayload;

  void check() const;
};

#endif
