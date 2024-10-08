#ifndef Geometry_ForwardGeometry_ZdcHardcodeGeometryLoader_H
#define Geometry_ForwardGeometry_ZdcHardcodeGeometryLoader_H 1

#include "Geometry/CaloGeometry/interface/CaloVGeometryLoader.h"
#include "Geometry/ForwardGeometry/interface/ZdcTopology.h"

class CaloCellGeometry;
class CaloSubdetectorGeometry;
class HcalZDCDetId;

/** \class ZdcHardcodeGeometryLoader
 *
 * \author E. Garcia - UIC
*/

class ZdcHardcodeGeometryLoader {
public:
  typedef CaloSubdetectorGeometry* ReturnType;

  explicit ZdcHardcodeGeometryLoader(const ZdcTopology& ht);
  ~ZdcHardcodeGeometryLoader() {}

  ReturnType load(DetId::Detector det, int subdet);
  ReturnType load();
  void setAddRPD(bool flag) { m_zdcAddRPD = flag; }

private:
  void init();
  void fill(HcalZDCDetId::Section section, CaloSubdetectorGeometry* cg);
  void makeCell(const HcalZDCDetId& detId, CaloSubdetectorGeometry* geom) const;

  //  ZdcTopology* theTopology;
  const ZdcTopology* extTopology;
  bool m_zdcAddRPD;
  float theEMSectiondX;
  float theEMSectiondY;
  float theEMSectiondZ;
  float theLUMSectiondX;
  float theLUMSectiondY;
  float theLUMSectiondZ;
  float theHADSectiondX;
  float theHADSectiondY;
  float theHADSectiondZ;
  float theRPDSectiondX;
  float theRPDSectiondY;
  float theRPDSectiondZ;
};

#endif
