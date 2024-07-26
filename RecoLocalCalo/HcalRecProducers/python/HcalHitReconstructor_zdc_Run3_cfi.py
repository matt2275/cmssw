import FWCore.ParameterSet.Config as cms

zdcreco = cms.EDProducer(
    "ZdcHitReconstructor_Run3",
    digiLabelQIE10ZDC = cms.InputTag("hcalDigis:ZDC"),
    Subdetector = cms.string('ZDC'),
    dropZSmarkedPassed = cms.bool(True),
    skipRPD = cms.bool(False),
    recoMethod = cms.int32(1),
    correctionMethodEM = cms.int32(1),
    correctionMethodHAD = cms.int32(1),
    correctionMethodRPD = cms.int32(0),
    ootpuRatioEM = cms.double(3.0),
    ootpuRatioHAD = cms.double(3.0),
    ootpuRatioRPD = cms.double(-1.0),
    ootpuFracEM = cms.double(.4),
    ootpuFracHAD = cms.double(.4),
    ootpuFracRPD = cms.double(1.0),
    chargeRatiosEM = cms.vdouble(1.0,.23157, .10477, .06312),
    chargeRatiosHAD = cms.vdouble(1.0,.23157, .10477, .06312),
    chargeRatiosRPD = cms.vdouble(1.0,.23157, .10477, .06312),
    bxTs = cms.vuint32(0,2,4),
    nTs = cms.int32(6),
    forceSOI = cms.bool(False),
    signalSOI = cms.vuint32(2),    
    noiseSOI = cms.vuint32(1),    

    setSaturationFlags = cms.bool(True),    
    saturationParameters=  cms.PSet(maxADCvalue=cms.int32(255))
    ) # zdcreco

zdcreco2023Trigger = cms.EDProducer(
    "ZdcHitReconstructor_Run3",
    digiLabelQIE10ZDC = cms.InputTag("hcalDigis:ZDC"),
    Subdetector = cms.string('ZDC'),
    dropZSmarkedPassed = cms.bool(True),
    skipRPD = cms.bool(True),
    recoMethod = cms.int32(1),
    correctionMethodEM = cms.int32(0),
    correctionMethodHAD = cms.int32(0),
    correctionMethodRPD = cms.int32(0),
    ootpuRatioEM = cms.double(-1.0),
    ootpuRatioHAD = cms.double(-1.0),
    ootpuRatioRPD = cms.double(-1.0),
    ootpuFracEM = cms.double(.3789), # fracion used was 97/256
    ootpuFracHAD = cms.double(.3798), 
    ootpuFracRPD = cms.double(1.0),
    chargeRatiosEM = cms.vdouble(1.0,.23157, .10477, .06312),
    chargeRatiosHAD = cms.vdouble(1.0,.23157, .10477, .06312),
    chargeRatiosRPD = cms.vdouble(1.0,.23157, .10477, .06312),
    bxTs = cms.vuint32(0,2,4),
    nTs = cms.int32(6),
    forceSOI = cms.bool(True),
    signalSOI = cms.vuint32(2),    
    noiseSOI = cms.vuint32(1),    

    setSaturationFlags = cms.bool(True),    
    saturationParameters=  cms.PSet(maxADCvalue=cms.int32(255))
    ) # zdcreco
