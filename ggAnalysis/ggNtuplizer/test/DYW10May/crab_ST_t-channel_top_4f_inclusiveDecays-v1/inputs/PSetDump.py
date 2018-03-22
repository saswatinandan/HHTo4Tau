import FWCore.ParameterSet.Config as cms

process = cms.Process("ggKit")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)
process.CondDB = cms.PSet(
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('')
)

process.CondDBCommon = cms.PSet(
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('')
)

process.METSignificanceParams = cms.PSet(
    dRMatch = cms.double(0.4),
    jetThreshold = cms.double(15),
    jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
    jpar = cms.vdouble(1.29, 1.19, 1.07, 1.13, 1.12),
    pjpar = cms.vdouble(-0.04, 0.6504)
)

process.METSignificance_params = cms.PSet(
    EB_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    EB_PhiResPar = cms.vdouble(0.00502),
    EE_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    EE_PhiResPar = cms.vdouble(0.02511),
    HB_EtResPar = cms.vdouble(0.0, 1.22, 0.05),
    HB_PhiResPar = cms.vdouble(0.02511),
    HE_EtResPar = cms.vdouble(0.0, 1.3, 0.05),
    HE_PhiResPar = cms.vdouble(0.02511),
    HF_EtResPar = cms.vdouble(0.0, 1.82, 0.09),
    HF_PhiResPar = cms.vdouble(0.05022),
    HO_EtResPar = cms.vdouble(0.0, 1.3, 0.005),
    HO_PhiResPar = cms.vdouble(0.02511),
    PF_EtResType1 = cms.vdouble(0.05, 0, 0),
    PF_EtResType2 = cms.vdouble(0.05, 0, 0),
    PF_EtResType3 = cms.vdouble(0.05, 0, 0),
    PF_EtResType4 = cms.vdouble(0.042, 0.1, 0.0),
    PF_EtResType5 = cms.vdouble(0.41, 0.52, 0.25),
    PF_EtResType6 = cms.vdouble(0.0, 1.22, 0.05),
    PF_EtResType7 = cms.vdouble(0.0, 1.22, 0.05),
    PF_PhiResType1 = cms.vdouble(0.002),
    PF_PhiResType2 = cms.vdouble(0.002),
    PF_PhiResType3 = cms.vdouble(0.002),
    PF_PhiResType4 = cms.vdouble(0.0028, 0.0, 0.0022),
    PF_PhiResType5 = cms.vdouble(0.1, 0.1, 0.13),
    PF_PhiResType6 = cms.vdouble(0.02511),
    PF_PhiResType7 = cms.vdouble(0.02511),
    jdphi0 = cms.vdouble(0.034, 0.034, 0.034, 0.034, 0.032, 
        0.031, 0.028, 0.027, 0.027, 0.027),
    jdphi1 = cms.vdouble(0.034, 0.035, 0.035, 0.035, 0.035, 
        0.034, 0.031, 0.03, 0.029, 0.027),
    jdphi2 = cms.vdouble(0.04, 0.04, 0.04, 0.04, 0.04, 
        0.038, 0.036, 0.035, 0.034, 0.033),
    jdphi3 = cms.vdouble(0.042, 0.043, 0.044, 0.043, 0.041, 
        0.039, 0.039, 0.036, 0.034, 0.031),
    jdphi4 = cms.vdouble(0.042, 0.042, 0.043, 0.042, 0.038, 
        0.036, 0.036, 0.033, 0.031, 0.031),
    jdphi5 = cms.vdouble(0.069, 0.069, 0.064, 0.058, 0.053, 
        0.049, 0.049, 0.043, 0.039, 0.04),
    jdphi6 = cms.vdouble(0.084, 0.08, 0.072, 0.065, 0.066, 
        0.06, 0.051, 0.049, 0.045, 0.045),
    jdphi7 = cms.vdouble(0.077, 0.072, 0.059, 0.05, 0.045, 
        0.042, 0.039, 0.039, 0.037, 0.031),
    jdphi8 = cms.vdouble(0.059, 0.057, 0.051, 0.044, 0.038, 
        0.035, 0.037, 0.032, 0.028, 0.028),
    jdphi9 = cms.vdouble(0.062, 0.059, 0.053, 0.047, 0.042, 
        0.045, 0.036, 0.032, 0.034, 0.044),
    jdpt0 = cms.vdouble(0.749, 0.829, 1.099, 1.355, 1.584, 
        1.807, 2.035, 2.217, 2.378, 2.591),
    jdpt1 = cms.vdouble(0.718, 0.813, 1.133, 1.384, 1.588, 
        1.841, 2.115, 2.379, 2.508, 2.772),
    jdpt2 = cms.vdouble(0.841, 0.937, 1.316, 1.605, 1.919, 
        2.295, 2.562, 2.722, 2.943, 3.293),
    jdpt3 = cms.vdouble(0.929, 1.04, 1.46, 1.74, 2.042, 
        2.289, 2.639, 2.837, 2.946, 2.971),
    jdpt4 = cms.vdouble(0.85, 0.961, 1.337, 1.593, 1.854, 
        2.005, 2.209, 2.533, 2.812, 3.047),
    jdpt5 = cms.vdouble(1.049, 1.149, 1.607, 1.869, 2.012, 
        2.219, 2.289, 2.412, 2.695, 2.865),
    jdpt6 = cms.vdouble(1.213, 1.298, 1.716, 2.015, 2.191, 
        2.612, 2.863, 2.879, 2.925, 2.902),
    jdpt7 = cms.vdouble(1.094, 1.139, 1.436, 1.672, 1.831, 
        2.05, 2.267, 2.549, 2.785, 2.86),
    jdpt8 = cms.vdouble(0.889, 0.939, 1.166, 1.365, 1.553, 
        1.805, 2.06, 2.22, 2.268, 2.247),
    jdpt9 = cms.vdouble(0.843, 0.885, 1.245, 1.665, 1.944, 
        1.981, 1.972, 2.875, 3.923, 7.51),
    ptresolthreshold = cms.double(10.0),
    resolutionsAlgo = cms.string('AK5PF'),
    resolutionsEra = cms.string('Spring10')
)

process.combinedSecondaryVertexCommon = cms.PSet(
    SoftLeptonFlip = cms.bool(False),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    trackFlip = cms.bool(False),
    trackMultiplicityMin = cms.uint32(2),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(False)
)

process.ghostTrackCommon = cms.PSet(
    charmCut = cms.double(1.5),
    minimumTrackWeight = cms.double(0.5),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig')
)

process.hggPhotonIDCuts = cms.PSet(
**dict(
    [
        ("cutsublead_drtotk_25_990" , cms.vdouble(0.26, 0.029, 0.0062, 0.0055) ),
        ("cutsublead_drtotk_25_991" , cms.vdouble(0.98, 0.029, 0.0109, 0.0111) ),
        ("cutsublead_drtotk_25_9910" , cms.vdouble(1, 1, 1, 1) ),
        ("cutsublead_drtotk_25_9911" , cms.vdouble(1, 1, 1, 1) ),
        ("cutsublead_drtotk_25_992" , cms.vdouble(1, 0.029, 0.021, 0.028) ),
        ("cutsublead_drtotk_25_993" , cms.vdouble(1, 0.062, 0.97, 0.97) ),
        ("cutsublead_drtotk_25_994" , cms.vdouble(1, 0.97, 1, 1) ),
        ("cutsublead_drtotk_25_995" , cms.vdouble(1, 1, 1, 1) ),
        ("cutsublead_drtotk_25_996" , cms.vdouble(1, 1, 1, 1) ),
        ("cutsublead_drtotk_25_996c0" , cms.vdouble(0.32, 0.99, 0.0097, 0.0064, 0.0048, 
        0.0061) ),
        ("cutsublead_drtotk_25_996c1" , cms.vdouble(0.98, 1, 0.012, 0.0111, 0.0086, 
        0.0122) ),
        ("cutsublead_drtotk_25_996c10" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsublead_drtotk_25_996c11" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsublead_drtotk_25_996c2" , cms.vdouble(1, 1, 0.0146, 0.021, 0.36, 
        0.026) ),
        ("cutsublead_drtotk_25_996c3" , cms.vdouble(1, 1, 0.022, 0.97, 0.98, 
        0.97) ),
        ("cutsublead_drtotk_25_996c4" , cms.vdouble(1, 1, 0.091, 1, 1, 
        1) ),
        ("cutsublead_drtotk_25_996c5" , cms.vdouble(1, 1, 0.97, 1, 1, 
        1) ),
        ("cutsublead_drtotk_25_996c6" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsublead_drtotk_25_996c7" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsublead_drtotk_25_996c8" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsublead_drtotk_25_996c9" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsublead_drtotk_25_997" , cms.vdouble(1, 1, 1, 1) ),
        ("cutsublead_drtotk_25_998" , cms.vdouble(1, 1, 1, 1) ),
        ("cutsublead_drtotk_25_999" , cms.vdouble(1, 1, 1, 1) ),
        ("cutsubleadhovere0" , cms.vdouble(0.09, 0.089, 0.101, 0.073) ),
        ("cutsubleadhovere1" , cms.vdouble(0.089, 0.079, 0.09, 0.061) ),
        ("cutsubleadhovere10" , cms.vdouble(0.0002, 0, 1.6e-05, 0) ),
        ("cutsubleadhovere11" , cms.vdouble(99.0, 99.0, 99.0, 99.0) ),
        ("cutsubleadhovere2" , cms.vdouble(0.087, 0.065, 0.087, 0.05) ),
        ("cutsubleadhovere3" , cms.vdouble(0.082, 0.062, 0.065, 0.048) ),
        ("cutsubleadhovere4" , cms.vdouble(0.076, 0.03, 0.047, 0.046) ),
        ("cutsubleadhovere5" , cms.vdouble(0.048, 0.0189, 0.032, 0.0085) ),
        ("cutsubleadhovere6" , cms.vdouble(0.042, 0.0173, 0.023, 0.0085) ),
        ("cutsubleadhovere6c0" , cms.vdouble(0.091, 0.079, 0.084, 0.103, 0.06, 
        0.099) ),
        ("cutsubleadhovere6c1" , cms.vdouble(0.089, 0.076, 0.08, 0.091, 0.051, 
        0.063) ),
        ("cutsubleadhovere6c10" , cms.vdouble(0.026, 0.0082, 3.2e-05, 3.2e-05, 0.0073, 
        0.0083) ),
        ("cutsubleadhovere6c11" , cms.vdouble(99.0, 99.0, 99.0, 99.0, 99.0, 
        99.0) ),
        ("cutsubleadhovere6c2" , cms.vdouble(0.083, 0.07, 0.061, 0.088, 0.05, 
        0.061) ),
        ("cutsubleadhovere6c3" , cms.vdouble(0.082, 0.06, 0.049, 0.077, 0.05, 
        0.061) ),
        ("cutsubleadhovere6c4" , cms.vdouble(0.068, 0.06, 0.027, 0.055, 0.047, 
        0.047) ),
        ("cutsubleadhovere6c5" , cms.vdouble(0.058, 0.023, 0.0174, 0.055, 0.047, 
        0.047) ),
        ("cutsubleadhovere6c6" , cms.vdouble(0.038, 0.022, 0.014, 0.043, 0.0116, 
        0.047) ),
        ("cutsubleadhovere6c7" , cms.vdouble(0.032, 0.0154, 0.0111, 0.023, 0.0107, 
        0.0138) ),
        ("cutsubleadhovere6c8" , cms.vdouble(0.026, 0.0137, 0.0111, 0.0192, 0.0073, 
        0.0083) ),
        ("cutsubleadhovere6c9" , cms.vdouble(0.026, 0.0137, 0.00114, 0.00068, 0.0073, 
        0.0083) ),
        ("cutsubleadhovere7" , cms.vdouble(0.037, 0.00049, 0.0198, 0.00024) ),
        ("cutsubleadhovere8" , cms.vdouble(0.028, 1.4e-05, 0.0198, 7e-06) ),
        ("cutsubleadhovere9" , cms.vdouble(0.0071, 0, 0.00055, 0) ),
        ("cutsubleadisosumoet0" , cms.vdouble(8.2, 4.1, 5.4, 2.6) ),
        ("cutsubleadisosumoet1" , cms.vdouble(6.4, 3.2, 3.4, 2.2) ),
        ("cutsubleadisosumoet10" , cms.vdouble(0.146, 0.058, -0.04, -1.16) ),
        ("cutsubleadisosumoet11" , cms.vdouble(7.0, 7.0, 4.0, 4.0) ),
        ("cutsubleadisosumoet2" , cms.vdouble(4.7, 2.8, 2.5, 1.46) ),
        ("cutsubleadisosumoet3" , cms.vdouble(3.8, 2.2, 1.77, 1.29) ),
        ("cutsubleadisosumoet4" , cms.vdouble(3.2, 1.76, 1.39, 1.18) ),
        ("cutsubleadisosumoet5" , cms.vdouble(2.6, 1.31, 1.33, 0.82) ),
        ("cutsubleadisosumoet6" , cms.vdouble(1.85, 0.96, 1.21, -0.029) ),
        ("cutsubleadisosumoet6c0" , cms.vdouble(14.1, 11.7, 9.8, 11, 9.2, 
        8.9) ),
        ("cutsubleadisosumoet6c1" , cms.vdouble(12.4, 10.2, 9.2, 10, 8.3, 
        8.7) ),
        ("cutsubleadisosumoet6c10" , cms.vdouble(5.2, 5.6, 4.8, 4.4, 4.6, 
        4.5) ),
        ("cutsubleadisosumoet6c11" , cms.vdouble(7.0, 7.0, 7.0, 4.0, 4.0, 
        4.0) ),
        ("cutsubleadisosumoet6c2" , cms.vdouble(11.2, 9.3, 9, 9.2, 7.9, 
        8.2) ),
        ("cutsubleadisosumoet6c3" , cms.vdouble(10, 8.8, 8.7, 8.5, 7.9, 
        7.5) ),
        ("cutsubleadisosumoet6c4" , cms.vdouble(9.3, 8.1, 7.5, 7.9, 6.3, 
        6.6) ),
        ("cutsubleadisosumoet6c5" , cms.vdouble(8.6, 6.6, 6.9, 7.1, 6, 
        6.3) ),
        ("cutsubleadisosumoet6c6" , cms.vdouble(8, 6.6, 6.5, 5.3, 4.9, 
        5.7) ),
        ("cutsubleadisosumoet6c7" , cms.vdouble(6.8, 6.1, 5, 4.8, 4.6, 
        4.9) ),
        ("cutsubleadisosumoet6c8" , cms.vdouble(6.5, 5.7, 4.8, 4.4, 4.6, 
        4.5) ),
        ("cutsubleadisosumoet6c9" , cms.vdouble(6.1, 5.6, 4.8, 4.4, 4.6, 
        4.5) ),
        ("cutsubleadisosumoet7" , cms.vdouble(1.31, 0.3, 1.15, -0.029) ),
        ("cutsubleadisosumoet8" , cms.vdouble(0.94, 0.112, 0.039, -0.029) ),
        ("cutsubleadisosumoet9" , cms.vdouble(0.3, 0.072, -0.04, -0.97) ),
        ("cutsubleadisosumoetbad0" , cms.vdouble(67, 69, 85, 7.2) ),
        ("cutsubleadisosumoetbad1" , cms.vdouble(64, 10.8, 13, 3.5) ),
        ("cutsubleadisosumoetbad10" , cms.vdouble(0.42, 0.41, -0.84, -1.77) ),
        ("cutsubleadisosumoetbad11" , cms.vdouble(100.0, 100.0, 100.0, 100.0) ),
        ("cutsubleadisosumoetbad2" , cms.vdouble(62, 5.2, 7.3, 2.5) ),
        ("cutsubleadisosumoetbad3" , cms.vdouble(11.7, 3.4, 3.9, 1.84) ),
        ("cutsubleadisosumoetbad4" , cms.vdouble(6.1, 2.7, 2.8, 0.66) ),
        ("cutsubleadisosumoetbad5" , cms.vdouble(5.1, 1.62, 1.38, -0.22) ),
        ("cutsubleadisosumoetbad6" , cms.vdouble(3.7, 0.97, 1.38, -0.88) ),
        ("cutsubleadisosumoetbad6c0" , cms.vdouble(73, 59, 85, 94, 74, 
        14.4) ),
        ("cutsubleadisosumoetbad6c1" , cms.vdouble(71, 50, 22, 24, 18.8, 
        12.7) ),
        ("cutsubleadisosumoetbad6c10" , cms.vdouble(7.7, 7.1, 6.1, 5.3, 5.2, 
        6.5) ),
        ("cutsubleadisosumoetbad6c11" , cms.vdouble(100.0, 100.0, 100.0, 100.0, 100.0, 
        100.0) ),
        ("cutsubleadisosumoetbad6c2" , cms.vdouble(71, 16.8, 13.7, 15.6, 13.2, 
        10.4) ),
        ("cutsubleadisosumoetbad6c3" , cms.vdouble(21, 14.4, 11.5, 12.9, 10.2, 
        9.7) ),
        ("cutsubleadisosumoetbad6c4" , cms.vdouble(15.4, 10.8, 10.2, 11.5, 8.5, 
        8.8) ),
        ("cutsubleadisosumoetbad6c5" , cms.vdouble(13.4, 10.4, 9.2, 9.3, 7, 
        7.3) ),
        ("cutsubleadisosumoetbad6c6" , cms.vdouble(11.4, 8.9, 8.3, 8.1, 7, 
        6.7) ),
        ("cutsubleadisosumoetbad6c7" , cms.vdouble(9.2, 8.1, 7.7, 7.4, 7, 
        6.7) ),
        ("cutsubleadisosumoetbad6c8" , cms.vdouble(8.3, 7.7, 6.4, 6.8, 5.3, 
        6.5) ),
        ("cutsubleadisosumoetbad6c9" , cms.vdouble(7.8, 7.4, 6.2, 5.4, 5.3, 
        6.5) ),
        ("cutsubleadisosumoetbad7" , cms.vdouble(1.72, 0.69, 1.14, -0.88) ),
        ("cutsubleadisosumoetbad8" , cms.vdouble(0.86, 0.45, -0.37, -0.88) ),
        ("cutsubleadisosumoetbad9" , cms.vdouble(0.59, 0.42, -0.84, -1.77) ),
        ("cutsubleadpf_drtotk_25_990" , cms.vdouble(0.99, 0.96, 0.0064, 0.0066, 0.0041, 
        0.0032) ),
        ("cutsubleadpf_drtotk_25_991" , cms.vdouble(1, 1, 0.0114, 0.0081, 0.0068, 
        0.0109) ),
        ("cutsubleadpf_drtotk_25_9910" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsubleadpf_drtotk_25_9911" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsubleadpf_drtotk_25_992" , cms.vdouble(1, 1, 0.0125, 0.0094, 0.02, 
        0.0145) ),
        ("cutsubleadpf_drtotk_25_993" , cms.vdouble(1, 1, 0.06, 0.0163, 0.057, 
        0.0168) ),
        ("cutsubleadpf_drtotk_25_994" , cms.vdouble(1, 1, 0.44, 0.36, 0.39, 
        0.0182) ),
        ("cutsubleadpf_drtotk_25_995" , cms.vdouble(1, 1, 0.99, 0.99, 0.96, 
        0.019) ),
        ("cutsubleadpf_drtotk_25_996" , cms.vdouble(1, 1, 1, 1, 1, 
        0.53) ),
        ("cutsubleadpf_drtotk_25_997" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsubleadpf_drtotk_25_998" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsubleadpf_drtotk_25_999" , cms.vdouble(1, 1, 1, 1, 1, 
        1) ),
        ("cutsubleadpfchisooet0" , cms.vdouble(6.3, 4.8, 4.4, 4.4, 7.1, 
        4.2) ),
        ("cutsubleadpfchisooet1" , cms.vdouble(5.6, 3.7, 3.5, 3.7, 4.5, 
        3.4) ),
        ("cutsubleadpfchisooet10" , cms.vdouble(0.96, 0.149, 0.0135, 7e-05, 0.159, 
        0) ),
        ("cutsubleadpfchisooet11" , cms.vdouble(6.0, 6.0, 6.0, 5.0, 5.0, 
        5.0) ),
        ("cutsubleadpfchisooet2" , cms.vdouble(4.7, 3.7, 2.5, 3.4, 4.5, 
        2.3) ),
        ("cutsubleadpfchisooet3" , cms.vdouble(4, 3.3, 2.4, 3.2, 4.3, 
        2) ),
        ("cutsubleadpfchisooet4" , cms.vdouble(3.4, 2.3, 2.4, 2.4, 1.9, 
        1.35) ),
        ("cutsubleadpfchisooet5" , cms.vdouble(3.2, 2.1, 2.1, 1.97, 1.33, 
        0.028) ),
        ("cutsubleadpfchisooet6" , cms.vdouble(2.7, 2.1, 1.68, 1.25, 1.1, 
        0.028) ),
        ("cutsubleadpfchisooet7" , cms.vdouble(1.32, 1.3, 1.45, 0.76, 0.97, 
        0.00047) ),
        ("cutsubleadpfchisooet8" , cms.vdouble(1.07, 1, 1.34, 0.69, 0.69, 
        5e-06) ),
        ("cutsubleadpfchisooet9" , cms.vdouble(0.98, 0.68, 1.34, 0.007, 0.52, 
        0) ),
        ("cutsubleadpfhovere0" , cms.vdouble(0.129, 0.143, 0.112, 0.146, 0.136, 
        0.147) ),
        ("cutsubleadpfhovere1" , cms.vdouble(0.129, 0.142, 0.11, 0.146, 0.136, 
        0.075) ),
        ("cutsubleadpfhovere10" , cms.vdouble(0.02, 0.0164, 0.021, 0.029, 0.018, 
        0.0178) ),
        ("cutsubleadpfhovere11" , cms.vdouble(99.0, 99.0, 99.0, 99.0, 99.0, 
        99.0) ),
        ("cutsubleadpfhovere2" , cms.vdouble(0.123, 0.142, 0.11, 0.146, 0.134, 
        0.061) ),
        ("cutsubleadpfhovere3" , cms.vdouble(0.123, 0.142, 0.11, 0.146, 0.108, 
        0.059) ),
        ("cutsubleadpfhovere4" , cms.vdouble(0.123, 0.142, 0.11, 0.146, 0.09, 
        0.059) ),
        ("cutsubleadpfhovere5" , cms.vdouble(0.123, 0.142, 0.11, 0.108, 0.068, 
        0.059) ),
        ("cutsubleadpfhovere6" , cms.vdouble(0.123, 0.122, 0.11, 0.108, 0.026, 
        0.059) ),
        ("cutsubleadpfhovere7" , cms.vdouble(0.09, 0.091, 0.11, 0.067, 0.0198, 
        0.053) ),
        ("cutsubleadpfhovere8" , cms.vdouble(0.086, 0.091, 0.059, 0.057, 0.0198, 
        0.036) ),
        ("cutsubleadpfhovere9" , cms.vdouble(0.039, 0.04, 0.059, 0.032, 0.0198, 
        0.0178) ),
        ("cutsubleadpfisosumoet0" , cms.vdouble(9.2, 7.4, 6, 9.3, 10, 
        7) ),
        ("cutsubleadpfisosumoet1" , cms.vdouble(8, 6.5, 5.1, 6.2, 6.7, 
        4.8) ),
        ("cutsubleadpfisosumoet10" , cms.vdouble(2.2, 2.1, 1.37, 2.1, 2.2, 
        1.7) ),
        ("cutsubleadpfisosumoet11" , cms.vdouble(8.0, 8.0, 8.0, 7.0, 7.0, 
        7.0) ),
        ("cutsubleadpfisosumoet2" , cms.vdouble(7.2, 6.2, 4.7, 5.7, 6.7, 
        4.3) ),
        ("cutsubleadpfisosumoet3" , cms.vdouble(6.1, 5.5, 4.7, 5.6, 6.3, 
        3.8) ),
        ("cutsubleadpfisosumoet4" , cms.vdouble(5.6, 4.8, 4.3, 5.6, 3.9, 
        3) ),
        ("cutsubleadpfisosumoet5" , cms.vdouble(5.2, 4.6, 3.6, 4, 3.2, 
        2.5) ),
        ("cutsubleadpfisosumoet6" , cms.vdouble(4.8, 4.3, 2.7, 2.7, 2.5, 
        2.4) ),
        ("cutsubleadpfisosumoet7" , cms.vdouble(2.8, 2.8, 2.4, 2.5, 2.4, 
        2.3) ),
        ("cutsubleadpfisosumoet8" , cms.vdouble(2.5, 2.4, 2.2, 2.3, 2.2, 
        1.98) ),
        ("cutsubleadpfisosumoet9" , cms.vdouble(2.4, 2.2, 1.88, 2.2, 2.2, 
        1.7) ),
        ("cutsubleadpfisosumoetbad0" , cms.vdouble(47, 55, 11.3, 19.3, 58, 
        7.7) ),
        ("cutsubleadpfisosumoetbad1" , cms.vdouble(44, 11.4, 8.2, 10.1, 9.5, 
        6.5) ),
        ("cutsubleadpfisosumoetbad10" , cms.vdouble(4.9, 3.6, 1.17, 2.1, 3.1, 
        1.55) ),
        ("cutsubleadpfisosumoetbad11" , cms.vdouble(100.0, 100.0, 100.0, 100.0, 100.0, 
        100.0) ),
        ("cutsubleadpfisosumoetbad2" , cms.vdouble(11.8, 8.2, 6.8, 8.4, 6.5, 
        5.6) ),
        ("cutsubleadpfisosumoetbad3" , cms.vdouble(10.9, 7.5, 6.2, 6.9, 5.6, 
        5.1) ),
        ("cutsubleadpfisosumoetbad4" , cms.vdouble(9.8, 6.7, 5.6, 5.3, 4.6, 
        4) ),
        ("cutsubleadpfisosumoetbad5" , cms.vdouble(9, 6.4, 4.6, 4.9, 4.2, 
        3.8) ),
        ("cutsubleadpfisosumoetbad6" , cms.vdouble(6.2, 6, 4.2, 4, 4.2, 
        3.4) ),
        ("cutsubleadpfisosumoetbad7" , cms.vdouble(6.2, 5.4, 3.5, 4, 4.2, 
        2.6) ),
        ("cutsubleadpfisosumoetbad8" , cms.vdouble(5.7, 4.7, 2.7, 4, 4.2, 
        2.1) ),
        ("cutsubleadpfisosumoetbad9" , cms.vdouble(5.6, 3.6, 2.1, 2.4, 4.2, 
        1.55) ),
        ("cutsubleadpfr90" , cms.vdouble(0.94, 0.9, 0.23, 0.94, 0.9, 
        0.23) ),
        ("cutsubleadpfr91" , cms.vdouble(0.94, 0.9, 0.23, 0.94, 0.9, 
        0.23) ),
        ("cutsubleadpfr910" , cms.vdouble(0.94, 0.9, 0.46, 0.95, 0.9, 
        0.76) ),
        ("cutsubleadpfr911" , cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0) ),
        ("cutsubleadpfr92" , cms.vdouble(0.94, 0.9, 0.24, 0.94, 0.9, 
        0.23) ),
        ("cutsubleadpfr93" , cms.vdouble(0.94, 0.9, 0.24, 0.94, 0.9, 
        0.23) ),
        ("cutsubleadpfr94" , cms.vdouble(0.94, 0.9, 0.38, 0.94, 0.9, 
        0.23) ),
        ("cutsubleadpfr95" , cms.vdouble(0.94, 0.9, 0.38, 0.94, 0.9, 
        0.23) ),
        ("cutsubleadpfr96" , cms.vdouble(0.94, 0.9, 0.38, 0.94, 0.9, 
        0.23) ),
        ("cutsubleadpfr97" , cms.vdouble(0.94, 0.9, 0.38, 0.94, 0.9, 
        0.32) ),
        ("cutsubleadpfr98" , cms.vdouble(0.94, 0.9, 0.42, 0.95, 0.9, 
        0.53) ),
        ("cutsubleadpfr99" , cms.vdouble(0.94, 0.9, 0.42, 0.95, 0.9, 
        0.67) ),
        ("cutsubleadpfsieie0" , cms.vdouble(0.0116, 0.0114, 0.0103, 0.032, 0.031, 
        0.03) ),
        ("cutsubleadpfsieie1" , cms.vdouble(0.0116, 0.0114, 0.0101, 0.029, 0.03, 
        0.029) ),
        ("cutsubleadpfsieie10" , cms.vdouble(0.0107, 0.0101, 0.009, 0.024, 0.023, 
        0.023) ),
        ("cutsubleadpfsieie11" , cms.vdouble(0.012, 0.012, 0.012, 0.03, 0.03, 
        0.03) ),
        ("cutsubleadpfsieie2" , cms.vdouble(0.0113, 0.0111, 0.0101, 0.028, 0.03, 
        0.028) ),
        ("cutsubleadpfsieie3" , cms.vdouble(0.0112, 0.0108, 0.01, 0.028, 0.029, 
        0.027) ),
        ("cutsubleadpfsieie4" , cms.vdouble(0.0108, 0.0107, 0.0099, 0.028, 0.028, 
        0.027) ),
        ("cutsubleadpfsieie5" , cms.vdouble(0.0107, 0.0107, 0.0098, 0.028, 0.027, 
        0.027) ),
        ("cutsubleadpfsieie6" , cms.vdouble(0.0107, 0.0105, 0.0098, 0.027, 0.027, 
        0.027) ),
        ("cutsubleadpfsieie7" , cms.vdouble(0.0107, 0.0105, 0.0098, 0.027, 0.027, 
        0.027) ),
        ("cutsubleadpfsieie8" , cms.vdouble(0.0107, 0.0104, 0.0097, 0.027, 0.024, 
        0.025) ),
        ("cutsubleadpfsieie9" , cms.vdouble(0.0107, 0.0103, 0.0097, 0.025, 0.023, 
        0.024) ),
        ("cutsubleadr90" , cms.vdouble(0.94, 0.31, 0.92, 0.29) ),
        ("cutsubleadr91" , cms.vdouble(0.94, 0.32, 0.94, 0.29) ),
        ("cutsubleadr910" , cms.vdouble(0.95, 0.69, 0.97, 0.85) ),
        ("cutsubleadr911" , cms.vdouble(0.0, 0.0, 0.0, 0.0) ),
        ("cutsubleadr92" , cms.vdouble(0.94, 0.34, 0.94, 0.29) ),
        ("cutsubleadr93" , cms.vdouble(0.94, 0.36, 0.94, 0.32) ),
        ("cutsubleadr94" , cms.vdouble(0.94, 0.41, 0.94, 0.34) ),
        ("cutsubleadr95" , cms.vdouble(0.94, 0.47, 0.94, 0.52) ),
        ("cutsubleadr96" , cms.vdouble(0.94, 0.69, 0.97, 0.52) ),
        ("cutsubleadr96c0" , cms.vdouble(0.94, 0.9, 0.26, 0.94, 0.9, 
        0.28) ),
        ("cutsubleadr96c1" , cms.vdouble(0.94, 0.9, 0.27, 0.94, 0.9, 
        0.28) ),
        ("cutsubleadr96c10" , cms.vdouble(0.94, 0.9, 0.52, 0.95, 0.92, 
        0.76) ),
        ("cutsubleadr96c11" , cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0) ),
        ("cutsubleadr96c2" , cms.vdouble(0.94, 0.9, 0.32, 0.94, 0.9, 
        0.28) ),
        ("cutsubleadr96c3" , cms.vdouble(0.94, 0.9, 0.32, 0.94, 0.9, 
        0.28) ),
        ("cutsubleadr96c4" , cms.vdouble(0.94, 0.9, 0.35, 0.94, 0.9, 
        0.3) ),
        ("cutsubleadr96c5" , cms.vdouble(0.94, 0.9, 0.37, 0.94, 0.9, 
        0.58) ),
        ("cutsubleadr96c6" , cms.vdouble(0.94, 0.9, 0.37, 0.95, 0.9, 
        0.58) ),
        ("cutsubleadr96c7" , cms.vdouble(0.94, 0.9, 0.39, 0.95, 0.9, 
        0.69) ),
        ("cutsubleadr96c8" , cms.vdouble(0.94, 0.9, 0.52, 0.95, 0.92, 
        0.76) ),
    ] +
    [
        ("cutsubleadr96c9" , cms.vdouble(0.94, 0.9, 0.52, 0.95, 0.92, 
        0.76) ),
        ("cutsubleadr97" , cms.vdouble(0.94, 0.69, 0.97, 0.73) ),
        ("cutsubleadr98" , cms.vdouble(0.94, 0.69, 0.97, 0.73) ),
        ("cutsubleadr99" , cms.vdouble(0.95, 0.69, 0.97, 0.84) ),
        ("cutsubleadsieie0" , cms.vdouble(0.0112, 0.0102, 0.029, 0.028) ),
        ("cutsubleadsieie1" , cms.vdouble(0.0109, 0.01, 0.029, 0.028) ),
        ("cutsubleadsieie10" , cms.vdouble(0.0094, 0.009, 0.024, 0.023) ),
        ("cutsubleadsieie11" , cms.vdouble(0.012, 0.012, 0.03, 0.03) ),
        ("cutsubleadsieie2" , cms.vdouble(0.0107, 0.0099, 0.028, 0.027) ),
        ("cutsubleadsieie3" , cms.vdouble(0.0106, 0.0097, 0.028, 0.027) ),
        ("cutsubleadsieie4" , cms.vdouble(0.0104, 0.0094, 0.028, 0.025) ),
        ("cutsubleadsieie5" , cms.vdouble(0.0101, 0.0093, 0.027, 0.023) ),
        ("cutsubleadsieie6" , cms.vdouble(0.0099, 0.0092, 0.027, 0.023) ),
        ("cutsubleadsieie6c0" , cms.vdouble(0.0114, 0.011, 0.0101, 0.029, 0.029, 
        0.028) ),
        ("cutsubleadsieie6c1" , cms.vdouble(0.011, 0.0107, 0.0098, 0.029, 0.029, 
        0.027) ),
        ("cutsubleadsieie6c10" , cms.vdouble(0.0098, 0.0091, 0.0086, 0.024, 0.023, 
        0.022) ),
        ("cutsubleadsieie6c11" , cms.vdouble(0.012, 0.012, 0.012, 0.03, 0.03, 
        0.03) ),
        ("cutsubleadsieie6c2" , cms.vdouble(0.0108, 0.0105, 0.0097, 0.028, 0.028, 
        0.027) ),
        ("cutsubleadsieie6c3" , cms.vdouble(0.0106, 0.0102, 0.0095, 0.028, 0.028, 
        0.026) ),
        ("cutsubleadsieie6c4" , cms.vdouble(0.0105, 0.01, 0.0094, 0.028, 0.027, 
        0.026) ),
        ("cutsubleadsieie6c5" , cms.vdouble(0.0104, 0.0098, 0.0094, 0.028, 0.027, 
        0.024) ),
        ("cutsubleadsieie6c6" , cms.vdouble(0.0101, 0.0097, 0.0093, 0.028, 0.025, 
        0.023) ),
        ("cutsubleadsieie6c7" , cms.vdouble(0.0099, 0.0092, 0.0092, 0.028, 0.025, 
        0.022) ),
        ("cutsubleadsieie6c8" , cms.vdouble(0.0098, 0.0092, 0.0088, 0.025, 0.025, 
        0.022) ),
        ("cutsubleadsieie6c9" , cms.vdouble(0.0098, 0.0091, 0.0086, 0.024, 0.023, 
        0.022) ),
        ("cutsubleadsieie7" , cms.vdouble(0.0098, 0.009, 0.026, 0.023) ),
        ("cutsubleadsieie8" , cms.vdouble(0.0097, 0.009, 0.026, 0.023) ),
        ("cutsubleadsieie9" , cms.vdouble(0.0094, 0.009, 0.024, 0.023) ),
        ("cutsubleadtrkisooetom0" , cms.vdouble(7.5, 4.5, 5.2, 2.5) ),
        ("cutsubleadtrkisooetom1" , cms.vdouble(6.4, 3.4, 3.8, 2.1) ),
        ("cutsubleadtrkisooetom10" , cms.vdouble(0.43, 0.0111, 0.0075, 0.0073) ),
        ("cutsubleadtrkisooetom11" , cms.vdouble(7.0, 7.0, 4.0, 4.0) ),
        ("cutsubleadtrkisooetom2" , cms.vdouble(4.7, 2.9, 3.8, 1.63) ),
        ("cutsubleadtrkisooetom3" , cms.vdouble(3.5, 2.2, 2.3, 1.45) ),
        ("cutsubleadtrkisooetom4" , cms.vdouble(3.4, 1.86, 1.67, 1.44) ),
        ("cutsubleadtrkisooetom5" , cms.vdouble(2.9, 1.6, 1.55, 1.44) ),
        ("cutsubleadtrkisooetom6" , cms.vdouble(1.93, 1.4, 1.48, 0.056) ),
        ("cutsubleadtrkisooetom6c0" , cms.vdouble(7.9, 5.6, 4.5, 5.9, 4.2, 
        3) ),
        ("cutsubleadtrkisooetom6c1" , cms.vdouble(6.6, 4.7, 4.4, 5.7, 2.4, 
        2.2) ),
        ("cutsubleadtrkisooetom6c10" , cms.vdouble(1.6, 0.033, 0.00082, 0.68, 0.48, 
        0.26) ),
        ("cutsubleadtrkisooetom6c11" , cms.vdouble(7.0, 7.0, 7.0, 4.0, 4.0, 
        4.0) ),
        ("cutsubleadtrkisooetom6c2" , cms.vdouble(5.8, 4, 2.9, 4.6, 2.3, 
        1.93) ),
        ("cutsubleadtrkisooetom6c3" , cms.vdouble(4, 3.4, 2.5, 1.89, 1.7, 
        1.67) ),
        ("cutsubleadtrkisooetom6c4" , cms.vdouble(3.6, 2.9, 1.96, 1.63, 1.43, 
        1.36) ),
        ("cutsubleadtrkisooetom6c5" , cms.vdouble(2.7, 2.9, 1.56, 1.56, 1.06, 
        1.36) ),
        ("cutsubleadtrkisooetom6c6" , cms.vdouble(2.3, 2.6, 1.54, 1.56, 1.06, 
        1.07) ),
        ("cutsubleadtrkisooetom6c7" , cms.vdouble(2.1, 1.79, 1.07, 1.56, 0.8, 
        0.26) ),
        ("cutsubleadtrkisooetom6c8" , cms.vdouble(2.1, 1.55, 0.25, 0.8, 0.48, 
        0.26) ),
        ("cutsubleadtrkisooetom6c9" , cms.vdouble(1.83, 1.2, 0.029, 0.68, 0.48, 
        0.26) ),
        ("cutsubleadtrkisooetom7" , cms.vdouble(1.42, 0.76, 1.48, 0.056) ),
        ("cutsubleadtrkisooetom8" , cms.vdouble(1.21, 0.51, 0.27, 0.056) ),
        ("cutsubleadtrkisooetom9" , cms.vdouble(0.43, 0.4, 0.0075, 0.056) ),
        ]
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.mvaEleID_PHYS14_PU20bx25_nonTrig_V1_producer_config = cms.PSet(
    mvaName = cms.string('ElectronMVAEstimatorRun2Phys14NonTrig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml')
)

process.mvaEleID_PHYS14_PU20bx25_nonTrig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.253, 0.081, -0.081, 0.965, 0.917, 
            0.683),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(False)
)

process.mvaEleID_PHYS14_PU20bx25_nonTrig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.483, -0.267, -0.323, 0.933, 0.825, 
            0.337),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(False)
)

process.mvaEleID_Spring15_25ns_Trig_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
)

process.mvaEleID_Spring15_25ns_Trig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
        mvaCuts = cms.vdouble(0.988153, 0.96791, 0.841729),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-Trig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_Trig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
        mvaCuts = cms.vdouble(0.972153, 0.922126, 0.610764),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-Trig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring15NonTrig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml')
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(0.287435, 0.221846, -0.303263, 0.967083, 0.929117, 
            0.726311),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.083313, -0.235222, -0.67099, 0.913286, 0.805013, 
            0.358969),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_wpLoose = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.265, -0.556, -0.551, -0.072, -0.286, 
            -0.267),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_50ns_Trig_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
    mvaTag = cms.string('50nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
)

process.mvaEleID_Spring15_50ns_Trig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Categories"),
        mvaCuts = cms.vdouble(0.981841, 0.946762, 0.79704),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-50ns-Trig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_50ns_Trig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Categories"),
        mvaCuts = cms.vdouble(0.953843, 0.849994, 0.514118),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-50ns-Trig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring16_GeneralPurpose_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring16GeneralPurpose'),
    mvaTag = cms.string('V1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml')
)

process.mvaEleID_Spring16_HZZ_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring16HZZ'),
    mvaTag = cms.string('V1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml')
)

process.mvaPhoID_Spring15_25ns_nonTrig_V2p1_producer_config = cms.PSet(
    esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
    full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
    full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
    full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
    full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
    full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
    mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
    mvaTag = cms.string('25nsV2p1'),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    useValueMaps = cms.bool(False),
    weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EB_V2.weights.xml', 
        'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EE_V2.weights.xml')
)

process.mvaPhoID_Spring15_25ns_nonTrig_V2p1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Categories"),
        mvaCuts = cms.vdouble(0.374, 0.336),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring15NonTrig25nsV2p1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-Spring15-25ns-nonTrig-V2p1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaPhoID_Spring15_50ns_nonTrig_V2p1_producer_config = cms.PSet(
    esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
    full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
    full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
    full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
    full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
    full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
    mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
    mvaTag = cms.string('50nsV2p1'),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    useValueMaps = cms.bool(False),
    weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EB_V2.weights.xml', 
        'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EE_V2.weights.xml')
)

process.mvaPhoID_Spring15_50ns_nonTrig_V2p1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring15NonTrig50nsV2p1Categories"),
        mvaCuts = cms.vdouble(0.29538, 0.45837),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring15NonTrig50nsV2p1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-Spring15-50ns-nonTrig-V2p1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaPhoID_Spring16_nonTrig_V1_producer_config = cms.PSet(
    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt'),
    mvaName = cms.string('PhotonMVAEstimatorRun2Spring16NonTrig'),
    mvaTag = cms.string('V1'),
    phoChargedIsolation = cms.InputTag("egmPhotonIsolation","h+-DR030-"),
    phoIsoCutoff = cms.double(2.5),
    phoIsoPtScalingCoeff = cms.vdouble(0.0053, 0.0034),
    phoPhotonIsolation = cms.InputTag("egmPhotonIsolation","gamma-DR030-"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationWithConeVeto"),
    rho = cms.InputTag("fixedGridRhoAll"),
    weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EB_V3.weights.xml', 
        'RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EE_V3.weights.xml')
)

process.mvaPhoID_Spring16_nonTrig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Categories"),
        mvaCuts = cms.vdouble(0.68, 0.6),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-Spring16-nonTrig-V1-wp80')
)

process.mvaPhoID_Spring16_nonTrig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('PhoMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Categories"),
        mvaCuts = cms.vdouble(0.2, 0.2),
        mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaPhoID-Spring16-nonTrig-V1-wp90')
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    allowUnscheduled = cms.untracked.bool(True)
)

process.pset = cms.PSet(
    etaMax = cms.double(-2.901376),
    etaMin = cms.double(-5.2),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('egammaHFMinus'),
    px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
    py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
    type = cms.int32(7),
    varType = cms.int32(0)
)

process.regressionModifier = cms.PSet(
    eOverP_ECALTRKThr = cms.double(0.025),
    ecalrechitsEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    ecalrechitsEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
    electron_config = cms.PSet(
        regressionKey_ecalonly = cms.vstring('electron_eb_ECALonly_lowpt', 
            'electron_eb_ECALonly', 
            'electron_ee_ECALonly_lowpt', 
            'electron_ee_ECALonly'),
        regressionKey_ecaltrk = cms.vstring('electron_eb_ECALTRK_lowpt', 
            'electron_eb_ECALTRK', 
            'electron_ee_ECALTRK_lowpt', 
            'electron_ee_ECALTRK'),
        uncertaintyKey_ecalonly = cms.vstring('electron_eb_ECALonly_lowpt_var', 
            'electron_eb_ECALonly_var', 
            'electron_ee_ECALonly_lowpt_var', 
            'electron_ee_ECALonly_var'),
        uncertaintyKey_ecaltrk = cms.vstring('electron_eb_ECALTRK_lowpt_var', 
            'electron_eb_ECALTRK_var', 
            'electron_ee_ECALTRK_lowpt_var', 
            'electron_ee_ECALTRK_var')
    ),
    epDiffSig_ECALTRKThr = cms.double(15.0),
    epSig_ECALTRKThr = cms.double(10.0),
    highEnergy_ECALTRKThr = cms.double(200.0),
    lowEnergy_ECALTRKThr = cms.double(50.0),
    lowEnergy_ECALonlyThr = cms.double(300.0),
    modifierName = cms.string('EGExtraInfoModifierFromDBUser'),
    photon_config = cms.PSet(
        regressionKey_ecalonly = cms.vstring('photon_eb_ECALonly_lowpt', 
            'photon_eb_ECALonly', 
            'photon_ee_ECALonly_lowpt', 
            'photon_ee_ECALonly'),
        uncertaintyKey_ecalonly = cms.vstring('photon_eb_ECALonly_lowpt_var', 
            'photon_eb_ECALonly_var', 
            'photon_ee_ECALonly_lowpt_var', 
            'photon_ee_ECALonly_var')
    ),
    rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
    useLocalFile = cms.bool(False)
)

process.softPFElectronCommon = cms.PSet(
    gbrForestLabel = cms.string('btag_SoftPFElectron_BDT'),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SoftLepton/data/SoftPFElectron_BDT.weights.xml.gz')
)

process.softPFMuonCommon = cms.PSet(
    gbrForestLabel = cms.string('btag_SoftPFMuon_BDT'),
    useAdaBoost = cms.bool(True),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SoftLepton/data/SoftPFMuon_BDT.weights.xml.gz')
)

process.trackPseudoSelectionBlock = cms.PSet(
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    )
)

process.trackSelectionBlock = cms.PSet(
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    )
)

process.variableJTAPars = cms.PSet(
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5)
)

process.IsoConeDefinitions = cms.VPSet(cms.PSet(
    coneSize = cms.double(0.3),
    isolateAgainst = cms.string('h+'),
    isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
    miniAODVertexCodes = cms.vuint32(2, 3),
    vertexIndex = cms.int32(0)
), 
    cms.PSet(
        coneSize = cms.double(0.3),
        isolateAgainst = cms.string('h0'),
        isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
        miniAODVertexCodes = cms.vuint32(2, 3),
        vertexIndex = cms.int32(0)
    ), 
    cms.PSet(
        coneSize = cms.double(0.3),
        isolateAgainst = cms.string('gamma'),
        isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
        miniAODVertexCodes = cms.vuint32(2, 3),
        vertexIndex = cms.int32(0)
    ))

process.c_vs_b_vars_vpset = cms.VPSet(cms.PSet(
    default = cms.double(-1),
    name = cms.string('vertexLeptonCategory'),
    taggingVarName = cms.string('vertexLeptonCategory')
), 
    cms.PSet(
        default = cms.double(-100),
        idx = cms.int32(0),
        name = cms.string('trackSip2dSig_0'),
        taggingVarName = cms.string('trackSip2dSig')
    ), 
    cms.PSet(
        default = cms.double(-100),
        idx = cms.int32(1),
        name = cms.string('trackSip2dSig_1'),
        taggingVarName = cms.string('trackSip2dSig')
    ), 
    cms.PSet(
        default = cms.double(-100),
        idx = cms.int32(0),
        name = cms.string('trackSip3dSig_0'),
        taggingVarName = cms.string('trackSip3dSig')
    ), 
    cms.PSet(
        default = cms.double(-100),
        idx = cms.int32(1),
        name = cms.string('trackSip3dSig_1'),
        taggingVarName = cms.string('trackSip3dSig')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('trackPtRel_0'),
        taggingVarName = cms.string('trackPtRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('trackPtRel_1'),
        taggingVarName = cms.string('trackPtRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('trackPPar_0'),
        taggingVarName = cms.string('trackPPar')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('trackPPar_1'),
        taggingVarName = cms.string('trackPPar')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('trackEtaRel_0'),
        taggingVarName = cms.string('trackEtaRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('trackEtaRel_1'),
        taggingVarName = cms.string('trackEtaRel')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('trackDeltaR_0'),
        taggingVarName = cms.string('trackDeltaR')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(1),
        name = cms.string('trackDeltaR_1'),
        taggingVarName = cms.string('trackDeltaR')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('trackPtRatio_0'),
        taggingVarName = cms.string('trackPtRatio')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(1),
        name = cms.string('trackPtRatio_1'),
        taggingVarName = cms.string('trackPtRatio')
    ), 
    cms.PSet(
        default = cms.double(1.1),
        idx = cms.int32(0),
        name = cms.string('trackPParRatio_0'),
        taggingVarName = cms.string('trackPParRatio')
    ), 
    cms.PSet(
        default = cms.double(1.1),
        idx = cms.int32(1),
        name = cms.string('trackPParRatio_1'),
        taggingVarName = cms.string('trackPParRatio')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('trackJetDist_0'),
        taggingVarName = cms.string('trackJetDist')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(1),
        name = cms.string('trackJetDist_1'),
        taggingVarName = cms.string('trackJetDist')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('trackDecayLenVal_0'),
        taggingVarName = cms.string('trackDecayLenVal')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(1),
        name = cms.string('trackDecayLenVal_1'),
        taggingVarName = cms.string('trackDecayLenVal')
    ), 
    cms.PSet(
        default = cms.double(0),
        name = cms.string('jetNSecondaryVertices'),
        taggingVarName = cms.string('jetNSecondaryVertices')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        name = cms.string('jetNTracks'),
        taggingVarName = cms.string('jetNTracks')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        name = cms.string('trackSumJetEtRatio'),
        taggingVarName = cms.string('trackSumJetEtRatio')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        name = cms.string('trackSumJetDeltaR'),
        taggingVarName = cms.string('trackSumJetDeltaR')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('vertexMass_0'),
        taggingVarName = cms.string('vertexMass')
    ), 
    cms.PSet(
        default = cms.double(-10),
        idx = cms.int32(0),
        name = cms.string('vertexEnergyRatio_0'),
        taggingVarName = cms.string('vertexEnergyRatio')
    ), 
    cms.PSet(
        default = cms.double(-999),
        idx = cms.int32(0),
        name = cms.string('trackSip2dSigAboveCharm_0'),
        taggingVarName = cms.string('trackSip2dSigAboveCharm')
    ), 
    cms.PSet(
        default = cms.double(-999),
        idx = cms.int32(0),
        name = cms.string('trackSip3dSigAboveCharm_0'),
        taggingVarName = cms.string('trackSip3dSigAboveCharm')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('flightDistance2dSig_0'),
        taggingVarName = cms.string('flightDistance2dSig')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('flightDistance3dSig_0'),
        taggingVarName = cms.string('flightDistance3dSig')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('vertexJetDeltaR_0'),
        taggingVarName = cms.string('vertexJetDeltaR')
    ), 
    cms.PSet(
        default = cms.double(0),
        idx = cms.int32(0),
        name = cms.string('vertexNTracks_0'),
        taggingVarName = cms.string('vertexNTracks')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('massVertexEnergyFraction_0'),
        taggingVarName = cms.string('massVertexEnergyFraction')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('vertexBoostOverSqrtJetPt_0'),
        taggingVarName = cms.string('vertexBoostOverSqrtJetPt')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('leptonPtRel_0'),
        taggingVarName = cms.string('leptonPtRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('leptonPtRel_1'),
        taggingVarName = cms.string('leptonPtRel')
    ), 
    cms.PSet(
        default = cms.double(-10000),
        idx = cms.int32(0),
        name = cms.string('leptonSip3d_0'),
        taggingVarName = cms.string('leptonSip3d')
    ), 
    cms.PSet(
        default = cms.double(-10000),
        idx = cms.int32(1),
        name = cms.string('leptonSip3d_1'),
        taggingVarName = cms.string('leptonSip3d')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('leptonDeltaR_0'),
        taggingVarName = cms.string('leptonDeltaR')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('leptonDeltaR_1'),
        taggingVarName = cms.string('leptonDeltaR')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('leptonRatioRel_0'),
        taggingVarName = cms.string('leptonRatioRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('leptonRatioRel_1'),
        taggingVarName = cms.string('leptonRatioRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('leptonEtaRel_0'),
        taggingVarName = cms.string('leptonEtaRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('leptonEtaRel_1'),
        taggingVarName = cms.string('leptonEtaRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('leptonRatio_0'),
        taggingVarName = cms.string('leptonRatio')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('leptonRatio_1'),
        taggingVarName = cms.string('leptonRatio')
    ))

process.c_vs_l_vars_vpset = cms.VPSet(cms.PSet(
    default = cms.double(-1),
    name = cms.string('vertexLeptonCategory'),
    taggingVarName = cms.string('vertexLeptonCategory')
), 
    cms.PSet(
        default = cms.double(-100),
        idx = cms.int32(0),
        name = cms.string('trackSip2dSig_0'),
        taggingVarName = cms.string('trackSip2dSig')
    ), 
    cms.PSet(
        default = cms.double(-100),
        idx = cms.int32(1),
        name = cms.string('trackSip2dSig_1'),
        taggingVarName = cms.string('trackSip2dSig')
    ), 
    cms.PSet(
        default = cms.double(-100),
        idx = cms.int32(0),
        name = cms.string('trackSip3dSig_0'),
        taggingVarName = cms.string('trackSip3dSig')
    ), 
    cms.PSet(
        default = cms.double(-100),
        idx = cms.int32(1),
        name = cms.string('trackSip3dSig_1'),
        taggingVarName = cms.string('trackSip3dSig')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('trackPtRel_0'),
        taggingVarName = cms.string('trackPtRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('trackPtRel_1'),
        taggingVarName = cms.string('trackPtRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('trackPPar_0'),
        taggingVarName = cms.string('trackPPar')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('trackPPar_1'),
        taggingVarName = cms.string('trackPPar')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('trackEtaRel_0'),
        taggingVarName = cms.string('trackEtaRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('trackEtaRel_1'),
        taggingVarName = cms.string('trackEtaRel')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('trackDeltaR_0'),
        taggingVarName = cms.string('trackDeltaR')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(1),
        name = cms.string('trackDeltaR_1'),
        taggingVarName = cms.string('trackDeltaR')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('trackPtRatio_0'),
        taggingVarName = cms.string('trackPtRatio')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(1),
        name = cms.string('trackPtRatio_1'),
        taggingVarName = cms.string('trackPtRatio')
    ), 
    cms.PSet(
        default = cms.double(1.1),
        idx = cms.int32(0),
        name = cms.string('trackPParRatio_0'),
        taggingVarName = cms.string('trackPParRatio')
    ), 
    cms.PSet(
        default = cms.double(1.1),
        idx = cms.int32(1),
        name = cms.string('trackPParRatio_1'),
        taggingVarName = cms.string('trackPParRatio')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('trackJetDist_0'),
        taggingVarName = cms.string('trackJetDist')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(1),
        name = cms.string('trackJetDist_1'),
        taggingVarName = cms.string('trackJetDist')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('trackDecayLenVal_0'),
        taggingVarName = cms.string('trackDecayLenVal')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(1),
        name = cms.string('trackDecayLenVal_1'),
        taggingVarName = cms.string('trackDecayLenVal')
    ), 
    cms.PSet(
        default = cms.double(0),
        name = cms.string('jetNSecondaryVertices'),
        taggingVarName = cms.string('jetNSecondaryVertices')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        name = cms.string('jetNTracks'),
        taggingVarName = cms.string('jetNTracks')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        name = cms.string('trackSumJetEtRatio'),
        taggingVarName = cms.string('trackSumJetEtRatio')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        name = cms.string('trackSumJetDeltaR'),
        taggingVarName = cms.string('trackSumJetDeltaR')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('vertexMass_0'),
        taggingVarName = cms.string('vertexMass')
    ), 
    cms.PSet(
        default = cms.double(-10),
        idx = cms.int32(0),
        name = cms.string('vertexEnergyRatio_0'),
        taggingVarName = cms.string('vertexEnergyRatio')
    ), 
    cms.PSet(
        default = cms.double(-999),
        idx = cms.int32(0),
        name = cms.string('trackSip2dSigAboveCharm_0'),
        taggingVarName = cms.string('trackSip2dSigAboveCharm')
    ), 
    cms.PSet(
        default = cms.double(-999),
        idx = cms.int32(0),
        name = cms.string('trackSip3dSigAboveCharm_0'),
        taggingVarName = cms.string('trackSip3dSigAboveCharm')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('flightDistance2dSig_0'),
        taggingVarName = cms.string('flightDistance2dSig')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('flightDistance3dSig_0'),
        taggingVarName = cms.string('flightDistance3dSig')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('vertexJetDeltaR_0'),
        taggingVarName = cms.string('vertexJetDeltaR')
    ), 
    cms.PSet(
        default = cms.double(0),
        idx = cms.int32(0),
        name = cms.string('vertexNTracks_0'),
        taggingVarName = cms.string('vertexNTracks')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('massVertexEnergyFraction_0'),
        taggingVarName = cms.string('massVertexEnergyFraction')
    ), 
    cms.PSet(
        default = cms.double(-0.1),
        idx = cms.int32(0),
        name = cms.string('vertexBoostOverSqrtJetPt_0'),
        taggingVarName = cms.string('vertexBoostOverSqrtJetPt')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('leptonPtRel_0'),
        taggingVarName = cms.string('leptonPtRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('leptonPtRel_1'),
        taggingVarName = cms.string('leptonPtRel')
    ), 
    cms.PSet(
        default = cms.double(-10000),
        idx = cms.int32(0),
        name = cms.string('leptonSip3d_0'),
        taggingVarName = cms.string('leptonSip3d')
    ), 
    cms.PSet(
        default = cms.double(-10000),
        idx = cms.int32(1),
        name = cms.string('leptonSip3d_1'),
        taggingVarName = cms.string('leptonSip3d')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('leptonDeltaR_0'),
        taggingVarName = cms.string('leptonDeltaR')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('leptonDeltaR_1'),
        taggingVarName = cms.string('leptonDeltaR')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('leptonRatioRel_0'),
        taggingVarName = cms.string('leptonRatioRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('leptonRatioRel_1'),
        taggingVarName = cms.string('leptonRatioRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('leptonEtaRel_0'),
        taggingVarName = cms.string('leptonEtaRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('leptonEtaRel_1'),
        taggingVarName = cms.string('leptonEtaRel')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(0),
        name = cms.string('leptonRatio_0'),
        taggingVarName = cms.string('leptonRatio')
    ), 
    cms.PSet(
        default = cms.double(-1),
        idx = cms.int32(1),
        name = cms.string('leptonRatio_1'),
        taggingVarName = cms.string('leptonRatio')
    ))

process.egamma_modifications = cms.VPSet(cms.PSet(
    eOverP_ECALTRKThr = cms.double(0.025),
    ecalrechitsEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    ecalrechitsEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
    electron_config = cms.PSet(
        regressionKey_ecalonly = cms.vstring('electron_eb_ECALonly_lowpt', 
            'electron_eb_ECALonly', 
            'electron_ee_ECALonly_lowpt', 
            'electron_ee_ECALonly'),
        regressionKey_ecaltrk = cms.vstring('electron_eb_ECALTRK_lowpt', 
            'electron_eb_ECALTRK', 
            'electron_ee_ECALTRK_lowpt', 
            'electron_ee_ECALTRK'),
        uncertaintyKey_ecalonly = cms.vstring('electron_eb_ECALonly_lowpt_var', 
            'electron_eb_ECALonly_var', 
            'electron_ee_ECALonly_lowpt_var', 
            'electron_ee_ECALonly_var'),
        uncertaintyKey_ecaltrk = cms.vstring('electron_eb_ECALTRK_lowpt_var', 
            'electron_eb_ECALTRK_var', 
            'electron_ee_ECALTRK_lowpt_var', 
            'electron_ee_ECALTRK_var')
    ),
    epDiffSig_ECALTRKThr = cms.double(15.0),
    epSig_ECALTRKThr = cms.double(10.0),
    highEnergy_ECALTRKThr = cms.double(200.0),
    lowEnergy_ECALTRKThr = cms.double(50.0),
    lowEnergy_ECALonlyThr = cms.double(300.0),
    modifierName = cms.string('EGExtraInfoModifierFromDBUser'),
    photon_config = cms.PSet(
        regressionKey_ecalonly = cms.vstring('photon_eb_ECALonly_lowpt', 
            'photon_eb_ECALonly', 
            'photon_ee_ECALonly_lowpt', 
            'photon_ee_ECALonly'),
        uncertaintyKey_ecalonly = cms.vstring('photon_eb_ECALonly_lowpt_var', 
            'photon_eb_ECALonly_var', 
            'photon_ee_ECALonly_lowpt_var', 
            'photon_ee_ECALonly_var')
    ),
    rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
    useLocalFile = cms.bool(False)
))

process.multPhiCorrParams_T0pcT1T2Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0pcT1T2Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0pcT1Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0pcT1Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0pcTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0pcTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0rtT1T2Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0rtT1T2Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0rtT1Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0rtT1Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0rtTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T0rtTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T1T2Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T1T2Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T1Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_T1Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.multPhiCorrParams_Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.mvaConfigsForEleProducer = cms.VPSet(cms.PSet(
    mvaName = cms.string('ElectronMVAEstimatorRun2Phys14NonTrig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml')
), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring15NonTrig'),
        mvaTag = cms.string('25nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
        mvaTag = cms.string('50nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
        mvaTag = cms.string('25nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring16HZZ'),
        mvaTag = cms.string('V1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring16GeneralPurpose'),
        mvaTag = cms.string('V1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml')
    ))

process.mvaConfigsForPhoProducer = cms.VPSet(cms.PSet(
    esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
    full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
    full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
    full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
    full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
    full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
    mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
    mvaTag = cms.string('50nsV2p1'),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    useValueMaps = cms.bool(False),
    weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EB_V2.weights.xml', 
        'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EE_V2.weights.xml')
), 
    cms.PSet(
        esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
        full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
        full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
        full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
        full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
        full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
        full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
        mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
        mvaTag = cms.string('25nsV2p1'),
        phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
        phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
        phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        useValueMaps = cms.bool(False),
        weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EB_V2.weights.xml', 
            'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EE_V2.weights.xml')
    ), 
    cms.PSet(
        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt'),
        mvaName = cms.string('PhotonMVAEstimatorRun2Spring16NonTrig'),
        mvaTag = cms.string('V1'),
        phoChargedIsolation = cms.InputTag("egmPhotonIsolation","h+-DR030-"),
        phoIsoCutoff = cms.double(2.5),
        phoIsoPtScalingCoeff = cms.vdouble(0.0053, 0.0034),
        phoPhotonIsolation = cms.InputTag("egmPhotonIsolation","gamma-DR030-"),
        phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationWithConeVeto"),
        rho = cms.InputTag("fixedGridRhoAll"),
        weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EB_V3.weights.xml', 
            'RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EE_V3.weights.xml')
    ))

process.patMultPhiCorrParams_T0pcT1SmearTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1SmearTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1T2SmearTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1T2SmearTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1T2Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1T2Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1SmearTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1SmearTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1T2SmearTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1T2SmearTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1T2Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1T2Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.puppiCentral = cms.VPSet(cms.PSet(
    algoId = cms.int32(5),
    applyLowPUCorr = cms.bool(False),
    combOpt = cms.int32(0),
    cone = cms.double(0.4),
    rmsPtMin = cms.double(0.1),
    rmsScaleFactor = cms.double(1.0),
    useCharged = cms.bool(True)
))

process.puppiForward = cms.VPSet(cms.PSet(
    algoId = cms.int32(5),
    applyLowPUCorr = cms.bool(False),
    combOpt = cms.int32(0),
    cone = cms.double(0.4),
    rmsPtMin = cms.double(0.5),
    rmsScaleFactor = cms.double(1.0),
    useCharged = cms.bool(False)
))

process.PFCandAssoMap = cms.EDProducer("PFCand_AssoMap",
    AssociationType = cms.InputTag("Both"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ConversionsCollection = cms.InputTag("allConversions"),
    FinalAssociation = cms.untracked.int32(1),
    GetCleanedCollections = cms.bool(False),
    MaxNumberOfAssociations = cms.int32(1),
    NIVertexCollection = cms.InputTag("particleFlowDisplacedVertex"),
    PFCandidateCollection = cms.InputTag("particleFlow"),
    V0KshortCollection = cms.InputTag("generalV0Candidates","Kshort"),
    V0LambdaCollection = cms.InputTag("generalV0Candidates","Lambda"),
    VertexCollection = cms.InputTag("offlinePrimaryVertices"),
    doReassociation = cms.bool(True),
    ignoreMissingCollection = cms.bool(True),
    nTrackWeight = cms.double(0.001)
)


process.ak4CaloL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector")
)


process.ak4CaloL1FastL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloL6SLBCorrector")
)


process.ak4CaloL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloResidualCorrector")
)


process.ak4CaloL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak4CaloL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector")
)


process.ak4CaloL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloResidualCorrector")
)


process.ak4CaloL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4CaloL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector")
)


process.ak4CaloL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloL6SLBCorrector")
)


process.ak4CaloL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloResidualCorrector")
)


process.ak4CaloL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2Relative')
)


process.ak4CaloL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L3Absolute')
)


process.ak4CaloL6SLBCorrector = cms.EDProducer("L6SLBCorrectorProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4CaloJetsSoftMuonTagInfos")
)


process.ak4CaloResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.ak4JPTL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4JPTL1FastjetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector")
)


process.ak4JPTL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4JPTL1FastjetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector", "ak4JPTResidualCorrector")
)


process.ak4JPTL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak4JPTL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector")
)


process.ak4JPTL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector", "ak4JPTResidualCorrector")
)


process.ak4JPTL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4JPTL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector")
)


process.ak4JPTL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector", "ak4JPTResidualCorrector")
)


process.ak4JPTL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2Relative')
)


process.ak4JPTL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L3Absolute')
)


process.ak4JPTResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2L3Residual')
)


process.ak4L1JPTOffsetCorrector = cms.EDProducer("L1JPTOffsetCorrectorProducer",
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.InputTag("ak4CaloL1OffsetCorrector")
)


process.ak4PFCHSL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1FastjetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector")
)


process.ak4PFCHSL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1FastjetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector", "ak4PFCHSResidualCorrector")
)


process.ak4PFCHSL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFCHSL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1OffsetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector")
)


process.ak4PFCHSL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1OffsetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector", "ak4PFCHSResidualCorrector")
)


process.ak4PFCHSL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4PFCHSL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector")
)


process.ak4PFCHSL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector", "ak4PFCHSResidualCorrector")
)


process.ak4PFCHSL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2Relative')
)


process.ak4PFCHSL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L3Absolute')
)


process.ak4PFCHSResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak4PFL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1FastjetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector")
)


process.ak4PFL1FastL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1FastjetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFL6SLBCorrector")
)


process.ak4PFL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1FastjetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFResidualCorrector")
)


process.ak4PFL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1OffsetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector")
)


process.ak4PFL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1OffsetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFResidualCorrector")
)


process.ak4PFL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4PFL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector")
)


process.ak4PFL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFL6SLBCorrector")
)


process.ak4PFL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFResidualCorrector")
)


process.ak4PFL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2Relative')
)


process.ak4PFL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L3Absolute')
)


process.ak4PFL6SLBCorrector = cms.EDProducer("L6SLBCorrectorProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4PFJetsSoftMuonTagInfos")
)


process.ak4PFPuppiL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1FastjetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector")
)


process.ak4PFPuppiL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1FastjetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector", "ak4PFPuppiResidualCorrector")
)


process.ak4PFPuppiL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFPuppiL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1OffsetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector")
)


process.ak4PFPuppiL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1OffsetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector", "ak4PFPuppiResidualCorrector")
)


process.ak4PFPuppiL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4PFPuppiL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector")
)


process.ak4PFPuppiL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector", "ak4PFPuppiResidualCorrector")
)


process.ak4PFPuppiL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L2Relative')
)


process.ak4PFPuppiL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L3Absolute')
)


process.ak4PFPuppiResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L2L3Residual')
)


process.ak4PFResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2L3Residual')
)


process.ak4TrackL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4TrackL2RelativeCorrector", "ak4TrackL3AbsoluteCorrector")
)


process.ak4TrackL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4TRK'),
    level = cms.string('L2Relative')
)


process.ak4TrackL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4TRK'),
    level = cms.string('L3Absolute')
)


process.basicJetsForMet = cms.EDProducer("PATJetCleanerForType1MET",
    jetCorrEtaMax = cms.double(9.9),
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("patJetsReapplyJEC"),
    type1JetPtThreshold = cms.double(15.0)
)


process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_23Jan_ele'),
    electrons = cms.InputTag("slimmedElectrons"),
    gbrForestName = cms.vstring('electron_eb_ECALTRK_lowpt', 
        'electron_eb_ECALTRK', 
        'electron_ee_ECALTRK_lowpt', 
        'electron_ee_ECALTRK', 
        'electron_eb_ECALTRK_lowpt_var', 
        'electron_eb_ECALTRK_var', 
        'electron_ee_ECALTRK_lowpt_var', 
        'electron_ee_ECALTRK_var'),
    isMC = cms.bool(True),
    isSynchronization = cms.bool(False)
)


process.calibratedPatPhotons = cms.EDProducer("CalibratedPatPhotonProducerRun2",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Moriond17_23Jan_ele'),
    isMC = cms.bool(True),
    isSynchronization = cms.bool(False),
    photons = cms.InputTag("slimmedPhotons")
)


process.caloMetT1 = cms.EDProducer("CorrectedCaloMETProducer",
    src = cms.InputTag("caloMetM"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrCaloMetType1","type1"))
)


process.caloMetT1T2 = cms.EDProducer("CorrectedCaloMETProducer",
    src = cms.InputTag("caloMetM"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrCaloMetType1","type1"), cms.InputTag("corrCaloMetType2"))
)


process.cleanedPatJets = cms.EDProducer("PATJetCleaner",
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("slimmedElectrons")
        ),
        muons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("slimmedMuons")
        )
    ),
    finalCut = cms.string(''),
    preselection = cms.string(''),
    src = cms.InputTag("jetSelectorForMet")
)


process.corrCaloMetType1 = cms.EDProducer("CaloJetMETcorrInputProducer",
    jetCorrEtaMax = cms.double(9.9),
    jetCorrLabel = cms.InputTag("ak4CaloL2L3Corrector"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    src = cms.InputTag("ak4CaloJets"),
    srcMET = cms.InputTag("caloMetM"),
    type1JetPtThreshold = cms.double(20.0)
)


process.corrCaloMetType2 = cms.EDProducer("Type2CorrectionProducer",
    srcUnclEnergySums = cms.VInputTag(cms.InputTag("corrCaloMetType1","type2"), cms.InputTag("muCaloMetCorr")),
    type2CorrFormula = cms.string('A + B*TMath::Exp(-C*x)'),
    type2CorrParameter = cms.PSet(
        A = cms.double(2.0),
        B = cms.double(1.3),
        C = cms.double(0.1)
    )
)


process.corrPfMetType0PfCand = cms.EDProducer("Type0PFMETcorrInputProducer",
    correction = cms.PSet(
        formula = cms.string('(x<35)?(-( [0]+x*[1]+pow(x, 2)*[2]+pow(x, 3)*[3] )):(-( [0]+35*[1]+pow(35, 2)*[2]+pow(35, 3)*[3] ))'),
        par0 = cms.double(-0.181414),
        par1 = cms.double(-0.476934),
        par2 = cms.double(0.00863564),
        par3 = cms.double(-4.94181e-05)
    ),
    minDz = cms.double(0.2),
    srcHardScatterVertex = cms.InputTag("selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0"),
    srcPFCandidateToVertexAssociations = cms.InputTag("pfCandidateToVertexAssociation")
)


process.corrPfMetType1 = cms.EDProducer("PFJetMETcorrInputProducer",
    jetCorrEtaMax = cms.double(9.9),
    jetCorrLabel = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrLabelRes = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector"),
    offsetCorrLabel = cms.InputTag("ak4PFCHSL1FastjetCorrector"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("ak4PFJetsCHS"),
    type1JetPtThreshold = cms.double(15.0)
)


process.corrPfMetType2 = cms.EDProducer("Type2CorrectionProducer",
    srcUnclEnergySums = cms.VInputTag(cms.InputTag("corrPfMetType1","type2"), cms.InputTag("corrPfMetType1","offset"), cms.InputTag("pfCandMETcorr")),
    type2CorrFormula = cms.string('A'),
    type2CorrParameter = cms.PSet(
        A = cms.double(1.4)
    )
)


process.egmGsfElectronIDs = cms.EDProducer("VersionedGsfElectronIdProducer",
    physicsObjectIDs = cms.VPSet(cms.PSet(
        idDefinition = cms.PSet(
            cutFlow = cms.VPSet(cms.PSet(
                cutName = cms.string('MinPtCut'),
                isIgnored = cms.bool(False),
                minPt = cms.double(5.0),
                needsAdditionalProducts = cms.bool(False)
            ), 
                cms.PSet(
                    allowedEtaRanges = cms.VPSet(cms.PSet(
                        maxEta = cms.double(1.479),
                        minEta = cms.double(0.0)
                    ), 
                        cms.PSet(
                            maxEta = cms.double(2.5),
                            minEta = cms.double(1.479)
                        )),
                    cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False),
                    useAbsEta = cms.bool(True)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleDEtaInSeedCut'),
                    dEtaInSeedCutValueEB = cms.double(0.00477),
                    dEtaInSeedCutValueEE = cms.double(0.00868),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleDPhiInCut'),
                    dPhiInCutValueEB = cms.double(0.222),
                    dPhiInCutValueEE = cms.double(0.213),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut'),
                    full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.011),
                    full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0314),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleHadronicOverEMCut'),
                    hadronicOverEMCutValueEB = cms.double(0.298),
                    hadronicOverEMCutValueEE = cms.double(0.101),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                    eInverseMinusPInverseCutValueEB = cms.double(0.241),
                    eInverseMinusPInverseCutValueEE = cms.double(0.14),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleEffAreaPFIsoCut'),
                    effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
                    isIgnored = cms.bool(False),
                    isRelativeIso = cms.bool(True),
                    isoCutEBHighPt = cms.double(0.0994),
                    isoCutEBLowPt = cms.double(0.0994),
                    isoCutEEHighPt = cms.double(0.107),
                    isoCutEELowPt = cms.double(0.107),
                    needsAdditionalProducts = cms.bool(True),
                    ptCutOff = cms.double(20.0),
                    rho = cms.InputTag("fixedGridRhoFastjetAll")
                ), 
                cms.PSet(
                    beamspotSrc = cms.InputTag("offlineBeamSpot"),
                    conversionSrc = cms.InputTag("allConversions"),
                    conversionSrcMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
                    cutName = cms.string('GsfEleConversionVetoCut'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('GsfEleMissingHitsCut'),
                    isIgnored = cms.bool(False),
                    maxMissingHitsEB = cms.uint32(1),
                    maxMissingHitsEE = cms.uint32(1),
                    needsAdditionalProducts = cms.bool(False)
                )),
            idName = cms.string('cutBasedElectronID-Summer16-80X-V1-loose')
        ),
        idMD5 = cms.string('c1c4c739f1ba0791d40168c123183475'),
        isPOGApproved = cms.untracked.bool(False)
    ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDEtaInSeedCut'),
                        dEtaInSeedCutValueEB = cms.double(0.00311),
                        dEtaInSeedCutValueEE = cms.double(0.00609),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDPhiInCut'),
                        dPhiInCutValueEB = cms.double(0.103),
                        dPhiInCutValueEE = cms.double(0.045),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut'),
                        full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.00998),
                        full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0298),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleHadronicOverEMCut'),
                        hadronicOverEMCutValueEB = cms.double(0.253),
                        hadronicOverEMCutValueEE = cms.double(0.0878),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                        eInverseMinusPInverseCutValueEB = cms.double(0.134),
                        eInverseMinusPInverseCutValueEE = cms.double(0.13),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEffAreaPFIsoCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
                        isIgnored = cms.bool(False),
                        isRelativeIso = cms.bool(True),
                        isoCutEBHighPt = cms.double(0.0695),
                        isoCutEBLowPt = cms.double(0.0695),
                        isoCutEEHighPt = cms.double(0.0821),
                        isoCutEELowPt = cms.double(0.0821),
                        needsAdditionalProducts = cms.bool(True),
                        ptCutOff = cms.double(20.0),
                        rho = cms.InputTag("fixedGridRhoFastjetAll")
                    ), 
                    cms.PSet(
                        beamspotSrc = cms.InputTag("offlineBeamSpot"),
                        conversionSrc = cms.InputTag("allConversions"),
                        conversionSrcMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
                        cutName = cms.string('GsfEleConversionVetoCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleMissingHitsCut'),
                        isIgnored = cms.bool(False),
                        maxMissingHitsEB = cms.uint32(1),
                        maxMissingHitsEE = cms.uint32(1),
                        needsAdditionalProducts = cms.bool(False)
                    )),
                idName = cms.string('cutBasedElectronID-Summer16-80X-V1-medium')
            ),
            idMD5 = cms.string('71b43f74a27d2fd3d27416afd22e8692'),
            isPOGApproved = cms.untracked.bool(False)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDEtaInSeedCut'),
                        dEtaInSeedCutValueEB = cms.double(0.00308),
                        dEtaInSeedCutValueEE = cms.double(0.00605),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDPhiInCut'),
                        dPhiInCutValueEB = cms.double(0.0816),
                        dPhiInCutValueEE = cms.double(0.0394),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut'),
                        full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.00998),
                        full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.0292),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleHadronicOverEMCut'),
                        hadronicOverEMCutValueEB = cms.double(0.0414),
                        hadronicOverEMCutValueEE = cms.double(0.0641),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                        eInverseMinusPInverseCutValueEB = cms.double(0.0129),
                        eInverseMinusPInverseCutValueEE = cms.double(0.0129),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEffAreaPFIsoCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
                        isIgnored = cms.bool(False),
                        isRelativeIso = cms.bool(True),
                        isoCutEBHighPt = cms.double(0.0588),
                        isoCutEBLowPt = cms.double(0.0588),
                        isoCutEEHighPt = cms.double(0.0571),
                        isoCutEELowPt = cms.double(0.0571),
                        needsAdditionalProducts = cms.bool(True),
                        ptCutOff = cms.double(20.0),
                        rho = cms.InputTag("fixedGridRhoFastjetAll")
                    ), 
                    cms.PSet(
                        beamspotSrc = cms.InputTag("offlineBeamSpot"),
                        conversionSrc = cms.InputTag("allConversions"),
                        conversionSrcMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
                        cutName = cms.string('GsfEleConversionVetoCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleMissingHitsCut'),
                        isIgnored = cms.bool(False),
                        maxMissingHitsEB = cms.uint32(1),
                        maxMissingHitsEE = cms.uint32(1),
                        needsAdditionalProducts = cms.bool(False)
                    )),
                idName = cms.string('cutBasedElectronID-Summer16-80X-V1-tight')
            ),
            idMD5 = cms.string('ca2a9db2976d80ba2c13f9bfccdc32f2'),
            isPOGApproved = cms.untracked.bool(False)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDEtaInSeedCut'),
                        dEtaInSeedCutValueEB = cms.double(0.00749),
                        dEtaInSeedCutValueEE = cms.double(0.00895),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDPhiInCut'),
                        dPhiInCutValueEB = cms.double(0.228),
                        dPhiInCutValueEE = cms.double(0.213),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut'),
                        full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.0115),
                        full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.037),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleHadronicOverEMCut'),
                        hadronicOverEMCutValueEB = cms.double(0.356),
                        hadronicOverEMCutValueEE = cms.double(0.211),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                        eInverseMinusPInverseCutValueEB = cms.double(0.299),
                        eInverseMinusPInverseCutValueEE = cms.double(0.15),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEffAreaPFIsoCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
                        isIgnored = cms.bool(False),
                        isRelativeIso = cms.bool(True),
                        isoCutEBHighPt = cms.double(0.175),
                        isoCutEBLowPt = cms.double(0.175),
                        isoCutEEHighPt = cms.double(0.159),
                        isoCutEELowPt = cms.double(0.159),
                        needsAdditionalProducts = cms.bool(True),
                        ptCutOff = cms.double(20.0),
                        rho = cms.InputTag("fixedGridRhoFastjetAll")
                    ), 
                    cms.PSet(
                        beamspotSrc = cms.InputTag("offlineBeamSpot"),
                        conversionSrc = cms.InputTag("allConversions"),
                        conversionSrcMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
                        cutName = cms.string('GsfEleConversionVetoCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleMissingHitsCut'),
                        isIgnored = cms.bool(False),
                        maxMissingHitsEB = cms.uint32(2),
                        maxMissingHitsEE = cms.uint32(3),
                        needsAdditionalProducts = cms.bool(False)
                    )),
                idName = cms.string('cutBasedElectronID-Summer16-80X-V1-veto')
            ),
            idMD5 = cms.string('0025c1841da1ab64a08d703ded72409b'),
            isPOGApproved = cms.untracked.bool(False)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaCut'),
                        full5x5SigmaIEtaIEtaCutValueEB = cms.double(0.011),
                        full5x5SigmaIEtaIEtaCutValueEE = cms.double(0.031),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDEtaInSeedCut'),
                        dEtaInSeedCutValueEB = cms.double(0.004),
                        dEtaInSeedCutValueEE = cms.double(1e+30),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDPhiInCut'),
                        dPhiInCutValueEB = cms.double(0.02),
                        dPhiInCutValueEE = cms.double(1e+30),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleHadronicOverEMCut'),
                        hadronicOverEMCutValueEB = cms.double(0.06),
                        hadronicOverEMCutValueEE = cms.double(0.065),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEInverseMinusPInverseCut'),
                        eInverseMinusPInverseCutValueEB = cms.double(0.013),
                        eInverseMinusPInverseCutValueEE = cms.double(0.013),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleCalPFClusterIsoCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_HLT_ecalPFClusterIso.txt'),
                        isIgnored = cms.bool(False),
                        isRelativeIso = cms.bool(True),
                        isoCutEBHighPt = cms.double(0.16),
                        isoCutEBLowPt = cms.double(0.16),
                        isoCutEEHighPt = cms.double(0.12),
                        isoCutEELowPt = cms.double(0.12),
                        isoType = cms.int32(0),
                        needsAdditionalProducts = cms.bool(True),
                        ptCutOff = cms.double(20.0),
                        rho = cms.InputTag("fixedGridRhoFastjetCentralCalo")
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleCalPFClusterIsoCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_HLT_hcalPFClusterIso.txt'),
                        isIgnored = cms.bool(False),
                        isRelativeIso = cms.bool(True),
                        isoCutEBHighPt = cms.double(0.12),
                        isoCutEBLowPt = cms.double(0.12),
                        isoCutEEHighPt = cms.double(0.12),
                        isoCutEELowPt = cms.double(0.12),
                        isoType = cms.int32(1),
                        needsAdditionalProducts = cms.bool(True),
                        ptCutOff = cms.double(20.0),
                        rho = cms.InputTag("fixedGridRhoFastjetCentralCalo")
                    ), 
                    cms.PSet(
                        constTermEB = cms.double(0.0),
                        constTermEE = cms.double(0.0),
                        cutName = cms.string('GsfEleTrkPtIsoCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        slopeStartEB = cms.double(0.0),
                        slopeStartEE = cms.double(0.0),
                        slopeTermEB = cms.double(0.08),
                        slopeTermEE = cms.double(0.08)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleNormalizedGsfChi2Cut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        normalizedGsfChi2CutValueEB = cms.double(1e+30),
                        normalizedGsfChi2CutValueEE = cms.double(3.0)
                    )),
                idName = cms.string('cutBasedElectronHLTPreselection-Summer16-V1')
            ),
            idMD5 = cms.string('aef5f00cc25a63b44192c6fc449f7dc0'),
            isPOGApproved = cms.untracked.bool(False)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(35.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.4442),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.566)
                            )),
                        cutName = cms.string('GsfEleSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDEtaInSeedCut'),
                        dEtaInSeedCutValueEB = cms.double(0.004),
                        dEtaInSeedCutValueEE = cms.double(0.006),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDPhiInCut'),
                        dPhiInCutValueEB = cms.double(0.06),
                        dPhiInCutValueEE = cms.double(0.06),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        cutName = cms.string('GsfEleFull5x5SigmaIEtaIEtaWithSatCut'),
                        isIgnored = cms.bool(False),
                        maxNrSatCrysIn5x5EB = cms.int32(0),
                        maxNrSatCrysIn5x5EE = cms.int32(0),
                        maxSigmaIEtaIEtaEB = cms.double(9999),
                        maxSigmaIEtaIEtaEE = cms.double(0.03),
                        needsAdditionalProducts = cms.bool(True),
                        nrSatCrysValueMap = cms.InputTag("heepIDVarValueMaps","eleNrSaturateIn5x5")
                    ), 
                    cms.PSet(
                        cutName = cms.string('GsfEleFull5x5E2x5OverE5x5WithSatCut'),
                        isIgnored = cms.bool(False),
                        maxNrSatCrysIn5x5EB = cms.int32(0),
                        maxNrSatCrysIn5x5EE = cms.int32(0),
                        minE1x5OverE5x5EB = cms.double(0.83),
                        minE1x5OverE5x5EE = cms.double(-1.0),
                        minE2x5OverE5x5EB = cms.double(0.94),
                        minE2x5OverE5x5EE = cms.double(-1.0),
                        needsAdditionalProducts = cms.bool(True),
                        nrSatCrysValueMap = cms.InputTag("heepIDVarValueMaps","eleNrSaturateIn5x5")
                    ), 
                    cms.PSet(
                        constTermEB = cms.double(1.0),
                        constTermEE = cms.double(5),
                        cutName = cms.string('GsfEleHadronicOverEMLinearCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        slopeStartEB = cms.double(0.0),
                        slopeStartEE = cms.double(0.0),
                        slopeTermEB = cms.double(0.05),
                        slopeTermEE = cms.double(0.05)
                    ), 
                    cms.PSet(
                        constTermEB = cms.double(5.0),
                        constTermEE = cms.double(5.0),
                        cutName = cms.string('GsfEleValueMapIsoRhoCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag(""),
                        rhoEAEB = cms.double(0.0),
                        rhoEAEE = cms.double(0.0),
                        rhoEtStartEB = cms.double(999999.0),
                        rhoEtStartEE = cms.double(999999.0),
                        slopeStartEB = cms.double(0.0),
                        slopeStartEE = cms.double(0.0),
                        slopeTermEB = cms.double(0.0),
                        slopeTermEE = cms.double(0.0),
                        value = cms.InputTag("heepIDVarValueMaps","eleTrkPtIso")
                    ), 
                    cms.PSet(
                        constTermEB = cms.double(2.0),
                        constTermEE = cms.double(2.5),
                        cutName = cms.string('GsfEleEmHadD1IsoRhoCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        rhoConstant = cms.double(0.28),
                        slopeStartEB = cms.double(0.0),
                        slopeStartEE = cms.double(50.0),
                        slopeTermEB = cms.double(0.03),
                        slopeTermEE = cms.double(0.03)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleDxyCut'),
                        dxyCutValueEB = cms.double(0.02),
                        dxyCutValueEE = cms.double(0.05),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        vertexSrc = cms.InputTag("offlinePrimaryVertices"),
                        vertexSrcMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleMissingHitsCut'),
                        isIgnored = cms.bool(False),
                        maxMissingHitsEB = cms.uint32(1),
                        maxMissingHitsEE = cms.uint32(1),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('GsfEleEcalDrivenCut'),
                        ecalDrivenEB = cms.int32(1),
                        ecalDrivenEE = cms.int32(1),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    )),
                idName = cms.string('heepElectronID-HEEPV70')
            ),
            idMD5 = cms.string('45d49ea5f46f3f13f579d208e5e3412b'),
            isPOGApproved = cms.untracked.bool(False)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
                    mvaCuts = cms.vdouble(0.940962684155, 0.899208843708, 0.758484721184),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Spring16-GeneralPurpose-V1-wp80')
            ),
            idMD5 = cms.string('b490bc0b0af2d5f3e9efea562370af2a'),
            isPOGApproved = cms.untracked.bool(False)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
                    mvaCuts = cms.vdouble(0.836695742607, 0.715337944031, 0.356799721718),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Spring16-GeneralPurpose-V1-wp90')
            ),
            idMD5 = cms.string('14c153aaf3c207deb3ad4932586647a7'),
            isPOGApproved = cms.untracked.bool(False)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16HZZV1Categories"),
                    mvaCuts = cms.vdouble(-0.211, -0.396, -0.215, -0.87, -0.838, 
                        -0.763),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16HZZV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Spring16-HZZ-V1-wpLoose')
            ),
            idMD5 = cms.string('1797cc03eb62387e10266fca72ea10cd'),
            isPOGApproved = cms.untracked.bool(False)
        )),
    physicsObjectSrc = cms.InputTag("slimmedElectrons")
)


process.egmPhotonIDs = cms.EDProducer("VersionedPhotonIdProducer",
    physicsObjectIDs = cms.VPSet(cms.PSet(
        idDefinition = cms.PSet(
            cutFlow = cms.VPSet(cms.PSet(
                cutName = cms.string('MinPtCut'),
                isIgnored = cms.bool(False),
                minPt = cms.double(5.0),
                needsAdditionalProducts = cms.bool(False)
            ), 
                cms.PSet(
                    allowedEtaRanges = cms.VPSet(cms.PSet(
                        maxEta = cms.double(1.479),
                        minEta = cms.double(0.0)
                    ), 
                        cms.PSet(
                            maxEta = cms.double(2.5),
                            minEta = cms.double(1.479)
                        )),
                    cutName = cms.string('PhoSCEtaMultiRangeCut'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False),
                    useAbsEta = cms.bool(True)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('PhoSingleTowerHadOverEmCut'),
                    hadronicOverEMCutValueEB = cms.double(0.0597),
                    hadronicOverEMCutValueEE = cms.double(0.0481),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('PhoFull5x5SigmaIEtaIEtaCut'),
                    cutValueEB = cms.double(0.01031),
                    cutValueEE = cms.double(0.03013),
                    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                cms.PSet(
                    C1_EB = cms.double(1.295),
                    C1_EE = cms.double(1.011),
                    C2_EB = cms.double(0),
                    C2_EE = cms.double(0.0),
                    anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('PhoAnyPFIsoWithEACut'),
                    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfChargedHadrons_90percentBased.txt'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True),
                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                    useRelativeIso = cms.bool(False)
                ), 
                cms.PSet(
                    C1_EB = cms.double(10.91),
                    C1_EE = cms.double(5.931),
                    C2_EB = cms.double(0.0148),
                    C2_EE = cms.double(0.0163),
                    C3_EB = cms.double(1.7e-05),
                    C3_EE = cms.double(1.4e-05),
                    anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('PhoAnyPFIsoWithEAAndQuadScalingCut'),
                    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased.txt'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True),
                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                    useRelativeIso = cms.bool(False)
                ), 
                cms.PSet(
                    C1_EB = cms.double(3.63),
                    C1_EE = cms.double(6.641),
                    C2_EB = cms.double(0.0047),
                    C2_EE = cms.double(0.0034),
                    anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
                    barrelCutOff = cms.double(1.479),
                    cutName = cms.string('PhoAnyPFIsoWithEACut'),
                    effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt'),
                    isIgnored = cms.bool(False),
                    needsAdditionalProducts = cms.bool(True),
                    rho = cms.InputTag("fixedGridRhoFastjetAll"),
                    useRelativeIso = cms.bool(False)
                )),
            idName = cms.string('cutBasedPhotonID-Spring16-V2p2-loose')
        ),
        idMD5 = cms.string('d6ce6a4f3476294bf0a3261e00170daf'),
        isPOGApproved = cms.untracked.bool(False)
    ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('PhoSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoSingleTowerHadOverEmCut'),
                        hadronicOverEMCutValueEB = cms.double(0.0396),
                        hadronicOverEMCutValueEE = cms.double(0.0219),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoFull5x5SigmaIEtaIEtaCut'),
                        cutValueEB = cms.double(0.01022),
                        cutValueEE = cms.double(0.03001),
                        full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(0.441),
                        C1_EE = cms.double(0.442),
                        C2_EB = cms.double(0.0),
                        C2_EE = cms.double(0.0),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEACut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfChargedHadrons_90percentBased.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(2.725),
                        C1_EE = cms.double(1.715),
                        C2_EB = cms.double(0.0148),
                        C2_EE = cms.double(0.0163),
                        C3_EB = cms.double(1.7e-05),
                        C3_EE = cms.double(1.4e-05),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEAAndQuadScalingCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(2.571),
                        C1_EE = cms.double(3.863),
                        C2_EB = cms.double(0.0047),
                        C2_EE = cms.double(0.0034),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEACut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    )),
                idName = cms.string('cutBasedPhotonID-Spring16-V2p2-medium')
            ),
            idMD5 = cms.string('c739cfd0b6287b8586da187c06d4053f'),
            isPOGApproved = cms.untracked.bool(False)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('MinPtCut'),
                    isIgnored = cms.bool(False),
                    minPt = cms.double(5.0),
                    needsAdditionalProducts = cms.bool(False)
                ), 
                    cms.PSet(
                        allowedEtaRanges = cms.VPSet(cms.PSet(
                            maxEta = cms.double(1.479),
                            minEta = cms.double(0.0)
                        ), 
                            cms.PSet(
                                maxEta = cms.double(2.5),
                                minEta = cms.double(1.479)
                            )),
                        cutName = cms.string('PhoSCEtaMultiRangeCut'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False),
                        useAbsEta = cms.bool(True)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoSingleTowerHadOverEmCut'),
                        hadronicOverEMCutValueEB = cms.double(0.0269),
                        hadronicOverEMCutValueEE = cms.double(0.0213),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoFull5x5SigmaIEtaIEtaCut'),
                        cutValueEB = cms.double(0.00994),
                        cutValueEE = cms.double(0.03),
                        full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(0.202),
                        C1_EE = cms.double(0.034),
                        C2_EB = cms.double(0.0),
                        C2_EE = cms.double(0.0),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEACut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfChargedHadrons_90percentBased.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(0.264),
                        C1_EE = cms.double(0.586),
                        C2_EB = cms.double(0.0148),
                        C2_EE = cms.double(0.0163),
                        C3_EB = cms.double(1.7e-05),
                        C3_EE = cms.double(1.4e-05),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEAAndQuadScalingCut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    ), 
                    cms.PSet(
                        C1_EB = cms.double(2.362),
                        C1_EE = cms.double(2.617),
                        C2_EB = cms.double(0.0047),
                        C2_EE = cms.double(0.0034),
                        anyPFIsoMap = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
                        barrelCutOff = cms.double(1.479),
                        cutName = cms.string('PhoAnyPFIsoWithEACut'),
                        effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt'),
                        isIgnored = cms.bool(False),
                        needsAdditionalProducts = cms.bool(True),
                        rho = cms.InputTag("fixedGridRhoFastjetAll"),
                        useRelativeIso = cms.bool(False)
                    )),
                idName = cms.string('cutBasedPhotonID-Spring16-V2p2-tight')
            ),
            idMD5 = cms.string('bdb623bdb1a15c13545020a919dd9530'),
            isPOGApproved = cms.untracked.bool(False)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('PhoMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Categories"),
                    mvaCuts = cms.vdouble(0.68, 0.6),
                    mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaPhoID-Spring16-nonTrig-V1-wp80')
            ),
            idMD5 = cms.string('beb95233f7d1e033ad9e20cf3d804ba0'),
            isPOGApproved = cms.untracked.bool(False)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('PhoMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Categories"),
                    mvaCuts = cms.vdouble(0.2, 0.2),
                    mvaValueMapName = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaPhoID-Spring16-nonTrig-V1-wp90')
            ),
            idMD5 = cms.string('36efe663348f95de0bc1cfa8dc7fa8fe'),
            isPOGApproved = cms.untracked.bool(False)
        )),
    physicsObjectSrc = cms.InputTag("slimmedPhotons")
)


process.egmPhotonIsolation = cms.EDProducer("CITKPFIsolationSumProducer",
    isolationConeDefinitions = cms.VPSet(cms.PSet(
        coneSize = cms.double(0.3),
        isolateAgainst = cms.string('h+'),
        isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
        miniAODVertexCodes = cms.vuint32(2, 3),
        vertexIndex = cms.int32(0)
    ), 
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('h0'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            vertexIndex = cms.int32(0)
        ), 
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('gamma'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            vertexIndex = cms.int32(0)
        )),
    srcForIsolationCone = cms.InputTag("packedPFCandidates"),
    srcToIsolate = cms.InputTag("slimmedPhotons")
)


process.egmPhotonIsolationMiniAOD = cms.EDProducer("CITKPFIsolationSumProducer",
    isolationConeDefinitions = cms.VPSet(cms.PSet(
        coneSize = cms.double(0.3),
        isolateAgainst = cms.string('h+'),
        isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
        miniAODVertexCodes = cms.vuint32(2, 3),
        vertexIndex = cms.int32(0)
    ), 
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('h0'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            vertexIndex = cms.int32(0)
        ), 
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('gamma'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            vertexIndex = cms.int32(0)
        )),
    srcForIsolationCone = cms.InputTag("packedPFCandidates"),
    srcToIsolate = cms.InputTag("slimmedPhotons")
)


process.egmPhotonIsolationMiniAODPUPPI = cms.EDProducer("CITKPFIsolationSumProducerForPUPPI",
    isolationConeDefinitions = cms.VPSet(cms.PSet(
        coneSize = cms.double(0.3),
        isolateAgainst = cms.string('h+'),
        isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
        miniAODVertexCodes = cms.vuint32(2, 3),
        vertexIndex = cms.int32(0)
    ), 
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('h0'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            vertexIndex = cms.int32(0)
        ), 
        cms.PSet(
            coneSize = cms.double(0.3),
            isolateAgainst = cms.string('gamma'),
            isolationAlgo = cms.string('PhotonPFIsolationWithMapBasedVeto'),
            miniAODVertexCodes = cms.vuint32(2, 3),
            vertexIndex = cms.int32(0)
        )),
    puppiValueMap = cms.InputTag("puppi"),
    srcForIsolationCone = cms.InputTag("packedPFCandidates"),
    srcToIsolate = cms.InputTag("slimmedPhotons")
)


process.elPFIsoDepositChargedAllPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositChargedAllPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedElectronsPFBRECO"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositChargedPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositChargedPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedElectronsPFBRECO"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositGammaPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFCandWithSuperClusterExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        MissHitSCMatch_Veto = cms.bool(True),
        SCMatch_Veto = cms.bool(False),
        inputCandView = cms.InputTag("pfAllPhotonsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositGammaPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFCandWithSuperClusterExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        MissHitSCMatch_Veto = cms.bool(True),
        SCMatch_Veto = cms.bool(False),
        inputCandView = cms.InputTag("pfAllPhotonsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedElectronsPFBRECO"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositNeutralPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositNeutralPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedElectronsPFBRECO"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositPUPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons"),
    trackType = cms.string('candidate')
)


process.elPFIsoDepositPUPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedElectronsPFBRECO"),
    trackType = cms.string('candidate')
)


process.elPFIsoValueCharged03NoPFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositCharged"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged03NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged03NoPFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged03PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositCharged"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged03PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged04NoPFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositCharged"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged04NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged04NoPFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged04PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositCharged"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueCharged04PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll03NoPFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAll"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll03NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll03NoPFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll03PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAll"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll03PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll04NoPFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAll"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll04NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll04NoPFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll04PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAll"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueChargedAll04PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma03NoPFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGamma"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma03NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma03NoPFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma03PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGamma"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma03PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma04NoPFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGamma"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma04NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma04NoPFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma04PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGamma"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueGamma04PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.08)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral03NoPFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutral"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral03NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral03NoPFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral03PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutral"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral03PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral04NoPFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutral"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral04NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral04NoPFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral04PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutral"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValueNeutral04PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU03NoPFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPU"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU03NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU03NoPFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU03PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPU"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU03PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU04NoPFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPU"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU04NoPFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU04NoPFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU04PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPU"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.elPFIsoValuePU04PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("elPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.015)'),
        weight = cms.string('1')
    ))
)


process.electronIDValueMapProducer = cms.EDProducer("ElectronIDValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    src = cms.InputTag("gedGsfElectrons"),
    srcMiniAOD = cms.InputTag("slimmedElectrons")
)


process.electronMVAValueMapProducer = cms.EDProducer("ElectronMVAValueMapProducer",
    mvaConfigurations = cms.VPSet(cms.PSet(
        mvaName = cms.string('ElectronMVAEstimatorRun2Phys14NonTrig'),
        mvaTag = cms.string('25nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml')
    ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring15NonTrig'),
            mvaTag = cms.string('25nsV1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
            mvaTag = cms.string('50nsV1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
            mvaTag = cms.string('25nsV1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring16HZZ'),
            mvaTag = cms.string('V1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring16GeneralPurpose'),
            mvaTag = cms.string('V1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml')
        )),
    src = cms.InputTag("gedGsfElectrons"),
    srcMiniAOD = cms.InputTag("slimmedElectrons")
)


process.electronMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(True),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(0.5),
    maxDeltaR = cms.double(0.5),
    mcPdgId = cms.vint32(11),
    mcStatus = cms.vint32(1),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("gedGsfElectrons")
)


process.electronRegressionValueMapProducer = cms.EDProducer("ElectronRegressionValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    src = cms.InputTag("gedGsfElectrons"),
    srcMiniAOD = cms.InputTag("slimmedElectrons","","@skipCurrentProcess"),
    useFull5x5 = cms.bool(False)
)


process.genMetExtractor = cms.EDProducer("GenMETExtractor",
    metSource = cms.InputTag("slimmedMETs","","@skipCurrentProcess")
)


process.heepIDVarValueMaps = cms.EDProducer("ElectronHEEPIDValueMapProducer",
    beamSpot = cms.InputTag("offlineBeamSpot"),
    candsAOD = cms.VInputTag("packedCandsForTkIso", "lostTracksForTkIso"),
    candsMiniAOD = cms.VInputTag("packedPFCandidates", "lostTracks"),
    dataFormat = cms.int32(0),
    ebRecHitsAOD = cms.InputTag("reducedEcalRecHitsEB"),
    ebRecHitsMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeRecHitsAOD = cms.InputTag("reducedEcalRecHitsEB"),
    eeRecHitsMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    elesAOD = cms.InputTag("gedGsfElectrons"),
    elesMiniAOD = cms.InputTag("slimmedElectrons"),
    trkIsoConfig = cms.PSet(
        barrelCuts = cms.PSet(
            algosToReject = cms.vstring(),
            allowedQualities = cms.vstring(),
            maxDPtPt = cms.double(0.1),
            maxDR = cms.double(0.3),
            maxDZ = cms.double(0.1),
            minDEta = cms.double(0.005),
            minDR = cms.double(0.0),
            minHits = cms.int32(8),
            minPixelHits = cms.int32(1),
            minPt = cms.double(1.0)
        ),
        endcapCuts = cms.PSet(
            algosToReject = cms.vstring(),
            allowedQualities = cms.vstring(),
            maxDPtPt = cms.double(0.1),
            maxDR = cms.double(0.3),
            maxDZ = cms.double(0.5),
            minDEta = cms.double(0.005),
            minDR = cms.double(0.0),
            minHits = cms.int32(8),
            minPixelHits = cms.int32(1),
            minPt = cms.double(1.0)
        )
    )
)


process.hpsPFTauChargedIsoPtSum = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawPUsumPt = cms.bool(False),
    storeRawSumPt = cms.bool(True),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByLoosePileupWeightedIsolation3Hits = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(True),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(True),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByMediumPileupWeightedIsolation3Hits = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(True),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(True),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(1.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(True),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(True),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawSumPt = cms.bool(True),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByRawPileupWeightedIsolation3Hits = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(True),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawSumPt = cms.bool(True),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauDiscriminationByTightPileupWeightedIsolation3Hits = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(True),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(True),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(True),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(0.8),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauFootprintCorrection = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawFootprintCorrection = cms.bool(True),
    storeRawPUsumPt = cms.bool(False),
    storeRawSumPt = cms.bool(False),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauNeutralIsoPtSum = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(True),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawPUsumPt = cms.bool(False),
    storeRawSumPt = cms.bool(True),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauNeutralIsoPtSumWeight = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(True),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(True),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawPUsumPt = cms.bool(False),
    storeRawSumPt = cms.bool(True),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauPUcorrPtSum = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(True),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawPUsumPt = cms.bool(True),
    storeRawSumPt = cms.bool(False),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.hpsPFTauPhotonPtSumOutsideSignalCone = cms.EDProducer("PFRecoTauDiscriminationByIsolation",
    ApplyDiscriminationByECALIsolation = cms.bool(False),
    ApplyDiscriminationByTrackerIsolation = cms.bool(False),
    ApplyDiscriminationByWeightedECALIsolation = cms.bool(False),
    PFTauProducer = cms.InputTag("hpsPFTauProducer"),
    Prediscriminants = cms.PSet(
        BooleanOperator = cms.string('and'),
        decayMode = cms.PSet(
            Producer = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
            cut = cms.double(0.5)
        )
    ),
    UseAllPFCandsForWeights = cms.bool(False),
    WeightECALIsolation = cms.double(1.0),
    applyDeltaBetaCorrection = cms.bool(False),
    applyFootprintCorrection = cms.bool(False),
    applyOccupancyCut = cms.bool(False),
    applyPhotonPtSumOutsideSignalConeCut = cms.bool(True),
    applyRelativeSumPtCut = cms.bool(False),
    applyRhoCorrection = cms.bool(False),
    applySumPtCut = cms.bool(False),
    customOuterCone = cms.double(0.5),
    deltaBetaFactor = cms.string('0.2000'),
    deltaBetaPUTrackPtCutOverride = cms.double(0.5),
    footprintCorrections = cms.VPSet(cms.PSet(
        offset = cms.string('0.0'),
        selection = cms.string('decayMode() = 0')
    ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 1 || decayMode() = 2')
        ), 
        cms.PSet(
            offset = cms.string('2.7'),
            selection = cms.string('decayMode() = 5')
        ), 
        cms.PSet(
            offset = cms.string('0.0'),
            selection = cms.string('decayMode() = 6')
        ), 
        cms.PSet(
            offset = cms.string('max(2.0, 0.22*pt() - 2.0)'),
            selection = cms.string('decayMode() = 10')
        )),
    isoConeSizeForDeltaBeta = cms.double(0.8),
    maxAbsPhotonSumPt_outsideSignalCone = cms.double(1000000000.0),
    maxRelPhotonSumPt_outsideSignalCone = cms.double(0.1),
    maximumOccupancy = cms.uint32(0),
    maximumSumPtCut = cms.double(2.5),
    minTauPtForNoIso = cms.double(-99.0),
    particleFlowSrc = cms.InputTag("particleFlow"),
    qualityCuts = cms.PSet(
        isolationQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.2),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.03),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        leadingTrkOrPFCandOption = cms.string('leadPFCand'),
        primaryVertexSrc = cms.InputTag("offlinePrimaryVertices"),
        pvFindingAlgo = cms.string('closestInDeltaZ'),
        recoverLeadingTrk = cms.bool(False),
        signalQualityCuts = cms.PSet(
            maxDeltaZ = cms.double(0.4),
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minNeutralHadronEt = cms.double(30.0),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        ),
        vertexTrackFiltering = cms.bool(False),
        vxAssocQualityCuts = cms.PSet(
            maxTrackChi2 = cms.double(100.0),
            maxTransverseImpactParameter = cms.double(0.1),
            minGammaEt = cms.double(0.5),
            minTrackHits = cms.uint32(3),
            minTrackPixelHits = cms.uint32(0),
            minTrackPt = cms.double(0.5),
            minTrackVertexWeight = cms.double(-1.0)
        )
    ),
    relativeSumPtCut = cms.double(0.0),
    relativeSumPtOffset = cms.double(0.0),
    rhoConeSize = cms.double(0.5),
    rhoProducer = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoUEOffsetCorrection = cms.double(1.0),
    storeRawPUsumPt = cms.bool(False),
    storeRawPhotonSumPt_outsideSignalCone = cms.bool(True),
    storeRawSumPt = cms.bool(False),
    verbosity = cms.int32(0),
    vertexSrc = cms.InputTag("offlinePrimaryVertices")
)


process.isoDeposits = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag(""),
    trackType = cms.string('candidate')
)


process.metrawCalo = cms.EDProducer("RecoMETExtractor",
    correctionLevel = cms.string('rawCalo'),
    metSource = cms.InputTag("slimmedMETs","","@skipCurrentProcess")
)


process.modifiedElectrons = cms.EDProducer("ModifiedElectronProducer",
    modifierConfig = cms.PSet(
        modifications = cms.VPSet()
    ),
    src = cms.InputTag("slimmedElectrons","","@skipCurrentProcess")
)


process.muCaloMetCorr = cms.EDProducer("MuonMETcorrInputProducer",
    src = cms.InputTag("muons"),
    srcMuonCorrections = cms.InputTag("muonMETValueMapProducer","muCorrData")
)


process.muPFIsoDepositChargedAllPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("muons"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositChargedAllPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedMuonsPFBRECO"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositChargedPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("muons"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositChargedPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedMuonsPFBRECO"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositGammaPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotonsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("muons"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositGammaPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotonsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedMuonsPFBRECO"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositNeutralPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("muons"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositNeutralPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedMuonsPFBRECO"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositPUPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("muons"),
    trackType = cms.string('candidate')
)


process.muPFIsoDepositPUPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedMuonsPFBRECO"),
    trackType = cms.string('candidate')
)


process.muPFIsoValueCharged03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositCharged"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueCharged03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueCharged03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueCharged04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositCharged"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueCharged04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueCharged04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueChargedAll03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAll"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueChargedAll03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueChargedAll03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueChargedAll04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAll"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueChargedAll04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueChargedAll04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGamma03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGamma03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGamma03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGamma04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGamma04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGamma04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGammaHighThreshold03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGammaHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGammaHighThreshold03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGammaHighThreshold04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGammaHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueGammaHighThreshold04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutral03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutral03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutral03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutral04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutral04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutral04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutralHighThreshold03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutralHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutralHighThreshold03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutralHighThreshold04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutralHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValueNeutralHighThreshold04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValuePU03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPU"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValuePU03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValuePU03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValuePU04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPU"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValuePU04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFIsoValuePU04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueCharged03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositCharged"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueCharged03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueCharged03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueCharged04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositCharged"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueCharged04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueCharged04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueChargedAll03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAll"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueChargedAll03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueChargedAll03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueChargedAll04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAll"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueChargedAll04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueChargedAll04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGamma03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGamma03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGamma03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGamma04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGamma04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGamma04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGammaHighThreshold03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGammaHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGammaHighThreshold03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGammaHighThreshold04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGammaHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueGammaHighThreshold04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutral03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutral03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutral03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutral04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutral04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutral04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutralHighThreshold03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutralHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutralHighThreshold03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutralHighThreshold04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutralHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValueNeutralHighThreshold04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValuePU03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPU"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValuePU03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValuePU03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValuePU04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPU"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValuePU04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFMeanDRIsoValuePU04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('meanDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueCharged03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositCharged"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueCharged03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueCharged03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueCharged04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositCharged"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueCharged04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueCharged04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueChargedAll03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAll"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueChargedAll03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueChargedAll03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueChargedAll04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAll"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueChargedAll04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueChargedAll04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring('0.0001', 
            'Threshold(0.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGamma03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGamma03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGamma03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGamma04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGamma04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGamma04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGammaHighThreshold03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGammaHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGammaHighThreshold03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGammaHighThreshold04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGamma"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGammaHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueGammaHighThreshold04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutral03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutral03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutral03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutral04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutral04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutral04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutralHighThreshold03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutralHighThreshold03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutralHighThreshold03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutralHighThreshold04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutral"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutralHighThreshold04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValueNeutralHighThreshold04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(1.0)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValuePU03 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPU"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValuePU03PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValuePU03PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.3),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValuePU04 = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPU"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValuePU04PAT = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPAT"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muPFSumDRIsoValuePU04PFBRECO = cms.EDProducer("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        deltaR = cms.double(0.4),
        mode = cms.string('sumDR'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("muPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring('0.01', 
            'Threshold(0.5)'),
        weight = cms.string('1')
    ))
)


process.muonMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(True),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(0.5),
    maxDeltaR = cms.double(0.5),
    mcPdgId = cms.vint32(13),
    mcStatus = cms.vint32(1),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("muons")
)


process.newV4pfBoostedDoubleSVAK8TagInfosUpdatedJECAK8 = cms.EDProducer("BoostedDoubleSVProducer",
    R0 = cms.double(0.8),
    beta = cms.double(1.0),
    maxSVDeltaRToJet = cms.double(0.7),
    svTagInfos = cms.InputTag("newV4pfInclusiveSecondaryVertexFinderAK8TagInfosUpdatedJECAK8"),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.8),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    )
)


process.newV4pfBoostedDoubleSecondaryVertexAK8BJetTagsUpdatedJECAK8 = cms.EDProducer("JetTagProducer",
    jetTagComputer = cms.string('candidateBoostedDoubleSecondaryVertexAK8Computer'),
    tagInfos = cms.VInputTag(cms.InputTag("newV4pfBoostedDoubleSVAK8TagInfosUpdatedJECAK8"))
)


process.newV4pfImpactParameterAK8TagInfosUpdatedJECAK8 = cms.EDProducer("CandIPProducer",
    candidates = cms.InputTag("packedPFCandidates"),
    computeGhostTrack = cms.bool(False),
    computeProbabilities = cms.bool(False),
    ghostTrackPriorDeltaR = cms.double(0.03),
    jetDirectionUsingGhostTrack = cms.bool(False),
    jetDirectionUsingTracks = cms.bool(False),
    jets = cms.InputTag("updatedPatJetsUpdatedJECAK8"),
    maxDeltaR = cms.double(0.8),
    maximumChiSquared = cms.double(5.0),
    maximumLongitudinalImpactParameter = cms.double(17.0),
    maximumTransverseImpactParameter = cms.double(0.2),
    minimumNumberOfHits = cms.int32(0),
    minimumNumberOfPixelHits = cms.int32(1),
    minimumTransverseMomentum = cms.double(1.0),
    primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    useTrackQuality = cms.bool(False)
)


process.newV4pfInclusiveSecondaryVertexFinderAK8TagInfosUpdatedJECAK8 = cms.EDProducer("CandSecondaryVertexProducer",
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    constraint = cms.string('BeamSpot'),
    extSVCollection = cms.InputTag("slimmedSecondaryVertices"),
    extSVDeltaRToJet = cms.double(0.8),
    minimumTrackWeight = cms.double(0.5),
    trackIPTagInfos = cms.InputTag("newV4pfImpactParameterAK8TagInfosUpdatedJECAK8"),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.8),
        maxDecayLen = cms.double(99999.9),
        maxDistToAxis = cms.double(0.2),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(1),
        ptMin = cms.double(1.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip3dSig'),
    useExternalSV = cms.bool(True),
    usePVError = cms.bool(True),
    vertexCuts = cms.PSet(
        distSig2dMax = cms.double(99999.9),
        distSig2dMin = cms.double(2.0),
        distSig3dMax = cms.double(99999.9),
        distSig3dMin = cms.double(-99999.9),
        distVal2dMax = cms.double(2.5),
        distVal2dMin = cms.double(0.01),
        distVal3dMax = cms.double(99999.9),
        distVal3dMin = cms.double(-99999.9),
        fracPV = cms.double(0.79),
        massMax = cms.double(6.5),
        maxDeltaRToJetAxis = cms.double(0.8),
        minimumTrackWeight = cms.double(0.5),
        multiplicityMin = cms.uint32(2),
        useTrackWeights = cms.bool(True),
        v0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        )
    ),
    vertexReco = cms.PSet(
        finder = cms.string('avr'),
        minweight = cms.double(0.5),
        primcut = cms.double(1.8),
        seccut = cms.double(6.0),
        smoothing = cms.bool(False),
        weightthreshold = cms.double(0.001)
    ),
    vertexSelection = cms.PSet(
        sortCriterium = cms.string('dist3dError')
    )
)


process.particleFlowDisplacedVertex = cms.EDProducer("PFDisplacedVertexProducer",
    avfParameters = cms.PSet(
        Tini = cms.double(256.0),
        ratio = cms.double(0.25),
        sigmacut = cms.double(6.0)
    ),
    debug = cms.untracked.bool(False),
    longSize = cms.double(5),
    mainVertexLabel = cms.InputTag("offlinePrimaryVertices"),
    minAdaptWeight = cms.double(0.5),
    offlineBeamSpotLabel = cms.InputTag("offlineBeamSpot"),
    primaryVertexCut = cms.double(1.8),
    switchOff2TrackVertex = cms.untracked.bool(True),
    tecCut = cms.double(220),
    tobCut = cms.double(100),
    tracksSelectorParameters = cms.PSet(
        bSelectTracks = cms.bool(True),
        dxy_min = cms.double(0.2),
        nChi2_max = cms.double(5.0),
        nChi2_min = cms.double(0.5),
        nHits_min = cms.int32(6),
        nOuterHits_max = cms.int32(9),
        pt_min = cms.double(0.2),
        quality = cms.string('HighPurity')
    ),
    transvSize = cms.double(1.0),
    verbose = cms.untracked.bool(False),
    vertexCandidatesLabel = cms.InputTag("particleFlowDisplacedVertexCandidate"),
    vertexIdentifierParameters = cms.PSet(
        angles = cms.vdouble(15, 15),
        bIdentifyVertices = cms.bool(True),
        logPrimSec_min = cms.double(0.0),
        looper_eta_max = cms.double(0.1),
        masses = cms.vdouble(0.05, 0.485, 0.515, 0.48, 0.52, 
            1.107, 1.125, 0.2),
        pt_kink_min = cms.double(3.0),
        pt_min = cms.double(0.5)
    )
)


process.particleFlowPtrs = cms.EDProducer("PFCandidateFwdPtrProducer",
    src = cms.InputTag("particleFlow")
)


process.patCaloMet = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetTrue"),
    metSource = cms.InputTag("metrawCalo"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.29, 1.19, 1.07, 1.13, 1.12),
        pjpar = cms.vdouble(-0.04, 0.6504)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatPhotons"),
    srcPFCands = cms.InputTag("particleFlow"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patElectrons = cms.EDProducer("PATElectronProducer",
    addEfficiencies = cms.bool(False),
    addElectronID = cms.bool(True),
    addGenMatch = cms.bool(True),
    addPFClusterIso = cms.bool(False),
    addResolutions = cms.bool(False),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    efficiencies = cms.PSet(

    ),
    electronIDSources = cms.PSet(
        cutBasedElectronHLTPreselection_Summer16_V1 = cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"),
        cutBasedElectronID_Summer16_80X_V1_loose = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"),
        cutBasedElectronID_Summer16_80X_V1_medium = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"),
        cutBasedElectronID_Summer16_80X_V1_tight = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"),
        cutBasedElectronID_Summer16_80X_V1_veto = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"),
        eidLoose = cms.InputTag("eidLoose"),
        eidRobustHighEnergy = cms.InputTag("eidRobustHighEnergy"),
        eidRobustLoose = cms.InputTag("eidRobustLoose"),
        eidRobustTight = cms.InputTag("eidRobustTight"),
        eidTight = cms.InputTag("eidTight"),
        heepElectronID_HEEPV70 = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"),
        mvaEleID_Spring16_GeneralPurpose_V1_wp80 = cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
        mvaEleID_Spring16_GeneralPurpose_V1_wp90 = cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
        mvaEleID_Spring16_HZZ_V1_wpLoose = cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose")
    ),
    electronSource = cms.InputTag("gedGsfElectrons"),
    embedBasicClusters = cms.bool(True),
    embedGenMatch = cms.bool(True),
    embedGsfElectronCore = cms.bool(True),
    embedGsfTrack = cms.bool(True),
    embedHighLevelSelection = cms.bool(True),
    embedPFCandidate = cms.bool(True),
    embedPflowBasicClusters = cms.bool(True),
    embedPflowPreshowerClusters = cms.bool(True),
    embedPflowSuperCluster = cms.bool(True),
    embedPreshowerClusters = cms.bool(True),
    embedRecHits = cms.bool(True),
    embedSeedCluster = cms.bool(True),
    embedSuperCluster = cms.bool(True),
    embedTrack = cms.bool(True),
    genParticleMatch = cms.InputTag("electronMatch"),
    isoDeposits = cms.PSet(
        pfChargedAll = cms.InputTag("elPFIsoDepositChargedAllPAT"),
        pfChargedHadrons = cms.InputTag("elPFIsoDepositChargedPAT"),
        pfNeutralHadrons = cms.InputTag("elPFIsoDepositNeutralPAT"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoDepositPUPAT"),
        pfPhotons = cms.InputTag("elPFIsoDepositGammaPAT")
    ),
    isolationValues = cms.PSet(
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll04PFIdPAT"),
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged04PFIdPAT"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral04PFIdPAT"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU04PFIdPAT"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma04PFIdPAT")
    ),
    isolationValuesNoPFId = cms.PSet(
        pfChargedAll = cms.InputTag("elPFIsoValueChargedAll04NoPFIdPAT"),
        pfChargedHadrons = cms.InputTag("elPFIsoValueCharged04NoPFIdPAT"),
        pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral04NoPFIdPAT"),
        pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU04NoPFIdPAT"),
        pfPhotons = cms.InputTag("elPFIsoValueGamma04NoPFIdPAT")
    ),
    pfCandidateMap = cms.InputTag("particleFlow","electrons"),
    pfElectronSource = cms.InputTag("particleFlow"),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    resolutions = cms.PSet(

    ),
    useParticleFlow = cms.bool(False),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag(cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"), 
                cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"), 
                cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"), cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"), 
                cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"), 
                cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"), cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"), 
                cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"), 
                cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"), cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"), 
                cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"), 
                cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"), cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose"), 
                cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"), 
                cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"), 
                cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"), cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"), 
                cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"), 
                cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"), cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"), 
                cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"), 
                cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"), cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"), 
                cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"), cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"), 
                cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"), cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp80"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-GeneralPurpose-V1-wp90"), cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring16-HZZ-V1-wpLoose"))
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    ),
    userIsolation = cms.PSet(

    )
)


process.patHemispheres = cms.EDProducer("PATHemisphereProducer",
    combinationMethod = cms.int32(3),
    maxElectronEta = cms.double(5),
    maxJetEta = cms.double(5),
    maxMuonEta = cms.double(5),
    maxPhotonEta = cms.double(5),
    maxTauEta = cms.double(-1),
    minElectronEt = cms.double(7),
    minJetEt = cms.double(30),
    minMuonEt = cms.double(7),
    minPhotonEt = cms.double(200000),
    minTauEt = cms.double(1000000),
    patElectrons = cms.InputTag("cleanLayer1Electrons"),
    patJets = cms.InputTag("cleanLayer1Jets"),
    patMets = cms.InputTag("layer1METs"),
    patMuons = cms.InputTag("cleanLayer1Muons"),
    patPhotons = cms.InputTag("cleanLayer1Photons"),
    patTaus = cms.InputTag("cleanLayer1Taus"),
    seedMethod = cms.int32(3)
)


process.patJetCharge = cms.EDProducer("JetChargeProducer",
    exp = cms.double(1.0),
    src = cms.InputTag("ak4JetTracksAssociatorAtVertexPF"),
    var = cms.string('Pt')
)


process.patJetCorrFactors = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring('L1FastJet', 
        'L2Relative', 
        'L3Absolute'),
    payload = cms.string('AK4PFchs'),
    primaryVertices = cms.InputTag("offlinePrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("ak4PFJetsCHS"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.patJetCorrFactorsReapplyJEC = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring('L1FastJet', 
        'L2Relative', 
        'L3Absolute'),
    payload = cms.string('AK4PFchs'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedJets"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.patJetCorrFactorsTransientCorrectedUpdatedJECAK8 = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring('L1FastJet', 
        'L2Relative', 
        'L3Absolute'),
    payload = cms.string('AK8PFchs'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("updatedPatJetsUpdatedJECAK8"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.patJetCorrFactorsUpdatedJEC = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring('L1FastJet', 
        'L2Relative', 
        'L3Absolute'),
    payload = cms.string('AK4PFchs'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedJets"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.patJetCorrFactorsUpdatedJECAK8 = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring(),
    payload = cms.string('AK4PFchs'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedJetsAK8"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.patJetFlavourAssociation = cms.EDProducer("JetFlavourClustering",
    bHadrons = cms.InputTag("patJetPartons","bHadrons"),
    cHadrons = cms.InputTag("patJetPartons","cHadrons"),
    ghostRescaling = cms.double(1e-18),
    hadronFlavourHasPriority = cms.bool(False),
    jetAlgorithm = cms.string('AntiKt'),
    jets = cms.InputTag("ak4PFJetsCHS"),
    leptons = cms.InputTag("patJetPartons","leptons"),
    partons = cms.InputTag("patJetPartons","algorithmicPartons"),
    rParam = cms.double(0.4)
)


process.patJetFlavourAssociationLegacy = cms.EDProducer("JetFlavourIdentifier",
    physicsDefinition = cms.bool(False),
    srcByReference = cms.InputTag("patJetPartonAssociationLegacy")
)


process.patJetGenJetMatch = cms.EDProducer("GenJetMatcher",
    checkCharge = cms.bool(False),
    matched = cms.InputTag("ak4GenJets"),
    maxDeltaR = cms.double(0.4),
    mcPdgId = cms.vint32(),
    mcStatus = cms.vint32(),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("ak4PFJetsCHS")
)


process.patJetPartonAssociationLegacy = cms.EDProducer("JetPartonMatcher",
    coneSizeToAssociate = cms.double(0.3),
    jets = cms.InputTag("ak4PFJetsCHS"),
    partons = cms.InputTag("patJetPartonsLegacy")
)


process.patJetPartonMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(False),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(3.0),
    maxDeltaR = cms.double(0.4),
    mcPdgId = cms.vint32(1, 2, 3, 4, 5, 
        21),
    mcStatus = cms.vint32(3, 23),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("ak4PFJetsCHS")
)


process.patJetPartons = cms.EDProducer("HadronAndPartonSelector",
    particles = cms.InputTag("genParticles"),
    partonMode = cms.string('Auto'),
    src = cms.InputTag("generator")
)


process.patJetPartonsLegacy = cms.EDProducer("PartonSelector",
    src = cms.InputTag("genParticles"),
    withLeptons = cms.bool(False)
)


process.patJets = cms.EDProducer("PATJetProducer",
    JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociation"),
    JetPartonMapSource = cms.InputTag("patJetFlavourAssociationLegacy"),
    addAssociatedTracks = cms.bool(True),
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addEfficiencies = cms.bool(False),
    addGenJetMatch = cms.bool(True),
    addGenPartonMatch = cms.bool(True),
    addJetCharge = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addJetFlavourInfo = cms.bool(True),
    addJetID = cms.bool(False),
    addPartonJetMatch = cms.bool(False),
    addResolutions = cms.bool(False),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(cms.InputTag("pfJetBProbabilityBJetTags"), cms.InputTag("pfJetProbabilityBJetTags"), cms.InputTag("pfTrackCountingHighEffBJetTags"), cms.InputTag("pfSimpleSecondaryVertexHighEffBJetTags"), cms.InputTag("pfSimpleInclusiveSecondaryVertexHighEffBJetTags"), 
        cms.InputTag("pfCombinedSecondaryVertexV2BJetTags"), cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"), cms.InputTag("softPFMuonBJetTags"), cms.InputTag("softPFElectronBJetTags"), cms.InputTag("pfCombinedMVAV2BJetTags"), 
        cms.InputTag("pfCombinedCvsLJetTags"), cms.InputTag("pfCombinedCvsBJetTags")),
    efficiencies = cms.PSet(

    ),
    embedGenJetMatch = cms.bool(True),
    embedGenPartonMatch = cms.bool(True),
    embedPFCandidates = cms.bool(False),
    genJetMatch = cms.InputTag("patJetGenJetMatch"),
    genPartonMatch = cms.InputTag("patJetPartonMatch"),
    getJetMCFlavour = cms.bool(True),
    jetChargeSource = cms.InputTag("patJetCharge"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactors")),
    jetIDMap = cms.InputTag("ak4JetID"),
    jetSource = cms.InputTag("ak4PFJetsCHS"),
    partonJetSource = cms.InputTag("NOT_IMPLEMENTED"),
    resolutions = cms.PSet(

    ),
    tagInfoSources = cms.VInputTag(),
    trackAssociationSource = cms.InputTag("ak4JetTracksAssociatorAtVertexPF"),
    useLegacyJetMCFlavour = cms.bool(False),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patJetsReapplyJEC = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC")),
    jetSource = cms.InputTag("slimmedJets"),
    tagInfoSources = cms.VInputTag(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patMETs = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(True),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetTrue"),
    metSource = cms.InputTag("pfMetT1"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.29, 1.19, 1.07, 1.13, 1.12),
        pjpar = cms.vdouble(-0.04, 0.6504)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatPhotons"),
    srcPFCands = cms.InputTag("particleFlow"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patMuons = cms.EDProducer("PATMuonProducer",
    addEfficiencies = cms.bool(False),
    addGenMatch = cms.bool(True),
    addResolutions = cms.bool(False),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    caloMETMuonCorrs = cms.InputTag("muonMETValueMapProducer","muCorrData"),
    efficiencies = cms.PSet(

    ),
    embedCaloMETMuonCorrs = cms.bool(True),
    embedCombinedMuon = cms.bool(True),
    embedDytMuon = cms.bool(True),
    embedGenMatch = cms.bool(True),
    embedHighLevelSelection = cms.bool(True),
    embedMuonBestTrack = cms.bool(True),
    embedPFCandidate = cms.bool(True),
    embedPfEcalEnergy = cms.bool(True),
    embedPickyMuon = cms.bool(True),
    embedStandAloneMuon = cms.bool(True),
    embedTcMETMuonCorrs = cms.bool(False),
    embedTpfmsMuon = cms.bool(True),
    embedTrack = cms.bool(False),
    embedTunePMuonBestTrack = cms.bool(True),
    forceBestTrackEmbedding = cms.bool(False),
    genParticleMatch = cms.InputTag("muonMatch"),
    isoDeposits = cms.PSet(
        pfChargedAll = cms.InputTag("muPFIsoDepositChargedAllPAT"),
        pfChargedHadrons = cms.InputTag("muPFIsoDepositChargedPAT"),
        pfNeutralHadrons = cms.InputTag("muPFIsoDepositNeutralPAT"),
        pfPUChargedHadrons = cms.InputTag("muPFIsoDepositPUPAT"),
        pfPhotons = cms.InputTag("muPFIsoDepositGammaPAT")
    ),
    isolationValues = cms.PSet(
        pfChargedAll = cms.InputTag("muPFIsoValueChargedAll04PAT"),
        pfChargedHadrons = cms.InputTag("muPFIsoValueCharged04PAT"),
        pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral04PAT"),
        pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU04PAT"),
        pfPhotons = cms.InputTag("muPFIsoValueGamma04PAT")
    ),
    muonSource = cms.InputTag("muons"),
    pfMuonSource = cms.InputTag("particleFlow"),
    pvSrc = cms.InputTag("offlinePrimaryVertices"),
    resolutions = cms.PSet(

    ),
    tcMETMuonCorrs = cms.InputTag("muonTCMETValueMapProducer","muCorrData"),
    useParticleFlow = cms.bool(False),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    ),
    userIsolation = cms.PSet(

    )
)


process.patPFMet = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(True),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(True),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetExtractor"),
    metSource = cms.InputTag("pfMet"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.29, 1.19, 1.07, 1.13, 1.12),
        pjpar = cms.vdouble(-0.04, 0.6504)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag(cms.InputTag("slimmedElectrons"), cms.InputTag("slimmedMuons"), cms.InputTag("slimmedPhotons")),
    srcPFCands = cms.InputTag("packedPFCandidates"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patPFMetT0Corr = cms.EDProducer("Type0PFMETcorrInputProducer",
    correction = cms.PSet(
        formula = cms.string('(x<35)?(-( [0]+x*[1]+pow(x, 2)*[2]+pow(x, 3)*[3] )):(-( [0]+35*[1]+pow(35, 2)*[2]+pow(35, 3)*[3] ))'),
        par0 = cms.double(-0.181414),
        par1 = cms.double(-0.476934),
        par2 = cms.double(0.00863564),
        par3 = cms.double(-4.94181e-05)
    ),
    minDz = cms.double(0.2),
    srcHardScatterVertex = cms.InputTag("selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0"),
    srcPFCandidateToVertexAssociations = cms.InputTag("pfCandidateToVertexAssociation")
)


process.patPFMetT0pcT1 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1T2 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1T2Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1T2Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT0pcT1T2TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT0pcT1Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT0pcT1TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"))
)


process.patPFMetT1ElectronEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnDown"))
)


process.patPFMetT1ElectronEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnUp"))
)


process.patPFMetT1JetEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnDown"))
)


process.patPFMetT1JetEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnUp"))
)


process.patPFMetT1JetResDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetResDown"))
)


process.patPFMetT1JetResUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetResUp"))
)


process.patPFMetT1MuonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnDown"))
)


process.patPFMetT1MuonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnUp"))
)


process.patPFMetT1PhotonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnDown"))
)


process.patPFMetT1PhotonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnUp"))
)


process.patPFMetT1Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"))
)


process.patPFMetT1SmearElectronEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnDown"))
)


process.patPFMetT1SmearElectronEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnUp"))
)


process.patPFMetT1SmearJetEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnDown"))
)


process.patPFMetT1SmearJetEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnUp"))
)


process.patPFMetT1SmearJetResDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrSmearedJetResDown"))
)


process.patPFMetT1SmearJetResUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrSmearedJetResUp"))
)


process.patPFMetT1SmearMuonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnDown"))
)


process.patPFMetT1SmearMuonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnUp"))
)


process.patPFMetT1SmearPhotonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnDown"))
)


process.patPFMetT1SmearPhotonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnUp"))
)


process.patPFMetT1SmearTauEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnDown"))
)


process.patPFMetT1SmearTauEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnUp"))
)


process.patPFMetT1SmearUnclusteredEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnDown"))
)


process.patPFMetT1SmearUnclusteredEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnUp"))
)


process.patPFMetT1T2 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"))
)


process.patPFMetT1T2Corr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("cleanedPatJets"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT1T2Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"))
)


process.patPFMetT1T2SmearCorr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("selectedPatJetsForMetT1T2SmearCorr"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT1T2Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1T2TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1TauEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnDown"))
)


process.patPFMetT1TauEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnUp"))
)


process.patPFMetT1Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1UnclusteredEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnDown"))
)


process.patPFMetT1UnclusteredEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnUp"))
)


process.patPFMetT2Corr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("cleanedPatJets"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT2SmearCorr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("selectedPatJetsForMetT2SmearCorr"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetTxy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetTxyCorr = cms.EDProducer("MultShiftMETcorrInputProducer",
    parameters = cms.VPSet(cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
        cms.PSet(
            etaMax = cms.double(0),
            etaMin = cms.double(-2.7),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hEtaMinus'),
            px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
            py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
            type = cms.int32(1),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(1.392),
            etaMin = cms.double(-1.392),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0Barrel'),
            px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
            py = cms.vdouble(0.00798098092474, -0.000103998219585),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(3),
            etaMin = cms.double(1.392),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0EndcapPlus'),
            px = cms.vdouble(-0.00305719113962, -0.00032676418359),
            py = cms.vdouble(-0.00345131507897, 0.000164816815994),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-1.392),
            etaMin = cms.double(-3.0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0EndcapMinus'),
            px = cms.vdouble(-0.000159031461755, 0.00012231873804),
            py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(1.479),
            etaMin = cms.double(-1.479),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaBarrel'),
            px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
            py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(3.0),
            etaMin = cms.double(1.479),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaEndcapPlus'),
            px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
            py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-1.479),
            etaMin = cms.double(-3.0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaEndcapMinus'),
            px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
            py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(5.2),
            etaMin = cms.double(2.901376),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hHFPlus'),
            px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
            py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
            type = cms.int32(6),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-2.901376),
            etaMin = cms.double(-5.2),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hHFMinus'),
            px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
            py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
            type = cms.int32(6),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(5.2),
            etaMin = cms.double(2.901376),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('egammaHFPlus'),
            px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
            py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
            type = cms.int32(7),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-2.901376),
            etaMin = cms.double(-5.2),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('egammaHFMinus'),
            px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
            py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
            type = cms.int32(7),
            varType = cms.int32(0)
        )),
    srcPFlow = cms.InputTag("packedPFCandidates"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.patPhotons = cms.EDProducer("PATPhotonProducer",
    addEfficiencies = cms.bool(False),
    addGenMatch = cms.bool(True),
    addPFClusterIso = cms.bool(False),
    addPhotonID = cms.bool(True),
    addResolutions = cms.bool(False),
    beamLineSrc = cms.InputTag("offlineBeamSpot"),
    efficiencies = cms.PSet(

    ),
    electronSource = cms.InputTag("gedGsfElectrons"),
    embedBasicClusters = cms.bool(True),
    embedGenMatch = cms.bool(True),
    embedPreshowerClusters = cms.bool(True),
    embedRecHits = cms.bool(True),
    embedSeedCluster = cms.bool(True),
    embedSuperCluster = cms.bool(True),
    genParticleMatch = cms.InputTag("photonMatch"),
    isoDeposits = cms.PSet(
        pfChargedAll = cms.InputTag("phPFIsoDepositChargedAllPAT"),
        pfChargedHadrons = cms.InputTag("phPFIsoDepositChargedPAT"),
        pfNeutralHadrons = cms.InputTag("phPFIsoDepositNeutralPAT"),
        pfPUChargedHadrons = cms.InputTag("phPFIsoDepositPUPAT"),
        pfPhotons = cms.InputTag("phPFIsoDepositGammaPAT")
    ),
    isolationValues = cms.PSet(
        pfChargedAll = cms.InputTag("phPFIsoValueChargedAll04PFIdPAT"),
        pfChargedHadrons = cms.InputTag("phPFIsoValueCharged04PFIdPAT"),
        pfNeutralHadrons = cms.InputTag("phPFIsoValueNeutral04PFIdPAT"),
        pfPUChargedHadrons = cms.InputTag("phPFIsoValuePU04PFIdPAT"),
        pfPhotons = cms.InputTag("phPFIsoValueGamma04PFIdPAT")
    ),
    photonIDSources = cms.PSet(
        PhotonCutBasedIDLoose = cms.InputTag("PhotonIDProdGED","PhotonCutBasedIDLoose"),
        PhotonCutBasedIDTight = cms.InputTag("PhotonIDProdGED","PhotonCutBasedIDTight"),
        cutBasedPhotonID_Spring16_V2p2_loose = cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"),
        cutBasedPhotonID_Spring16_V2p2_medium = cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"),
        cutBasedPhotonID_Spring16_V2p2_tight = cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"),
        mvaPhoID_Spring16_nonTrig_V1_wp80 = cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"),
        mvaPhoID_Spring16_nonTrig_V1_wp90 = cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90")
    ),
    photonSource = cms.InputTag("gedPhotons"),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    resolutions = cms.PSet(

    ),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag(cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90"), 
                cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90"), 
                cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90"), 
                cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90"), 
                cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90"), 
                cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90"), 
                cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90"), 
                cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90"), 
                cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90"), 
                cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"), cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp80"), cms.InputTag("egmPhotonIDs","mvaPhoID-Spring16-nonTrig-V1-wp90"))
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    ),
    userIsolation = cms.PSet(

    )
)


process.patSmearedJets = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    variation = cms.int32(0)
)


process.patTaus = cms.EDProducer("PATTauProducer",
    addEfficiencies = cms.bool(False),
    addGenJetMatch = cms.bool(True),
    addGenMatch = cms.bool(True),
    addResolutions = cms.bool(False),
    addTauID = cms.bool(True),
    addTauJetCorrFactors = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    embedGenJetMatch = cms.bool(True),
    embedGenMatch = cms.bool(True),
    embedIsolationPFCands = cms.bool(False),
    embedIsolationPFChargedHadrCands = cms.bool(False),
    embedIsolationPFGammaCands = cms.bool(False),
    embedIsolationPFNeutralHadrCands = cms.bool(False),
    embedIsolationTracks = cms.bool(False),
    embedLeadPFCand = cms.bool(False),
    embedLeadPFChargedHadrCand = cms.bool(False),
    embedLeadPFNeutralCand = cms.bool(False),
    embedLeadTrack = cms.bool(False),
    embedSignalPFCands = cms.bool(False),
    embedSignalPFChargedHadrCands = cms.bool(False),
    embedSignalPFGammaCands = cms.bool(False),
    embedSignalPFNeutralHadrCands = cms.bool(False),
    embedSignalTracks = cms.bool(False),
    genJetMatch = cms.InputTag("tauGenJetMatch"),
    genParticleMatch = cms.InputTag("tauMatch"),
    isoDeposits = cms.PSet(

    ),
    resolutions = cms.PSet(

    ),
    tauIDSources = cms.PSet(
        againstElectronLooseMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6LooseElectronRejection"),
        againstElectronMVA6Raw = cms.InputTag("hpsPFTauDiscriminationByMVA6rawElectronRejection"),
        againstElectronMVA6category = cms.InputTag("hpsPFTauDiscriminationByMVA6rawElectronRejection","category"),
        againstElectronMediumMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6MediumElectronRejection"),
        againstElectronTightMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6TightElectronRejection"),
        againstElectronVLooseMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6VLooseElectronRejection"),
        againstElectronVTightMVA6 = cms.InputTag("hpsPFTauDiscriminationByMVA6VTightElectronRejection"),
        againstMuonLoose3 = cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection3"),
        againstMuonTight3 = cms.InputTag("hpsPFTauDiscriminationByTightMuonRejection3"),
        byCombinedIsolationDeltaBetaCorrRaw3Hits = cms.InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits"),
        byIsolationMVArun2v1DBdR03oldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw"),
        byIsolationMVArun2v1DBnewDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw"),
        byIsolationMVArun2v1DBoldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw"),
        byIsolationMVArun2v1PWdR03oldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw"),
        byIsolationMVArun2v1PWnewDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw"),
        byIsolationMVArun2v1PWoldDMwLTraw = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw"),
        byLooseCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits"),
        byLooseIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT"),
        byLooseIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT"),
        byLooseIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT"),
        byLooseIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT"),
        byLooseIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT"),
        byLooseIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT"),
        byMediumCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits"),
        byMediumIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT"),
        byMediumIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT"),
        byMediumIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT"),
        byMediumIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT"),
        byMediumIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT"),
        byMediumIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT"),
        byPhotonPtSumOutsideSignalCone = cms.InputTag("hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone"),
        byTightCombinedIsolationDeltaBetaCorr3Hits = cms.InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits"),
        byTightIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT"),
        byTightIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT"),
        byTightIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT"),
        byTightIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT"),
        byTightIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT"),
        byTightIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT"),
        byVLooseIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT"),
        byVLooseIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT"),
        byVLooseIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT"),
        byVLooseIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT"),
        byVLooseIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT"),
        byVLooseIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT"),
        byVTightIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT"),
        byVTightIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT"),
        byVTightIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT"),
        byVTightIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT"),
        byVTightIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT"),
        byVTightIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT"),
        byVVTightIsolationMVArun2v1DBdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT"),
        byVVTightIsolationMVArun2v1DBnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT"),
        byVVTightIsolationMVArun2v1DBoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT"),
        byVVTightIsolationMVArun2v1PWdR03oldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT"),
        byVVTightIsolationMVArun2v1PWnewDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT"),
        byVVTightIsolationMVArun2v1PWoldDMwLT = cms.InputTag("hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT"),
        chargedIsoPtSum = cms.InputTag("hpsPFTauChargedIsoPtSum"),
        chargedIsoPtSumdR03 = cms.InputTag("hpsPFTauChargedIsoPtSumdR03"),
        decayModeFinding = cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
        decayModeFindingNewDMs = cms.InputTag("hpsPFTauDiscriminationByDecayModeFindingNewDMs"),
        footprintCorrection = cms.InputTag("hpsPFTauFootprintCorrection"),
        footprintCorrectiondR03 = cms.InputTag("hpsPFTauFootprintCorrectiondR03"),
        neutralIsoPtSum = cms.InputTag("hpsPFTauNeutralIsoPtSum"),
        neutralIsoPtSumWeight = cms.InputTag("hpsPFTauNeutralIsoPtSumWeight"),
        neutralIsoPtSumWeightdR03 = cms.InputTag("hpsPFTauNeutralIsoPtSumWeightdR03"),
        neutralIsoPtSumdR03 = cms.InputTag("hpsPFTauNeutralIsoPtSumdR03"),
        photonPtSumOutsideSignalCone = cms.InputTag("hpsPFTauPhotonPtSumOutsideSignalCone"),
        photonPtSumOutsideSignalConedR03 = cms.InputTag("hpsPFTauPhotonPtSumOutsideSignalConedR03"),
        puCorrPtSum = cms.InputTag("hpsPFTauPUcorrPtSum")
    ),
    tauJetCorrFactorsSource = cms.VInputTag(cms.InputTag("patTauJetCorrFactors")),
    tauSource = cms.InputTag("hpsPFTauProducer"),
    tauTransverseImpactParameterSource = cms.InputTag("hpsPFTauTransverseImpactParameters"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    ),
    userIsolation = cms.PSet(

    )
)


process.pfCandMETcorr = cms.EDProducer("PFCandMETcorrInputProducer",
    src = cms.InputTag("pfCandsNotInJetsForMetCorr")
)


process.pfCandidateToVertexAssociation = cms.EDProducer("PFCand_AssoMap",
    AssociationType = cms.InputTag("Both"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ConversionsCollection = cms.InputTag("allConversions"),
    FinalAssociation = cms.untracked.int32(1),
    GetCleanedCollections = cms.bool(False),
    MaxNumberOfAssociations = cms.int32(1),
    NIVertexCollection = cms.InputTag("particleFlowDisplacedVertex"),
    PFCandidateCollection = cms.InputTag("particleFlow"),
    UseBeamSpotCompatibility = cms.untracked.bool(True),
    V0KshortCollection = cms.InputTag("generalV0Candidates","Kshort"),
    V0LambdaCollection = cms.InputTag("generalV0Candidates","Lambda"),
    VertexCollection = cms.InputTag("offlinePrimaryVertices"),
    doReassociation = cms.bool(True),
    ignoreMissingCollection = cms.bool(True),
    nTrackWeight = cms.double(0.001)
)


process.pfCandsForUnclusteredUnc = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEleNoMuNoTau"),
    veto = cms.InputTag("slimmedPhotons")
)


process.pfCandsNoJets = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("packedPFCandidates"),
    veto = cms.InputTag("cleanedPatJets")
)


process.pfCandsNoJetsNoEle = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJets"),
    veto = cms.InputTag("slimmedElectrons")
)


process.pfCandsNoJetsNoEleNoMu = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEle"),
    veto = cms.InputTag("slimmedMuons")
)


process.pfCandsNoJetsNoEleNoMuNoTau = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEleNoMu"),
    veto = cms.InputTag("slimmedTaus")
)


process.pfCandsNotInJetsForMetCorr = cms.EDProducer("PFCandidateFromFwdPtrProducer",
    src = cms.InputTag("pfCandsNotInJetsPtrForMetCorr")
)


process.pfCandsNotInJetsPtrForMetCorr = cms.EDProducer("TPPFJetsOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlowPtrs"),
    enable = cms.bool(True),
    name = cms.untracked.string('noJet'),
    topCollection = cms.InputTag("pfJetsPtrForMetCorr"),
    verbose = cms.untracked.bool(False)
)


process.pfJetsPtrForMetCorr = cms.EDProducer("PFJetFwdPtrProducer",
    src = cms.InputTag("ak4PFJets")
)


process.pfMETcorrType0 = cms.EDProducer("Type0PFMETcorrInputProducer",
    correction = cms.PSet(
        formula = cms.string('(x<35)?(-( [0]+x*[1]+pow(x, 2)*[2]+pow(x, 3)*[3] )):(-( [0]+35*[1]+pow(35, 2)*[2]+pow(35, 3)*[3] ))'),
        par0 = cms.double(-0.181414),
        par1 = cms.double(-0.476934),
        par2 = cms.double(0.00863564),
        par3 = cms.double(-4.94181e-05)
    ),
    minDz = cms.double(0.2),
    srcHardScatterVertex = cms.InputTag("selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0"),
    srcPFCandidateToVertexAssociations = cms.InputTag("pfCandidateToVertexAssociation")
)


process.pfMEtMultShiftCorr = cms.EDProducer("MultShiftMETcorrInputProducer",
    parameters = cms.VPSet(cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
        cms.PSet(
            etaMax = cms.double(0),
            etaMin = cms.double(-2.7),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hEtaMinus'),
            px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
            py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
            type = cms.int32(1),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(1.392),
            etaMin = cms.double(-1.392),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0Barrel'),
            px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
            py = cms.vdouble(0.00798098092474, -0.000103998219585),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(3),
            etaMin = cms.double(1.392),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0EndcapPlus'),
            px = cms.vdouble(-0.00305719113962, -0.00032676418359),
            py = cms.vdouble(-0.00345131507897, 0.000164816815994),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-1.392),
            etaMin = cms.double(-3.0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0EndcapMinus'),
            px = cms.vdouble(-0.000159031461755, 0.00012231873804),
            py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(1.479),
            etaMin = cms.double(-1.479),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaBarrel'),
            px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
            py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(3.0),
            etaMin = cms.double(1.479),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaEndcapPlus'),
            px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
            py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-1.479),
            etaMin = cms.double(-3.0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaEndcapMinus'),
            px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
            py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(5.2),
            etaMin = cms.double(2.901376),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hHFPlus'),
            px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
            py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
            type = cms.int32(6),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-2.901376),
            etaMin = cms.double(-5.2),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hHFMinus'),
            px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
            py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
            type = cms.int32(6),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(5.2),
            etaMin = cms.double(2.901376),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('egammaHFPlus'),
            px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
            py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
            type = cms.int32(7),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-2.901376),
            etaMin = cms.double(-5.2),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('egammaHFMinus'),
            px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
            py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
            type = cms.int32(7),
            varType = cms.int32(0)
        )),
    srcPFlow = cms.InputTag("particleFlow"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.pfMet = cms.EDProducer("RecoMETExtractor",
    correctionLevel = cms.string('raw'),
    metSource = cms.InputTag("slimmedMETs","","@skipCurrentProcess")
)


process.pfMetT0pc = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0PfCand"))
)


process.pfMetT0pcT1 = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0PfCand"), cms.InputTag("corrPfMetType1","type1"))
)


process.pfMetT0pcT1T2Txy = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0PfCand"), cms.InputTag("corrPfMetType1","type1"), cms.InputTag("corrPfMetType2"), cms.InputTag("corrPfMetXYMult"))
)


process.pfMetT0pcT1Txy = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0PfCand"), cms.InputTag("corrPfMetType1","type1"), cms.InputTag("corrPfMetXYMult"))
)


process.pfMetT0pcTxy = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0PfCand"), cms.InputTag("corrPfMetXYMult"))
)


process.pfMetT0rt = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0RecoTrack"))
)


process.pfMetT0rtT1 = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0RecoTrack"), cms.InputTag("corrPfMetType1","type1"))
)


process.pfMetT0rtT1T2 = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0RecoTrackForType2"), cms.InputTag("corrPfMetType1","type1"), cms.InputTag("corrPfMetType2"))
)


process.pfMetT0rtT1T2Txy = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0RecoTrackForType2"), cms.InputTag("corrPfMetType1","type1"), cms.InputTag("corrPfMetType2"), cms.InputTag("corrPfMetXYMult"))
)


process.pfMetT0rtT1Txy = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0RecoTrack"), cms.InputTag("corrPfMetType1","type1"), cms.InputTag("corrPfMetXYMult"))
)


process.pfMetT0rtT2 = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0RecoTrackForType2"), cms.InputTag("corrPfMetType2"))
)


process.pfMetT0rtTxy = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType0RecoTrack"), cms.InputTag("corrPfMetXYMult"))
)


process.pfMetT1 = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType1","type1"))
)


process.pfMetT1T2 = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType1","type1"), cms.InputTag("corrPfMetType2"))
)


process.pfMetT1T2Txy = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType1","type1"), cms.InputTag("corrPfMetType2"), cms.InputTag("corrPfMetXYMult"))
)


process.pfMetT1Txy = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetType1","type1"), cms.InputTag("corrPfMetXYMult"))
)


process.pfMetTxy = cms.EDProducer("CorrectedPFMETProducer",
    src = cms.InputTag("pfMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("corrPfMetXYMult"))
)


process.pfNoJet = cms.EDProducer("TPPFJetsOnPFCandidates",
    bottomCollection = cms.InputTag("pfNoElectronJME"),
    enable = cms.bool(True),
    name = cms.untracked.string('noJet'),
    topCollection = cms.InputTag("pfJetsPtrs"),
    verbose = cms.untracked.bool(False)
)


process.pfNoPileUpIsoPFBRECO = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlowPtrs"),
    enable = cms.bool(True),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    topCollection = cms.InputTag("pfPileUpIsoPFBRECO"),
    verbose = cms.untracked.bool(False)
)


process.pfNoPileUpPFBRECO = cms.EDProducer("TPPFCandidatesOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlowPtrs"),
    enable = cms.bool(True),
    name = cms.untracked.string('pileUpOnPFCandidates'),
    topCollection = cms.InputTag("pfPileUpPFBRECO"),
    verbose = cms.untracked.bool(False)
)


process.pfPileUpIsoPFBRECO = cms.EDProducer("PFPileUp",
    Enable = cms.bool(True),
    PFCandidates = cms.InputTag("particleFlowPtrs"),
    Vertices = cms.InputTag("offlinePrimaryVertices"),
    checkClosestZVertex = cms.bool(True),
    verbose = cms.untracked.bool(False)
)


process.pfPileUpPFBRECO = cms.EDProducer("PFPileUp",
    Enable = cms.bool(True),
    PFCandidates = cms.InputTag("particleFlowPtrs"),
    Vertices = cms.InputTag("offlinePrimaryVertices"),
    checkClosestZVertex = cms.bool(True),
    verbose = cms.untracked.bool(False)
)


process.phPFIsoDepositChargedAllPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedPhotons"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositChargedAllPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedPhotonsPFBRECO"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositChargedPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedPhotons"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositChargedPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedPhotonsPFBRECO"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositGammaPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFCandWithSuperClusterExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        MissHitSCMatch_Veto = cms.bool(False),
        SCMatch_Veto = cms.bool(True),
        inputCandView = cms.InputTag("pfAllPhotonsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedPhotons"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositGammaPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFCandWithSuperClusterExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        MissHitSCMatch_Veto = cms.bool(False),
        SCMatch_Veto = cms.bool(True),
        inputCandView = cms.InputTag("pfAllPhotonsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedPhotonsPFBRECO"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositNeutralPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedPhotons"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositNeutralPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadronsPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedPhotonsPFBRECO"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositPUPAT = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("gedPhotons"),
    trackType = cms.string('candidate')
)


process.phPFIsoDepositPUPFBRECO = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(99999.99),
        Diff_z = cms.double(99999.99),
        inputCandView = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("pfSelectedPhotonsPFBRECO"),
    trackType = cms.string('candidate')
)


process.phPFIsoValueCharged03PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositCharged"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueCharged03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueCharged03PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueCharged04PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositCharged"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueCharged04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueCharged04PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueChargedAll03PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedAll"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueChargedAll03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueChargedAll03PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueChargedAll04PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedAll"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueChargedAll04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedAllPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueChargedAll04PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositChargedAllPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueGamma03PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositGamma"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.05)'),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueGamma03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.05)'),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueGamma03PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.05)'),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueGamma04PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositGamma"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.05)'),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueGamma04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositGammaPAT"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.05)'),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueGamma04PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositGammaPFBRECO"),
        vetos = cms.vstring('EcalEndcaps:ConeVeto(0.05)'),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueNeutral03PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositNeutral"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueNeutral03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueNeutral03PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueNeutral04PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositNeutral"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueNeutral04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositNeutralPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValueNeutral04PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositNeutralPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValuePU03PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositPU"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValuePU03PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositPUPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValuePU03PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.3),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValuePU04PFId = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositPU"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValuePU04PFIdPAT = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositPUPAT"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.phPFIsoValuePU04PFIdPFBRECO = cms.EDProducer("PFCandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        PivotCoordinatesForEBEE = cms.bool(True),
        deltaR = cms.double(0.4),
        mode = cms.string('sum'),
        skipDefaultVeto = cms.bool(True),
        src = cms.InputTag("phPFIsoDepositPUPFBRECO"),
        vetos = cms.vstring(),
        weight = cms.string('1')
    ))
)


process.photonIDValueMapProducer = cms.EDProducer("PhotonIDValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
    pfCandidates = cms.InputTag("particleFlow"),
    pfCandidatesMiniAOD = cms.InputTag("packedPFCandidates"),
    src = cms.InputTag("gedPhotons"),
    srcMiniAOD = cms.InputTag("slimmedPhotons"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    verticesMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.photonMVAValueMapProducer = cms.EDProducer("PhotonMVAValueMapProducer",
    mvaConfigurations = cms.VPSet(cms.PSet(
        esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
        full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
        full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
        full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
        full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
        full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
        full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
        mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
        mvaTag = cms.string('50nsV2p1'),
        phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
        phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
        phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        useValueMaps = cms.bool(False),
        weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EB_V2.weights.xml', 
            'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_50ns_EE_V2.weights.xml')
    ), 
        cms.PSet(
            esEffSigmaRRMap = cms.InputTag("photonIDValueMapProducer","phoESEffSigmaRR"),
            full5x5E1x3Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E1x3"),
            full5x5E2x2Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x2"),
            full5x5E2x5MaxMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5E2x5Max"),
            full5x5E5x5Map = cms.InputTag("photonIDValueMapProducer","phoFull5x5E5x5"),
            full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
            full5x5SigmaIEtaIPhiMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIPhi"),
            mvaName = cms.string('PhotonMVAEstimatorRun2Spring15NonTrig'),
            mvaTag = cms.string('25nsV2p1'),
            phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
            phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
            phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
            rho = cms.InputTag("fixedGridRhoFastjetAll"),
            useValueMaps = cms.bool(False),
            weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EB_V2.weights.xml', 
                'RecoEgamma/PhotonIdentification/data/Spring15/photon_general_MVA_Spring15_25ns_EE_V2.weights.xml')
        ), 
        cms.PSet(
            effAreasConfigFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Spring16/effAreaPhotons_cone03_pfPhotons_90percentBased.txt'),
            mvaName = cms.string('PhotonMVAEstimatorRun2Spring16NonTrig'),
            mvaTag = cms.string('V1'),
            phoChargedIsolation = cms.InputTag("egmPhotonIsolation","h+-DR030-"),
            phoIsoCutoff = cms.double(2.5),
            phoIsoPtScalingCoeff = cms.vdouble(0.0053, 0.0034),
            phoPhotonIsolation = cms.InputTag("egmPhotonIsolation","gamma-DR030-"),
            phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolationWithConeVeto"),
            rho = cms.InputTag("fixedGridRhoAll"),
            weightFileNames = cms.vstring('RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EB_V3.weights.xml', 
                'RecoEgamma/PhotonIdentification/data/Spring16/photon_general_MVA_Spring16_EE_V3.weights.xml')
        )),
    src = cms.InputTag("gedPhotons"),
    srcMiniAOD = cms.InputTag("slimmedPhotons")
)


process.photonMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(True),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(1.0),
    maxDeltaR = cms.double(0.2),
    mcPdgId = cms.vint32(22),
    mcStatus = cms.vint32(1),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("gedPhotons")
)


process.photonRegressionValueMapProducer = cms.EDProducer("PhotonRegressionValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    src = cms.InputTag("gedPhotons"),
    srcMiniAOD = cms.InputTag("slimmedPhotons","","@skipCurrentProcess"),
    useFull5x5 = cms.bool(False)
)


process.puppi = cms.EDProducer("PuppiProducer",
    DeltaZCut = cms.double(0.3),
    MinPuppiWeight = cms.double(0.01),
    UseDeltaZCut = cms.bool(True),
    algos = cms.VPSet(cms.PSet(
        EtaMaxExtrap = cms.double(2.0),
        MedEtaSF = cms.vdouble(1.0),
        MinNeutralPt = cms.vdouble(0.2),
        MinNeutralPtSlope = cms.vdouble(0.015),
        RMSEtaSF = cms.vdouble(1.0),
        etaMax = cms.vdouble(2.5),
        etaMin = cms.vdouble(0.0),
        ptMin = cms.vdouble(0.0),
        puppiAlgos = cms.VPSet(cms.PSet(
            algoId = cms.int32(5),
            applyLowPUCorr = cms.bool(False),
            combOpt = cms.int32(0),
            cone = cms.double(0.4),
            rmsPtMin = cms.double(0.1),
            rmsScaleFactor = cms.double(1.0),
            useCharged = cms.bool(True)
        ))
    ), 
        cms.PSet(
            EtaMaxExtrap = cms.double(2.0),
            MedEtaSF = cms.vdouble(0.9, 0.75),
            MinNeutralPt = cms.vdouble(1.7, 2.0),
            MinNeutralPtSlope = cms.vdouble(0.08, 0.08),
            RMSEtaSF = cms.vdouble(1.2, 0.95),
            etaMax = cms.vdouble(3.0, 10.0),
            etaMin = cms.vdouble(2.5, 3.0),
            ptMin = cms.vdouble(0.0, 0.0),
            puppiAlgos = cms.VPSet(cms.PSet(
                algoId = cms.int32(5),
                applyLowPUCorr = cms.bool(False),
                combOpt = cms.int32(0),
                cone = cms.double(0.4),
                rmsPtMin = cms.double(0.5),
                rmsScaleFactor = cms.double(1.0),
                useCharged = cms.bool(False)
            ))
        )),
    applyCHS = cms.bool(True),
    candName = cms.InputTag("packedPFCandidates"),
    clonePackedCands = cms.bool(False),
    invertPuppi = cms.bool(False),
    puppiDiagnostics = cms.bool(False),
    puppiForLeptons = cms.bool(False),
    useExistingWeights = cms.bool(False),
    useExp = cms.bool(False),
    useWeightsNoLep = cms.bool(False),
    vertexName = cms.InputTag("offlineSlimmedPrimaryVertices"),
    vtxNdofCut = cms.int32(4),
    vtxZCut = cms.double(24)
)


process.shiftedPatElectronEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("slimmedElectrons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.006+0*x):(0.015+0*x))')
)


process.shiftedPatElectronEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("slimmedElectrons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.006+0*x):(0.015+0*x))')
)


process.shiftedPatJetEnDown = cms.EDProducer("ShiftedPATJetProducer",
    addResidualJES = cms.bool(True),
    jetCorrLabelUpToL3 = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector"),
    jetCorrPayloadName = cms.string('AK4PFchs'),
    jetCorrUncertaintyTag = cms.string('Uncertainty'),
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("cleanedPatJets")
)


process.shiftedPatJetEnUp = cms.EDProducer("ShiftedPATJetProducer",
    addResidualJES = cms.bool(True),
    jetCorrLabelUpToL3 = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector"),
    jetCorrPayloadName = cms.string('AK4PFchs'),
    jetCorrUncertaintyTag = cms.string('Uncertainty'),
    shiftBy = cms.double(1.0),
    src = cms.InputTag("cleanedPatJets")
)


process.shiftedPatJetResDown = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    variation = cms.int32(-101)
)


process.shiftedPatJetResUp = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    variation = cms.int32(101)
)


process.shiftedPatMETCorrElectronEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedElectrons"),
    srcShifted = cms.InputTag("shiftedPatElectronEnDown")
)


process.shiftedPatMETCorrElectronEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedElectrons"),
    srcShifted = cms.InputTag("shiftedPatElectronEnUp")
)


process.shiftedPatMETCorrJetEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetEnDown")
)


process.shiftedPatMETCorrJetEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetEnUp")
)


process.shiftedPatMETCorrJetResDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetResDown")
)


process.shiftedPatMETCorrJetResUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetResUp")
)


process.shiftedPatMETCorrMuonEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedMuons"),
    srcShifted = cms.InputTag("shiftedPatMuonEnDown")
)


process.shiftedPatMETCorrMuonEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedMuons"),
    srcShifted = cms.InputTag("shiftedPatMuonEnUp")
)


process.shiftedPatMETCorrPhotonEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedPhotons"),
    srcShifted = cms.InputTag("shiftedPatPhotonEnDown")
)


process.shiftedPatMETCorrPhotonEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedPhotons"),
    srcShifted = cms.InputTag("shiftedPatPhotonEnUp")
)


process.shiftedPatMETCorrSmearedJetResDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("selectedPatJetsForMetT1T2SmearCorr"),
    srcShifted = cms.InputTag("shiftedPatSmearedJetResDown")
)


process.shiftedPatMETCorrSmearedJetResUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("selectedPatJetsForMetT1T2SmearCorr"),
    srcShifted = cms.InputTag("shiftedPatSmearedJetResUp")
)


process.shiftedPatMETCorrTauEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedTaus"),
    srcShifted = cms.InputTag("shiftedPatTauEnDown")
)


process.shiftedPatMETCorrTauEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedTaus"),
    srcShifted = cms.InputTag("shiftedPatTauEnUp")
)


process.shiftedPatMETCorrUnclusteredEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfCandsForUnclusteredUnc"),
    srcShifted = cms.InputTag("shiftedPatUnclusteredEnDown")
)


process.shiftedPatMETCorrUnclusteredEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfCandsForUnclusteredUnc"),
    srcShifted = cms.InputTag("shiftedPatUnclusteredEnUp")
)


process.shiftedPatMuonEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("slimmedMuons"),
    uncertainty = cms.string('((x<100)?(0.002+0*y):(0.05+0*y))')
)


process.shiftedPatMuonEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("slimmedMuons"),
    uncertainty = cms.string('((x<100)?(0.002+0*y):(0.05+0*y))')
)


process.shiftedPatPhotonEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("slimmedPhotons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.01+0*x):(0.025+0*x))')
)


process.shiftedPatPhotonEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("slimmedPhotons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.01+0*x):(0.025+0*x))')
)


process.shiftedPatSmearedJetResDown = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    variation = cms.int32(-1)
)


process.shiftedPatSmearedJetResUp = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    variation = cms.int32(1)
)


process.shiftedPatTauEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("slimmedTaus"),
    uncertainty = cms.string('0.03+0*x*y')
)


process.shiftedPatTauEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("slimmedTaus"),
    uncertainty = cms.string('0.03+0*x*y')
)


process.shiftedPatUnclusteredEnDown = cms.EDProducer("ShiftedParticleProducer",
    binning = cms.VPSet(cms.PSet(
        binSelection = cms.string('charge!=0'),
        binUncertainty = cms.string('sqrt(pow(0.00009*x,2)+pow(0.0085/sqrt(sin(2*atan(exp(-y)))),2))')
    ), 
        cms.PSet(
            binSelection = cms.string('pdgId==130'),
            binUncertainty = cms.string('((abs(y)<1.3)?(min(0.25,sqrt(0.64/x+0.0025))):(min(0.30,sqrt(1.0/x+0.0016))))'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==22'),
            binUncertainty = cms.string('sqrt(0.0009/x+0.000001)+0*y'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==1 || pdgId==2'),
            binUncertainty = cms.string('sqrt(1./x+0.0025)+0*y'),
            energyDependency = cms.bool(True)
        )),
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfCandsForUnclusteredUnc")
)


process.shiftedPatUnclusteredEnUp = cms.EDProducer("ShiftedParticleProducer",
    binning = cms.VPSet(cms.PSet(
        binSelection = cms.string('charge!=0'),
        binUncertainty = cms.string('sqrt(pow(0.00009*x,2)+pow(0.0085/sqrt(sin(2*atan(exp(-y)))),2))')
    ), 
        cms.PSet(
            binSelection = cms.string('pdgId==130'),
            binUncertainty = cms.string('((abs(y)<1.3)?(min(0.25,sqrt(0.64/x+0.0025))):(min(0.30,sqrt(1.0/x+0.0016))))'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==22'),
            binUncertainty = cms.string('sqrt(0.0009/x+0.000001)+0*y'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==1 || pdgId==2'),
            binUncertainty = cms.string('sqrt(1./x+0.0025)+0*y'),
            energyDependency = cms.bool(True)
        )),
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfCandsForUnclusteredUnc")
)


process.slimmedElectrons = cms.EDProducer("ModifiedElectronProducer",
    modifierConfig = cms.PSet(
        modifications = cms.VPSet(cms.PSet(
            eOverP_ECALTRKThr = cms.double(0.025),
            ecalrechitsEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
            ecalrechitsEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
            electron_config = cms.PSet(
                regressionKey_ecalonly = cms.vstring('electron_eb_ECALonly_lowpt', 
                    'electron_eb_ECALonly', 
                    'electron_ee_ECALonly_lowpt', 
                    'electron_ee_ECALonly'),
                regressionKey_ecaltrk = cms.vstring('electron_eb_ECALTRK_lowpt', 
                    'electron_eb_ECALTRK', 
                    'electron_ee_ECALTRK_lowpt', 
                    'electron_ee_ECALTRK'),
                uncertaintyKey_ecalonly = cms.vstring('electron_eb_ECALonly_lowpt_var', 
                    'electron_eb_ECALonly_var', 
                    'electron_ee_ECALonly_lowpt_var', 
                    'electron_ee_ECALonly_var'),
                uncertaintyKey_ecaltrk = cms.vstring('electron_eb_ECALTRK_lowpt_var', 
                    'electron_eb_ECALTRK_var', 
                    'electron_ee_ECALTRK_lowpt_var', 
                    'electron_ee_ECALTRK_var')
            ),
            epDiffSig_ECALTRKThr = cms.double(15.0),
            epSig_ECALTRKThr = cms.double(10.0),
            highEnergy_ECALTRKThr = cms.double(200.0),
            lowEnergy_ECALTRKThr = cms.double(50.0),
            lowEnergy_ECALonlyThr = cms.double(300.0),
            modifierName = cms.string('EGExtraInfoModifierFromDBUser'),
            photon_config = cms.PSet(
                regressionKey_ecalonly = cms.vstring('photon_eb_ECALonly_lowpt', 
                    'photon_eb_ECALonly', 
                    'photon_ee_ECALonly_lowpt', 
                    'photon_ee_ECALonly'),
                uncertaintyKey_ecalonly = cms.vstring('photon_eb_ECALonly_lowpt_var', 
                    'photon_eb_ECALonly_var', 
                    'photon_ee_ECALonly_lowpt_var', 
                    'photon_ee_ECALonly_var')
            ),
            rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
            useLocalFile = cms.bool(False)
        ))
    ),
    src = cms.InputTag("slimmedElectrons","","@skipCurrentProcess")
)


process.slimmedMETs = cms.EDProducer("PATMETSlimmer",
    caloMET = cms.InputTag("patCaloMet"),
    rawVariation = cms.InputTag("patPFMet"),
    runningOnMiniAOD = cms.bool(True),
    src = cms.InputTag("patPFMetT1"),
    t01Variation = cms.InputTag("slimmedMETs","","@skipCurrentProcess"),
    t1SmearedVarsAndUncs = cms.InputTag("patPFMetT1Smear%s"),
    t1Uncertainties = cms.InputTag("patPFMetT1%s"),
    tXYUncForT1 = cms.InputTag("patPFMetT1Txy")
)


process.tauGenJetMatch = cms.EDProducer("GenJetMatcher",
    checkCharge = cms.bool(False),
    matched = cms.InputTag("tauGenJetsSelectorAllHadrons"),
    maxDPtRel = cms.double(3.0),
    maxDeltaR = cms.double(0.1),
    mcPdgId = cms.vint32(),
    mcStatus = cms.vint32(),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer")
)


process.tauGenJets = cms.EDProducer("TauGenJetProducer",
    GenParticles = cms.InputTag("genParticles"),
    includeNeutrinos = cms.bool(False),
    verbose = cms.untracked.bool(False)
)


process.tauIsoDepositPFCandidates = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(10000.0),
        Diff_z = cms.double(10000.0),
        candidateSource = cms.InputTag("particleFlow"),
        dRmatchPFTau = cms.double(0.1),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("hpsPFTauProducer")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer"),
    trackType = cms.string('candidate')
)


process.tauIsoDepositPFChargedHadrons = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(0.1),
        Diff_z = cms.double(0.2),
        candidateSource = cms.InputTag("pfAllChargedHadronsPFBRECO"),
        dRmatchPFTau = cms.double(0.1),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("hpsPFTauProducer")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer"),
    trackType = cms.string('candidate')
)


process.tauIsoDepositPFGammas = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(10000.0),
        Diff_z = cms.double(10000.0),
        candidateSource = cms.InputTag("pfAllPhotonsPFBRECO"),
        dRmatchPFTau = cms.double(0.1),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("hpsPFTauProducer")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer"),
    trackType = cms.string('candidate')
)


process.tauIsoDepositPFNeutralHadrons = cms.EDProducer("CandIsoDepositProducer",
    ExtractorPSet = cms.PSet(
        ComponentName = cms.string('PFTauExtractor'),
        DR_Max = cms.double(1.0),
        DR_Veto = cms.double(0.0),
        DepositLabel = cms.untracked.string(''),
        Diff_r = cms.double(10000.0),
        Diff_z = cms.double(10000.0),
        candidateSource = cms.InputTag("pfAllNeutralHadronsPFBRECO"),
        dRmatchPFTau = cms.double(0.1),
        dRvetoPFTauSignalConeConstituents = cms.double(0.01),
        tauSource = cms.InputTag("hpsPFTauProducer")
    ),
    MultipleDepositsFlag = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer"),
    trackType = cms.string('candidate')
)


process.tauMatch = cms.EDProducer("MCMatcher",
    checkCharge = cms.bool(True),
    matched = cms.InputTag("genParticles"),
    maxDPtRel = cms.double(999.9),
    maxDeltaR = cms.double(999.9),
    mcPdgId = cms.vint32(15),
    mcStatus = cms.vint32(2),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(False),
    src = cms.InputTag("hpsPFTauProducer")
)


process.updatedPatJetsTransientCorrectedUpdatedJECAK8 = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(cms.InputTag("newV4pfBoostedDoubleSecondaryVertexAK8BJetTagsUpdatedJECAK8")),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsTransientCorrectedUpdatedJECAK8")),
    jetSource = cms.InputTag("updatedPatJetsUpdatedJECAK8"),
    tagInfoSources = cms.VInputTag(cms.InputTag("newV4pfBoostedDoubleSVAK8TagInfosUpdatedJECAK8"), cms.InputTag("newV4pfImpactParameterAK8TagInfosUpdatedJECAK8"), cms.InputTag("newV4pfInclusiveSecondaryVertexFinderAK8TagInfosUpdatedJECAK8")),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.updatedPatJetsUpdatedJEC = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(False),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsUpdatedJEC")),
    jetSource = cms.InputTag("slimmedJets"),
    tagInfoSources = cms.VInputTag(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.updatedPatJetsUpdatedJECAK8 = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsUpdatedJECAK8")),
    jetSource = cms.InputTag("slimmedJetsAK8"),
    tagInfoSources = cms.VInputTag(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.BadChargedCandidateFilter = cms.EDFilter("BadChargedCandidateFilter",
    PFCandidates = cms.InputTag("packedPFCandidates"),
    debug = cms.bool(False),
    innerTrackRelErr = cms.double(1.0),
    maxDR = cms.double(1e-05),
    minMuonPt = cms.double(100.0),
    minMuonTrackRelErr = cms.double(2.0),
    minPtDiffRel = cms.double(1e-05),
    muons = cms.InputTag("slimmedMuons"),
    segmentCompatibility = cms.double(0.3),
    taggingMode = cms.bool(False)
)


process.BadPFMuonFilter = cms.EDFilter("BadPFMuonFilter",
    PFCandidates = cms.InputTag("packedPFCandidates"),
    algo = cms.int32(14),
    debug = cms.bool(False),
    innerTrackRelErr = cms.double(1.0),
    minDZ = cms.double(0.1),
    minMuPt = cms.double(100),
    minPtError = cms.double(2.0),
    muons = cms.InputTag("slimmedMuons"),
    segmentCompatibility = cms.double(0.3),
    taggingMode = cms.bool(False)
)


process.jetSelectorForMet = cms.EDFilter("PATJetSelector",
    cut = cms.string('pt>15 && abs(eta)<9.9'),
    src = cms.InputTag("basicJetsForMet")
)


process.pfAllChargedHadronsPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(211, -211, 321, -321, 999211, 
        2212, -2212),
    src = cms.InputTag("pfNoPileUpIsoPFBRECO")
)


process.pfAllChargedParticlesPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(211, -211, 321, -321, 999211, 
        2212, -2212, 11, -11, 13, 
        -13),
    src = cms.InputTag("pfNoPileUpIsoPFBRECO")
)


process.pfAllNeutralHadronsAndPhotonsPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(22, 111, 130, 310, 2112),
    src = cms.InputTag("pfNoPileUpIsoPFBRECO")
)


process.pfAllNeutralHadronsPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(111, 130, 310, 2112),
    src = cms.InputTag("pfNoPileUpIsoPFBRECO")
)


process.pfAllPhotonsPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(22),
    src = cms.InputTag("pfNoPileUpIsoPFBRECO")
)


process.pfPileUpAllChargedParticlesPFBRECO = cms.EDFilter("PFCandidateFwdPtrCollectionPdgIdFilter",
    makeClones = cms.bool(True),
    pdgId = cms.vint32(211, -211, 321, -321, 999211, 
        2212, -2212, 11, -11, 13, 
        -13),
    src = cms.InputTag("pfPileUpIsoPFBRECO")
)


process.selectedPatElectrons = cms.EDFilter("PATElectronSelector",
    cut = cms.string(''),
    src = cms.InputTag("patElectrons")
)


process.selectedPatJets = cms.EDFilter("PATJetSelector",
    cut = cms.string(''),
    src = cms.InputTag("patJets")
)


process.selectedPatJetsForMetT1T2Corr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) < 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patJets")
)


process.selectedPatJetsForMetT1T2SmearCorr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) < 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patSmearedJets")
)


process.selectedPatJetsForMetT2Corr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) > 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patJets")
)


process.selectedPatJetsForMetT2SmearCorr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) > 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patSmearedJets")
)


process.selectedPatMuons = cms.EDFilter("PATMuonSelector",
    cut = cms.string(''),
    src = cms.InputTag("patMuons")
)


process.selectedPatPhotons = cms.EDFilter("PATPhotonSelector",
    cut = cms.string(''),
    src = cms.InputTag("patPhotons")
)


process.selectedPatTaus = cms.EDFilter("PATTauSelector",
    cut = cms.string(''),
    src = cms.InputTag("patTaus")
)


process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0 = cms.EDFilter("PATSingleVertexSelector",
    filter = cms.bool(False),
    mode = cms.string('firstVertex'),
    vertices = cms.InputTag("selectedVerticesForPFMEtCorrType0")
)


process.selectedUpdatedPatJetsUpdatedJEC = cms.EDFilter("PATJetSelector",
    cut = cms.string(''),
    src = cms.InputTag("updatedPatJetsUpdatedJEC")
)


process.selectedUpdatedPatJetsUpdatedJECAK8 = cms.EDFilter("PATJetSelector",
    cut = cms.string(''),
    src = cms.InputTag("updatedPatJetsTransientCorrectedUpdatedJECAK8")
)


process.selectedVerticesForPFMEtCorrType0 = cms.EDFilter("VertexSelector",
    cut = cms.string('isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2.'),
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices")
)


process.tauGenJetsSelectorAllHadrons = cms.EDFilter("TauGenJetDecayModeSelector",
    filter = cms.bool(False),
    select = cms.vstring('oneProng0Pi0', 
        'oneProng1Pi0', 
        'oneProng2Pi0', 
        'oneProngOther', 
        'threeProng0Pi0', 
        'threeProng1Pi0', 
        'threeProngOther', 
        'rare'),
    src = cms.InputTag("tauGenJets")
)


process.ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
    BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter"),
    BadPFMuonFilter = cms.InputTag("BadPFMuonFilter"),
    LHEEventLabel = cms.InputTag("externalLHEProducer"),
    PFAllCandidates = cms.InputTag("particleFlow"),
    TrackLabel = cms.InputTag("generalTracks"),
    VtxBSLabel = cms.InputTag("offlinePrimaryVerticesWithBS"),
    VtxLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
    addFilterInfoAOD = cms.bool(False),
    addFilterInfoMINIAOD = cms.bool(True),
    ak4JetSrc = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
    ak8JetSrc = cms.InputTag("selectedUpdatedPatJetsUpdatedJECAK8"),
    calibelectronSrc = cms.InputTag("calibratedPatElectrons"),
    calibphotonSrc = cms.InputTag("calibratedPatPhotons"),
    development = cms.bool(False),
    doGenParticles = cms.bool(True),
    doNoHFMET = cms.bool(False),
    dumpJets = cms.bool(True),
    dumpPDFSystWeight = cms.bool(False),
    dumpPhotons = cms.bool(True),
    dumpSoftDrop = cms.bool(False),
    dumpSubJets = cms.bool(False),
    dumpTaus = cms.bool(True),
    ebReducedRecHitCollection = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEgamma","reducedEERecHits"),
    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"),
    eleHLTIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronHLTPreselection-Summer16-V1"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-loose"),
    eleMVAHZZValuesMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16HZZV1Values"),
    eleMVAValuesMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-medium"),
    elePFClusEcalIsoProducer = cms.InputTag("electronEcalPFClusterIsolationProducer"),
    elePFClusHcalIsoProducer = cms.InputTag("electronHcalPFClusterIsolationProducer"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-tight"),
    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Summer16-80X-V1-veto"),
    electronSrc = cms.InputTag("slimmedElectrons"),
    esReducedRecHitCollection = cms.InputTag("reducedEgamma","reducedESRecHits"),
    genParticleSrc = cms.InputTag("prunedGenParticles"),
    generatorLabel = cms.InputTag("generator"),
    gsfElectronLabel = cms.InputTag("gsfElectrons"),
    gsfTrackSrc = cms.InputTag("electronGsfTracks"),
    hggPhotonIDConfiguration = cms.PSet(
    **dict(
        [
            ("cutsublead_drtotk_25_990" , cms.vdouble(0.26, 0.029, 0.0062, 0.0055) ),
            ("cutsublead_drtotk_25_991" , cms.vdouble(0.98, 0.029, 0.0109, 0.0111) ),
            ("cutsublead_drtotk_25_9910" , cms.vdouble(1, 1, 1, 1) ),
            ("cutsublead_drtotk_25_9911" , cms.vdouble(1, 1, 1, 1) ),
            ("cutsublead_drtotk_25_992" , cms.vdouble(1, 0.029, 0.021, 0.028) ),
            ("cutsublead_drtotk_25_993" , cms.vdouble(1, 0.062, 0.97, 0.97) ),
            ("cutsublead_drtotk_25_994" , cms.vdouble(1, 0.97, 1, 1) ),
            ("cutsublead_drtotk_25_995" , cms.vdouble(1, 1, 1, 1) ),
            ("cutsublead_drtotk_25_996" , cms.vdouble(1, 1, 1, 1) ),
            ("cutsublead_drtotk_25_996c0" , cms.vdouble(0.32, 0.99, 0.0097, 0.0064, 0.0048, 
            0.0061) ),
            ("cutsublead_drtotk_25_996c1" , cms.vdouble(0.98, 1, 0.012, 0.0111, 0.0086, 
            0.0122) ),
            ("cutsublead_drtotk_25_996c10" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsublead_drtotk_25_996c11" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsublead_drtotk_25_996c2" , cms.vdouble(1, 1, 0.0146, 0.021, 0.36, 
            0.026) ),
            ("cutsublead_drtotk_25_996c3" , cms.vdouble(1, 1, 0.022, 0.97, 0.98, 
            0.97) ),
            ("cutsublead_drtotk_25_996c4" , cms.vdouble(1, 1, 0.091, 1, 1, 
            1) ),
            ("cutsublead_drtotk_25_996c5" , cms.vdouble(1, 1, 0.97, 1, 1, 
            1) ),
            ("cutsublead_drtotk_25_996c6" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsublead_drtotk_25_996c7" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsublead_drtotk_25_996c8" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsublead_drtotk_25_996c9" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsublead_drtotk_25_997" , cms.vdouble(1, 1, 1, 1) ),
            ("cutsublead_drtotk_25_998" , cms.vdouble(1, 1, 1, 1) ),
            ("cutsublead_drtotk_25_999" , cms.vdouble(1, 1, 1, 1) ),
            ("cutsubleadhovere0" , cms.vdouble(0.09, 0.089, 0.101, 0.073) ),
            ("cutsubleadhovere1" , cms.vdouble(0.089, 0.079, 0.09, 0.061) ),
            ("cutsubleadhovere10" , cms.vdouble(0.0002, 0, 1.6e-05, 0) ),
            ("cutsubleadhovere11" , cms.vdouble(99.0, 99.0, 99.0, 99.0) ),
            ("cutsubleadhovere2" , cms.vdouble(0.087, 0.065, 0.087, 0.05) ),
            ("cutsubleadhovere3" , cms.vdouble(0.082, 0.062, 0.065, 0.048) ),
            ("cutsubleadhovere4" , cms.vdouble(0.076, 0.03, 0.047, 0.046) ),
            ("cutsubleadhovere5" , cms.vdouble(0.048, 0.0189, 0.032, 0.0085) ),
            ("cutsubleadhovere6" , cms.vdouble(0.042, 0.0173, 0.023, 0.0085) ),
            ("cutsubleadhovere6c0" , cms.vdouble(0.091, 0.079, 0.084, 0.103, 0.06, 
            0.099) ),
            ("cutsubleadhovere6c1" , cms.vdouble(0.089, 0.076, 0.08, 0.091, 0.051, 
            0.063) ),
            ("cutsubleadhovere6c10" , cms.vdouble(0.026, 0.0082, 3.2e-05, 3.2e-05, 0.0073, 
            0.0083) ),
            ("cutsubleadhovere6c11" , cms.vdouble(99.0, 99.0, 99.0, 99.0, 99.0, 
            99.0) ),
            ("cutsubleadhovere6c2" , cms.vdouble(0.083, 0.07, 0.061, 0.088, 0.05, 
            0.061) ),
            ("cutsubleadhovere6c3" , cms.vdouble(0.082, 0.06, 0.049, 0.077, 0.05, 
            0.061) ),
            ("cutsubleadhovere6c4" , cms.vdouble(0.068, 0.06, 0.027, 0.055, 0.047, 
            0.047) ),
            ("cutsubleadhovere6c5" , cms.vdouble(0.058, 0.023, 0.0174, 0.055, 0.047, 
            0.047) ),
            ("cutsubleadhovere6c6" , cms.vdouble(0.038, 0.022, 0.014, 0.043, 0.0116, 
            0.047) ),
            ("cutsubleadhovere6c7" , cms.vdouble(0.032, 0.0154, 0.0111, 0.023, 0.0107, 
            0.0138) ),
            ("cutsubleadhovere6c8" , cms.vdouble(0.026, 0.0137, 0.0111, 0.0192, 0.0073, 
            0.0083) ),
            ("cutsubleadhovere6c9" , cms.vdouble(0.026, 0.0137, 0.00114, 0.00068, 0.0073, 
            0.0083) ),
            ("cutsubleadhovere7" , cms.vdouble(0.037, 0.00049, 0.0198, 0.00024) ),
            ("cutsubleadhovere8" , cms.vdouble(0.028, 1.4e-05, 0.0198, 7e-06) ),
            ("cutsubleadhovere9" , cms.vdouble(0.0071, 0, 0.00055, 0) ),
            ("cutsubleadisosumoet0" , cms.vdouble(8.2, 4.1, 5.4, 2.6) ),
            ("cutsubleadisosumoet1" , cms.vdouble(6.4, 3.2, 3.4, 2.2) ),
            ("cutsubleadisosumoet10" , cms.vdouble(0.146, 0.058, -0.04, -1.16) ),
            ("cutsubleadisosumoet11" , cms.vdouble(7.0, 7.0, 4.0, 4.0) ),
            ("cutsubleadisosumoet2" , cms.vdouble(4.7, 2.8, 2.5, 1.46) ),
            ("cutsubleadisosumoet3" , cms.vdouble(3.8, 2.2, 1.77, 1.29) ),
            ("cutsubleadisosumoet4" , cms.vdouble(3.2, 1.76, 1.39, 1.18) ),
            ("cutsubleadisosumoet5" , cms.vdouble(2.6, 1.31, 1.33, 0.82) ),
            ("cutsubleadisosumoet6" , cms.vdouble(1.85, 0.96, 1.21, -0.029) ),
            ("cutsubleadisosumoet6c0" , cms.vdouble(14.1, 11.7, 9.8, 11, 9.2, 
            8.9) ),
            ("cutsubleadisosumoet6c1" , cms.vdouble(12.4, 10.2, 9.2, 10, 8.3, 
            8.7) ),
            ("cutsubleadisosumoet6c10" , cms.vdouble(5.2, 5.6, 4.8, 4.4, 4.6, 
            4.5) ),
            ("cutsubleadisosumoet6c11" , cms.vdouble(7.0, 7.0, 7.0, 4.0, 4.0, 
            4.0) ),
            ("cutsubleadisosumoet6c2" , cms.vdouble(11.2, 9.3, 9, 9.2, 7.9, 
            8.2) ),
            ("cutsubleadisosumoet6c3" , cms.vdouble(10, 8.8, 8.7, 8.5, 7.9, 
            7.5) ),
            ("cutsubleadisosumoet6c4" , cms.vdouble(9.3, 8.1, 7.5, 7.9, 6.3, 
            6.6) ),
            ("cutsubleadisosumoet6c5" , cms.vdouble(8.6, 6.6, 6.9, 7.1, 6, 
            6.3) ),
            ("cutsubleadisosumoet6c6" , cms.vdouble(8, 6.6, 6.5, 5.3, 4.9, 
            5.7) ),
            ("cutsubleadisosumoet6c7" , cms.vdouble(6.8, 6.1, 5, 4.8, 4.6, 
            4.9) ),
            ("cutsubleadisosumoet6c8" , cms.vdouble(6.5, 5.7, 4.8, 4.4, 4.6, 
            4.5) ),
            ("cutsubleadisosumoet6c9" , cms.vdouble(6.1, 5.6, 4.8, 4.4, 4.6, 
            4.5) ),
            ("cutsubleadisosumoet7" , cms.vdouble(1.31, 0.3, 1.15, -0.029) ),
            ("cutsubleadisosumoet8" , cms.vdouble(0.94, 0.112, 0.039, -0.029) ),
            ("cutsubleadisosumoet9" , cms.vdouble(0.3, 0.072, -0.04, -0.97) ),
            ("cutsubleadisosumoetbad0" , cms.vdouble(67, 69, 85, 7.2) ),
            ("cutsubleadisosumoetbad1" , cms.vdouble(64, 10.8, 13, 3.5) ),
            ("cutsubleadisosumoetbad10" , cms.vdouble(0.42, 0.41, -0.84, -1.77) ),
            ("cutsubleadisosumoetbad11" , cms.vdouble(100.0, 100.0, 100.0, 100.0) ),
            ("cutsubleadisosumoetbad2" , cms.vdouble(62, 5.2, 7.3, 2.5) ),
            ("cutsubleadisosumoetbad3" , cms.vdouble(11.7, 3.4, 3.9, 1.84) ),
            ("cutsubleadisosumoetbad4" , cms.vdouble(6.1, 2.7, 2.8, 0.66) ),
            ("cutsubleadisosumoetbad5" , cms.vdouble(5.1, 1.62, 1.38, -0.22) ),
            ("cutsubleadisosumoetbad6" , cms.vdouble(3.7, 0.97, 1.38, -0.88) ),
            ("cutsubleadisosumoetbad6c0" , cms.vdouble(73, 59, 85, 94, 74, 
            14.4) ),
            ("cutsubleadisosumoetbad6c1" , cms.vdouble(71, 50, 22, 24, 18.8, 
            12.7) ),
            ("cutsubleadisosumoetbad6c10" , cms.vdouble(7.7, 7.1, 6.1, 5.3, 5.2, 
            6.5) ),
            ("cutsubleadisosumoetbad6c11" , cms.vdouble(100.0, 100.0, 100.0, 100.0, 100.0, 
            100.0) ),
            ("cutsubleadisosumoetbad6c2" , cms.vdouble(71, 16.8, 13.7, 15.6, 13.2, 
            10.4) ),
            ("cutsubleadisosumoetbad6c3" , cms.vdouble(21, 14.4, 11.5, 12.9, 10.2, 
            9.7) ),
            ("cutsubleadisosumoetbad6c4" , cms.vdouble(15.4, 10.8, 10.2, 11.5, 8.5, 
            8.8) ),
            ("cutsubleadisosumoetbad6c5" , cms.vdouble(13.4, 10.4, 9.2, 9.3, 7, 
            7.3) ),
            ("cutsubleadisosumoetbad6c6" , cms.vdouble(11.4, 8.9, 8.3, 8.1, 7, 
            6.7) ),
            ("cutsubleadisosumoetbad6c7" , cms.vdouble(9.2, 8.1, 7.7, 7.4, 7, 
            6.7) ),
            ("cutsubleadisosumoetbad6c8" , cms.vdouble(8.3, 7.7, 6.4, 6.8, 5.3, 
            6.5) ),
            ("cutsubleadisosumoetbad6c9" , cms.vdouble(7.8, 7.4, 6.2, 5.4, 5.3, 
            6.5) ),
            ("cutsubleadisosumoetbad7" , cms.vdouble(1.72, 0.69, 1.14, -0.88) ),
            ("cutsubleadisosumoetbad8" , cms.vdouble(0.86, 0.45, -0.37, -0.88) ),
            ("cutsubleadisosumoetbad9" , cms.vdouble(0.59, 0.42, -0.84, -1.77) ),
            ("cutsubleadpf_drtotk_25_990" , cms.vdouble(0.99, 0.96, 0.0064, 0.0066, 0.0041, 
            0.0032) ),
            ("cutsubleadpf_drtotk_25_991" , cms.vdouble(1, 1, 0.0114, 0.0081, 0.0068, 
            0.0109) ),
            ("cutsubleadpf_drtotk_25_9910" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsubleadpf_drtotk_25_9911" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsubleadpf_drtotk_25_992" , cms.vdouble(1, 1, 0.0125, 0.0094, 0.02, 
            0.0145) ),
            ("cutsubleadpf_drtotk_25_993" , cms.vdouble(1, 1, 0.06, 0.0163, 0.057, 
            0.0168) ),
            ("cutsubleadpf_drtotk_25_994" , cms.vdouble(1, 1, 0.44, 0.36, 0.39, 
            0.0182) ),
            ("cutsubleadpf_drtotk_25_995" , cms.vdouble(1, 1, 0.99, 0.99, 0.96, 
            0.019) ),
            ("cutsubleadpf_drtotk_25_996" , cms.vdouble(1, 1, 1, 1, 1, 
            0.53) ),
            ("cutsubleadpf_drtotk_25_997" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsubleadpf_drtotk_25_998" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsubleadpf_drtotk_25_999" , cms.vdouble(1, 1, 1, 1, 1, 
            1) ),
            ("cutsubleadpfchisooet0" , cms.vdouble(6.3, 4.8, 4.4, 4.4, 7.1, 
            4.2) ),
            ("cutsubleadpfchisooet1" , cms.vdouble(5.6, 3.7, 3.5, 3.7, 4.5, 
            3.4) ),
            ("cutsubleadpfchisooet10" , cms.vdouble(0.96, 0.149, 0.0135, 7e-05, 0.159, 
            0) ),
            ("cutsubleadpfchisooet11" , cms.vdouble(6.0, 6.0, 6.0, 5.0, 5.0, 
            5.0) ),
            ("cutsubleadpfchisooet2" , cms.vdouble(4.7, 3.7, 2.5, 3.4, 4.5, 
            2.3) ),
            ("cutsubleadpfchisooet3" , cms.vdouble(4, 3.3, 2.4, 3.2, 4.3, 
            2) ),
            ("cutsubleadpfchisooet4" , cms.vdouble(3.4, 2.3, 2.4, 2.4, 1.9, 
            1.35) ),
            ("cutsubleadpfchisooet5" , cms.vdouble(3.2, 2.1, 2.1, 1.97, 1.33, 
            0.028) ),
            ("cutsubleadpfchisooet6" , cms.vdouble(2.7, 2.1, 1.68, 1.25, 1.1, 
            0.028) ),
            ("cutsubleadpfchisooet7" , cms.vdouble(1.32, 1.3, 1.45, 0.76, 0.97, 
            0.00047) ),
            ("cutsubleadpfchisooet8" , cms.vdouble(1.07, 1, 1.34, 0.69, 0.69, 
            5e-06) ),
            ("cutsubleadpfchisooet9" , cms.vdouble(0.98, 0.68, 1.34, 0.007, 0.52, 
            0) ),
            ("cutsubleadpfhovere0" , cms.vdouble(0.129, 0.143, 0.112, 0.146, 0.136, 
            0.147) ),
            ("cutsubleadpfhovere1" , cms.vdouble(0.129, 0.142, 0.11, 0.146, 0.136, 
            0.075) ),
            ("cutsubleadpfhovere10" , cms.vdouble(0.02, 0.0164, 0.021, 0.029, 0.018, 
            0.0178) ),
            ("cutsubleadpfhovere11" , cms.vdouble(99.0, 99.0, 99.0, 99.0, 99.0, 
            99.0) ),
            ("cutsubleadpfhovere2" , cms.vdouble(0.123, 0.142, 0.11, 0.146, 0.134, 
            0.061) ),
            ("cutsubleadpfhovere3" , cms.vdouble(0.123, 0.142, 0.11, 0.146, 0.108, 
            0.059) ),
            ("cutsubleadpfhovere4" , cms.vdouble(0.123, 0.142, 0.11, 0.146, 0.09, 
            0.059) ),
            ("cutsubleadpfhovere5" , cms.vdouble(0.123, 0.142, 0.11, 0.108, 0.068, 
            0.059) ),
            ("cutsubleadpfhovere6" , cms.vdouble(0.123, 0.122, 0.11, 0.108, 0.026, 
            0.059) ),
            ("cutsubleadpfhovere7" , cms.vdouble(0.09, 0.091, 0.11, 0.067, 0.0198, 
            0.053) ),
            ("cutsubleadpfhovere8" , cms.vdouble(0.086, 0.091, 0.059, 0.057, 0.0198, 
            0.036) ),
            ("cutsubleadpfhovere9" , cms.vdouble(0.039, 0.04, 0.059, 0.032, 0.0198, 
            0.0178) ),
            ("cutsubleadpfisosumoet0" , cms.vdouble(9.2, 7.4, 6, 9.3, 10, 
            7) ),
            ("cutsubleadpfisosumoet1" , cms.vdouble(8, 6.5, 5.1, 6.2, 6.7, 
            4.8) ),
            ("cutsubleadpfisosumoet10" , cms.vdouble(2.2, 2.1, 1.37, 2.1, 2.2, 
            1.7) ),
            ("cutsubleadpfisosumoet11" , cms.vdouble(8.0, 8.0, 8.0, 7.0, 7.0, 
            7.0) ),
            ("cutsubleadpfisosumoet2" , cms.vdouble(7.2, 6.2, 4.7, 5.7, 6.7, 
            4.3) ),
            ("cutsubleadpfisosumoet3" , cms.vdouble(6.1, 5.5, 4.7, 5.6, 6.3, 
            3.8) ),
            ("cutsubleadpfisosumoet4" , cms.vdouble(5.6, 4.8, 4.3, 5.6, 3.9, 
            3) ),
            ("cutsubleadpfisosumoet5" , cms.vdouble(5.2, 4.6, 3.6, 4, 3.2, 
            2.5) ),
            ("cutsubleadpfisosumoet6" , cms.vdouble(4.8, 4.3, 2.7, 2.7, 2.5, 
            2.4) ),
            ("cutsubleadpfisosumoet7" , cms.vdouble(2.8, 2.8, 2.4, 2.5, 2.4, 
            2.3) ),
            ("cutsubleadpfisosumoet8" , cms.vdouble(2.5, 2.4, 2.2, 2.3, 2.2, 
            1.98) ),
            ("cutsubleadpfisosumoet9" , cms.vdouble(2.4, 2.2, 1.88, 2.2, 2.2, 
            1.7) ),
            ("cutsubleadpfisosumoetbad0" , cms.vdouble(47, 55, 11.3, 19.3, 58, 
            7.7) ),
            ("cutsubleadpfisosumoetbad1" , cms.vdouble(44, 11.4, 8.2, 10.1, 9.5, 
            6.5) ),
            ("cutsubleadpfisosumoetbad10" , cms.vdouble(4.9, 3.6, 1.17, 2.1, 3.1, 
            1.55) ),
            ("cutsubleadpfisosumoetbad11" , cms.vdouble(100.0, 100.0, 100.0, 100.0, 100.0, 
            100.0) ),
            ("cutsubleadpfisosumoetbad2" , cms.vdouble(11.8, 8.2, 6.8, 8.4, 6.5, 
            5.6) ),
            ("cutsubleadpfisosumoetbad3" , cms.vdouble(10.9, 7.5, 6.2, 6.9, 5.6, 
            5.1) ),
            ("cutsubleadpfisosumoetbad4" , cms.vdouble(9.8, 6.7, 5.6, 5.3, 4.6, 
            4) ),
            ("cutsubleadpfisosumoetbad5" , cms.vdouble(9, 6.4, 4.6, 4.9, 4.2, 
            3.8) ),
            ("cutsubleadpfisosumoetbad6" , cms.vdouble(6.2, 6, 4.2, 4, 4.2, 
            3.4) ),
            ("cutsubleadpfisosumoetbad7" , cms.vdouble(6.2, 5.4, 3.5, 4, 4.2, 
            2.6) ),
            ("cutsubleadpfisosumoetbad8" , cms.vdouble(5.7, 4.7, 2.7, 4, 4.2, 
            2.1) ),
            ("cutsubleadpfisosumoetbad9" , cms.vdouble(5.6, 3.6, 2.1, 2.4, 4.2, 
            1.55) ),
            ("cutsubleadpfr90" , cms.vdouble(0.94, 0.9, 0.23, 0.94, 0.9, 
            0.23) ),
            ("cutsubleadpfr91" , cms.vdouble(0.94, 0.9, 0.23, 0.94, 0.9, 
            0.23) ),
            ("cutsubleadpfr910" , cms.vdouble(0.94, 0.9, 0.46, 0.95, 0.9, 
            0.76) ),
            ("cutsubleadpfr911" , cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0) ),
            ("cutsubleadpfr92" , cms.vdouble(0.94, 0.9, 0.24, 0.94, 0.9, 
            0.23) ),
            ("cutsubleadpfr93" , cms.vdouble(0.94, 0.9, 0.24, 0.94, 0.9, 
            0.23) ),
            ("cutsubleadpfr94" , cms.vdouble(0.94, 0.9, 0.38, 0.94, 0.9, 
            0.23) ),
            ("cutsubleadpfr95" , cms.vdouble(0.94, 0.9, 0.38, 0.94, 0.9, 
            0.23) ),
            ("cutsubleadpfr96" , cms.vdouble(0.94, 0.9, 0.38, 0.94, 0.9, 
            0.23) ),
            ("cutsubleadpfr97" , cms.vdouble(0.94, 0.9, 0.38, 0.94, 0.9, 
            0.32) ),
            ("cutsubleadpfr98" , cms.vdouble(0.94, 0.9, 0.42, 0.95, 0.9, 
            0.53) ),
            ("cutsubleadpfr99" , cms.vdouble(0.94, 0.9, 0.42, 0.95, 0.9, 
            0.67) ),
            ("cutsubleadpfsieie0" , cms.vdouble(0.0116, 0.0114, 0.0103, 0.032, 0.031, 
            0.03) ),
            ("cutsubleadpfsieie1" , cms.vdouble(0.0116, 0.0114, 0.0101, 0.029, 0.03, 
            0.029) ),
            ("cutsubleadpfsieie10" , cms.vdouble(0.0107, 0.0101, 0.009, 0.024, 0.023, 
            0.023) ),
            ("cutsubleadpfsieie11" , cms.vdouble(0.012, 0.012, 0.012, 0.03, 0.03, 
            0.03) ),
            ("cutsubleadpfsieie2" , cms.vdouble(0.0113, 0.0111, 0.0101, 0.028, 0.03, 
            0.028) ),
            ("cutsubleadpfsieie3" , cms.vdouble(0.0112, 0.0108, 0.01, 0.028, 0.029, 
            0.027) ),
            ("cutsubleadpfsieie4" , cms.vdouble(0.0108, 0.0107, 0.0099, 0.028, 0.028, 
            0.027) ),
            ("cutsubleadpfsieie5" , cms.vdouble(0.0107, 0.0107, 0.0098, 0.028, 0.027, 
            0.027) ),
            ("cutsubleadpfsieie6" , cms.vdouble(0.0107, 0.0105, 0.0098, 0.027, 0.027, 
            0.027) ),
            ("cutsubleadpfsieie7" , cms.vdouble(0.0107, 0.0105, 0.0098, 0.027, 0.027, 
            0.027) ),
            ("cutsubleadpfsieie8" , cms.vdouble(0.0107, 0.0104, 0.0097, 0.027, 0.024, 
            0.025) ),
            ("cutsubleadpfsieie9" , cms.vdouble(0.0107, 0.0103, 0.0097, 0.025, 0.023, 
            0.024) ),
            ("cutsubleadr90" , cms.vdouble(0.94, 0.31, 0.92, 0.29) ),
            ("cutsubleadr91" , cms.vdouble(0.94, 0.32, 0.94, 0.29) ),
            ("cutsubleadr910" , cms.vdouble(0.95, 0.69, 0.97, 0.85) ),
            ("cutsubleadr911" , cms.vdouble(0.0, 0.0, 0.0, 0.0) ),
            ("cutsubleadr92" , cms.vdouble(0.94, 0.34, 0.94, 0.29) ),
            ("cutsubleadr93" , cms.vdouble(0.94, 0.36, 0.94, 0.32) ),
            ("cutsubleadr94" , cms.vdouble(0.94, 0.41, 0.94, 0.34) ),
            ("cutsubleadr95" , cms.vdouble(0.94, 0.47, 0.94, 0.52) ),
            ("cutsubleadr96" , cms.vdouble(0.94, 0.69, 0.97, 0.52) ),
            ("cutsubleadr96c0" , cms.vdouble(0.94, 0.9, 0.26, 0.94, 0.9, 
            0.28) ),
            ("cutsubleadr96c1" , cms.vdouble(0.94, 0.9, 0.27, 0.94, 0.9, 
            0.28) ),
            ("cutsubleadr96c10" , cms.vdouble(0.94, 0.9, 0.52, 0.95, 0.92, 
            0.76) ),
            ("cutsubleadr96c11" , cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0) ),
            ("cutsubleadr96c2" , cms.vdouble(0.94, 0.9, 0.32, 0.94, 0.9, 
            0.28) ),
            ("cutsubleadr96c3" , cms.vdouble(0.94, 0.9, 0.32, 0.94, 0.9, 
            0.28) ),
            ("cutsubleadr96c4" , cms.vdouble(0.94, 0.9, 0.35, 0.94, 0.9, 
            0.3) ),
            ("cutsubleadr96c5" , cms.vdouble(0.94, 0.9, 0.37, 0.94, 0.9, 
            0.58) ),
            ("cutsubleadr96c6" , cms.vdouble(0.94, 0.9, 0.37, 0.95, 0.9, 
            0.58) ),
            ("cutsubleadr96c7" , cms.vdouble(0.94, 0.9, 0.39, 0.95, 0.9, 
            0.69) ),
            ("cutsubleadr96c8" , cms.vdouble(0.94, 0.9, 0.52, 0.95, 0.92, 
            0.76) ),
        ] +
        [
            ("cutsubleadr96c9" , cms.vdouble(0.94, 0.9, 0.52, 0.95, 0.92, 
            0.76) ),
            ("cutsubleadr97" , cms.vdouble(0.94, 0.69, 0.97, 0.73) ),
            ("cutsubleadr98" , cms.vdouble(0.94, 0.69, 0.97, 0.73) ),
            ("cutsubleadr99" , cms.vdouble(0.95, 0.69, 0.97, 0.84) ),
            ("cutsubleadsieie0" , cms.vdouble(0.0112, 0.0102, 0.029, 0.028) ),
            ("cutsubleadsieie1" , cms.vdouble(0.0109, 0.01, 0.029, 0.028) ),
            ("cutsubleadsieie10" , cms.vdouble(0.0094, 0.009, 0.024, 0.023) ),
            ("cutsubleadsieie11" , cms.vdouble(0.012, 0.012, 0.03, 0.03) ),
            ("cutsubleadsieie2" , cms.vdouble(0.0107, 0.0099, 0.028, 0.027) ),
            ("cutsubleadsieie3" , cms.vdouble(0.0106, 0.0097, 0.028, 0.027) ),
            ("cutsubleadsieie4" , cms.vdouble(0.0104, 0.0094, 0.028, 0.025) ),
            ("cutsubleadsieie5" , cms.vdouble(0.0101, 0.0093, 0.027, 0.023) ),
            ("cutsubleadsieie6" , cms.vdouble(0.0099, 0.0092, 0.027, 0.023) ),
            ("cutsubleadsieie6c0" , cms.vdouble(0.0114, 0.011, 0.0101, 0.029, 0.029, 
            0.028) ),
            ("cutsubleadsieie6c1" , cms.vdouble(0.011, 0.0107, 0.0098, 0.029, 0.029, 
            0.027) ),
            ("cutsubleadsieie6c10" , cms.vdouble(0.0098, 0.0091, 0.0086, 0.024, 0.023, 
            0.022) ),
            ("cutsubleadsieie6c11" , cms.vdouble(0.012, 0.012, 0.012, 0.03, 0.03, 
            0.03) ),
            ("cutsubleadsieie6c2" , cms.vdouble(0.0108, 0.0105, 0.0097, 0.028, 0.028, 
            0.027) ),
            ("cutsubleadsieie6c3" , cms.vdouble(0.0106, 0.0102, 0.0095, 0.028, 0.028, 
            0.026) ),
            ("cutsubleadsieie6c4" , cms.vdouble(0.0105, 0.01, 0.0094, 0.028, 0.027, 
            0.026) ),
            ("cutsubleadsieie6c5" , cms.vdouble(0.0104, 0.0098, 0.0094, 0.028, 0.027, 
            0.024) ),
            ("cutsubleadsieie6c6" , cms.vdouble(0.0101, 0.0097, 0.0093, 0.028, 0.025, 
            0.023) ),
            ("cutsubleadsieie6c7" , cms.vdouble(0.0099, 0.0092, 0.0092, 0.028, 0.025, 
            0.022) ),
            ("cutsubleadsieie6c8" , cms.vdouble(0.0098, 0.0092, 0.0088, 0.025, 0.025, 
            0.022) ),
            ("cutsubleadsieie6c9" , cms.vdouble(0.0098, 0.0091, 0.0086, 0.024, 0.023, 
            0.022) ),
            ("cutsubleadsieie7" , cms.vdouble(0.0098, 0.009, 0.026, 0.023) ),
            ("cutsubleadsieie8" , cms.vdouble(0.0097, 0.009, 0.026, 0.023) ),
            ("cutsubleadsieie9" , cms.vdouble(0.0094, 0.009, 0.024, 0.023) ),
            ("cutsubleadtrkisooetom0" , cms.vdouble(7.5, 4.5, 5.2, 2.5) ),
            ("cutsubleadtrkisooetom1" , cms.vdouble(6.4, 3.4, 3.8, 2.1) ),
            ("cutsubleadtrkisooetom10" , cms.vdouble(0.43, 0.0111, 0.0075, 0.0073) ),
            ("cutsubleadtrkisooetom11" , cms.vdouble(7.0, 7.0, 4.0, 4.0) ),
            ("cutsubleadtrkisooetom2" , cms.vdouble(4.7, 2.9, 3.8, 1.63) ),
            ("cutsubleadtrkisooetom3" , cms.vdouble(3.5, 2.2, 2.3, 1.45) ),
            ("cutsubleadtrkisooetom4" , cms.vdouble(3.4, 1.86, 1.67, 1.44) ),
            ("cutsubleadtrkisooetom5" , cms.vdouble(2.9, 1.6, 1.55, 1.44) ),
            ("cutsubleadtrkisooetom6" , cms.vdouble(1.93, 1.4, 1.48, 0.056) ),
            ("cutsubleadtrkisooetom6c0" , cms.vdouble(7.9, 5.6, 4.5, 5.9, 4.2, 
            3) ),
            ("cutsubleadtrkisooetom6c1" , cms.vdouble(6.6, 4.7, 4.4, 5.7, 2.4, 
            2.2) ),
            ("cutsubleadtrkisooetom6c10" , cms.vdouble(1.6, 0.033, 0.00082, 0.68, 0.48, 
            0.26) ),
            ("cutsubleadtrkisooetom6c11" , cms.vdouble(7.0, 7.0, 7.0, 4.0, 4.0, 
            4.0) ),
            ("cutsubleadtrkisooetom6c2" , cms.vdouble(5.8, 4, 2.9, 4.6, 2.3, 
            1.93) ),
            ("cutsubleadtrkisooetom6c3" , cms.vdouble(4, 3.4, 2.5, 1.89, 1.7, 
            1.67) ),
            ("cutsubleadtrkisooetom6c4" , cms.vdouble(3.6, 2.9, 1.96, 1.63, 1.43, 
            1.36) ),
            ("cutsubleadtrkisooetom6c5" , cms.vdouble(2.7, 2.9, 1.56, 1.56, 1.06, 
            1.36) ),
            ("cutsubleadtrkisooetom6c6" , cms.vdouble(2.3, 2.6, 1.54, 1.56, 1.06, 
            1.07) ),
            ("cutsubleadtrkisooetom6c7" , cms.vdouble(2.1, 1.79, 1.07, 1.56, 0.8, 
            0.26) ),
            ("cutsubleadtrkisooetom6c8" , cms.vdouble(2.1, 1.55, 0.25, 0.8, 0.48, 
            0.26) ),
            ("cutsubleadtrkisooetom6c9" , cms.vdouble(1.83, 1.2, 0.029, 0.68, 0.48, 
            0.26) ),
            ("cutsubleadtrkisooetom7" , cms.vdouble(1.42, 0.76, 1.48, 0.056) ),
            ("cutsubleadtrkisooetom8" , cms.vdouble(1.21, 0.51, 0.27, 0.056) ),
            ("cutsubleadtrkisooetom9" , cms.vdouble(0.43, 0.4, 0.0075, 0.056) ),
            ]
        )
    ),
    isAOD = cms.bool(False),
    jecAK8PayloadNames = cms.vstring('Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt', 
        'Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt'),
    muonSrc = cms.InputTag("slimmedMuons"),
    newParticles = cms.vint32(4000011, 4000013, 1000006, 1000021, 1000022, 
        1000024, 1000025, 35),
    packedPFCands = cms.InputTag("packedPFCandidates"),
    patTriggerResults = cms.InputTag("TriggerResults","","PAT"),
    pfMETLabel = cms.InputTag("slimmedMETs"),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoChargedIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD","h+-DR030-"),
    phoChargedIsolation_PUPPI = cms.InputTag("egmPhotonIsolationMiniAODPUPPI","h+-DR030-"),
    phoLooseIdMap = cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-loose"),
    phoMVAValuesMap = cms.InputTag("photonMVAValueMapProducer","PhotonMVAEstimatorRun2Spring16NonTrigV1Values"),
    phoMediumIdMap = cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-medium"),
    phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
    phoNeutralHadronIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD","h0-DR030-"),
    phoNeutralHadronIsolation_PUPPI = cms.InputTag("egmPhotonIsolationMiniAODPUPPI","h0-DR030-"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    phoPhotonIsolation_CITK = cms.InputTag("egmPhotonIsolationMiniAOD","gamma-DR030-"),
    phoPhotonIsolation_PUPPI = cms.InputTag("egmPhotonIsolationMiniAODPUPPI","gamma-DR030-"),
    phoTightIdMap = cms.InputTag("egmPhotonIDs","cutBasedPhotonID-Spring16-V2p2-tight"),
    phoWorstChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoWorstChargedIsolation"),
    photonSrc = cms.InputTag("slimmedPhotons"),
    pileupCollection = cms.InputTag("slimmedAddPileupInfo"),
    recoPhotonSrc = cms.InputTag("reducedEgamma","reducedGedPhotonCores"),
    rhoCentralLabel = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    rhoLabel = cms.InputTag("fixedGridRhoFastjetAll"),
    runHFElectrons = cms.bool(False),
    runOnParticleGun = cms.bool(False),
    runOnSherpa = cms.bool(False),
    tauSrc = cms.InputTag("slimmedTaus"),
    trgFilterDeltaPtCut = cms.double(0.5),
    trgFilterDeltaRCut = cms.double(0.3),
    triggerEvent = cms.InputTag("selectedPatTrigger"),
    triggerResults = cms.InputTag("TriggerResults","","HLT")
)


process.ggNtuplizerGen = cms.EDAnalyzer("GenParticleInfo")


process.patCandidateSummary = cms.EDAnalyzer("CandidateSummaryTable",
    candidates = cms.VInputTag(cms.InputTag("patElectrons"), cms.InputTag("patMuons"), cms.InputTag("patTaus"), cms.InputTag("patPhotons"), cms.InputTag("patJets"), 
        cms.InputTag("patMETs")),
    logName = cms.untracked.string('patCandidates|PATSummaryTables')
)


process.selectedPatCandidateSummary = cms.EDAnalyzer("CandidateSummaryTable",
    candidates = cms.VInputTag(cms.InputTag("selectedPatElectrons"), cms.InputTag("selectedPatMuons"), cms.InputTag("selectedPatTaus"), cms.InputTag("selectedPatPhotons"), cms.InputTag("selectedPatJets")),
    logName = cms.untracked.string('selectedPatCanddiates|PATSummaryTables')
)


process.type0PFMEtCorrectionPFCandToVertexAssociationForValidation = cms.Sequence(cms.ignore(process.selectedVerticesForPFMEtCorrType0)+cms.ignore(process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0)+process.particleFlowDisplacedVertex+process.pfCandidateToVertexAssociation)


process.pfMEtSysShiftCorrSequence = cms.Sequence(process.pfMEtMultShiftCorr)


process.ak4L1JPTOffsetCorrectorChain = cms.Sequence(process.ak4CaloL1OffsetCorrector+process.ak4L1JPTOffsetCorrector)


process.patPFMetT1T2CorrSequence = cms.Sequence(process.patPFMetT1T2Corr)


process.patPFMetSmearCorrSequence = cms.Sequence(process.patSmearedJets+process.selectedPatJetsForMetT1T2SmearCorr+process.patPFMetT1T2SmearCorr)


process.ak4PFCHSL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFCHSL1FastjetCorrector+process.ak4PFCHSL2RelativeCorrector+process.ak4PFCHSL3AbsoluteCorrector+process.ak4PFCHSResidualCorrector+process.ak4PFCHSL1FastL2L3ResidualCorrector)


process.patShrinkingConePFTauDiscrimination = cms.Sequence()


process.ak4PFCHSL1FastL2L3CorrectorChain = cms.Sequence(process.ak4PFCHSL1FastjetCorrector+process.ak4PFCHSL2RelativeCorrector+process.ak4PFCHSL3AbsoluteCorrector+process.ak4PFCHSL1FastL2L3Corrector)


process.ak4CaloL1L2L3CorrectorChain = cms.Sequence(process.ak4CaloL1OffsetCorrector+process.ak4CaloL2RelativeCorrector+process.ak4CaloL3AbsoluteCorrector+process.ak4CaloL1L2L3Corrector)


process.ak4PFL1FastL2L3CorrectorChain = cms.Sequence(process.ak4PFL1FastjetCorrector+process.ak4PFL2RelativeCorrector+process.ak4PFL3AbsoluteCorrector+process.ak4PFL1FastL2L3Corrector)


process.patMetModuleSequence = cms.Sequence(process.pfMet+process.genMetExtractor+process.patJetCorrFactorsReapplyJEC+process.patJetsReapplyJEC+process.basicJetsForMet+process.jetSelectorForMet+process.cleanedPatJets+process.metrawCalo+process.patPFMet)


process.pfSortByTypePFBRECOSequence = cms.Sequence(process.pfAllNeutralHadronsPFBRECO+process.pfAllChargedHadronsPFBRECO+process.pfAllPhotonsPFBRECO+process.pfAllChargedParticlesPFBRECO+process.pfPileUpAllChargedParticlesPFBRECO+process.pfAllNeutralHadronsAndPhotonsPFBRECO)


process.ak4TrackL2L3CorrectorChain = cms.Sequence(process.ak4TrackL2RelativeCorrector+process.ak4TrackL3AbsoluteCorrector+process.ak4TrackL2L3Corrector)


process.patPFMetT1SmearpatShiftedModuleSequence = cms.Sequence(process.patPFMetT1SmearJetResDown+process.patPFMetT1SmearJetResUp+process.patPFMetT1SmearMuonEnUp+process.patPFMetT1SmearMuonEnDown+process.patPFMetT1SmearJetEnUp+process.patPFMetT1SmearJetEnDown+process.patPFMetT1SmearTauEnUp+process.patPFMetT1SmearTauEnDown+process.patPFMetT1SmearPhotonEnUp+process.patPFMetT1SmearPhotonEnDown+process.patPFMetT1SmearElectronEnDown+process.patPFMetT1SmearElectronEnUp+process.patPFMetT1SmearUnclusteredEnUp+process.patPFMetT1SmearUnclusteredEnDown)


process.electronPFIsolationDepositsPFBRECOSequence = cms.Sequence(process.elPFIsoDepositChargedPFBRECO+process.elPFIsoDepositChargedAllPFBRECO+process.elPFIsoDepositGammaPFBRECO+process.elPFIsoDepositNeutralPFBRECO+process.elPFIsoDepositPUPFBRECO)


process.ak4PFL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFL1FastjetCorrector+process.ak4PFL2RelativeCorrector+process.ak4PFL3AbsoluteCorrector+process.ak4PFResidualCorrector+process.ak4PFL1FastL2L3ResidualCorrector)


process.ak4JPTL1L2L3CorrectorChain = cms.Sequence(process.ak4L1JPTOffsetCorrectorChain+process.ak4JPTL2RelativeCorrector+process.ak4JPTL3AbsoluteCorrector+process.ak4JPTL1L2L3Corrector)


process.ak4PFPuppiL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFPuppiL1FastjetCorrector+process.ak4PFPuppiL2RelativeCorrector+process.ak4PFPuppiL3AbsoluteCorrector+process.ak4PFPuppiResidualCorrector+process.ak4PFPuppiL1FastL2L3ResidualCorrector)


process.patJetFlavourIdLegacy = cms.Sequence(process.patJetPartonsLegacy+process.patJetPartonAssociationLegacy+process.patJetFlavourAssociationLegacy)


process.pfNoPileUpIsoPFBRECOSequence = cms.Sequence(process.pfPileUpIsoPFBRECO+process.pfNoPileUpIsoPFBRECO)


process.ak4PFCHSL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFCHSL2RelativeCorrector+process.ak4PFCHSL3AbsoluteCorrector+process.ak4PFCHSResidualCorrector+process.ak4PFCHSL2L3ResidualCorrector)


process.patPFCandidateIsoDepositSelection = cms.Sequence(process.pfNoPileUpIsoPFBRECOSequence+process.pfSortByTypePFBRECOSequence)


process.updateHPSPFTaus = cms.Sequence()


process.selectedPatCandidates = cms.Sequence(process.selectedPatElectrons+process.selectedPatMuons+process.selectedPatTaus+process.selectedPatPhotons+process.selectedPatJets+process.selectedPatCandidateSummary)


process.ak4JPTL1FastL2L3CorrectorChain = cms.Sequence(process.ak4JPTL1FastjetCorrector+process.ak4JPTL2RelativeCorrector+process.ak4JPTL3AbsoluteCorrector+process.ak4JPTL1FastL2L3Corrector)


process.ak4PFCHSL2L3CorrectorChain = cms.Sequence(process.ak4PFCHSL2RelativeCorrector+process.ak4PFCHSL3AbsoluteCorrector+process.ak4PFCHSL2L3Corrector)


process.ak4PFL1L2L3CorrectorChain = cms.Sequence(process.ak4PFL1OffsetCorrector+process.ak4PFL2RelativeCorrector+process.ak4PFL3AbsoluteCorrector+process.ak4PFL1L2L3Corrector)


process.patCaloTauDiscrimination = cms.Sequence()


process.photonPFIsolationDepositsPFBRECOSequence = cms.Sequence(process.phPFIsoDepositChargedPFBRECO+process.phPFIsoDepositChargedAllPFBRECO+process.phPFIsoDepositGammaPFBRECO+process.phPFIsoDepositNeutralPFBRECO+process.phPFIsoDepositPUPFBRECO)


process.ak4CaloL1FastL2L3CorrectorChain = cms.Sequence(process.ak4CaloL1FastjetCorrector+process.ak4CaloL2RelativeCorrector+process.ak4CaloL3AbsoluteCorrector+process.ak4CaloL1FastL2L3Corrector)


process.patPFMetT2CorrSequence = cms.Sequence(process.patPFMetT2Corr)


process.ak4PFL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFL1OffsetCorrector+process.ak4PFL2RelativeCorrector+process.ak4PFL3AbsoluteCorrector+process.ak4PFResidualCorrector+process.ak4PFL1L2L3ResidualCorrector)


process.electronPFIsolationDepositsPATSequence = cms.Sequence(process.elPFIsoDepositChargedPAT+process.elPFIsoDepositChargedAllPAT+process.elPFIsoDepositGammaPAT+process.elPFIsoDepositNeutralPAT+process.elPFIsoDepositPUPAT)


process.ak4JPTL2L3ResidualCorrectorChain = cms.Sequence(process.ak4L1JPTOffsetCorrectorChain+process.ak4JPTL2RelativeCorrector+process.ak4JPTL3AbsoluteCorrector+process.ak4JPTResidualCorrector+process.ak4JPTL2L3ResidualCorrector)


process.ak4PFL1FastL2L3L6CorrectorChain = cms.Sequence(process.ak4PFL1FastjetCorrector+process.ak4PFL2RelativeCorrector+process.ak4PFL3AbsoluteCorrector+process.ak4PFL6SLBCorrector+process.ak4PFL1FastL2L3L6Corrector)


process.photonPFIsolationDepositsPATSequence = cms.Sequence(process.phPFIsoDepositChargedPAT+process.phPFIsoDepositChargedAllPAT+process.phPFIsoDepositGammaPAT+process.phPFIsoDepositNeutralPAT+process.phPFIsoDepositPUPAT)


process.patPFMetT1SmearpatMetUncertaintySequence = cms.Sequence(process.shiftedPatSmearedJetResDown+process.shiftedPatMETCorrSmearedJetResDown+process.shiftedPatSmearedJetResUp+process.shiftedPatMETCorrSmearedJetResUp)


process.patJetFlavourId = cms.Sequence(process.patJetPartons+process.patJetFlavourAssociation)


process.ak4PFPuppiL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFPuppiL2RelativeCorrector+process.ak4PFPuppiL3AbsoluteCorrector+process.ak4PFPuppiResidualCorrector+process.ak4PFPuppiL2L3ResidualCorrector)


process.ak4JPTL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4JPTL1FastjetCorrector+process.ak4JPTL2RelativeCorrector+process.ak4JPTL3AbsoluteCorrector+process.ak4JPTResidualCorrector+process.ak4JPTL1FastL2L3ResidualCorrector)


process.ak4PFL2L3L6CorrectorChain = cms.Sequence(process.ak4PFL2RelativeCorrector+process.ak4PFL3AbsoluteCorrector+process.ak4PFL6SLBCorrector+process.ak4PFL2L3L6Corrector)


process.ak4PFPuppiL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFPuppiL1OffsetCorrector+process.ak4PFPuppiL2RelativeCorrector+process.ak4PFPuppiL3AbsoluteCorrector+process.ak4PFPuppiResidualCorrector+process.ak4PFPuppiL1L2L3ResidualCorrector)


process.patPFTauIsolation = cms.Sequence(process.tauIsoDepositPFCandidates+process.tauIsoDepositPFChargedHadrons+process.tauIsoDepositPFNeutralHadrons+process.tauIsoDepositPFGammas)


process.pfNoPileUpPFBRECOSequence = cms.Sequence(process.pfPileUpPFBRECO+process.pfNoPileUpPFBRECO)


process.type0PFMEtCorrectionPFCandToVertexAssociation = cms.Sequence(process.selectedVerticesForPFMEtCorrType0+process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0+process.particleFlowDisplacedVertex+process.pfCandidateToVertexAssociation)


process.ak4CaloL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4CaloL1FastjetCorrector+process.ak4CaloL2RelativeCorrector+process.ak4CaloL3AbsoluteCorrector+process.ak4CaloResidualCorrector+process.ak4CaloL1FastL2L3ResidualCorrector)


process.regressionApplication = cms.Sequence(process.slimmedElectrons)


process.ak4CaloL2L3L6CorrectorChain = cms.Sequence(process.ak4CaloL2RelativeCorrector+process.ak4CaloL3AbsoluteCorrector+process.ak4CaloL6SLBCorrector+process.ak4CaloL2L3L6Corrector)


process.ak4PFPuppiL1FastL2L3CorrectorChain = cms.Sequence(process.ak4PFPuppiL1FastjetCorrector+process.ak4PFPuppiL2RelativeCorrector+process.ak4PFPuppiL3AbsoluteCorrector+process.ak4PFPuppiL1FastL2L3Corrector)


process.patPFMetT2SmearCorrSequence = cms.Sequence(process.patSmearedJets+process.selectedPatJetsForMetT1T2SmearCorr+process.selectedPatJetsForMetT2SmearCorr+process.patPFMetT1T2SmearCorr+process.patPFMetT2SmearCorr)


process.ak4JPTL2L3CorrectorChain = cms.Sequence(process.ak4L1JPTOffsetCorrectorChain+process.ak4JPTL2RelativeCorrector+process.ak4JPTL3AbsoluteCorrector+process.ak4JPTL2L3Corrector)


process.ggMETFiltersSequence = cms.Sequence(process.BadPFMuonFilter+process.BadChargedCandidateFilter)


process.ak4CaloL2L3CorrectorChain = cms.Sequence(process.ak4CaloL2RelativeCorrector+process.ak4CaloL3AbsoluteCorrector+process.ak4CaloL2L3Corrector)


process.ak4JPTL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4L1JPTOffsetCorrectorChain+process.ak4JPTL2RelativeCorrector+process.ak4JPTL3AbsoluteCorrector+process.ak4JPTResidualCorrector+process.ak4JPTL1L2L3ResidualCorrector)


process.ak4PFCHSL1L2L3CorrectorChain = cms.Sequence(process.ak4PFCHSL1OffsetCorrector+process.ak4PFCHSL2RelativeCorrector+process.ak4PFCHSL3AbsoluteCorrector+process.ak4PFCHSL1L2L3Corrector)


process.patPFMetTxyCorrSequence = cms.Sequence(process.patPFMetTxyCorr)


process.producePatPFMETCorrectionsUnc = cms.Sequence(process.patPFMet+process.pfCandsNotInJetsForMetCorr+process.selectedPatJetsForMetT1T2Corr+process.selectedPatJetsForMetT2Corr+process.patPFMetT1T2Corr+process.patPFMetT2Corr+process.type0PFMEtCorrectionPFCandToVertexAssociation+process.patPFMetT0Corr+process.pfCandMETcorr)


process.pfPhotonIsolationPATSequence = cms.Sequence(process.photonPFIsolationDepositsPATSequence+process.phPFIsoValueCharged03PFIdPAT+process.phPFIsoValueChargedAll03PFIdPAT+process.phPFIsoValueGamma03PFIdPAT+process.phPFIsoValueNeutral03PFIdPAT+process.phPFIsoValuePU03PFIdPAT+process.phPFIsoValueCharged04PFIdPAT+process.phPFIsoValueChargedAll04PFIdPAT+process.phPFIsoValueGamma04PFIdPAT+process.phPFIsoValueNeutral04PFIdPAT+process.phPFIsoValuePU04PFIdPAT)


process.pfParticleSelectionPFBRECOSequence = cms.Sequence(process.pfNoPileUpIsoPFBRECOSequence+process.pfNoPileUpPFBRECOSequence+process.pfSortByTypePFBRECOSequence)


process.ggPhotonIso_CITK_PUPPI_Sequence = cms.Sequence(process.puppi+process.egmPhotonIsolationMiniAOD+process.egmPhotonIsolationMiniAODPUPPI)


process.producePatPFMETCorrections = cms.Sequence(process.patPFMet+process.pfCandsNotInJetsForMetCorr+process.selectedPatJetsForMetT1T2Corr+process.selectedPatJetsForMetT2Corr+process.patPFMetT1T2Corr+process.patPFMetT2Corr+process.type0PFMEtCorrectionPFCandToVertexAssociation+process.patPFMetT0Corr+process.pfCandMETcorr+process.patPFMetT1+process.patPFMetT1T2+process.patPFMetT0pcT1+process.patPFMetT0pcT1T2)


process.patJetCorrections = cms.Sequence(process.patJetCorrFactors)


process.muonPFIsolationDepositsPATSequence = cms.Sequence(process.muPFIsoDepositChargedPAT+process.muPFIsoDepositChargedAllPAT+process.muPFIsoDepositGammaPAT+process.muPFIsoDepositNeutralPAT+process.muPFIsoDepositPUPAT)


process.patFixedConePFTauDiscrimination = cms.Sequence()


process.type0PFMEtCorrection = cms.Sequence(process.type0PFMEtCorrectionPFCandToVertexAssociation+process.pfMETcorrType0)


process.muonPFIsolationPATSequence = cms.Sequence(process.muonPFIsolationDepositsPATSequence+process.muPFIsoValueCharged03PAT+process.muPFMeanDRIsoValueCharged03PAT+process.muPFSumDRIsoValueCharged03PAT+process.muPFIsoValueChargedAll03PAT+process.muPFMeanDRIsoValueChargedAll03PAT+process.muPFSumDRIsoValueChargedAll03PAT+process.muPFIsoValueGamma03PAT+process.muPFMeanDRIsoValueGamma03PAT+process.muPFSumDRIsoValueGamma03PAT+process.muPFIsoValueNeutral03PAT+process.muPFMeanDRIsoValueNeutral03PAT+process.muPFSumDRIsoValueNeutral03PAT+process.muPFIsoValueGammaHighThreshold03PAT+process.muPFMeanDRIsoValueGammaHighThreshold03PAT+process.muPFSumDRIsoValueGammaHighThreshold03PAT+process.muPFIsoValueNeutralHighThreshold03PAT+process.muPFMeanDRIsoValueNeutralHighThreshold03PAT+process.muPFSumDRIsoValueNeutralHighThreshold03PAT+process.muPFIsoValuePU03PAT+process.muPFMeanDRIsoValuePU03PAT+process.muPFSumDRIsoValuePU03PAT+process.muPFIsoValueCharged04PAT+process.muPFMeanDRIsoValueCharged04PAT+process.muPFSumDRIsoValueCharged04PAT+process.muPFIsoValueChargedAll04PAT+process.muPFMeanDRIsoValueChargedAll04PAT+process.muPFSumDRIsoValueChargedAll04PAT+process.muPFIsoValueGamma04PAT+process.muPFMeanDRIsoValueGamma04PAT+process.muPFSumDRIsoValueGamma04PAT+process.muPFIsoValueNeutral04PAT+process.muPFMeanDRIsoValueNeutral04PAT+process.muPFSumDRIsoValueNeutral04PAT+process.muPFIsoValueGammaHighThreshold04PAT+process.muPFMeanDRIsoValueGammaHighThreshold04PAT+process.muPFSumDRIsoValueGammaHighThreshold04PAT+process.muPFIsoValueNeutralHighThreshold04PAT+process.muPFMeanDRIsoValueNeutralHighThreshold04PAT+process.muPFSumDRIsoValueNeutralHighThreshold04PAT+process.muPFIsoValuePU04PAT+process.muPFMeanDRIsoValuePU04PAT+process.muPFSumDRIsoValuePU04PAT)


process.egmGsfElectronIDSequence = cms.Sequence(process.electronMVAValueMapProducer+process.heepIDVarValueMaps+process.heepIDVarValueMaps+process.heepIDVarValueMaps+process.heepIDVarValueMaps+process.heepIDVarValueMaps+process.heepIDVarValueMaps+process.heepIDVarValueMaps+process.heepIDVarValueMaps+process.heepIDVarValueMaps+process.heepIDVarValueMaps+process.egmGsfElectronIDs+process.electronRegressionValueMapProducer)


process.pfElectronIsolationPATSequence = cms.Sequence(process.electronPFIsolationDepositsPATSequence+process.elPFIsoValueCharged03PFIdPAT+process.elPFIsoValueChargedAll03PFIdPAT+process.elPFIsoValueGamma03PFIdPAT+process.elPFIsoValueNeutral03PFIdPAT+process.elPFIsoValuePU03PFIdPAT+process.elPFIsoValueCharged04PFIdPAT+process.elPFIsoValueChargedAll04PFIdPAT+process.elPFIsoValueGamma04PFIdPAT+process.elPFIsoValueNeutral04PFIdPAT+process.elPFIsoValuePU04PFIdPAT+process.elPFIsoValueCharged03NoPFIdPAT+process.elPFIsoValueChargedAll03NoPFIdPAT+process.elPFIsoValueGamma03NoPFIdPAT+process.elPFIsoValueNeutral03NoPFIdPAT+process.elPFIsoValuePU03NoPFIdPAT+process.elPFIsoValueCharged04NoPFIdPAT+process.elPFIsoValueChargedAll04NoPFIdPAT+process.elPFIsoValueGamma04NoPFIdPAT+process.elPFIsoValueNeutral04NoPFIdPAT+process.elPFIsoValuePU04NoPFIdPAT)


process.correctionTermsPfMetType0PFCandidateForValidation = cms.Sequence(process.type0PFMEtCorrectionPFCandToVertexAssociationForValidation+process.corrPfMetType0PfCand)


process.muonPFIsolationDepositsPFBRECOSequence = cms.Sequence(process.muPFIsoDepositChargedPFBRECO+process.muPFIsoDepositChargedAllPFBRECO+process.muPFIsoDepositGammaPFBRECO+process.muPFIsoDepositNeutralPFBRECO+process.muPFIsoDepositPUPFBRECO)


process.egmPhotonIsolationMiniAODSequence = cms.Sequence(process.egmPhotonIsolation)


process.ak4CaloL1FastL2L3L6CorrectorChain = cms.Sequence(process.ak4CaloL1FastjetCorrector+process.ak4CaloL2RelativeCorrector+process.ak4CaloL3AbsoluteCorrector+process.ak4CaloL6SLBCorrector+process.ak4CaloL1FastL2L3L6Corrector)


process.type0PFMEtCorrectionPFCandToVertexAssociationForValidationMiniAOD = cms.Sequence(cms.ignore(process.selectedVerticesForPFMEtCorrType0)+cms.ignore(process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0)+process.particleFlowDisplacedVertex+process.pfCandidateToVertexAssociation)


process.ak4PFL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFL2RelativeCorrector+process.ak4PFL3AbsoluteCorrector+process.ak4PFResidualCorrector+process.ak4PFL2L3ResidualCorrector)


process.ak4PFPuppiL2L3CorrectorChain = cms.Sequence(process.ak4PFPuppiL2RelativeCorrector+process.ak4PFPuppiL3AbsoluteCorrector+process.ak4PFPuppiL2L3Corrector)


process.ak4CaloL2L3ResidualCorrectorChain = cms.Sequence(process.ak4CaloL2RelativeCorrector+process.ak4CaloL3AbsoluteCorrector+process.ak4CaloResidualCorrector+process.ak4CaloL2L3ResidualCorrector)


process.ak4PFCHSL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFCHSL1OffsetCorrector+process.ak4PFCHSL2RelativeCorrector+process.ak4PFCHSL3AbsoluteCorrector+process.ak4PFCHSResidualCorrector+process.ak4PFCHSL1L2L3ResidualCorrector)


process.ak4PFL2L3CorrectorChain = cms.Sequence(process.ak4PFL2RelativeCorrector+process.ak4PFL3AbsoluteCorrector+process.ak4PFL2L3Corrector)


process.ak4CaloL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4CaloL1OffsetCorrector+process.ak4CaloL2RelativeCorrector+process.ak4CaloL3AbsoluteCorrector+process.ak4CaloResidualCorrector+process.ak4CaloL1L2L3ResidualCorrector)


process.ggPhotonIso_CITK_Sequence = cms.Sequence(process.egmPhotonIsolationMiniAOD)


process.ak4PFPuppiL1L2L3CorrectorChain = cms.Sequence(process.ak4PFPuppiL1OffsetCorrector+process.ak4PFPuppiL2RelativeCorrector+process.ak4PFPuppiL3AbsoluteCorrector+process.ak4PFPuppiL1L2L3Corrector)


process.patMetUncertaintySequence = cms.Sequence(process.ak4PFCHSL1FastL2L3CorrectorChain+process.ak4PFCHSL1FastL2L3ResidualCorrectorChain+(process.shiftedPatJetResDown+process.shiftedPatMETCorrJetResDown+process.shiftedPatJetResUp+process.shiftedPatMETCorrJetResUp+process.pfCandsNoJets+process.pfCandsNoJetsNoEle+process.pfCandsNoJetsNoEleNoMu+process.pfCandsNoJetsNoEleNoMuNoTau+process.pfCandsForUnclusteredUnc+process.shiftedPatMuonEnDown+process.shiftedPatMETCorrMuonEnDown+process.shiftedPatMuonEnUp+process.shiftedPatMETCorrMuonEnUp+process.shiftedPatJetEnDown+process.shiftedPatMETCorrJetEnDown+process.shiftedPatJetEnUp+process.shiftedPatMETCorrJetEnUp+process.shiftedPatTauEnDown+process.shiftedPatMETCorrTauEnDown+process.shiftedPatTauEnUp+process.shiftedPatMETCorrTauEnUp+process.shiftedPatPhotonEnDown+process.shiftedPatMETCorrPhotonEnDown+process.shiftedPatPhotonEnUp+process.shiftedPatMETCorrPhotonEnUp+process.shiftedPatElectronEnDown+process.shiftedPatMETCorrElectronEnDown+process.shiftedPatElectronEnUp+process.shiftedPatMETCorrElectronEnUp+process.shiftedPatUnclusteredEnDown+process.shiftedPatMETCorrUnclusteredEnDown+process.shiftedPatUnclusteredEnUp+process.shiftedPatMETCorrUnclusteredEnUp)+process.patPFMetT1SmearpatMetUncertaintySequence)


process.patHPSPFTauDiscrimination = cms.Sequence(process.updateHPSPFTaus)


process.correctionTermsPfMetType0PFCandidate = cms.Sequence(process.type0PFMEtCorrectionPFCandToVertexAssociation+process.corrPfMetType0PfCand)


process.correctionTermsPfMetType1Type2 = cms.Sequence(process.pfJetsPtrForMetCorr+process.particleFlowPtrs+process.pfCandsNotInJetsPtrForMetCorr+process.pfCandsNotInJetsForMetCorr+process.pfCandMETcorr+process.ak4PFCHSL1FastL2L3ResidualCorrectorChain+process.ak4PFCHSL1FastL2L3Corrector+process.corrPfMetType1+process.corrPfMetType2)


process.correctionTermsCaloMet = cms.Sequence(process.ak4CaloL2L3CorrectorChain+process.corrCaloMetType1+process.muCaloMetCorr+process.corrCaloMetType2)


process.makePatTaus = cms.Sequence(process.patHPSPFTauDiscrimination+process.patPFCandidateIsoDepositSelection+process.patPFTauIsolation+process.tauMatch+process.tauGenJets+process.tauGenJetsSelectorAllHadrons+process.tauGenJetMatch+process.patTaus)


process.egmPhotonIDSequence = cms.Sequence(process.egmPhotonIsolationMiniAODSequence+process.photonIDValueMapProducer+process.photonMVAValueMapProducer+process.egmPhotonIDs+process.photonRegressionValueMapProducer)


process.pfParticleSelectionForIsoSequence = cms.Sequence(process.particleFlowPtrs+process.pfParticleSelectionPFBRECOSequence)


process.patShiftedModuleSequence = cms.Sequence(process.patPFMetT1JetResUp+process.patPFMetT1JetResDown+process.patPFMetT1MuonEnUp+process.patPFMetT1MuonEnDown+process.patPFMetT1JetEnUp+process.patPFMetT1JetEnDown+process.patPFMetT1TauEnDown+process.patPFMetT1TauEnUp+process.patPFMetT1PhotonEnUp+process.patPFMetT1PhotonEnDown+process.patPFMetT1ElectronEnUp+process.patPFMetT1ElectronEnDown+process.patPFMetT1UnclusteredEnUp+process.patPFMetT1UnclusteredEnDown+process.patPFMetT1SmearpatShiftedModuleSequence)


process.patPFMetT0CorrSequence = cms.Sequence(process.type0PFMEtCorrectionPFCandToVertexAssociation+process.patPFMetT0Corr)


process.patMetCorrectionSequence = cms.Sequence(process.patPFMetT1T2CorrSequence+process.patPFMetT1+process.patPFMetTxyCorrSequence+process.patPFMetT1Txy+process.patPFMetTxy+process.patPFMetSmearCorrSequence+process.patPFMetT1Smear)


process.patMETCorrections = cms.Sequence(process.correctionTermsCaloMet+process.caloMetT1+process.caloMetT1T2+process.correctionTermsPfMetType1Type2+process.pfMetT1+process.pfMetT1T2)


process.makePatJets = cms.Sequence(process.patJetCorrections+process.patJetCharge+process.patJetPartonMatch+process.patJetGenJetMatch+process.patJetFlavourIdLegacy+process.patJetFlavourId+process.patJets)


process.makePatMuons = cms.Sequence(process.pfParticleSelectionForIsoSequence+process.muonPFIsolationPATSequence+process.muonMatch+process.patMuons)


process.makePatElectrons = cms.Sequence(process.pfParticleSelectionForIsoSequence+process.pfElectronIsolationPATSequence+process.electronMatch+process.patElectrons)


process.makePatMETs = cms.Sequence(process.patMETCorrections+process.patMETs)


process.makePatPhotons = cms.Sequence(process.pfParticleSelectionForIsoSequence+process.pfPhotonIsolationPATSequence+process.photonMatch+process.patPhotons)


process.fullPatMetSequence = cms.Sequence(process.patMetModuleSequence+process.patMetCorrectionSequence+process.patMetUncertaintySequence+process.patShiftedModuleSequence+process.patCaloMet+process.slimmedMETs)


process.patCandidates = cms.Sequence(process.makePatElectrons+process.makePatMuons+process.makePatTaus+process.makePatPhotons+process.makePatJets+process.makePatMETs+process.patCandidateSummary)


process.p = cms.Path(process.regressionApplication+process.ggMETFiltersSequence+process.calibratedPatElectrons+process.egmGsfElectronIDSequence+process.ggNtuplizer)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring('FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1000)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring('warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatElectrons = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(12345)
    ),
    calibratedPatPhotons = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(12345)
    ),
    ggNtuplizer = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(67890)
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('ggtree_mc.root')
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring('HCAL', 
        'ZDC', 
        'CASTOR', 
        'EcalBarrel', 
        'EcalEndcap', 
        'EcalPreshower', 
        'TOWER')
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False),
    hcalTopologyConstants = cms.PSet(
        maxDepthHB = cms.int32(2),
        maxDepthHE = cms.int32(3),
        mode = cms.string('HcalTopologyMode::LHC')
    )
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer")


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True),
    hcalTopologyConstants = cms.PSet(
        maxDepthHB = cms.int32(2),
        maxDepthHE = cms.int32(3),
        mode = cms.string('HcalTopologyMode::LHC')
    )
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.ParabolicParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
    label = cms.untracked.string('ParabolicMf'),
    valueOverride = cms.int32(-1),
    version = cms.string('Parabolic')
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    compatibiltyWith11 = cms.untracked.bool(True),
    useDDD = cms.untracked.bool(False)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer")


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducerFromDB",
    debugBuilder = cms.untracked.bool(False),
    label = cms.untracked.string(''),
    valueOverride = cms.int32(-1)
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.ak10PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak10PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFCHSL1Fastjet', 
        'ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute', 
        'ak10PFCHSResidual')
)


process.ak10PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFCHSL1Offset', 
        'ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute', 
        'ak10PFCHSResidual')
)


process.ak10PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak10PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute')
)


process.ak10PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute', 
        'ak10PFCHSResidual')
)


process.ak10PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L2Relative')
)


process.ak10PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L3Absolute')
)


process.ak10PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak10PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak10PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFL1Fastjet', 
        'ak10PFL2Relative', 
        'ak10PFL3Absolute', 
        'ak10PFResidual')
)


process.ak10PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFL1Offset', 
        'ak10PFL2Relative', 
        'ak10PFL3Absolute', 
        'ak10PFResidual')
)


process.ak10PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak10PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFL2Relative', 
        'ak10PFL3Absolute')
)


process.ak10PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFL2Relative', 
        'ak10PFL3Absolute', 
        'ak10PFResidual')
)


process.ak10PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L2Relative')
)


process.ak10PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L3Absolute')
)


process.ak10PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L2L3Residual')
)


process.ak1PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak1PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFCHSL1Fastjet', 
        'ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute', 
        'ak1PFCHSResidual')
)


process.ak1PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFCHSL1Offset', 
        'ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute', 
        'ak1PFCHSResidual')
)


process.ak1PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak1PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute')
)


process.ak1PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute', 
        'ak1PFCHSResidual')
)


process.ak1PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L2Relative')
)


process.ak1PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L3Absolute')
)


process.ak1PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak1PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak1PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFL1Fastjet', 
        'ak1PFL2Relative', 
        'ak1PFL3Absolute', 
        'ak1PFResidual')
)


process.ak1PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFL1Offset', 
        'ak1PFL2Relative', 
        'ak1PFL3Absolute', 
        'ak1PFResidual')
)


process.ak1PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak1PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFL2Relative', 
        'ak1PFL3Absolute')
)


process.ak1PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFL2Relative', 
        'ak1PFL3Absolute', 
        'ak1PFResidual')
)


process.ak1PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L2Relative')
)


process.ak1PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L3Absolute')
)


process.ak1PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L2L3Residual')
)


process.ak2PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak2PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFCHSL1Fastjet', 
        'ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute', 
        'ak2PFCHSResidual')
)


process.ak2PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFCHSL1Offset', 
        'ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute', 
        'ak2PFCHSResidual')
)


process.ak2PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak2PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute')
)


process.ak2PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute', 
        'ak2PFCHSResidual')
)


process.ak2PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L2Relative')
)


process.ak2PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L3Absolute')
)


process.ak2PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak2PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak2PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFL1Fastjet', 
        'ak2PFL2Relative', 
        'ak2PFL3Absolute', 
        'ak2PFResidual')
)


process.ak2PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFL1Offset', 
        'ak2PFL2Relative', 
        'ak2PFL3Absolute', 
        'ak2PFResidual')
)


process.ak2PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak2PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFL2Relative', 
        'ak2PFL3Absolute')
)


process.ak2PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFL2Relative', 
        'ak2PFL3Absolute', 
        'ak2PFResidual')
)


process.ak2PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L2Relative')
)


process.ak2PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L3Absolute')
)


process.ak2PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L2L3Residual')
)


process.ak3PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak3PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFCHSL1Fastjet', 
        'ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute', 
        'ak3PFCHSResidual')
)


process.ak3PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFCHSL1Offset', 
        'ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute', 
        'ak3PFCHSResidual')
)


process.ak3PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak3PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute')
)


process.ak3PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute', 
        'ak3PFCHSResidual')
)


process.ak3PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L2Relative')
)


process.ak3PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L3Absolute')
)


process.ak3PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak3PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak3PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFL1Fastjet', 
        'ak3PFL2Relative', 
        'ak3PFL3Absolute', 
        'ak3PFResidual')
)


process.ak3PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFL1Offset', 
        'ak3PFL2Relative', 
        'ak3PFL3Absolute', 
        'ak3PFResidual')
)


process.ak3PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak3PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFL2Relative', 
        'ak3PFL3Absolute')
)


process.ak3PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFL2Relative', 
        'ak3PFL3Absolute', 
        'ak3PFResidual')
)


process.ak3PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L2Relative')
)


process.ak3PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L3Absolute')
)


process.ak3PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L2L3Residual')
)


process.ak4CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute')
)


process.ak4CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloL6SLB')
)


process.ak4CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloResidual')
)


process.ak4CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak4CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Offset', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute')
)


process.ak4CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Offset', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloResidual')
)


process.ak4CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL2Relative', 
        'ak4CaloL3Absolute')
)


process.ak4CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloL6SLB')
)


process.ak4CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloResidual')
)


process.ak4CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2Relative')
)


process.ak4CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L3Absolute')
)


process.ak4CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4CaloJetsSoftMuonTagInfos")
)


process.ak4CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.ak4JPTL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4JPTL1Fastjet', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute')
)


process.ak4JPTL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4JPTL1Fastjet', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute', 
        'ak4JPTResidual')
)


process.ak4JPTL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak4JPTL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute')
)


process.ak4JPTL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute', 
        'ak4JPTResidual')
)


process.ak4JPTL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4JPTL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute')
)


process.ak4JPTL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute', 
        'ak4JPTResidual')
)


process.ak4JPTL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L2Relative')
)


process.ak4JPTL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L3Absolute')
)


process.ak4JPTResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L2L3Residual')
)


process.ak4L1JPTOffset = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK5JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak4CaloL1Offset')
)


process.ak4PFCHSL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL1Fastjet', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute')
)


process.ak4PFCHSL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL1Fastjet', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute', 
        'ak4PFCHSResidual')
)


process.ak4PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFCHSL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL1Offset', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute')
)


process.ak4PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL1Offset', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute', 
        'ak4PFCHSResidual')
)


process.ak4PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute')
)


process.ak4PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute', 
        'ak4PFCHSResidual')
)


process.ak4PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2Relative')
)


process.ak4PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L3Absolute')
)


process.ak4PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak4PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute')
)


process.ak4PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFL6SLB')
)


process.ak4PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFResidual')
)


process.ak4PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Offset', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute')
)


process.ak4PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Offset', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFResidual')
)


process.ak4PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL2Relative', 
        'ak4PFL3Absolute')
)


process.ak4PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFL6SLB')
)


process.ak4PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFResidual')
)


process.ak4PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2Relative')
)


process.ak4PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L3Absolute')
)


process.ak4PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4PFJetsSoftMuonTagInfos")
)


process.ak4PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2L3Residual')
)


process.ak4TrackL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ak4TrackL2Relative', 
        'ak4TrackL3Absolute')
)


process.ak4TrackL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4TrackL2Relative', 
        'ak4TrackL3Absolute')
)


process.ak4TrackL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5TRK'),
    level = cms.string('L2Relative')
)


process.ak4TrackL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5TRK'),
    level = cms.string('L3Absolute')
)


process.ak5PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak5PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFCHSL1Fastjet', 
        'ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute', 
        'ak5PFCHSResidual')
)


process.ak5PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFCHSL1Offset', 
        'ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute', 
        'ak5PFCHSResidual')
)


process.ak5PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak5PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute')
)


process.ak5PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute', 
        'ak5PFCHSResidual')
)


process.ak5PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L2Relative')
)


process.ak5PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L3Absolute')
)


process.ak5PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak5PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak5PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual')
)


process.ak5PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFL1Offset', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual')
)


process.ak5PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak5PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFL2Relative', 
        'ak5PFL3Absolute')
)


process.ak5PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual')
)


process.ak5PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2Relative')
)


process.ak5PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L3Absolute')
)


process.ak5PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2L3Residual')
)


process.ak6PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak6PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFCHSL1Fastjet', 
        'ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute', 
        'ak6PFCHSResidual')
)


process.ak6PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFCHSL1Offset', 
        'ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute', 
        'ak6PFCHSResidual')
)


process.ak6PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak6PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute')
)


process.ak6PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute', 
        'ak6PFCHSResidual')
)


process.ak6PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L2Relative')
)


process.ak6PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L3Absolute')
)


process.ak6PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak6PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak6PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFL1Fastjet', 
        'ak6PFL2Relative', 
        'ak6PFL3Absolute', 
        'ak6PFResidual')
)


process.ak6PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFL1Offset', 
        'ak6PFL2Relative', 
        'ak6PFL3Absolute', 
        'ak6PFResidual')
)


process.ak6PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak6PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFL2Relative', 
        'ak6PFL3Absolute')
)


process.ak6PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFL2Relative', 
        'ak6PFL3Absolute', 
        'ak6PFResidual')
)


process.ak6PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L2Relative')
)


process.ak6PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L3Absolute')
)


process.ak6PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L2L3Residual')
)


process.ak7CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute')
)


process.ak7CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloL6SLB')
)


process.ak7CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL1Fastjet', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual')
)


process.ak7CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak7CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute')
)


process.ak7CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual')
)


process.ak7CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL2Relative', 
        'ak7CaloL3Absolute')
)


process.ak7CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloL6SLB')
)


process.ak7CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual')
)


process.ak7CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L2Relative')
)


process.ak7CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L3Absolute')
)


process.ak7CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak7CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak7CaloJetsSoftMuonTagInfos")
)


process.ak7CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L2L3Residual')
)


process.ak7JPTL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7JPTL1Fastjet', 
        'ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute', 
        'ak7JPTResidual')
)


process.ak7JPTL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7JPT'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak7JPTL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7JPTL1Offset', 
        'ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute')
)


process.ak7JPTL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7JPTL1Offset', 
        'ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute', 
        'ak7JPTResidual')
)


process.ak7JPTL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7JPT'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7JPTL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute')
)


process.ak7L1JPTOffset = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK7JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak4CaloL1Offset')
)


process.ak7PFCHSL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL1Fastjet', 
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute')
)


process.ak7PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak7PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFCHSL1Fastjet', 
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute', 
        'ak7PFCHSResidual')
)


process.ak7PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFCHSL1Offset', 
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute', 
        'ak7PFCHSResidual')
)


process.ak7PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute')
)


process.ak7PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute', 
        'ak7PFCHSResidual')
)


process.ak7PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L2Relative')
)


process.ak7PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L3Absolute')
)


process.ak7PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak7PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute')
)


process.ak7PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFL6SLB')
)


process.ak7PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak7PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual')
)


process.ak7PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL1Offset', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute')
)


process.ak7PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL1Offset', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual')
)


process.ak7PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL2Relative', 
        'ak7PFL3Absolute')
)


process.ak7PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFL6SLB')
)


process.ak7PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual')
)


process.ak7PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L2Relative')
)


process.ak7PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L3Absolute')
)


process.ak7PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak7PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak7PFJetsSoftMuonTagInfos")
)


process.ak7PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L2L3Residual')
)


process.ak8PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak8PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFCHSL1Fastjet', 
        'ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute', 
        'ak8PFCHSResidual')
)


process.ak8PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFCHSL1Offset', 
        'ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute', 
        'ak8PFCHSResidual')
)


process.ak8PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak8PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute')
)


process.ak8PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute', 
        'ak8PFCHSResidual')
)


process.ak8PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L2Relative')
)


process.ak8PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L3Absolute')
)


process.ak8PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak8PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak8PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFL1Fastjet', 
        'ak8PFL2Relative', 
        'ak8PFL3Absolute', 
        'ak8PFResidual')
)


process.ak8PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFL1Offset', 
        'ak8PFL2Relative', 
        'ak8PFL3Absolute', 
        'ak8PFResidual')
)


process.ak8PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak8PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFL2Relative', 
        'ak8PFL3Absolute')
)


process.ak8PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFL2Relative', 
        'ak8PFL3Absolute', 
        'ak8PFResidual')
)


process.ak8PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L2Relative')
)


process.ak8PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L3Absolute')
)


process.ak8PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L2L3Residual')
)


process.ak9PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak9PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFCHSL1Fastjet', 
        'ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute', 
        'ak9PFCHSResidual')
)


process.ak9PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFCHSL1Offset', 
        'ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute', 
        'ak9PFCHSResidual')
)


process.ak9PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak9PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute')
)


process.ak9PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute', 
        'ak9PFCHSResidual')
)


process.ak9PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L2Relative')
)


process.ak9PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L3Absolute')
)


process.ak9PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak9PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak9PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFL1Fastjet', 
        'ak9PFL2Relative', 
        'ak9PFL3Absolute', 
        'ak9PFResidual')
)


process.ak9PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFL1Offset', 
        'ak9PFL2Relative', 
        'ak9PFL3Absolute', 
        'ak9PFResidual')
)


process.ak9PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak9PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFL2Relative', 
        'ak9PFL3Absolute')
)


process.ak9PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFL2Relative', 
        'ak9PFL3Absolute', 
        'ak9PFResidual')
)


process.ak9PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L2Relative')
)


process.ak9PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L3Absolute')
)


process.ak9PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L2L3Residual')
)


process.candidateBoostedDoubleSecondaryVertexAK8Computer = cms.ESProducer("CandidateBoostedDoubleSecondaryVertexESProducer",
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SecondaryVertex/data/BoostedDoubleSV_AK8_BDT_v4.weights.xml.gz')
)


process.candidateBoostedDoubleSecondaryVertexCA15Computer = cms.ESProducer("CandidateBoostedDoubleSecondaryVertexESProducer",
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SecondaryVertex/data/BoostedDoubleSV_CA15_BDT_v3.weights.xml.gz')
)


process.candidateChargeBTagComputer = cms.ESProducer("CandidateChargeBTagESProducer",
    gbrForestLabel = cms.string(''),
    jetChargeExp = cms.double(0.8),
    svChargeExp = cms.double(0.5),
    useAdaBoost = cms.bool(True),
    useCondDB = cms.bool(False),
    weightFile = cms.FileInPath('RecoBTag/Combined/data/ChargeBTag_16_1_2016.weights.xml.gz')
)


process.candidateCombinedMVAV2Computer = cms.ESProducer("CombinedMVAV2JetTagESProducer",
    gbrForestLabel = cms.string('btag_CombinedMVAv2_BDT'),
    jetTagComputers = cms.vstring('candidateJetProbabilityComputer', 
        'candidateJetBProbabilityComputer', 
        'candidateCombinedSecondaryVertexV2Computer', 
        'softPFMuonComputer', 
        'softPFElectronComputer'),
    mvaName = cms.string('bdt'),
    spectators = cms.vstring(),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    variables = cms.vstring('Jet_CSV', 
        'Jet_CSVIVF', 
        'Jet_JP', 
        'Jet_JBP', 
        'Jet_SoftMu', 
        'Jet_SoftEl'),
    weightFile = cms.FileInPath('RecoBTag/Combined/data/CombinedMVAV2_13_07_2015.weights.xml.gz')
)


process.candidateCombinedSecondaryVertexSoftLeptonComputer = cms.ESProducer("CandidateCombinedSecondaryVertexSoftLeptonESProducer",
    SoftLeptonFlip = cms.bool(False),
    calibrationRecords = cms.vstring('CombinedSVRecoVertexNoSoftLepton', 
        'CombinedSVPseudoVertexNoSoftLepton', 
        'CombinedSVNoVertexNoSoftLepton', 
        'CombinedSVRecoVertexSoftMuon', 
        'CombinedSVPseudoVertexSoftMuon', 
        'CombinedSVNoVertexSoftMuon', 
        'CombinedSVRecoVertexSoftElectron', 
        'CombinedSVPseudoVertexSoftElectron', 
        'CombinedSVNoVertexSoftElectron'),
    categoryVariableName = cms.string('vertexLeptonCategory'),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    recordLabel = cms.string(''),
    trackFlip = cms.bool(False),
    trackMultiplicityMin = cms.uint32(2),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(False)
)


process.candidateCombinedSecondaryVertexSoftLeptonCvsLComputer = cms.ESProducer("CandidateCombinedSecondaryVertexSoftLeptonCvsLESProducer",
    SoftLeptonFlip = cms.bool(False),
    calibrationRecords = cms.vstring('CombinedSVRecoVertexNoSoftLeptonCvsL', 
        'CombinedSVPseudoVertexNoSoftLeptonCvsL', 
        'CombinedSVNoVertexNoSoftLeptonCvsL', 
        'CombinedSVRecoVertexSoftMuonCvsL', 
        'CombinedSVPseudoVertexSoftMuonCvsL', 
        'CombinedSVNoVertexSoftMuonCvsL', 
        'CombinedSVRecoVertexSoftElectronCvsL', 
        'CombinedSVPseudoVertexSoftElectronCvsL', 
        'CombinedSVNoVertexSoftElectronCvsL'),
    categoryVariableName = cms.string('vertexLeptonCategory'),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    recordLabel = cms.string(''),
    trackFlip = cms.bool(False),
    trackMultiplicityMin = cms.uint32(2),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(False)
)


process.candidateCombinedSecondaryVertexV2Computer = cms.ESProducer("CandidateCombinedSecondaryVertexESProducer",
    SoftLeptonFlip = cms.bool(False),
    calibrationRecords = cms.vstring('CombinedSVIVFV2RecoVertex', 
        'CombinedSVIVFV2PseudoVertex', 
        'CombinedSVIVFV2NoVertex'),
    categoryVariableName = cms.string('vertexCategory'),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    recordLabel = cms.string(''),
    trackFlip = cms.bool(False),
    trackMultiplicityMin = cms.uint32(2),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(False)
)


process.candidateJetBProbabilityComputer = cms.ESProducer("CandidateJetBProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    numberOfBTracks = cms.uint32(4),
    trackIpSign = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.candidateJetProbabilityComputer = cms.ESProducer("CandidateJetProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(0.3),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    trackIpSign = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.candidateNegativeCombinedMVAV2Computer = cms.ESProducer("CombinedMVAV2JetTagESProducer",
    gbrForestLabel = cms.string('btag_CombinedMVAv2_BDT'),
    jetTagComputers = cms.vstring('candidateNegativeOnlyJetProbabilityComputer', 
        'candidateNegativeOnlyJetBProbabilityComputer', 
        'candidateNegativeCombinedSecondaryVertexV2Computer', 
        'negativeSoftPFMuonComputer', 
        'negativeSoftPFElectronComputer'),
    mvaName = cms.string('bdt'),
    spectators = cms.vstring(),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    variables = cms.vstring('Jet_CSV', 
        'Jet_CSVIVF', 
        'Jet_JP', 
        'Jet_JBP', 
        'Jet_SoftMu', 
        'Jet_SoftEl'),
    weightFile = cms.FileInPath('RecoBTag/Combined/data/CombinedMVAV2_13_07_2015.weights.xml.gz')
)


process.candidateNegativeCombinedSecondaryVertexV2Computer = cms.ESProducer("CandidateCombinedSecondaryVertexESProducer",
    SoftLeptonFlip = cms.bool(False),
    calibrationRecords = cms.vstring('CombinedSVIVFV2RecoVertex', 
        'CombinedSVIVFV2PseudoVertex', 
        'CombinedSVIVFV2NoVertex'),
    categoryVariableName = cms.string('vertexCategory'),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    recordLabel = cms.string(''),
    trackFlip = cms.bool(True),
    trackMultiplicityMin = cms.uint32(2),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(-2.0),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(True)
)


process.candidateNegativeOnlyJetBProbabilityComputer = cms.ESProducer("CandidateJetBProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    numberOfBTracks = cms.uint32(4),
    trackIpSign = cms.int32(-1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.candidateNegativeOnlyJetProbabilityComputer = cms.ESProducer("CandidateJetProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(0.3),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    trackIpSign = cms.int32(-1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.candidateNegativeTrackCounting3D2ndComputer = cms.ESProducer("CandidateNegativeTrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(-1),
    nthTrack = cms.int32(2),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True),
    useVariableJTA = cms.bool(False)
)


process.candidateNegativeTrackCounting3D3rdComputer = cms.ESProducer("CandidateNegativeTrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(-1),
    nthTrack = cms.int32(3),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True),
    useVariableJTA = cms.bool(False)
)


process.candidatePositiveCombinedMVAV2Computer = cms.ESProducer("CombinedMVAV2JetTagESProducer",
    gbrForestLabel = cms.string('btag_CombinedMVAv2_BDT'),
    jetTagComputers = cms.vstring('candidatePositiveOnlyJetProbabilityComputer', 
        'candidatePositiveOnlyJetBProbabilityComputer', 
        'candidatePositiveCombinedSecondaryVertexV2Computer', 
        'negativeSoftPFMuonComputer', 
        'negativeSoftPFElectronComputer'),
    mvaName = cms.string('bdt'),
    spectators = cms.vstring(),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    variables = cms.vstring('Jet_CSV', 
        'Jet_CSVIVF', 
        'Jet_JP', 
        'Jet_JBP', 
        'Jet_SoftMu', 
        'Jet_SoftEl'),
    weightFile = cms.FileInPath('RecoBTag/Combined/data/CombinedMVAV2_13_07_2015.weights.xml.gz')
)


process.candidatePositiveCombinedSecondaryVertexV2Computer = cms.ESProducer("CandidateCombinedSecondaryVertexESProducer",
    SoftLeptonFlip = cms.bool(False),
    calibrationRecords = cms.vstring('CombinedSVIVFV2RecoVertex', 
        'CombinedSVIVFV2PseudoVertex', 
        'CombinedSVIVFV2NoVertex'),
    categoryVariableName = cms.string('vertexCategory'),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    recordLabel = cms.string(''),
    trackFlip = cms.bool(False),
    trackMultiplicityMin = cms.uint32(2),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(0),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(0),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(False)
)


process.candidatePositiveOnlyJetBProbabilityComputer = cms.ESProducer("CandidateJetBProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    numberOfBTracks = cms.uint32(4),
    trackIpSign = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.candidatePositiveOnlyJetProbabilityComputer = cms.ESProducer("CandidateJetProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(0.3),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    trackIpSign = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.candidateSimpleSecondaryVertex2TrkComputer = cms.ESProducer("CandidateSimpleSecondaryVertexESProducer",
    minTracks = cms.uint32(2),
    unBoost = cms.bool(False),
    use3d = cms.bool(True),
    useSignificance = cms.bool(True)
)


process.candidateSimpleSecondaryVertex3TrkComputer = cms.ESProducer("CandidateSimpleSecondaryVertexESProducer",
    minTracks = cms.uint32(3),
    unBoost = cms.bool(False),
    use3d = cms.bool(True),
    useSignificance = cms.bool(True)
)


process.candidateTrackCounting3D2ndComputer = cms.ESProducer("CandidateTrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(-1),
    nthTrack = cms.int32(2),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True),
    useVariableJTA = cms.bool(False)
)


process.candidateTrackCounting3D3rdComputer = cms.ESProducer("CandidateTrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(-1),
    nthTrack = cms.int32(3),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True),
    useVariableJTA = cms.bool(False)
)


process.charmTagsComputerCvsB = cms.ESProducer("CharmTaggerESProducer",
    computer = cms.ESInputTag("combinedSecondaryVertexSoftLeptonComputer"),
    gbrForestLabel = cms.string(''),
    mvaName = cms.string('BDT'),
    slComputerCfg = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        calibrationRecords = cms.vstring('CombinedSVRecoVertexNoSoftLepton', 
            'CombinedSVPseudoVertexNoSoftLepton', 
            'CombinedSVNoVertexNoSoftLepton', 
            'CombinedSVRecoVertexSoftMuon', 
            'CombinedSVPseudoVertexSoftMuon', 
            'CombinedSVNoVertexSoftMuon', 
            'CombinedSVRecoVertexSoftElectron', 
            'CombinedSVPseudoVertexSoftElectron', 
            'CombinedSVNoVertexSoftElectron'),
        categoryVariableName = cms.string('vertexLeptonCategory'),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(False),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        recordLabel = cms.string(''),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useCategories = cms.bool(True),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    ),
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfos"), cms.InputTag("pfInclusiveSecondaryVertexFinderCvsLTagInfos"), cms.InputTag("softPFMuonsTagInfos"), cms.InputTag("softPFElectronsTagInfos")),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    variables = cms.VPSet(cms.PSet(
        default = cms.double(-1),
        name = cms.string('vertexLeptonCategory'),
        taggingVarName = cms.string('vertexLeptonCategory')
    ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSig_0'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip2dSig_1'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSig_0'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip3dSig_1'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPtRel_0'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPtRel_1'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPPar_0'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPPar_1'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackEtaRel_0'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackEtaRel_1'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDeltaR_0'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDeltaR_1'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackPtRatio_0'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackPtRatio_1'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(0),
            name = cms.string('trackPParRatio_0'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(1),
            name = cms.string('trackPParRatio_1'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackJetDist_0'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackJetDist_1'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDecayLenVal_0'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDecayLenVal_1'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(0),
            name = cms.string('jetNSecondaryVertices'),
            taggingVarName = cms.string('jetNSecondaryVertices')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('jetNTracks'),
            taggingVarName = cms.string('jetNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetEtRatio'),
            taggingVarName = cms.string('trackSumJetEtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetDeltaR'),
            taggingVarName = cms.string('trackSumJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexMass_0'),
            taggingVarName = cms.string('vertexMass')
        ), 
        cms.PSet(
            default = cms.double(-10),
            idx = cms.int32(0),
            name = cms.string('vertexEnergyRatio_0'),
            taggingVarName = cms.string('vertexEnergyRatio')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip2dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip3dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance2dSig_0'),
            taggingVarName = cms.string('flightDistance2dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance3dSig_0'),
            taggingVarName = cms.string('flightDistance3dSig')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexJetDeltaR_0'),
            taggingVarName = cms.string('vertexJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(0),
            idx = cms.int32(0),
            name = cms.string('vertexNTracks_0'),
            taggingVarName = cms.string('vertexNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('massVertexEnergyFraction_0'),
            taggingVarName = cms.string('massVertexEnergyFraction')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexBoostOverSqrtJetPt_0'),
            taggingVarName = cms.string('vertexBoostOverSqrtJetPt')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonPtRel_0'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonPtRel_1'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(0),
            name = cms.string('leptonSip3d_0'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(1),
            name = cms.string('leptonSip3d_1'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonDeltaR_0'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonDeltaR_1'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatioRel_0'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatioRel_1'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonEtaRel_0'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonEtaRel_1'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatio_0'),
            taggingVarName = cms.string('leptonRatio')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatio_1'),
            taggingVarName = cms.string('leptonRatio')
        )),
    weightFile = cms.FileInPath('RecoBTag/CTagging/data/c_vs_b_sklearn.weight.xml')
)


process.charmTagsComputerCvsL = cms.ESProducer("CharmTaggerESProducer",
    computer = cms.ESInputTag("combinedSecondaryVertexSoftLeptonComputer"),
    gbrForestLabel = cms.string(''),
    mvaName = cms.string('BDT'),
    slComputerCfg = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        calibrationRecords = cms.vstring('CombinedSVRecoVertexNoSoftLepton', 
            'CombinedSVPseudoVertexNoSoftLepton', 
            'CombinedSVNoVertexNoSoftLepton', 
            'CombinedSVRecoVertexSoftMuon', 
            'CombinedSVPseudoVertexSoftMuon', 
            'CombinedSVNoVertexSoftMuon', 
            'CombinedSVRecoVertexSoftElectron', 
            'CombinedSVPseudoVertexSoftElectron', 
            'CombinedSVNoVertexSoftElectron'),
        categoryVariableName = cms.string('vertexLeptonCategory'),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(False),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        recordLabel = cms.string(''),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useCategories = cms.bool(True),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    ),
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfos"), cms.InputTag("pfInclusiveSecondaryVertexFinderCvsLTagInfos"), cms.InputTag("softPFMuonsTagInfos"), cms.InputTag("softPFElectronsTagInfos")),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    variables = cms.VPSet(cms.PSet(
        default = cms.double(-1),
        name = cms.string('vertexLeptonCategory'),
        taggingVarName = cms.string('vertexLeptonCategory')
    ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSig_0'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip2dSig_1'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSig_0'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip3dSig_1'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPtRel_0'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPtRel_1'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPPar_0'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPPar_1'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackEtaRel_0'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackEtaRel_1'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDeltaR_0'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDeltaR_1'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackPtRatio_0'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackPtRatio_1'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(0),
            name = cms.string('trackPParRatio_0'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(1),
            name = cms.string('trackPParRatio_1'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackJetDist_0'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackJetDist_1'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDecayLenVal_0'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDecayLenVal_1'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(0),
            name = cms.string('jetNSecondaryVertices'),
            taggingVarName = cms.string('jetNSecondaryVertices')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('jetNTracks'),
            taggingVarName = cms.string('jetNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetEtRatio'),
            taggingVarName = cms.string('trackSumJetEtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetDeltaR'),
            taggingVarName = cms.string('trackSumJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexMass_0'),
            taggingVarName = cms.string('vertexMass')
        ), 
        cms.PSet(
            default = cms.double(-10),
            idx = cms.int32(0),
            name = cms.string('vertexEnergyRatio_0'),
            taggingVarName = cms.string('vertexEnergyRatio')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip2dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip3dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance2dSig_0'),
            taggingVarName = cms.string('flightDistance2dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance3dSig_0'),
            taggingVarName = cms.string('flightDistance3dSig')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexJetDeltaR_0'),
            taggingVarName = cms.string('vertexJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(0),
            idx = cms.int32(0),
            name = cms.string('vertexNTracks_0'),
            taggingVarName = cms.string('vertexNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('massVertexEnergyFraction_0'),
            taggingVarName = cms.string('massVertexEnergyFraction')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexBoostOverSqrtJetPt_0'),
            taggingVarName = cms.string('vertexBoostOverSqrtJetPt')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonPtRel_0'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonPtRel_1'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(0),
            name = cms.string('leptonSip3d_0'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(1),
            name = cms.string('leptonSip3d_1'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonDeltaR_0'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonDeltaR_1'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatioRel_0'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatioRel_1'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonEtaRel_0'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonEtaRel_1'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatio_0'),
            taggingVarName = cms.string('leptonRatio')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatio_1'),
            taggingVarName = cms.string('leptonRatio')
        )),
    weightFile = cms.FileInPath('RecoBTag/CTagging/data/c_vs_udsg_sklearn.weight.xml')
)


process.charmTagsNegativeComputerCvsB = cms.ESProducer("CharmTaggerESProducer",
    computer = cms.ESInputTag("combinedSecondaryVertexSoftLeptonComputer"),
    gbrForestLabel = cms.string(''),
    mvaName = cms.string('BDT'),
    slComputerCfg = cms.PSet(
        SoftLeptonFlip = cms.bool(True),
        calibrationRecords = cms.vstring('CombinedSVRecoVertexNoSoftLepton', 
            'CombinedSVPseudoVertexNoSoftLepton', 
            'CombinedSVNoVertexNoSoftLepton', 
            'CombinedSVRecoVertexSoftMuon', 
            'CombinedSVPseudoVertexSoftMuon', 
            'CombinedSVNoVertexSoftMuon', 
            'CombinedSVRecoVertexSoftElectron', 
            'CombinedSVPseudoVertexSoftElectron', 
            'CombinedSVNoVertexSoftElectron'),
        categoryVariableName = cms.string('vertexLeptonCategory'),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(False),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        recordLabel = cms.string(''),
        trackFlip = cms.bool(True),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(-2.0),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(0),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(0),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useCategories = cms.bool(True),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(True)
    ),
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfos"), cms.InputTag("pfInclusiveSecondaryVertexFinderCvsLTagInfos"), cms.InputTag("softPFMuonsTagInfos"), cms.InputTag("softPFElectronsTagInfos")),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    variables = cms.VPSet(cms.PSet(
        default = cms.double(-1),
        name = cms.string('vertexLeptonCategory'),
        taggingVarName = cms.string('vertexLeptonCategory')
    ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSig_0'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip2dSig_1'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSig_0'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip3dSig_1'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPtRel_0'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPtRel_1'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPPar_0'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPPar_1'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackEtaRel_0'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackEtaRel_1'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDeltaR_0'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDeltaR_1'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackPtRatio_0'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackPtRatio_1'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(0),
            name = cms.string('trackPParRatio_0'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(1),
            name = cms.string('trackPParRatio_1'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackJetDist_0'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackJetDist_1'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDecayLenVal_0'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDecayLenVal_1'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(0),
            name = cms.string('jetNSecondaryVertices'),
            taggingVarName = cms.string('jetNSecondaryVertices')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('jetNTracks'),
            taggingVarName = cms.string('jetNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetEtRatio'),
            taggingVarName = cms.string('trackSumJetEtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetDeltaR'),
            taggingVarName = cms.string('trackSumJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexMass_0'),
            taggingVarName = cms.string('vertexMass')
        ), 
        cms.PSet(
            default = cms.double(-10),
            idx = cms.int32(0),
            name = cms.string('vertexEnergyRatio_0'),
            taggingVarName = cms.string('vertexEnergyRatio')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip2dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip3dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance2dSig_0'),
            taggingVarName = cms.string('flightDistance2dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance3dSig_0'),
            taggingVarName = cms.string('flightDistance3dSig')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexJetDeltaR_0'),
            taggingVarName = cms.string('vertexJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(0),
            idx = cms.int32(0),
            name = cms.string('vertexNTracks_0'),
            taggingVarName = cms.string('vertexNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('massVertexEnergyFraction_0'),
            taggingVarName = cms.string('massVertexEnergyFraction')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexBoostOverSqrtJetPt_0'),
            taggingVarName = cms.string('vertexBoostOverSqrtJetPt')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonPtRel_0'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonPtRel_1'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(0),
            name = cms.string('leptonSip3d_0'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(1),
            name = cms.string('leptonSip3d_1'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonDeltaR_0'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonDeltaR_1'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatioRel_0'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatioRel_1'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonEtaRel_0'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonEtaRel_1'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatio_0'),
            taggingVarName = cms.string('leptonRatio')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatio_1'),
            taggingVarName = cms.string('leptonRatio')
        )),
    weightFile = cms.FileInPath('RecoBTag/CTagging/data/c_vs_b_sklearn.weight.xml')
)


process.charmTagsNegativeComputerCvsL = cms.ESProducer("CharmTaggerESProducer",
    computer = cms.ESInputTag("combinedSecondaryVertexSoftLeptonComputer"),
    gbrForestLabel = cms.string(''),
    mvaName = cms.string('BDT'),
    slComputerCfg = cms.PSet(
        SoftLeptonFlip = cms.bool(True),
        calibrationRecords = cms.vstring('CombinedSVRecoVertexNoSoftLepton', 
            'CombinedSVPseudoVertexNoSoftLepton', 
            'CombinedSVNoVertexNoSoftLepton', 
            'CombinedSVRecoVertexSoftMuon', 
            'CombinedSVPseudoVertexSoftMuon', 
            'CombinedSVNoVertexSoftMuon', 
            'CombinedSVRecoVertexSoftElectron', 
            'CombinedSVPseudoVertexSoftElectron', 
            'CombinedSVNoVertexSoftElectron'),
        categoryVariableName = cms.string('vertexLeptonCategory'),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(False),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        recordLabel = cms.string(''),
        trackFlip = cms.bool(True),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(-2.0),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(0),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(0),
            sip3dSigMin = cms.double(-99999.9),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useCategories = cms.bool(True),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(True)
    ),
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfos"), cms.InputTag("pfInclusiveSecondaryVertexFinderCvsLTagInfos"), cms.InputTag("softPFMuonsTagInfos"), cms.InputTag("softPFElectronsTagInfos")),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    variables = cms.VPSet(cms.PSet(
        default = cms.double(-1),
        name = cms.string('vertexLeptonCategory'),
        taggingVarName = cms.string('vertexLeptonCategory')
    ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSig_0'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip2dSig_1'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSig_0'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip3dSig_1'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPtRel_0'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPtRel_1'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPPar_0'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPPar_1'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackEtaRel_0'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackEtaRel_1'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDeltaR_0'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDeltaR_1'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackPtRatio_0'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackPtRatio_1'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(0),
            name = cms.string('trackPParRatio_0'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(1),
            name = cms.string('trackPParRatio_1'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackJetDist_0'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackJetDist_1'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDecayLenVal_0'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDecayLenVal_1'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(0),
            name = cms.string('jetNSecondaryVertices'),
            taggingVarName = cms.string('jetNSecondaryVertices')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('jetNTracks'),
            taggingVarName = cms.string('jetNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetEtRatio'),
            taggingVarName = cms.string('trackSumJetEtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetDeltaR'),
            taggingVarName = cms.string('trackSumJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexMass_0'),
            taggingVarName = cms.string('vertexMass')
        ), 
        cms.PSet(
            default = cms.double(-10),
            idx = cms.int32(0),
            name = cms.string('vertexEnergyRatio_0'),
            taggingVarName = cms.string('vertexEnergyRatio')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip2dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip3dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance2dSig_0'),
            taggingVarName = cms.string('flightDistance2dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance3dSig_0'),
            taggingVarName = cms.string('flightDistance3dSig')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexJetDeltaR_0'),
            taggingVarName = cms.string('vertexJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(0),
            idx = cms.int32(0),
            name = cms.string('vertexNTracks_0'),
            taggingVarName = cms.string('vertexNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('massVertexEnergyFraction_0'),
            taggingVarName = cms.string('massVertexEnergyFraction')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexBoostOverSqrtJetPt_0'),
            taggingVarName = cms.string('vertexBoostOverSqrtJetPt')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonPtRel_0'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonPtRel_1'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(0),
            name = cms.string('leptonSip3d_0'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(1),
            name = cms.string('leptonSip3d_1'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonDeltaR_0'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonDeltaR_1'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatioRel_0'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatioRel_1'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonEtaRel_0'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonEtaRel_1'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatio_0'),
            taggingVarName = cms.string('leptonRatio')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatio_1'),
            taggingVarName = cms.string('leptonRatio')
        )),
    weightFile = cms.FileInPath('RecoBTag/CTagging/data/c_vs_udsg_sklearn.weight.xml')
)


process.charmTagsPositiveComputerCvsB = cms.ESProducer("CharmTaggerESProducer",
    computer = cms.ESInputTag("combinedSecondaryVertexSoftLeptonComputer"),
    gbrForestLabel = cms.string(''),
    mvaName = cms.string('BDT'),
    slComputerCfg = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        calibrationRecords = cms.vstring('CombinedSVRecoVertexNoSoftLepton', 
            'CombinedSVPseudoVertexNoSoftLepton', 
            'CombinedSVNoVertexNoSoftLepton', 
            'CombinedSVRecoVertexSoftMuon', 
            'CombinedSVPseudoVertexSoftMuon', 
            'CombinedSVNoVertexSoftMuon', 
            'CombinedSVRecoVertexSoftElectron', 
            'CombinedSVPseudoVertexSoftElectron', 
            'CombinedSVNoVertexSoftElectron'),
        categoryVariableName = cms.string('vertexLeptonCategory'),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(False),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        recordLabel = cms.string(''),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(0),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(0),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useCategories = cms.bool(True),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    ),
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfos"), cms.InputTag("pfInclusiveSecondaryVertexFinderCvsLTagInfos"), cms.InputTag("softPFMuonsTagInfos"), cms.InputTag("softPFElectronsTagInfos")),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    variables = cms.VPSet(cms.PSet(
        default = cms.double(-1),
        name = cms.string('vertexLeptonCategory'),
        taggingVarName = cms.string('vertexLeptonCategory')
    ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSig_0'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip2dSig_1'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSig_0'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip3dSig_1'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPtRel_0'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPtRel_1'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPPar_0'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPPar_1'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackEtaRel_0'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackEtaRel_1'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDeltaR_0'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDeltaR_1'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackPtRatio_0'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackPtRatio_1'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(0),
            name = cms.string('trackPParRatio_0'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(1),
            name = cms.string('trackPParRatio_1'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackJetDist_0'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackJetDist_1'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDecayLenVal_0'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDecayLenVal_1'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(0),
            name = cms.string('jetNSecondaryVertices'),
            taggingVarName = cms.string('jetNSecondaryVertices')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('jetNTracks'),
            taggingVarName = cms.string('jetNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetEtRatio'),
            taggingVarName = cms.string('trackSumJetEtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetDeltaR'),
            taggingVarName = cms.string('trackSumJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexMass_0'),
            taggingVarName = cms.string('vertexMass')
        ), 
        cms.PSet(
            default = cms.double(-10),
            idx = cms.int32(0),
            name = cms.string('vertexEnergyRatio_0'),
            taggingVarName = cms.string('vertexEnergyRatio')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip2dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip3dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance2dSig_0'),
            taggingVarName = cms.string('flightDistance2dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance3dSig_0'),
            taggingVarName = cms.string('flightDistance3dSig')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexJetDeltaR_0'),
            taggingVarName = cms.string('vertexJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(0),
            idx = cms.int32(0),
            name = cms.string('vertexNTracks_0'),
            taggingVarName = cms.string('vertexNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('massVertexEnergyFraction_0'),
            taggingVarName = cms.string('massVertexEnergyFraction')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexBoostOverSqrtJetPt_0'),
            taggingVarName = cms.string('vertexBoostOverSqrtJetPt')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonPtRel_0'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonPtRel_1'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(0),
            name = cms.string('leptonSip3d_0'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(1),
            name = cms.string('leptonSip3d_1'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonDeltaR_0'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonDeltaR_1'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatioRel_0'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatioRel_1'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonEtaRel_0'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonEtaRel_1'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatio_0'),
            taggingVarName = cms.string('leptonRatio')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatio_1'),
            taggingVarName = cms.string('leptonRatio')
        )),
    weightFile = cms.FileInPath('RecoBTag/CTagging/data/c_vs_b_sklearn.weight.xml')
)


process.charmTagsPositiveComputerCvsL = cms.ESProducer("CharmTaggerESProducer",
    computer = cms.ESInputTag("combinedSecondaryVertexSoftLeptonComputer"),
    gbrForestLabel = cms.string(''),
    mvaName = cms.string('BDT'),
    slComputerCfg = cms.PSet(
        SoftLeptonFlip = cms.bool(False),
        calibrationRecords = cms.vstring('CombinedSVRecoVertexNoSoftLepton', 
            'CombinedSVPseudoVertexNoSoftLepton', 
            'CombinedSVNoVertexNoSoftLepton', 
            'CombinedSVRecoVertexSoftMuon', 
            'CombinedSVPseudoVertexSoftMuon', 
            'CombinedSVNoVertexSoftMuon', 
            'CombinedSVRecoVertexSoftElectron', 
            'CombinedSVPseudoVertexSoftElectron', 
            'CombinedSVNoVertexSoftElectron'),
        categoryVariableName = cms.string('vertexLeptonCategory'),
        charmCut = cms.double(1.5),
        correctVertexMass = cms.bool(False),
        minimumTrackWeight = cms.double(0.5),
        pseudoMultiplicityMin = cms.uint32(2),
        pseudoVertexV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.05)
        ),
        recordLabel = cms.string(''),
        trackFlip = cms.bool(False),
        trackMultiplicityMin = cms.uint32(2),
        trackPairV0Filter = cms.PSet(
            k0sMassWindow = cms.double(0.03)
        ),
        trackPseudoSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(2.0),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(0),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSelection = cms.PSet(
            a_dR = cms.double(-0.001053),
            a_pT = cms.double(0.005263),
            b_dR = cms.double(0.6263),
            b_pT = cms.double(0.3684),
            jetDeltaRMax = cms.double(0.3),
            maxDecayLen = cms.double(5),
            maxDistToAxis = cms.double(0.07),
            max_pT = cms.double(500),
            max_pT_dRcut = cms.double(0.1),
            max_pT_trackPTcut = cms.double(3),
            min_pT = cms.double(120),
            min_pT_dRcut = cms.double(0.5),
            normChi2Max = cms.double(99999.9),
            pixelHitsMin = cms.uint32(0),
            ptMin = cms.double(0.0),
            qualityClass = cms.string('any'),
            sip2dSigMax = cms.double(99999.9),
            sip2dSigMin = cms.double(-99999.9),
            sip2dValMax = cms.double(99999.9),
            sip2dValMin = cms.double(-99999.9),
            sip3dSigMax = cms.double(99999.9),
            sip3dSigMin = cms.double(0),
            sip3dValMax = cms.double(99999.9),
            sip3dValMin = cms.double(-99999.9),
            totalHitsMin = cms.uint32(0),
            useVariableJTA = cms.bool(False)
        ),
        trackSort = cms.string('sip2dSig'),
        useCategories = cms.bool(True),
        useTrackWeights = cms.bool(True),
        vertexFlip = cms.bool(False)
    ),
    tagInfos = cms.VInputTag(cms.InputTag("pfImpactParameterTagInfos"), cms.InputTag("pfInclusiveSecondaryVertexFinderCvsLTagInfos"), cms.InputTag("softPFMuonsTagInfos"), cms.InputTag("softPFElectronsTagInfos")),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    variables = cms.VPSet(cms.PSet(
        default = cms.double(-1),
        name = cms.string('vertexLeptonCategory'),
        taggingVarName = cms.string('vertexLeptonCategory')
    ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSig_0'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip2dSig_1'),
            taggingVarName = cms.string('trackSip2dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSig_0'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-100),
            idx = cms.int32(1),
            name = cms.string('trackSip3dSig_1'),
            taggingVarName = cms.string('trackSip3dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPtRel_0'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPtRel_1'),
            taggingVarName = cms.string('trackPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackPPar_0'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackPPar_1'),
            taggingVarName = cms.string('trackPPar')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('trackEtaRel_0'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('trackEtaRel_1'),
            taggingVarName = cms.string('trackEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDeltaR_0'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDeltaR_1'),
            taggingVarName = cms.string('trackDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackPtRatio_0'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackPtRatio_1'),
            taggingVarName = cms.string('trackPtRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(0),
            name = cms.string('trackPParRatio_0'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(1.1),
            idx = cms.int32(1),
            name = cms.string('trackPParRatio_1'),
            taggingVarName = cms.string('trackPParRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackJetDist_0'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackJetDist_1'),
            taggingVarName = cms.string('trackJetDist')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('trackDecayLenVal_0'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(1),
            name = cms.string('trackDecayLenVal_1'),
            taggingVarName = cms.string('trackDecayLenVal')
        ), 
        cms.PSet(
            default = cms.double(0),
            name = cms.string('jetNSecondaryVertices'),
            taggingVarName = cms.string('jetNSecondaryVertices')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('jetNTracks'),
            taggingVarName = cms.string('jetNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetEtRatio'),
            taggingVarName = cms.string('trackSumJetEtRatio')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            name = cms.string('trackSumJetDeltaR'),
            taggingVarName = cms.string('trackSumJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexMass_0'),
            taggingVarName = cms.string('vertexMass')
        ), 
        cms.PSet(
            default = cms.double(-10),
            idx = cms.int32(0),
            name = cms.string('vertexEnergyRatio_0'),
            taggingVarName = cms.string('vertexEnergyRatio')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip2dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip2dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-999),
            idx = cms.int32(0),
            name = cms.string('trackSip3dSigAboveCharm_0'),
            taggingVarName = cms.string('trackSip3dSigAboveCharm')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance2dSig_0'),
            taggingVarName = cms.string('flightDistance2dSig')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('flightDistance3dSig_0'),
            taggingVarName = cms.string('flightDistance3dSig')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexJetDeltaR_0'),
            taggingVarName = cms.string('vertexJetDeltaR')
        ), 
        cms.PSet(
            default = cms.double(0),
            idx = cms.int32(0),
            name = cms.string('vertexNTracks_0'),
            taggingVarName = cms.string('vertexNTracks')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('massVertexEnergyFraction_0'),
            taggingVarName = cms.string('massVertexEnergyFraction')
        ), 
        cms.PSet(
            default = cms.double(-0.1),
            idx = cms.int32(0),
            name = cms.string('vertexBoostOverSqrtJetPt_0'),
            taggingVarName = cms.string('vertexBoostOverSqrtJetPt')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonPtRel_0'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonPtRel_1'),
            taggingVarName = cms.string('leptonPtRel')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(0),
            name = cms.string('leptonSip3d_0'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-10000),
            idx = cms.int32(1),
            name = cms.string('leptonSip3d_1'),
            taggingVarName = cms.string('leptonSip3d')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonDeltaR_0'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonDeltaR_1'),
            taggingVarName = cms.string('leptonDeltaR')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatioRel_0'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatioRel_1'),
            taggingVarName = cms.string('leptonRatioRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonEtaRel_0'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonEtaRel_1'),
            taggingVarName = cms.string('leptonEtaRel')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(0),
            name = cms.string('leptonRatio_0'),
            taggingVarName = cms.string('leptonRatio')
        ), 
        cms.PSet(
            default = cms.double(-1),
            idx = cms.int32(1),
            name = cms.string('leptonRatio_1'),
            taggingVarName = cms.string('leptonRatio')
        )),
    weightFile = cms.FileInPath('RecoBTag/CTagging/data/c_vs_udsg_sklearn.weight.xml')
)


process.combinedMVAV2Computer = cms.ESProducer("CombinedMVAV2JetTagESProducer",
    jetTagComputers = cms.vstring('jetProbabilityComputer', 
        'jetBProbabilityComputer', 
        'combinedSecondaryVertexV2Computer', 
        'softPFMuonComputer', 
        'softPFElectronComputer'),
    mvaName = cms.string('bdt'),
    spectators = cms.vstring(),
    useAdaBoost = cms.bool(True),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    variables = cms.vstring('Jet_CSV', 
        'Jet_CSVIVF', 
        'Jet_JP', 
        'Jet_JBP', 
        'Jet_SoftMu', 
        'Jet_SoftEl'),
    weightFile = cms.FileInPath('RecoBTag/Combined/data/CombinedMVAV2_13_07_2015.weights.xml.gz')
)


process.combinedSecondaryVertexV2Computer = cms.ESProducer("CombinedSecondaryVertexESProducer",
    SoftLeptonFlip = cms.bool(False),
    calibrationRecords = cms.vstring('CombinedSVIVFV2RecoVertex', 
        'CombinedSVIVFV2PseudoVertex', 
        'CombinedSVIVFV2NoVertex'),
    categoryVariableName = cms.string('vertexCategory'),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    recordLabel = cms.string(''),
    trackFlip = cms.bool(False),
    trackMultiplicityMin = cms.uint32(2),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(False)
)


process.doubleVertex2TrkComputer = cms.ESProducer("SimpleSecondaryVertexESProducer",
    minTracks = cms.uint32(2),
    minVertices = cms.uint32(2),
    unBoost = cms.bool(False),
    use3d = cms.bool(True),
    useSignificance = cms.bool(True)
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.ghostTrackComputer = cms.ESProducer("GhostTrackESProducer",
    calibrationRecords = cms.vstring('GhostTrackRecoVertex', 
        'GhostTrackPseudoVertex', 
        'GhostTrackNoVertex'),
    categoryVariableName = cms.string('vertexCategory'),
    charmCut = cms.double(1.5),
    minimumTrackWeight = cms.double(0.5),
    recordLabel = cms.string(''),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True)
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.ic5CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute')
)


process.ic5CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloL6SLB')
)


process.ic5CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL1Fastjet', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual')
)


process.ic5CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ic5CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute')
)


process.ic5CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual')
)


process.ic5CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ic5CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL2Relative', 
        'ic5CaloL3Absolute')
)


process.ic5CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloL6SLB')
)


process.ic5CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual')
)


process.ic5CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L2Relative')
)


process.ic5CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L3Absolute')
)


process.ic5CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ic5CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ic5CaloJetsSoftMuonTagInfos")
)


process.ic5CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L2L3Residual')
)


process.ic5PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute')
)


process.ic5PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFL6SLB')
)


process.ic5PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual')
)


process.ic5PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ic5PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL1Offset', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute')
)


process.ic5PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL1Offset', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual')
)


process.ic5PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ic5PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL2Relative', 
        'ic5PFL3Absolute')
)


process.ic5PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFL6SLB')
)


process.ic5PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual')
)


process.ic5PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L2Relative')
)


process.ic5PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L3Absolute')
)


process.ic5PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ic5PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ic5PFJetsSoftMuonTagInfos")
)


process.ic5PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L2L3Residual')
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.impactParameterMVAComputer = cms.ESProducer("GenericMVAJetTagESProducer",
    calibrationRecord = cms.string('ImpactParameterMVA'),
    recordLabel = cms.string(''),
    useCategories = cms.bool(False)
)


process.jetBProbabilityComputer = cms.ESProducer("JetBProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    numberOfBTracks = cms.uint32(4),
    trackIpSign = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.jetProbabilityComputer = cms.ESProducer("JetProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(0.3),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    trackIpSign = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.kt4CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute')
)


process.kt4CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloL6SLB')
)


process.kt4CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL1Fastjet', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual')
)


process.kt4CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.kt4CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute')
)


process.kt4CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual')
)


process.kt4CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt4CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL2Relative', 
        'kt4CaloL3Absolute')
)


process.kt4CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloL6SLB')
)


process.kt4CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual')
)


process.kt4CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L2Relative')
)


process.kt4CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L3Absolute')
)


process.kt4CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt4CaloJetsSoftMuonTagInfos")
)


process.kt4CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L2L3Residual')
)


process.kt4PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute')
)


process.kt4PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFL6SLB')
)


process.kt4PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual')
)


process.kt4PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.kt4PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL1Offset', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute')
)


process.kt4PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL1Offset', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual')
)


process.kt4PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt4PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL2Relative', 
        'kt4PFL3Absolute')
)


process.kt4PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFL6SLB')
)


process.kt4PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual')
)


process.kt4PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L2Relative')
)


process.kt4PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L3Absolute')
)


process.kt4PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt4PFJetsSoftMuonTagInfos")
)


process.kt4PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L2L3Residual')
)


process.kt6CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute')
)


process.kt6CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloL6SLB')
)


process.kt6CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL1Fastjet', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual')
)


process.kt6CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.kt6CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute')
)


process.kt6CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual')
)


process.kt6CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt6CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL2Relative', 
        'kt6CaloL3Absolute')
)


process.kt6CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloL6SLB')
)


process.kt6CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual')
)


process.kt6CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L2Relative')
)


process.kt6CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L3Absolute')
)


process.kt6CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt6CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt6CaloJetsSoftMuonTagInfos")
)


process.kt6CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L2L3Residual')
)


process.kt6PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute')
)


process.kt6PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFL6SLB')
)


process.kt6PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual')
)


process.kt6PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.kt6PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL1Offset', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute')
)


process.kt6PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL1Offset', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual')
)


process.kt6PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt6PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL2Relative', 
        'kt6PFL3Absolute')
)


process.kt6PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFL6SLB')
)


process.kt6PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual')
)


process.kt6PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L2Relative')
)


process.kt6PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L3Absolute')
)


process.kt6PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt6PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt6PFJetsSoftMuonTagInfos")
)


process.kt6PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L2L3Residual')
)


process.negativeCombinedMVAV2Computer = cms.ESProducer("CombinedMVAV2JetTagESProducer",
    jetTagComputers = cms.vstring('negativeOnlyJetProbabilityComputer', 
        'negativeOnlyJetBProbabilityComputer', 
        'negativeCombinedSecondaryVertexV2Computer', 
        'negativeSoftPFMuonComputer', 
        'negativeSoftPFElectronComputer'),
    mvaName = cms.string('bdt'),
    spectators = cms.vstring(),
    useAdaBoost = cms.bool(True),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    variables = cms.vstring('Jet_CSV', 
        'Jet_CSVIVF', 
        'Jet_JP', 
        'Jet_JBP', 
        'Jet_SoftMu', 
        'Jet_SoftEl'),
    weightFile = cms.FileInPath('RecoBTag/Combined/data/CombinedMVAV2_13_07_2015.weights.xml.gz')
)


process.negativeCombinedSecondaryVertexV2Computer = cms.ESProducer("CombinedSecondaryVertexESProducer",
    SoftLeptonFlip = cms.bool(False),
    calibrationRecords = cms.vstring('CombinedSVIVFV2RecoVertex', 
        'CombinedSVIVFV2PseudoVertex', 
        'CombinedSVIVFV2NoVertex'),
    categoryVariableName = cms.string('vertexCategory'),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    recordLabel = cms.string(''),
    trackFlip = cms.bool(True),
    trackMultiplicityMin = cms.uint32(2),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(-2.0),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(0),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(True)
)


process.negativeOnlyJetBProbabilityComputer = cms.ESProducer("JetBProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    numberOfBTracks = cms.uint32(4),
    trackIpSign = cms.int32(-1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.negativeOnlyJetProbabilityComputer = cms.ESProducer("JetProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(0.3),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    trackIpSign = cms.int32(-1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.negativeSoftPFElectronByIP2dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('negative'),
    use3d = cms.bool(False)
)


process.negativeSoftPFElectronByIP3dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('negative'),
    use3d = cms.bool(True)
)


process.negativeSoftPFElectronByPtComputer = cms.ESProducer("LeptonTaggerByPtESProducer",
    ipSign = cms.string('negative')
)


process.negativeSoftPFElectronComputer = cms.ESProducer("ElectronTaggerESProducer",
    gbrForestLabel = cms.string('btag_SoftPFElectron_BDT'),
    ipSign = cms.string('negative'),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SoftLepton/data/SoftPFElectron_BDT.weights.xml.gz')
)


process.negativeSoftPFMuonByIP2dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('negative'),
    use3d = cms.bool(False)
)


process.negativeSoftPFMuonByIP3dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('negative'),
    use3d = cms.bool(True)
)


process.negativeSoftPFMuonByPtComputer = cms.ESProducer("LeptonTaggerByPtESProducer",
    ipSign = cms.string('negative')
)


process.negativeSoftPFMuonComputer = cms.ESProducer("MuonTaggerESProducer",
    gbrForestLabel = cms.string('btag_SoftPFMuon_BDT'),
    ipSign = cms.string('negative'),
    useAdaBoost = cms.bool(True),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SoftLepton/data/SoftPFMuon_BDT.weights.xml.gz')
)


process.negativeTrackCounting3D2ndComputer = cms.ESProducer("NegativeTrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(-1),
    nthTrack = cms.int32(2),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True),
    useVariableJTA = cms.bool(False)
)


process.negativeTrackCounting3D3rdComputer = cms.ESProducer("NegativeTrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(-1),
    nthTrack = cms.int32(3),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True),
    useVariableJTA = cms.bool(False)
)


process.positiveCombinedMVAV2Computer = cms.ESProducer("CombinedMVAV2JetTagESProducer",
    jetTagComputers = cms.vstring('positiveOnlyJetProbabilityComputer', 
        'positiveOnlyJetBProbabilityComputer', 
        'positiveCombinedSecondaryVertexV2Computer', 
        'positiveSoftPFMuonComputer', 
        'positiveSoftPFElectronComputer'),
    mvaName = cms.string('bdt'),
    spectators = cms.vstring(),
    useAdaBoost = cms.bool(True),
    useCondDB = cms.bool(False),
    useGBRForest = cms.bool(True),
    variables = cms.vstring('Jet_CSV', 
        'Jet_CSVIVF', 
        'Jet_JP', 
        'Jet_JBP', 
        'Jet_SoftMu', 
        'Jet_SoftEl'),
    weightFile = cms.FileInPath('RecoBTag/Combined/data/CombinedMVAV2_13_07_2015.weights.xml.gz')
)


process.positiveCombinedSecondaryVertexV2Computer = cms.ESProducer("CombinedSecondaryVertexESProducer",
    SoftLeptonFlip = cms.bool(False),
    calibrationRecords = cms.vstring('CombinedSVIVFV2RecoVertex', 
        'CombinedSVIVFV2PseudoVertex', 
        'CombinedSVIVFV2NoVertex'),
    categoryVariableName = cms.string('vertexCategory'),
    charmCut = cms.double(1.5),
    correctVertexMass = cms.bool(True),
    minimumTrackWeight = cms.double(0.5),
    pseudoMultiplicityMin = cms.uint32(2),
    pseudoVertexV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.05)
    ),
    recordLabel = cms.string(''),
    trackFlip = cms.bool(False),
    trackMultiplicityMin = cms.uint32(2),
    trackPairV0Filter = cms.PSet(
        k0sMassWindow = cms.double(0.03)
    ),
    trackPseudoSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(2.0),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(0),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
        maxDistToAxis = cms.double(0.07),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(0),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    ),
    trackSort = cms.string('sip2dSig'),
    useCategories = cms.bool(True),
    useTrackWeights = cms.bool(True),
    vertexFlip = cms.bool(False)
)


process.positiveOnlyJetBProbabilityComputer = cms.ESProducer("JetBProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    numberOfBTracks = cms.uint32(4),
    trackIpSign = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.positiveOnlyJetProbabilityComputer = cms.ESProducer("JetProbabilityESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(0.3),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumProbability = cms.double(0.005),
    trackIpSign = cms.int32(1),
    trackQualityClass = cms.string('any'),
    useVariableJTA = cms.bool(False)
)


process.positiveSoftPFElectronByIP2dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('positive'),
    use3d = cms.bool(False)
)


process.positiveSoftPFElectronByIP3dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('positive'),
    use3d = cms.bool(True)
)


process.positiveSoftPFElectronByPtComputer = cms.ESProducer("LeptonTaggerByPtESProducer",
    ipSign = cms.string('positive')
)


process.positiveSoftPFElectronComputer = cms.ESProducer("ElectronTaggerESProducer",
    gbrForestLabel = cms.string('btag_SoftPFElectron_BDT'),
    ipSign = cms.string('positive'),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SoftLepton/data/SoftPFElectron_BDT.weights.xml.gz')
)


process.positiveSoftPFMuonByIP2dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('positive'),
    use3d = cms.bool(False)
)


process.positiveSoftPFMuonByIP3dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('positive'),
    use3d = cms.bool(True)
)


process.positiveSoftPFMuonByPtComputer = cms.ESProducer("LeptonTaggerByPtESProducer",
    ipSign = cms.string('positive')
)


process.positiveSoftPFMuonComputer = cms.ESProducer("MuonTaggerESProducer",
    gbrForestLabel = cms.string('btag_SoftPFMuon_BDT'),
    ipSign = cms.string('positive'),
    useAdaBoost = cms.bool(True),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SoftLepton/data/SoftPFMuon_BDT.weights.xml.gz')
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiPixelQualityFromDbRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        ))
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(cms.PSet(
        Label = cms.untracked.string(''),
        NormalizationFactor = cms.untracked.double(1.0),
        Record = cms.string('SiStripApvGainRcd')
    ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiStripDetVOffRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.simpleSecondaryVertex2TrkComputer = cms.ESProducer("SimpleSecondaryVertexESProducer",
    minTracks = cms.uint32(2),
    unBoost = cms.bool(False),
    use3d = cms.bool(True),
    useSignificance = cms.bool(True)
)


process.simpleSecondaryVertex3TrkComputer = cms.ESProducer("SimpleSecondaryVertexESProducer",
    minTracks = cms.uint32(3),
    unBoost = cms.bool(False),
    use3d = cms.bool(True),
    useSignificance = cms.bool(True)
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.softPFElectronByIP2dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('any'),
    use3d = cms.bool(False)
)


process.softPFElectronByIP3dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('any'),
    use3d = cms.bool(True)
)


process.softPFElectronByPtComputer = cms.ESProducer("LeptonTaggerByPtESProducer",
    ipSign = cms.string('any')
)


process.softPFElectronComputer = cms.ESProducer("ElectronTaggerESProducer",
    gbrForestLabel = cms.string('btag_SoftPFElectron_BDT'),
    ipSign = cms.string('any'),
    useAdaBoost = cms.bool(False),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SoftLepton/data/SoftPFElectron_BDT.weights.xml.gz')
)


process.softPFMuonByIP2dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('any'),
    use3d = cms.bool(False)
)


process.softPFMuonByIP3dComputer = cms.ESProducer("LeptonTaggerByIPESProducer",
    ipSign = cms.string('any'),
    use3d = cms.bool(True)
)


process.softPFMuonByPtComputer = cms.ESProducer("LeptonTaggerByPtESProducer",
    ipSign = cms.string('any')
)


process.softPFMuonComputer = cms.ESProducer("MuonTaggerESProducer",
    gbrForestLabel = cms.string('btag_SoftPFMuon_BDT'),
    ipSign = cms.string('any'),
    useAdaBoost = cms.bool(True),
    useCondDB = cms.bool(True),
    useGBRForest = cms.bool(True),
    weightFile = cms.FileInPath('RecoBTag/SoftLepton/data/SoftPFMuon_BDT.weights.xml.gz')
)


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.trackCounting3D2ndComputer = cms.ESProducer("TrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(-1),
    nthTrack = cms.int32(2),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True),
    useVariableJTA = cms.bool(False)
)


process.trackCounting3D3rdComputer = cms.ESProducer("TrackCountingESProducer",
    a_dR = cms.double(-0.001053),
    a_pT = cms.double(0.005263),
    b_dR = cms.double(0.6263),
    b_pT = cms.double(0.3684),
    deltaR = cms.double(-1.0),
    impactParameterType = cms.int32(0),
    max_pT = cms.double(500),
    max_pT_dRcut = cms.double(0.1),
    max_pT_trackPTcut = cms.double(3),
    maximumDecayLength = cms.double(5.0),
    maximumDistanceToJetAxis = cms.double(0.07),
    min_pT = cms.double(120),
    min_pT_dRcut = cms.double(0.5),
    minimumImpactParameter = cms.double(-1),
    nthTrack = cms.int32(3),
    trackQualityClass = cms.string('any'),
    useSignedImpactParameterSig = cms.bool(True),
    useVariableJTA = cms.bool(False)
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.BTagRecord = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('JetTagComputerRecord')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('80X_mcRun2_asymptotic_2016_TrancheIV_v8'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet(cms.PSet(
        connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
        label = cms.untracked.string('electron_eb_ECALonly'),
        record = cms.string('GBRDWrapperRcd'),
        tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
    ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ))
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HERecalibration = cms.bool(False),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalibration = cms.bool(False),
    iLumi = cms.double(-1.0),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths')
)


process.jec = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('sqlite:Summer16_23Sep2016V4_MC.db'),
    toGet = cms.VPSet(cms.PSet(
        label = cms.untracked.string('AK4PFchs'),
        record = cms.string('JetCorrectionsRecord'),
        tag = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs')
    ), 
        cms.PSet(
            label = cms.untracked.string('AK8PFchs'),
            record = cms.string('JetCorrectionsRecord'),
            tag = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK8PFchs')
        ))
)


process.prefer("es_hardcode")

process.prefer("jec")

