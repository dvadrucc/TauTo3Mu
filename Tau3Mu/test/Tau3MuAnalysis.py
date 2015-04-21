import FWCore.ParameterSet.Config as cms

process = cms.Process('GenLevel1')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.MessageLogger = cms.Service("MessageLogger",
     cout = cms.untracked.PSet(
         threshold = cms.untracked.string('WARNING')
     ),
     destinations = cms.untracked.vstring('cout')
)

# Source
process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring(

#"file:/afs/cern.ch/work/d/dvadrucc/private/CMSSW_5_3_26/src/MyAnalysis/Tau3Mu/test/MinBias_TuneZ2_DsTau3Mu_cff_py_GEN_FASTSIM_HLT_PU.root"

 #   "file:TauDisplay3Tight.root"
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/752/F4AAA3E8-6CD9-E111-B00F-BCAEC53296F8.root",
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/752/F059F588-66D9-E111-825D-003048F1C832.root",
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/752/ECD1BC26-65D9-E111-895A-BCAEC5364C62.root",
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/752/EC7B95A0-6DD9-E111-9FEF-E0CB4E4408E3.root",
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/752/EA7D6EB8-6DD9-E111-98B8-003048D37524.root",
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/752/E6CFA1C1-6FD9-E111-9C95-BCAEC53296F4.root",
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/752/E0150FFC-73D9-E111-BFEC-003048D2BB90.root",
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/752/DA3D63DA-6AD9-E111-97F6-E0CB4E55365C.root",
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/752/C88DF7DA-6AD9-E111-A6A0-BCAEC518FF89.root",
#

#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/752/7C7CC5EF-7BD9-E111-A1EF-00237DDC5C24.root",

    
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/754/FAA3B6A2-B4D9-E111-9D07-001D09F24763.root",
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/754/F8AEE309-9FD9-E111-A44F-001D09F2906A.root",
#"/store/data/Run2012C/MuOnia/AOD/PromptReco-v2/000/199/754/F877A5BF-AED9-E111-8640-003048F117EA.root",    
"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_100_1_3Hy.root",
"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_67_1_tpw.root",
"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_70_1_Fxb.root",
"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_71_1_ZTD.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_72_1_W0T.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_73_1_Z91.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_74_1_yXY.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_75_1_SJq.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_76_1_Yx8.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_77_1_zec.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_78_1_RYY.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_80_1_uAL.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_81_1_cqF.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_82_1_Gt2.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_90_1_MDw.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_91_1_F2k.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_92_1_yYJ.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_93_1_3Zb.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_94_1_AME.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_95_1_QAs.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_96_1_Da5.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_97_1_50a.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_98_1_k4s.root",
#"/store/caf/user/cerminar/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/MinBias_TuneZ2_DsPhiTo2MuPhi_v2/c2a4db5a6695d8b884ae28dbddc41704/MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff_py_GEN_FASTSIM_HLT_PU_99_1_Npc.root"

#    
#
)
)

process.GlobalTag.globaltag ='START53_V7A::All'
#
#"START52_V9::All"
#'START53_V7::All'
#'START53_V7A::All'

#'GR_P_V40_AN1::All'
#"START52_V9::All"
#"START52_V11::All"

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.ana = cms.EDAnalyzer('Tau3MuAnalysis_V2',

OutFileName=cms.string("OUT_DATA.root"),

Debug=cms.bool(False),

IsMC=cms.bool(True),     #Be sure you run on MC otherwise you get a crash!

isSignal=cms.bool(False), #If true the gen matching is done on Signal else on Norm. sample
isBackground=cms.bool(False),

DiMuMassMin= cms.double(0.6),  # 
DiMuMassMax= cms.double(1.7),#

DiMuLxyMin= cms.double(-100),
DiMuLxyMax= cms.double(500),
DiMuLxySigMin= cms.double(1.),

DiMuVtxChi2Max= cms.double(20),
DiMuVprobMin=cms.double(0.1),

DiMuCosPointMin=cms.double(0.8),
DiMuTrackCosPointMin=cms.double(0.98),

GuessForTrackMass=cms.double(0.1), #Guess for the mass of the track | 0.1396 pion | 0.1057 muon |

DiMuTrackMassMin= cms.double(1.65),
DiMuTrackMassMax= cms.double(2.1),

DiMuTrackLxyMin= cms.double(-100),
DiMuTrackLxyMax= cms.double(700),
DiMuTrackLxySigMin= cms.double(1.), 

DiMuTrackVtxChi2Max= cms.double(20),
DiMuTrackVprobMin=cms.double(0.1),

MuPTCut=cms.double(1.0),
TrackPTCut=cms.double(.5),

HLT_paths = cms.vstring( # noting means passtrough
"HLT_Tau2Mu_ItTrack"
),

HLT_process = cms.string("HLT")

)

process.analysisPath = cms.Path(process.ana)
