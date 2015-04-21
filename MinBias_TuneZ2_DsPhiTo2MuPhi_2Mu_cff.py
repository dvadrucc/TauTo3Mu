import FWCore.ParameterSet.Config as cms

source = cms.Source("EmptySource")

generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(7000.0),
    crossSection = cms.untracked.double(71260000000.0),
    maxEventsToPrint = cms.untracked.int32(0),
   
    PythiaParameters = cms.PSet(
        pythiaUESettings = cms.vstring('MSTU(21)=1     ! Check on possible errors during program execution', 
            'MSTJ(22)=2     ! Decay those unstable particles', 
            'PARJ(71)=10 .  ! for which ctau  10 mm', 
            'MSTP(33)=0     ! no K factors in hard cross sections', 
            'MSTP(2)=1      ! which order running alphaS', 
            'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
            'MSTP(52)=2     ! work with LHAPDF', 
            'PARP(82)=1.832 ! pt cutoff for multiparton interactions', 
            'PARP(89)=1800. ! sqrts for which PARP82 is set', 
            'PARP(90)=0.275 ! Multiple interactions: rescaling power', 
            'MSTP(95)=6     ! CR (color reconnection parameters)', 
            'PARP(77)=1.016 ! CR', 
            'PARP(78)=0.538 ! CR', 
            'PARP(80)=0.1   ! Prob. colored parton from BBR', 
            'PARP(83)=0.356 ! Multiple interactions: matter distribution parameter', 
            'PARP(84)=0.651 ! Multiple interactions: matter distribution parameter', 
            'PARP(62)=1.025 ! ISR cutoff', 
            'MSTP(91)=1     ! Gaussian primordial kT', 
            'PARP(93)=10.0  ! primordial kT-max', 
            'MSTP(81)=21    ! multiple parton interactions 1 is Pythia default', 
            'MSTP(82)=4     ! Defines the multi-parton model'),
        processParameters = cms.vstring('MSEL=1'
	),
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    ),
	 ExternalDecays = cms.PSet(
        EvtGen = cms.untracked.PSet(
          operates_on_particles = cms.vint32(0), # 0=all
          use_default_decay = cms.untracked.bool(False),
          decay_table = cms.FileInPath('GeneratorInterface/ExternalDecays/data/DECAY_NOLONGLIFE.DEC'),
          particle_property_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/evt.pdl'),
          #user_decay_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/Ds_phipi_mumupi.dec'),
          user_decay_file = cms.FileInPath('MyAnalysis/Tau3Mu/data/Ds_phipi_mumupi.dec'),
          list_forced_decays = cms.vstring('MyD_s+','MyD_s-','MyPhi')
          ),
        parameterSets = cms.vstring('EvtGen')
        )			 
			      
)

# filter to select events with a Ds
DsFilter = cms.EDFilter("PythiaFilter",
       Status = cms.untracked.int32(2),
       MaxEta = cms.untracked.double(3),
       MinEta = cms.untracked.double(-3),
       MinPt = cms.untracked.double(5),
       ParticleID = cms.untracked.int32(431)  #D_s 
     
   )

#optionally, insert next filter to check all events have phi from Ds
phifromDfilter = cms.EDFilter("PythiaFilter",
       Status = cms.untracked.int32(2),
       MaxEta = cms.untracked.double(1000.0),
       MinEta = cms.untracked.double(-1000.0),
       MinPt = cms.untracked.double(0.0),
       ParticleID = cms.untracked.int32(333),  #phi  
       MotherID = cms.untracked.int32(431)    #D_s
   )

# mu from phi: filter used but basically no acceptance cuts
mufromphifilter = cms.EDFilter("PythiaFilter",
       Status = cms.untracked.int32(1),
       MaxEta = cms.untracked.double(1000),
       MinEta = cms.untracked.double(-1000),
       MinPt = cms.untracked.double(0.0),
       ParticleID = cms.untracked.int32(13),  #mu  
       MotherID = cms.untracked.int32(333)    #phi
   )

genParticlesForFilter = cms.EDProducer("GenParticleProducer",
				       saveBarCodes = cms.untracked.bool(True),
				       src = cms.InputTag("generator"),
				       abortOnUnknownPDGCode = cms.untracked.bool(False)
)


muonParticlesInAcc = cms.EDFilter("GenParticleSelector",
				  filter = cms.bool(False),
				  src = cms.InputTag("genParticlesForFilter"),
				  cut = cms.string('pt > 1. && abs(pdgId) == 13 && abs(eta) < 2.4'),
				  stableOnly = cms.bool(True)
				  )


twoMuonFilter = cms.EDFilter("CandViewCountFilter",
			       src = cms.InputTag("muonParticlesInAcc"),
			       minNumber = cms.uint32(2))

# Production Info
configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.3 $'),
    annotation = cms.untracked.string('MinBias_TuneZ2_DsPhiTo2MuPhi_2Mu_cff.py nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
    )



ProductionFilterSequence = cms.Sequence(generator *
					DsFilter *
					phifromDfilter *
					mufromphifilter *
					genParticlesForFilter *
					muonParticlesInAcc *
					twoMuonFilter)
