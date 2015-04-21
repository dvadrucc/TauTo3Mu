#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Candidate/interface/Particle.h"

#include "MyAnalysis/Tau3Mu/src/Histograms.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
//#include "DataFormats/MuonReco/interface/MuonEnergy.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h>

#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "Math/SMatrix.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include <TMatrixD.h>
#include <iostream>
#include <map>
#include <set>
#include <TROOT.h>

//#include "SMatrix.h"

using namespace std;
using namespace reco;
using namespace edm;

// USEFUL DOCUMENTS:
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate#GenPCand
// http://pdg.lbl.gov/mc_particle_id_contents.html


//-----------------------------------------------------------------


//-----------------------------------------------------------------

inline bool sortByPt(const reco::Candidate *part1, const reco::Candidate *part2) {

  return part1->pt() > part2->pt();
}

inline bool sortMuByPt(const reco::Muon mu1, const reco::Muon mu2) {

  return mu1.pt() > mu2.pt();
}

//-----------------------------------------------------------------
class ThreeMuFinder : public edm::EDAnalyzer {
public:
  explicit ThreeMuFinder(const edm::ParameterSet&);
  ~ThreeMuFinder();

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual bool isTight(const reco::Muon*);
  virtual void Initialize_TreeVars();
  virtual bool TriggerDecision(const edm::Event&);

  TFile* thefile;

  TTree* ExTree;

  TLorentzVector* _triMu4Mom, *_m1,*_m2,*_m3;
  int q1,q2,q3;
  int _Run,_Evt,_Lum,_nmuons;
  double Vprob;
  TVector3* SV,*SVe;

  std::vector<string> HLT_paths;
  std::string HLT_process;

  bool _TrigBit[10];
};

//
// constants, enums and typedefs
//

//////////////////////////////////////////////////////////////////
// generically maximum
template <class T> const T& max ( const T& a, const T& b ) {
  return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}

//

ThreeMuFinder::ThreeMuFinder(const edm::ParameterSet& cfg) {

  HLT_paths = cfg.getParameter<std::vector<string> > ("HLT_paths");
  HLT_process = cfg.getParameter<std::string> ("HLT_process");
}

ThreeMuFinder::~ThreeMuFinder() {}

//
// member functions
//

bool ThreeMuFinder::isTight(const Muon* recoMu){

  bool isTight=false;

  if (!recoMu->isGlobalMuon()) return false;
  if (!recoMu->isPFMuon()) return false;
  if (!recoMu->isTrackerMuon()) return false;
  if (recoMu->globalTrack()->normalizedChi2() < 5 && 
      recoMu->globalTrack()->hitPattern().numberOfValidMuonHits() > 5 &&
      recoMu->numberOfMatchedStations() > 1 &&
      recoMu->innerTrack()->hitPattern().numberOfValidPixelHits() > 1 &&
      recoMu->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) isTight=true;

  return isTight;
}

bool  ThreeMuFinder::TriggerDecision(const edm::Event& ev){

  //if (debug) cout << "Reading Trigger decision" << endl;

  bool passed=false;

  // check fired HLT paths
  edm::Handle<edm::TriggerResults> hltresults;
  edm::InputTag trigResultsTag("TriggerResults","",HLT_process);
  ev.getByLabel(trigResultsTag,hltresults);

  if (HLT_paths.size()==0){
    //if (debug) cout << "WARNING:No HLT Path Selected, the event will pass!!!" << endl;
    passed=true;
    return passed;
  }

  if (HLT_paths.size()>10){
    cout << "WARNING:Only the first 10 paths will be considered!!!" << endl;
  }


  if (hltresults.isValid()) {
    const edm::TriggerNames TrigNames_ = ev.triggerNames(*hltresults);
    const int ntrigs = hltresults->size();
    for (int itr=0; itr<ntrigs; itr++){
      if (!hltresults->accept(itr)) continue;
      string trigName=TrigNames_.triggerName(itr);
      //if (debug) cout<<"Found HLT path "<< trigName<<endl;
      int Tsize=HLT_paths.size();
      for (int i=0; i<min(Tsize,10); ++i){
	//if (debug) cout << "accepted " << trigName << endl; 
	if (trigName.find(HLT_paths[i])!=string::npos) {
	  passed=true;
	  _TrigBit[i]=true;
	}
      }
    }
  }
  else
    { 
      cout<<"Trigger results not found"<<endl;
    }

  ///if (passed && (debug)) cout << "Passed!!!!!!" << endl;

  return passed;

}



// ------------ method called to for each event  ------------
void ThreeMuFinder::analyze(const edm::Event& ev, const edm::EventSetup& iSetup) {

  Initialize_TreeVars();
  TriggerDecision(ev);

  _Run=ev.id().run();
  _Evt=ev.id().event();
  _Lum=ev.id().luminosityBlock();
      
  edm::ESHandle<TransientTrackBuilder> Builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", Builder); 
    
  MuonCollection muSkim; // to be used to find the best dimuon
    
  edm::Handle<MuonCollection> muons;
  ev.getByLabel("muons",muons);

  // check the validity of the collection
  if(muons.isValid()){
    for (MuonCollection::const_iterator recoMu = muons->begin(); recoMu!=muons->end(); ++recoMu){ // loop over all muons
      if (!isTight(&*recoMu) || muSkim.size() == 3) {_nmuons++; continue;}
      muSkim.push_back(*recoMu);
    }
  }
  
  if (muSkim.size()!=3) return;
  if (muSkim[0].pt() < 3 || muSkim[1].pt() < 3 || muSkim[2].pt() < 3) return;

  KalmanVertexFitter avf;

  TransientTrack tt1=Builder->build(muSkim[0].innerTrack());
  TransientTrack tt2=Builder->build(muSkim[1].innerTrack());
  TransientTrack tt3=Builder->build(muSkim[2].innerTrack());
  
  std::vector<TransientTrack> tt;
  
  tt.push_back(tt1);
  tt.push_back(tt2);
  tt.push_back(tt3);

  TransientVertex  tmpvtx=avf.vertex(tt);
  
  double vChi2 = tmpvtx.totalChiSquared();
  double vNDF = tmpvtx.degreesOfFreedom();
  
  double vProb(TMath::Prob(vChi2,(int)vNDF));
  
  if (vProb > 0.8 && vChi2/vNDF < 2) {

    TLorentzVector mu1=TLorentzVector(muSkim[0].px(),muSkim[0].py(),muSkim[0].pz(),muSkim[0].energy());
    TLorentzVector mu2=TLorentzVector(muSkim[1].px(),muSkim[1].py(),muSkim[1].pz(),muSkim[1].energy());
    TLorentzVector mu3=TLorentzVector(muSkim[2].px(),muSkim[2].py(),muSkim[2].pz(),muSkim[2].energy());

    TLorentzVector tot=mu1+mu2+mu3;

    _triMu4Mom->SetPxPyPzE(tot.Px(),tot.Py(),tot.Pz(),tot.E());

    _m1->SetPxPyPzE(muSkim[0].px(),muSkim[0].py(),muSkim[0].pz(),muSkim[0].energy());
    _m2->SetPxPyPzE(muSkim[1].px(),muSkim[1].py(),muSkim[1].pz(),muSkim[1].energy());
    _m3->SetPxPyPzE(muSkim[2].px(),muSkim[2].py(),muSkim[2].pz(),muSkim[2].energy());

    q1=muSkim[0].charge();
    q2=muSkim[1].charge();
    q3=muSkim[2].charge();

    SV->SetXYZ(tmpvtx.position().x(),tmpvtx.position().y(),tmpvtx.position().z());
    SVe->SetXYZ(tmpvtx.positionError().cxx(),tmpvtx.positionError().cyy(),tmpvtx.positionError().czz());

    Vprob=vProb;

    ExTree->Fill();
  }
}

void ThreeMuFinder::Initialize_TreeVars(){

  _triMu4Mom->SetPtEtaPhiM(0.,0.,0.,0.);

  _m1->SetPtEtaPhiM(0.,0.,0.,0.);
  _m2->SetPtEtaPhiM(0.,0.,0.,0.);
  _m3->SetPtEtaPhiM(0.,0.,0.,0.);

 for (int k=0; k<10; k++){
    _TrigBit[k]=false;
  }

  q1=0;q2=0;q3=0;  

  SV->SetXYZ(0.,0.,0.);
  SVe->SetXYZ(0.,0.,0.);
  Vprob=-1;

  _Run=0;
  _Lum=0;
  _Evt=0;
  _nmuons=0;
}

// ------------ method called once each job just before starting event loop  ------------
void ThreeMuFinder::beginJob() {

  thefile = new TFile ("OUT3Mu.root", "RECREATE" );
  thefile->cd();

  ExTree = new TTree("tree","tree");

  _triMu4Mom= new TLorentzVector(0.,0.,0.,0.);

  _m1= new TLorentzVector(0.,0.,0.,0.);
  _m2= new TLorentzVector(0.,0.,0.,0.);
  _m3= new TLorentzVector(0.,0.,0.,0.);

  SV=new TVector3(0.,0.,0.);
  SVe=new TVector3(0.,0.,0.);
  int Tsize=HLT_paths.size();
  for (int i=0; i< min(Tsize,10); ++i){
    ExTree->Branch(HLT_paths[i].c_str(), &_TrigBit[i],"_TrigBit/O");
  }

  ExTree->Branch("Run",&_Run   , "_Run/I");
  ExTree->Branch("Lumi",&_Lum   , "_Lum/I");
  ExTree->Branch("Event",&_Evt   , "_Evt/I");

  ExTree->Branch("triMu4Mom","TLorentzVector",&_triMu4Mom); 

  ExTree->Branch("Mu1","TLorentzVector",&_m1); 
  ExTree->Branch("Mu2","TLorentzVector",&_m2); 
  ExTree->Branch("Mu3","TLorentzVector",&_m3); 

  ExTree->Branch("NumberOfAdditionalTightMuons",&_nmuons   , "_nmuons/I");

  ExTree->Branch("q1",&q1,"q1/I");
  ExTree->Branch("q2",&q2,"q2/I");
  ExTree->Branch("q3",&q3,"q3/I");

  ExTree->Branch("Vprob",&Vprob,"Vprob/D");

  ExTree->Branch("SV","TVector3",&SV);
  ExTree->Branch("SVerr","TVector3",&SVe);

}

void 
ThreeMuFinder::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

void 
ThreeMuFinder::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
}


// ------------ method called nce each job just after ending the event loop  ------------
void 
ThreeMuFinder::endJob() {

  thefile->cd();
  thefile->Write();

  delete _triMu4Mom;
  delete SV;
  delete SVe;
  delete ExTree;
  
  thefile->Close();
  delete thefile;

}

#include "FWCore/Framework/interface/MakerMacros.h"  
DEFINE_FWK_MODULE( ThreeMuFinder );
