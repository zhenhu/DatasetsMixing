// -*- C++ -*-
//
// Package:    MixingRootupler
// Class:      MixingRootupler
// 
// Description: Merge MuOnia and ZeroBias datasets. Dump fourMuon decays
//
// Author:  Zhen Hu
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"

//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/Common/interface/View.h>
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

/* For event MIXING */
#include "FWCore/Sources/interface/VectorInputSource.h"
#include "FWCore/Sources/interface/VectorInputSourceDescription.h"
#include "FWCore/Sources/interface/VectorInputSourceFactory.h"
#include "FWCore/Version/interface/GetReleaseVersion.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/GetPassID.h"
#include "FWCore/Framework/src/SignallingProductRegistry.h"
#include "FWCore/Framework/interface/EventPrincipal.h"
#include "boost/cstdint.hpp"
#include "boost/shared_ptr.hpp"
#include "DataFormats/Provenance/interface/BranchIDListHelper.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/Provenance/interface/ModuleDescription.h"
#include "DataFormats/Provenance/interface/ThinnedAssociationsHelper.h"
/* For event MIXING */
namespace edm {
	class EventID;
	class ParameterSet;
}

std::vector<std::vector<pat::Muon>> muons_previousEvent;
std::vector<pat::Muon> muons_previousEvent_bestYMass;
std::vector<std::vector<pat::Muon>> zeroBiasMuons_previousEvent;
const ParticleMass muonMass(0.1056583);
float muonSigma = muonMass*1E-6;

//
// class declaration
//

class MixingRootupler:public edm::EDAnalyzer {
	public:
		explicit MixingRootupler(const edm::ParameterSet &);
		~MixingRootupler();

		static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

	private:
		UInt_t getTriggerBits(const edm::Event &);
		bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
		const  reco::Candidate* GetAncestor(const reco::Candidate *);
		int   tightMuon(edm::View<pat::Muon>::const_iterator rmu, reco::Vertex vertex);
		int   tightMuon(std::vector<pat::Muon>::iterator rmu, reco::Vertex vertex);
		int   mediumMuon(edm::View<pat::Muon>::const_iterator rmu); 
		int   mediumMuon(std::vector<pat::Muon>::iterator rmu);
		void fillUpsilonBestVertex(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs);
		void fillUpsilonBestMass(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs);
		void fourMuonFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);
		int fourMuonMixFit_ZeroBias(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, const std::vector<pat::Muon>* muons_ZeroBias, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);
		int fourMuonMixFit_ZeroBiasTwoEvts(pat::CompositeCandidate dimuonCand,std::vector<pat::Muon> muons_previousA_ZeroBias, std::vector<pat::Muon> muons_previousB_ZeroBias,  const std::vector<pat::Muon>* muons_ZeroBias, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV);
		virtual void beginJob();
		virtual void analyze(const edm::Event &, const edm::EventSetup &);
		virtual void endJob(const edm::Event &);

		virtual void beginRun(edm::Run const &, edm::EventSetup const &);
		virtual void endRun(edm::Run const &, edm::EventSetup const &);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);
		virtual void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &);

		// ----------member data ---------------------------
		std::string file_name;
		edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
		edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
		edm::EDGetTokenT<reco::BeamSpot> bs_Label;
		edm::EDGetTokenT<edm::View<pat::Muon>> muon_Label;
		edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;

		int  pdgid_;
		std::vector<double> OniaMassCuts_;
		bool isMC_;
		bool OnlyBest_;
		bool OnlyGen_;
		double upsilon_mass_;
		uint32_t triggerCuts_;

		/* For event MIXING */
		edm::InputTag mixMuonTag_;
		edm::InputTag mixPVTag_;
		typedef edm::VectorInputSource Source;
		std::unique_ptr<Source> constructSource(const edm::ParameterSet& sourceConfig);
		std::shared_ptr<edm::ProductRegistry> productRegistry_;
		std::unique_ptr<edm::VectorInputSource> const source_;
		std::unique_ptr<edm::ProcessConfiguration> processConfiguration_;
		std::unique_ptr<edm::EventPrincipal> eventPrincipal_;
		template <class T> static const T* getProduct(const edm::EventPrincipal& event, const edm::InputTag& tag);
		void addNextEventToMap(const edm::EventPrincipal& nextEvent);

		UInt_t run;
		UInt_t lumi;
		UInt_t event;
		Int_t  irank;
		UInt_t trigger;
		Int_t numPrimaryVertices;
		Int_t numPrimaryVertices_zerobias;
		Float_t pv_x;
		Float_t pv_y;
		Float_t pv_z;

		TLorentzVector dimuon_p4;
		TLorentzVector mu1_p4;
		TLorentzVector mu2_p4;
		Int_t mu1Charge;
		Int_t mu2Charge;
		Float_t mu1_d0;
		Float_t mu2_d0;
		Float_t mu1_d0err;
		Float_t mu2_d0err;
		Float_t mu1_dz;
		Float_t mu2_dz;
		Float_t mu1_dzerr;
		Float_t mu2_dzerr;
		Float_t mu1_vz;
		Float_t mu2_vz;
		Float_t mumufit_Mass;
		Float_t mumufit_MassErr;
		Float_t mumufit_VtxCL;
		Float_t mumufit_VtxCL2;
		Float_t mumufit_DecayVtxX;
		Float_t mumufit_DecayVtxY;
		Float_t mumufit_DecayVtxZ;
		Float_t mumufit_DecayVtxXE;
		Float_t mumufit_DecayVtxYE;
		Float_t mumufit_DecayVtxZE;
		TLorentzVector mumufit_p4;

		std::vector<Float_t> fourMuFit_Mass_allComb;
		std::vector<Float_t> fourMuFit_Mass;
		std::vector<Float_t> fourMuFit_MassErr;
		std::vector<Float_t> fourMuFit_Pt;
		std::vector<Float_t> fourMuFit_Eta;
		std::vector<Float_t> fourMuFit_Phi;
		std::vector<Float_t> fourMuFit_VtxX;
		std::vector<Float_t> fourMuFit_VtxY;
		std::vector<Float_t> fourMuFit_VtxZ;
		std::vector<Float_t> fourMuFit_VtxProb;
		std::vector<Float_t> fourMuFit_Chi2;
		std::vector<Int_t> fourMuFit_ndof;
		std::vector<Float_t> fourMuFit_mu1Pt;
		std::vector<Float_t> fourMuFit_mu1Eta;
		std::vector<Float_t> fourMuFit_mu1Phi;
		std::vector<Float_t> fourMuFit_mu1E;
		std::vector<Float_t> fourMuFit_mu2Pt;
		std::vector<Float_t> fourMuFit_mu2Eta;
		std::vector<Float_t> fourMuFit_mu2Phi;
		std::vector<Float_t> fourMuFit_mu2E;
		std::vector<Float_t> fourMuFit_mu3Pt;
		std::vector<Float_t> fourMuFit_mu3Eta;
		std::vector<Float_t> fourMuFit_mu3Phi;
		std::vector<Float_t> fourMuFit_mu3E;
		std::vector<Float_t> fourMuFit_mu4Pt;
		std::vector<Float_t> fourMuFit_mu4Eta;
		std::vector<Float_t> fourMuFit_mu4Phi;
		std::vector<Float_t> fourMuFit_mu4E;
		std::vector<Float_t> mu3_Pt;
		std::vector<Float_t> mu3_Eta;
		std::vector<Float_t> mu3_Phi;
		std::vector<Float_t> mu3_E;
		std::vector<Float_t> mu4_Pt;
		std::vector<Float_t> mu4_Eta;
		std::vector<Float_t> mu4_Phi;
		std::vector<Float_t> mu4_E;
		std::vector<Int_t> mu3Charge;
		std::vector<Int_t> mu4Charge;
		std::vector<Float_t> mu3_d0;
		std::vector<Float_t> mu4_d0;
		std::vector<Float_t> mu3_d0err;
		std::vector<Float_t> mu4_d0err;
		std::vector<Float_t> mu3_dz;
		std::vector<Float_t> mu4_dz;
		std::vector<Float_t> mu3_dzerr;
		std::vector<Float_t> mu4_dzerr;
		std::vector<Int_t> mu1_Tight;
		std::vector<Int_t> mu2_Tight;
		std::vector<Int_t> mu3_Tight;
		std::vector<Int_t> mu4_Tight;
		std::vector<Int_t> mu1_Medium;
		std::vector<Int_t> mu2_Medium;
		std::vector<Int_t> mu3_Medium;
		std::vector<Int_t> mu4_Medium;
		std::vector<Int_t> mu3_pdgID;
		std::vector<Int_t> mu4_pdgID;


		TLorentzVector dimuon_p4_bestYMass;
		TLorentzVector mu1_p4_bestYMass;
		TLorentzVector mu2_p4_bestYMass;
		Int_t mu1Charge_bestYMass;
		Int_t mu2Charge_bestYMass;
		Float_t mu1_d0_bestYMass;
		Float_t mu2_d0_bestYMass;
		Float_t mu1_d0err_bestYMass;
		Float_t mu2_d0err_bestYMass;
		Float_t mu1_dz_bestYMass;
		Float_t mu2_dz_bestYMass;
		Float_t mu1_dzerr_bestYMass;
		Float_t mu2_dzerr_bestYMass;
		Float_t mumufit_Mass_bestYMass;
		Float_t mumufit_MassErr_bestYMass;
		Float_t mumufit_VtxCL_bestYMass;
		Float_t mumufit_VtxCL2_bestYMass;
		Float_t mumufit_DecayVtxX_bestYMass;
		Float_t mumufit_DecayVtxY_bestYMass;
		Float_t mumufit_DecayVtxZ_bestYMass;
		Float_t mumufit_DecayVtxXE_bestYMass;
		Float_t mumufit_DecayVtxYE_bestYMass;
		Float_t mumufit_DecayVtxZE_bestYMass;
		TLorentzVector mumufit_p4_bestYMass;
		Int_t bestVertex_and_bestYMass;

		//std::vector<Float_t> fourMuFit_Mass_allComb_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_Mass_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_MassErr_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_Pt_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_Eta_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_Phi_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_VtxX_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_VtxY_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_VtxZ_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_VtxProb_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_Chi2_mix_ZeroBias;
		std::vector<Int_t> fourMuFit_ndof_mix_ZeroBias;
		std::vector<Int_t> fourMuFit_3plus1_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu1Pt_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu1Eta_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu1Phi_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu1E_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu2Pt_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu2Eta_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu2Phi_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu2E_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu3Pt_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu3Eta_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu3Phi_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu3E_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu4Pt_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu4Eta_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu4Phi_mix_ZeroBias;
		std::vector<Float_t> fourMuFit_mu4E_mix_ZeroBias;
		std::vector<Float_t> mu3_Pt_mix_ZeroBias;
		std::vector<Float_t> mu3_Eta_mix_ZeroBias;
		std::vector<Float_t> mu3_Phi_mix_ZeroBias;
		std::vector<Float_t> mu3_E_mix_ZeroBias;
		std::vector<Float_t> mu4_Pt_mix_ZeroBias;
		std::vector<Float_t> mu4_Eta_mix_ZeroBias;
		std::vector<Float_t> mu4_Phi_mix_ZeroBias;
		std::vector<Float_t> mu4_E_mix_ZeroBias;
		std::vector<Int_t> mu3Charge_mix_ZeroBias;
		std::vector<Int_t> mu4Charge_mix_ZeroBias;
		std::vector<Float_t> mu3_d0_mix_ZeroBias;
		std::vector<Float_t> mu4_d0_mix_ZeroBias;
		std::vector<Float_t> mu3_d0err_mix_ZeroBias;
		std::vector<Float_t> mu4_d0err_mix_ZeroBias;
		std::vector<Float_t> mu3_dz_mix_ZeroBias;
		std::vector<Float_t> mu4_dz_mix_ZeroBias;
		std::vector<Float_t> mu3_dzerr_mix_ZeroBias;
		std::vector<Float_t> mu4_dzerr_mix_ZeroBias;
		std::vector<Int_t> mu1_Tight_mix_ZeroBias;
		std::vector<Int_t> mu2_Tight_mix_ZeroBias;
		std::vector<Int_t> mu3_Tight_mix_ZeroBias;
		std::vector<Int_t> mu4_Tight_mix_ZeroBias;
		std::vector<Int_t> mu1_Medium_mix_ZeroBias;
		std::vector<Int_t> mu2_Medium_mix_ZeroBias;
		std::vector<Int_t> mu3_Medium_mix_ZeroBias;
		std::vector<Int_t> mu4_Medium_mix_ZeroBias;


		std::vector<Float_t> fourMuFit_Mass_mix3evts_ZeroBias;
		std::vector<Float_t> fourMuFit_VtxProb_mix3evts_ZeroBias;
		std::vector<Float_t> fourMuFit_Pt_mix3evts_ZeroBias;
		std::vector<Float_t> fourMuFit_Eta_mix3evts_ZeroBias;
		std::vector<Float_t> fourMuFit_Phi_mix3evts_ZeroBias;

		TTree *onia_tree;

		Int_t mother_pdgId;
		Int_t dimuon_pdgId;
		TLorentzVector gen_dimuon_p4;
		TLorentzVector gen_mu1_p4;
		TLorentzVector gen_mu2_p4;

		edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
		edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;
};

//
// constructors and destructor
//

MixingRootupler::MixingRootupler(const edm::ParameterSet & iConfig):
	dimuon_Label(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuons"))),
	primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"))),
	bs_Label(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBeamSpot"))),
	muon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter< edm::InputTag>("muons"))),
	triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
	pdgid_(iConfig.getParameter<uint32_t>("onia_pdgid")),
	OniaMassCuts_(iConfig.getParameter<std::vector<double>>("onia_mass_cuts")),
	isMC_(iConfig.getParameter<bool>("isMC")),
	OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
	OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
	upsilon_mass_(iConfig.getParameter<double>("upsilon_mass")),
	triggerCuts_(iConfig.getParameter<uint32_t>("triggerCuts")),
	/* For event MIXING */
	mixMuonTag_(iConfig.getParameter< edm::InputTag>("muons")),
	mixPVTag_(iConfig.getParameter< edm::InputTag>("primaryVertices")),
	productRegistry_(new edm::SignallingProductRegistry),
	source_(constructSource(iConfig.getParameter<edm::ParameterSet>("SecondSource"))),
	processConfiguration_(new edm::ProcessConfiguration(std::string("@MIXING"), edm::getReleaseVersion(), edm::getPassID())),
	eventPrincipal_()
{
	edm::Service < TFileService > fs;
	onia_tree = fs->make < TTree > ("oniaTree", "Tree of MuMuGamma");

	/* For event MIXING */
	// Use the empty parameter set for the parameter set ID of our "@MIXING" process.
	processConfiguration_->setParameterSetID(edm::ParameterSet::emptyParameterSetID());
	productRegistry_->setFrozen();
	eventPrincipal_.reset(new edm::EventPrincipal(source_->productRegistry(),
				std::make_shared<edm::BranchIDListHelper>(),
				std::make_shared<edm::ThinnedAssociationsHelper>(),
				*processConfiguration_,
				nullptr));

	if (!OnlyGen_) {
		onia_tree->Branch("run",     &run,     "run/I");
		onia_tree->Branch("lumi",     &lumi,     "lumi/I");
		onia_tree->Branch("event",   &event,   "event/I");
		onia_tree->Branch("irank",   &irank,   "irank/I");
		onia_tree->Branch("trigger", &trigger, "trigger/I");
		onia_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/I");
		onia_tree->Branch("numPrimaryVertices_zerobias", &numPrimaryVertices_zerobias, "numPrimaryVertices_zerobias/I");
		onia_tree->Branch("pv_x",     &pv_x,     "pv_x/F");
		onia_tree->Branch("pv_y",     &pv_y,     "pv_y/F");
		onia_tree->Branch("pv_z",     &pv_z,     "pv_z/F");

		onia_tree->Branch("mu1_p4",  "TLorentzVector", &mu1_p4);
		onia_tree->Branch("mu2_p4",  "TLorentzVector", &mu2_p4);
		onia_tree->Branch("mu1Charge",   &mu1Charge,    "mu1Charge/I");
		onia_tree->Branch("mu2Charge",   &mu2Charge,    "mu2Charge/I");
		onia_tree->Branch("mu1_d0",   &mu1_d0,    "mu1_d0/F");
		onia_tree->Branch("mu1_d0err",   &mu1_d0err,    "mu1_d0err/F");
		onia_tree->Branch("mu2_d0",   &mu2_d0,    "mu2_d0/F");
		onia_tree->Branch("mu2_d0err",   &mu2_d0err,    "mu2_d0err/F");
		onia_tree->Branch("mu1_dz",   &mu1_dz,    "mu1_dz/F");
		onia_tree->Branch("mu1_dzerr",   &mu1_dzerr,    "mu1_dzerr/F");
		onia_tree->Branch("mu2_dz",   &mu2_dz,    "mu2_dz/F");
		onia_tree->Branch("mu2_dzerr",   &mu2_dzerr,    "mu2_dzerr/F");
		onia_tree->Branch("mu1_vz",   &mu1_vz,    "mu1_vz/F");
		onia_tree->Branch("mu2_vz",   &mu2_vz,    "mu2_vz/F");
		onia_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
		onia_tree->Branch("mumufit_Mass",&mumufit_Mass,"mumufit_Mass/F");
		onia_tree->Branch("mumufit_MassErr",&mumufit_MassErr,"mumufit_MassErr/F");
		onia_tree->Branch("mumufit_VtxCL",&mumufit_VtxCL,"mumufit_VtxCL/F");
		onia_tree->Branch("mumufit_VtxCL2",&mumufit_VtxCL2,"mumufit_VtxCL2/F");
		onia_tree->Branch("mumufit_DecayVtxX",&mumufit_DecayVtxX,"mumufit_DecayVtxX/F");
		onia_tree->Branch("mumufit_DecayVtxY",&mumufit_DecayVtxY,"mumufit_DecayVtxY/F");
		onia_tree->Branch("mumufit_DecayVtxZ",&mumufit_DecayVtxZ,"mumufit_DecayVtxZ/F");
		onia_tree->Branch("mumufit_DecayVtxXE",&mumufit_DecayVtxXE,"mumufit_DecayVtxXE/F");
		onia_tree->Branch("mumufit_DecayVtxYE",&mumufit_DecayVtxYE,"mumufit_DecayVtxYE/F");
		onia_tree->Branch("mumufit_DecayVtxZE",&mumufit_DecayVtxZE,"mumufit_DecayVtxZE/F");
		onia_tree->Branch("mumufit_p4",  "TLorentzVector", &mumufit_p4);


		onia_tree->Branch("fourMuFit_Mass_allComb",&fourMuFit_Mass_allComb);
		onia_tree->Branch("fourMuFit_Mass",&fourMuFit_Mass);
		onia_tree->Branch("fourMuFit_MassErr",&fourMuFit_MassErr);
		onia_tree->Branch("fourMuFit_Pt",&fourMuFit_Pt);
		onia_tree->Branch("fourMuFit_Eta",&fourMuFit_Eta);
		onia_tree->Branch("fourMuFit_Phi",&fourMuFit_Phi);
		onia_tree->Branch("fourMuFit_VtxX",&fourMuFit_VtxX);
		onia_tree->Branch("fourMuFit_VtxY",&fourMuFit_VtxY);
		onia_tree->Branch("fourMuFit_VtxZ",&fourMuFit_VtxZ);
		onia_tree->Branch("fourMuFit_VtxProb",&fourMuFit_VtxProb);
		onia_tree->Branch("fourMuFit_Chi2",&fourMuFit_Chi2);
		onia_tree->Branch("fourMuFit_ndof",&fourMuFit_ndof);
		onia_tree->Branch("fourMuFit_mu1Pt",&fourMuFit_mu1Pt);
		onia_tree->Branch("fourMuFit_mu1Eta",&fourMuFit_mu1Eta);
		onia_tree->Branch("fourMuFit_mu1Phi",&fourMuFit_mu1Phi);
		onia_tree->Branch("fourMuFit_mu1E",&fourMuFit_mu1E);
		onia_tree->Branch("fourMuFit_mu2Pt",&fourMuFit_mu2Pt);
		onia_tree->Branch("fourMuFit_mu2Eta",&fourMuFit_mu2Eta);
		onia_tree->Branch("fourMuFit_mu2Phi",&fourMuFit_mu2Phi);
		onia_tree->Branch("fourMuFit_mu2E",&fourMuFit_mu2E);
		onia_tree->Branch("fourMuFit_mu3Pt",&fourMuFit_mu3Pt);
		onia_tree->Branch("fourMuFit_mu3Eta",&fourMuFit_mu3Eta);
		onia_tree->Branch("fourMuFit_mu3Phi",&fourMuFit_mu3Phi);
		onia_tree->Branch("fourMuFit_mu3E",&fourMuFit_mu3E);
		onia_tree->Branch("fourMuFit_mu4Pt",&fourMuFit_mu4Pt);
		onia_tree->Branch("fourMuFit_mu4Eta",&fourMuFit_mu4Eta);
		onia_tree->Branch("fourMuFit_mu4Phi",&fourMuFit_mu4Phi);
		onia_tree->Branch("fourMuFit_mu4E",&fourMuFit_mu4E);
		onia_tree->Branch("mu3_Pt",   &mu3_Pt);
		onia_tree->Branch("mu3_Eta",   &mu3_Eta);
		onia_tree->Branch("mu3_Phi",   &mu3_Phi);
		onia_tree->Branch("mu3_E",   &mu3_E);
		onia_tree->Branch("mu4_Pt",   &mu4_Pt);
		onia_tree->Branch("mu4_Eta",   &mu4_Eta);
		onia_tree->Branch("mu4_Phi",   &mu4_Phi);
		onia_tree->Branch("mu4_E",   &mu4_E);
		onia_tree->Branch("mu3Charge",   &mu3Charge);
		onia_tree->Branch("mu4Charge",   &mu4Charge);
		onia_tree->Branch("mu3_d0",   &mu3_d0);
		onia_tree->Branch("mu3_d0err",   &mu3_d0err);
		onia_tree->Branch("mu4_d0",   &mu4_d0);
		onia_tree->Branch("mu4_d0err",   &mu4_d0err);
		onia_tree->Branch("mu3_dz",   &mu3_dz);
		onia_tree->Branch("mu3_dzerr",   &mu3_dzerr);
		onia_tree->Branch("mu4_dz",   &mu4_dz);
		onia_tree->Branch("mu4_dzerr",   &mu4_dzerr);
		onia_tree->Branch("mu1_Tight",   &mu1_Tight);
		onia_tree->Branch("mu2_Tight",   &mu2_Tight);
		onia_tree->Branch("mu3_Tight",   &mu3_Tight);
		onia_tree->Branch("mu4_Tight",   &mu4_Tight);
		onia_tree->Branch("mu1_Medium",   &mu1_Medium);
		onia_tree->Branch("mu2_Medium",   &mu2_Medium);
		onia_tree->Branch("mu3_Medium",   &mu3_Medium);
		onia_tree->Branch("mu4_Medium",   &mu4_Medium);
		onia_tree->Branch("mu3_pdgID",   &mu3_pdgID);
		onia_tree->Branch("mu4_pdgID",   &mu4_pdgID);

		onia_tree->Branch("mu1_p4_bestYMass",  "TLorentzVector", &mu1_p4_bestYMass);
		onia_tree->Branch("mu2_p4_bestYMass",  "TLorentzVector", &mu2_p4_bestYMass);
		onia_tree->Branch("mu1Charge_bestYMass",   &mu1Charge_bestYMass,    "mu1Charge_bestYMass/I");
		onia_tree->Branch("mu2Charge_bestYMass",   &mu2Charge_bestYMass,    "mu2Charge_bestYMass/I");
		onia_tree->Branch("mu1_d0_bestYMass",   &mu1_d0_bestYMass,    "mu1_d0_bestYMass/F");
		onia_tree->Branch("mu1_d0err_bestYMass",   &mu1_d0err_bestYMass,    "mu1_d0err_bestYMass/F");
		onia_tree->Branch("mu2_d0_bestYMass",   &mu2_d0_bestYMass,    "mu2_d0_bestYMass/F");
		onia_tree->Branch("mu2_d0err_bestYMass",   &mu2_d0err_bestYMass,    "mu2_d0err_bestYMass/F");
		onia_tree->Branch("mu1_dz_bestYMass",   &mu1_dz_bestYMass,    "mu1_dz_bestYMass/F");
		onia_tree->Branch("mu1_dzerr_bestYMass",   &mu1_dzerr_bestYMass,    "mu1_dzerr_bestYMass/F");
		onia_tree->Branch("mu2_dz_bestYMass",   &mu2_dz_bestYMass,    "mu2_dz_bestYMass/F");
		onia_tree->Branch("mu2_dzerr_bestYMass",   &mu2_dzerr_bestYMass,    "mu2_dzerr_bestYMass/F");
		onia_tree->Branch("dimuon_p4_bestYMass", "TLorentzVector", &dimuon_p4_bestYMass);
		onia_tree->Branch("mumufit_Mass_bestYMass",&mumufit_Mass_bestYMass,"mumufit_Mass_bestYMass/F");
		onia_tree->Branch("mumufit_MassErr_bestYMass",&mumufit_MassErr_bestYMass,"mumufit_MassErr_bestYMass/F");
		onia_tree->Branch("mumufit_VtxCL_bestYMass",&mumufit_VtxCL_bestYMass,"mumufit_VtxCL_bestYMass/F");
		onia_tree->Branch("mumufit_VtxCL2_bestYMass",&mumufit_VtxCL2_bestYMass,"mumufit_VtxCL2_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxX_bestYMass",&mumufit_DecayVtxX_bestYMass,"mumufit_DecayVtxX_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxY_bestYMass",&mumufit_DecayVtxY_bestYMass,"mumufit_DecayVtxY_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxZ_bestYMass",&mumufit_DecayVtxZ_bestYMass,"mumufit_DecayVtxZ_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxXE_bestYMass",&mumufit_DecayVtxXE_bestYMass,"mumufit_DecayVtxXE_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxYE_bestYMass",&mumufit_DecayVtxYE_bestYMass,"mumufit_DecayVtxYE_bestYMass/F");
		onia_tree->Branch("mumufit_DecayVtxZE_bestYMass",&mumufit_DecayVtxZE_bestYMass,"mumufit_DecayVtxZE_bestYMass/F");
		onia_tree->Branch("mumufit_p4_bestYMass",  "TLorentzVector", &mumufit_p4_bestYMass);
		onia_tree->Branch("bestVertex_and_bestYMass", &bestVertex_and_bestYMass,"bestVertex_and_bestYMass/I");

		//onia_tree->Branch("fourMuFit_Mass_allComb_mix_ZeroBias",&fourMuFit_Mass_allComb_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_Mass_mix_ZeroBias",&fourMuFit_Mass_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_MassErr_mix_ZeroBias",&fourMuFit_MassErr_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_Pt_mix_ZeroBias",&fourMuFit_Pt_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_Eta_mix_ZeroBias",&fourMuFit_Eta_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_Phi_mix_ZeroBias",&fourMuFit_Phi_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_VtxX_mix_ZeroBias",&fourMuFit_VtxX_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_VtxY_mix_ZeroBias",&fourMuFit_VtxY_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_VtxZ_mix_ZeroBias",&fourMuFit_VtxZ_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_VtxProb_mix_ZeroBias",&fourMuFit_VtxProb_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_Chi2_mix_ZeroBias",&fourMuFit_Chi2_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_ndof_mix_ZeroBias",&fourMuFit_ndof_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_3plus1_mix_ZeroBias",  &fourMuFit_3plus1_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu1Pt_mix_ZeroBias",&fourMuFit_mu1Pt_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu1Eta_mix_ZeroBias",&fourMuFit_mu1Eta_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu1Phi_mix_ZeroBias",&fourMuFit_mu1Phi_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu1E_mix_ZeroBias",&fourMuFit_mu1E_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu2Pt_mix_ZeroBias",&fourMuFit_mu2Pt_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu2Eta_mix_ZeroBias",&fourMuFit_mu2Eta_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu2Phi_mix_ZeroBias",&fourMuFit_mu2Phi_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu2E_mix_ZeroBias",&fourMuFit_mu2E_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu3Pt_mix_ZeroBias",&fourMuFit_mu3Pt_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu3Eta_mix_ZeroBias",&fourMuFit_mu3Eta_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu3Phi_mix_ZeroBias",&fourMuFit_mu3Phi_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu3E_mix_ZeroBias",&fourMuFit_mu3E_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu4Pt_mix_ZeroBias",&fourMuFit_mu4Pt_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu4Eta_mix_ZeroBias",&fourMuFit_mu4Eta_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu4Phi_mix_ZeroBias",&fourMuFit_mu4Phi_mix_ZeroBias);
		onia_tree->Branch("fourMuFit_mu4E_mix_ZeroBias",&fourMuFit_mu4E_mix_ZeroBias);
		onia_tree->Branch("mu3_Pt_mix_ZeroBias",   &mu3_Pt_mix_ZeroBias);
		onia_tree->Branch("mu3_Eta_mix_ZeroBias",   &mu3_Eta_mix_ZeroBias);
		onia_tree->Branch("mu3_Phi_mix_ZeroBias",   &mu3_Phi_mix_ZeroBias);
		onia_tree->Branch("mu3_E_mix_ZeroBias",   &mu3_E_mix_ZeroBias);
		onia_tree->Branch("mu4_Pt_mix_ZeroBias",   &mu4_Pt_mix_ZeroBias);
		onia_tree->Branch("mu4_Eta_mix_ZeroBias",   &mu4_Eta_mix_ZeroBias);
		onia_tree->Branch("mu4_Phi_mix_ZeroBias",   &mu4_Phi_mix_ZeroBias);
		onia_tree->Branch("mu4_E_mix_ZeroBias",   &mu4_E_mix_ZeroBias);
		onia_tree->Branch("mu3Charge_mix_ZeroBias",   &mu3Charge_mix_ZeroBias);
		onia_tree->Branch("mu4Charge_mix_ZeroBias",   &mu4Charge_mix_ZeroBias);
		onia_tree->Branch("mu3_d0_mix_ZeroBias",   &mu3_d0_mix_ZeroBias);
		onia_tree->Branch("mu3_d0err_mix_ZeroBias",   &mu3_d0err_mix_ZeroBias);
		onia_tree->Branch("mu4_d0_mix_ZeroBias",   &mu4_d0_mix_ZeroBias);
		onia_tree->Branch("mu4_d0err_mix_ZeroBias",   &mu4_d0err_mix_ZeroBias);
		onia_tree->Branch("mu3_dz_mix_ZeroBias",   &mu3_dz_mix_ZeroBias);
		onia_tree->Branch("mu3_dzerr_mix_ZeroBias",   &mu3_dzerr_mix_ZeroBias);
		onia_tree->Branch("mu4_dz_mix_ZeroBias",   &mu4_dz_mix_ZeroBias);
		onia_tree->Branch("mu4_dzerr_mix_ZeroBias",   &mu4_dzerr_mix_ZeroBias);
		onia_tree->Branch("mu1_Tight_mix_ZeroBias",   &mu1_Tight_mix_ZeroBias);
		onia_tree->Branch("mu2_Tight_mix_ZeroBias",   &mu2_Tight_mix_ZeroBias);
		onia_tree->Branch("mu3_Tight_mix_ZeroBias",   &mu3_Tight_mix_ZeroBias);
		onia_tree->Branch("mu4_Tight_mix_ZeroBias",   &mu4_Tight_mix_ZeroBias);
		onia_tree->Branch("mu1_Medium_mix_ZeroBias",   &mu1_Medium_mix_ZeroBias);
		onia_tree->Branch("mu2_Medium_mix_ZeroBias",   &mu2_Medium_mix_ZeroBias);
		onia_tree->Branch("mu3_Medium_mix_ZeroBias",   &mu3_Medium_mix_ZeroBias);
		onia_tree->Branch("mu4_Medium_mix_ZeroBias",   &mu4_Medium_mix_ZeroBias);


		onia_tree->Branch("fourMuFit_Mass_mix3evts_ZeroBias",  &fourMuFit_Mass_mix3evts_ZeroBias);
		onia_tree->Branch("fourMuFit_Pt_mix3evts_ZeroBias",  &fourMuFit_Pt_mix3evts_ZeroBias);
		onia_tree->Branch("fourMuFit_Eta_mix3evts_ZeroBias",  &fourMuFit_Eta_mix3evts_ZeroBias);
		onia_tree->Branch("fourMuFit_Phi_mix3evts_ZeroBias",  &fourMuFit_Phi_mix3evts_ZeroBias);
		onia_tree->Branch("fourMuFit_VtxProb_mix3evts_ZeroBias",  &fourMuFit_VtxProb_mix3evts_ZeroBias);

	}

	if (isMC_ || OnlyGen_) {
		onia_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
		onia_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
		onia_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
		onia_tree->Branch("gen_mu1_p4",  "TLorentzVector",  &gen_mu1_p4);
		onia_tree->Branch("gen_mu2_p4",  "TLorentzVector",  &gen_mu2_p4);
	}
	genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
	packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");

}

MixingRootupler::~MixingRootupler() {}

//
// member functions
//

/* For event MIXING */
std::unique_ptr<MixingRootupler::Source> MixingRootupler::constructSource(const edm::ParameterSet& sourceConfig)
{
	const edm::VectorInputSourceFactory* sourceFactory = edm::VectorInputSourceFactory::get();
	edm::VectorInputSourceDescription description(productRegistry_, edm::PreallocationConfiguration());
	return sourceFactory->makeVectorInputSource(sourceConfig, description);
}
/* For event MIXING */
template <class T> const T* MixingRootupler::getProduct(const edm::EventPrincipal& event, const edm::InputTag& tag)
{
	//std::cout << "Retrieving product " << tag <<std::endl;
	// Note: The third argument to getProductByTag can be a nullptr
	// as long as unscheduled execution of an EDProducer cannot occur
	// as a result of this function call (and with the current implementation
	// of SpyEventMatcher unscheduled execution never happens).
	const std::shared_ptr< const edm::Wrapper<T> > productWrapper = edm::getProductByTag<T>(event,tag,nullptr);
	if (productWrapper) {
		return productWrapper->product();
	} else {
		return nullptr;
	}
}
/* For event MIXING */
void MixingRootupler::addNextEventToMap(const edm::EventPrincipal& nextEvent) {
	edm::EventID mixEventId = nextEvent.id();
	std::cout<<"Event ID for the second source: "<<mixEventId.event()<<std::endl;
}


const reco::Candidate* MixingRootupler::GetAncestor(const reco::Candidate* p) {
	if (p->numberOfMothers()) {
		if  ((p->mother(0))->pdgId() == p->pdgId()) return GetAncestor(p->mother(0));
		else return p->mother(0);
	}
	//std::cout << "GetAncestor: Inconsistet ancestor, particle does not have a mother " << std::endl;
	return p;
}

//Check recursively if any ancestor of particle is the given one
bool MixingRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
	if (ancestor == particle ) return true;
	for (size_t i=0; i< particle->numberOfMothers(); i++) {
		if (isAncestor(ancestor, particle->mother(i))) return true;
	}
	return false;
}

/* Grab Trigger information. Save it in variable trigger, trigger is an int between 0 and 127, in binary it is:
	(pass 2)(pass 1)(pass 0)
	ex. 7 = pass 0, 1 and 2
	ex. 6 = pass 2, 3
	ex. 1 = pass 0
	*/

UInt_t MixingRootupler::getTriggerBits(const edm::Event& iEvent ) {
	UInt_t itrigger = 0;
	edm::Handle<edm::TriggerResults> triggerResults_handle;
	iEvent.getByToken(triggerResults_Label, triggerResults_handle);
	//   if ( triggerResults_handle.isValid() ) { 
	//     std::string testTriggerName;
	//     const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
	//     for(unsigned int trigIndex = 0; trigIndex < TheTriggerNames.size(); trigIndex++){
	//     testTriggerName = TheTriggerNames.triggerName(trigIndex);
	//     std::cout<<testTriggerName.c_str()<<std::endl;
	//     }
	//   }
	if ( triggerResults_handle.isValid() ) { 
		const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
		std::vector <unsigned int> bits_0, bits_1, bits_2, bits_3, bits_4, bits_5, bits_6, bits_7, bits_8, bits_9;
		for ( int version = 1; version<20; version ++ ) {
			std::stringstream ss0,ss1,ss2,ss3,ss4,ss5,ss6,ss7,ss8,ss9;
			ss0<<"HLT_Dimuon25_Jpsi_v"<<version;
			//ss0<<"HLT_Dimuon16_Jpsi_v"<<version;
			bits_0.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss0.str()).label().c_str()));

			ss1<<"HLT_DiMuon0_Upsilon_Muon_v"<<version;
			//ss1<<"HLT_Dimuon13_PsiPrime_v"<<version;
			bits_1.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss1.str()).label().c_str()));

			ss2<<"HLT_Dimuon13_Upsilon_v"<<version;
			bits_2.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss2.str()).label().c_str()));

			ss3<<"HLT_DiMuon0_Jpsi_Muon_v"<<version;
			//ss3<<"HLT_Dimuon10_Jpsi_Barrel_v"<<version;
			bits_3.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss3.str()).label().c_str()));

			ss4<<"HLT_TripleMu5_v"<<version; 
			//ss4<<"HLT_Dimuon8_PsiPrime_Barrel_v"<<version;
			bits_4.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss4.str()).label().c_str()));

			ss5<<"HLT_Dimuon12_Upsilon_eta1p5_v"<<version;
			//ss5<<"HLT_Dimuon8_Upsilon_Barrel_v"<<version;
			bits_5.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss5.str()).label().c_str()));

			ss6<<"HLT_Dimuon20_Jpsi_v"<<version;
			bits_6.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss6.str()).label().c_str()));

			ss7<<"HLT_Trimuon2_Upsilon5_Muon_v"<<version;
			//ss7<<"HLT_Dimuon14_Phi_Barrel_Seagulls_v"<<version;
			//ss7<<"HLT_Dimuon0_Phi_Barrel_v"<<version;
			bits_7.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss7.str()).label().c_str()));

			ss8<<"HLT_Trimuon5_3p5_2_Upsilon_Muon_v"<<version;
			//ss8<<"HLT_DoubleMu4_3_Bs_v"<<version;
			bits_8.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss8.str()).label().c_str()));

			ss9<<"HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon_v"<<version;
			//ss9<<"HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"<<version;
			bits_9.push_back(TheTriggerNames.triggerIndex( edm::InputTag(ss9.str()).label().c_str()));
		}
		for (unsigned int i=0; i<bits_0.size(); i++) {
			unsigned int bit = bits_0[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 1;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_1.size(); i++) {
			unsigned int bit = bits_1[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 2;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_2.size(); i++) {
			unsigned int bit = bits_2[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 4;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_3.size(); i++) {
			unsigned int bit = bits_3[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 8;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_4.size(); i++) {
			unsigned int bit = bits_4[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 16;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_5.size(); i++) {
			unsigned int bit = bits_5[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 32;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_6.size(); i++) {
			unsigned int bit = bits_6[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 64;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_7.size(); i++) {
			unsigned int bit = bits_7[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) {
					itrigger += 128;
					break;
				}
			}
		}
		for (unsigned int i=0; i<bits_8.size(); i++) {
			unsigned int bit = bits_8[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) { 
					itrigger += 256;
					break;
				}   
			}   
		}   
		for (unsigned int i=0; i<bits_9.size(); i++) {
			unsigned int bit = bits_9[i];
			if ( bit < triggerResults_handle->size() ){
				if ( triggerResults_handle->accept( bit ) && !triggerResults_handle->error( bit ) ) { 
					itrigger += 512;
					break;
				}   
			}   
		}   
	}
	return itrigger;
}

// ------------ method called for each event  ------------
void MixingRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

	edm::Handle<pat::CompositeCandidateCollection> dimuons;
	iEvent.getByToken(dimuon_Label,dimuons);

	edm::Handle< edm::View<pat::Muon> > muons;
	iEvent.getByToken(muon_Label, muons);

	size_t fileNameHash = 0U;
	source_->loopOverEvents(*eventPrincipal_,fileNameHash,1,[this](auto const& iE, auto const&){ this->addNextEventToMap(iE); },nullptr,nullptr,false);
	const std::vector<pat::Muon> * zeroBiasMuons = getProduct< std::vector<pat::Muon> >(*eventPrincipal_,mixMuonTag_);

	//prepare muons from previous event for 3 events mixing or 4 events mixing 
	std::vector<pat::Muon> saveAsPreviousEvent;
	for (std::vector<pat::Muon>::const_iterator mu = zeroBiasMuons->begin(), muend = zeroBiasMuons->end(); mu != muend; ++mu){
		saveAsPreviousEvent.push_back(*mu);
	}
	zeroBiasMuons_previousEvent.push_back(saveAsPreviousEvent);

	const reco::VertexCollection * primaryVertices_handle_zerobias = getProduct< reco::VertexCollection >(*eventPrincipal_,mixPVTag_);
	numPrimaryVertices_zerobias=-1;	
	if (primaryVertices_handle_zerobias->size() > 0)  numPrimaryVertices_zerobias = (int) primaryVertices_handle_zerobias->size();

	edm::Handle<reco::VertexCollection> primaryVertices_handle;
	iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);
	reco::Vertex thePrimaryV;

	edm::Handle<reco::BeamSpot> theBeamSpot;
	iEvent.getByToken(bs_Label,theBeamSpot);
	reco::BeamSpot bs = *theBeamSpot;

	if ( primaryVertices_handle->begin() != primaryVertices_handle->end() ) {  
		thePrimaryV = reco::Vertex(*(primaryVertices_handle->begin()));
	}
	else {
		thePrimaryV = reco::Vertex(bs.position(), bs.covariance3D());
	}
	pv_x = thePrimaryV.x();
	pv_y = thePrimaryV.y();
	pv_z = thePrimaryV.z();

	if (!OnlyGen_) {
		numPrimaryVertices = -1;
		if (primaryVertices_handle.isValid()) 	numPrimaryVertices = (int) primaryVertices_handle->size();

		trigger = getTriggerBits(iEvent);
		run     = iEvent.id().run();
		lumi    = iEvent.id().luminosityBlock();
		event   = iEvent.id().event();
		//std::cout<<"trigger:"<<trigger<<std::endl;
	}

	dimuon_pdgId = 0;
	mother_pdgId = 0;
	irank = 0;

	dimuon_p4.SetPtEtaPhiM(0,0,0,0);
	mu1_p4.SetPtEtaPhiM(0,0,0,0);
	mu2_p4.SetPtEtaPhiM(0,0,0,0);
	mu1Charge = -10;
	mu2Charge = -10;
	mu1_d0 = -10;
	mu1_d0err = -10;
	mu2_d0 = -10;
	mu2_d0err = -10;
	mu1_dz = -1000;
	mu1_dzerr = -1000;
	mu2_dz = -1000;
	mu2_dzerr = -1000;
	mu1_vz = -1000;
	mu2_vz = -1000;
	mumufit_Mass = -10;
	mumufit_MassErr = -10;
	mumufit_VtxCL = -10;
	mumufit_VtxCL2 = -10;
	mumufit_DecayVtxX = -10;
	mumufit_DecayVtxY = -10;
	mumufit_DecayVtxZ = -10;
	mumufit_DecayVtxXE = -10;
	mumufit_DecayVtxYE = -10;
	mumufit_DecayVtxZE = -10;
	mumufit_p4.SetPtEtaPhiM(0,0,0,0);

	fourMuFit_Mass_allComb.clear();
	fourMuFit_Mass.clear();
	fourMuFit_MassErr.clear();
	fourMuFit_Pt.clear();
	fourMuFit_Eta.clear();
	fourMuFit_Phi.clear();
	fourMuFit_VtxX.clear();
	fourMuFit_VtxY.clear();
	fourMuFit_VtxZ.clear();
	fourMuFit_VtxProb.clear();
	fourMuFit_Chi2.clear();
	fourMuFit_ndof.clear();
	fourMuFit_mu1Pt.clear();
	fourMuFit_mu1Eta.clear();
	fourMuFit_mu1Phi.clear();
	fourMuFit_mu1E.clear();
	fourMuFit_mu2Pt.clear();
	fourMuFit_mu2Eta.clear();
	fourMuFit_mu2Phi.clear();
	fourMuFit_mu2E.clear();
	fourMuFit_mu3Pt.clear();
	fourMuFit_mu3Eta.clear();
	fourMuFit_mu3Phi.clear();
	fourMuFit_mu3E.clear();
	fourMuFit_mu4Pt.clear();
	fourMuFit_mu4Eta.clear();
	fourMuFit_mu4Phi.clear();
	fourMuFit_mu4E.clear();
	mu3_Pt.clear();
	mu3_Eta.clear();
	mu3_Phi.clear();
	mu3_E.clear();
	mu4_Pt.clear();
	mu4_Eta.clear();
	mu4_Phi.clear();
	mu4_E.clear();
	mu3_d0.clear();
	mu3_d0err.clear();
	mu4_d0.clear();
	mu4_d0err.clear();
	mu3_dz.clear();
	mu3_dzerr.clear();
	mu4_dz.clear();
	mu4_dzerr.clear();
	mu3Charge.clear(); 
	mu4Charge.clear(); 
	mu1_Tight.clear();
	mu2_Tight.clear();
	mu3_Tight.clear();
	mu4_Tight.clear();
	mu1_Medium.clear();
	mu2_Medium.clear();
	mu3_Medium.clear();
	mu4_Medium.clear();
	mu3_pdgID.clear();
	mu4_pdgID.clear();

	gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_mu1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_mu2_p4.SetPtEtaPhiM(0.,0.,0.,0.);


	dimuon_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu1_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu2_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	mu1Charge_bestYMass = -10; 
	mu2Charge_bestYMass = -10; 
	mu1_d0_bestYMass = -10; 
	mu1_d0err_bestYMass = -10; 
	mu2_d0_bestYMass = -10; 
	mu2_d0err_bestYMass = -10; 
	mu1_dz_bestYMass = -1000;
	mu1_dzerr_bestYMass = -1000;
	mu2_dz_bestYMass = -1000;
	mu2_dzerr_bestYMass = -1000;
	mumufit_Mass_bestYMass = -10; 
	mumufit_MassErr_bestYMass = -10; 
	mumufit_VtxCL_bestYMass = -10; 
	mumufit_VtxCL2_bestYMass = -10; 
	mumufit_DecayVtxX_bestYMass = -10; 
	mumufit_DecayVtxY_bestYMass = -10; 
	mumufit_DecayVtxZ_bestYMass = -10; 
	mumufit_DecayVtxXE_bestYMass = -10; 
	mumufit_DecayVtxYE_bestYMass = -10; 
	mumufit_DecayVtxZE_bestYMass = -10; 
	mumufit_p4_bestYMass.SetPtEtaPhiM(0,0,0,0);
	bestVertex_and_bestYMass = -1;

	//fourMuFit_Mass_allComb_mix_ZeroBias.clear();
	fourMuFit_Mass_mix_ZeroBias.clear();
	fourMuFit_MassErr_mix_ZeroBias.clear();
	fourMuFit_Pt_mix_ZeroBias.clear();
	fourMuFit_Eta_mix_ZeroBias.clear();
	fourMuFit_Phi_mix_ZeroBias.clear();
	fourMuFit_VtxX_mix_ZeroBias.clear();
	fourMuFit_VtxY_mix_ZeroBias.clear();
	fourMuFit_VtxZ_mix_ZeroBias.clear();
	fourMuFit_VtxProb_mix_ZeroBias.clear();
	fourMuFit_Chi2_mix_ZeroBias.clear();
	fourMuFit_ndof_mix_ZeroBias.clear();
	fourMuFit_3plus1_mix_ZeroBias.clear();
	fourMuFit_mu1Pt_mix_ZeroBias.clear();
	fourMuFit_mu1Eta_mix_ZeroBias.clear();
	fourMuFit_mu1Phi_mix_ZeroBias.clear();
	fourMuFit_mu1E_mix_ZeroBias.clear();
	fourMuFit_mu2Pt_mix_ZeroBias.clear();
	fourMuFit_mu2Eta_mix_ZeroBias.clear();
	fourMuFit_mu2Phi_mix_ZeroBias.clear();
	fourMuFit_mu2E_mix_ZeroBias.clear();
	fourMuFit_mu3Pt_mix_ZeroBias.clear();
	fourMuFit_mu3Eta_mix_ZeroBias.clear();
	fourMuFit_mu3Phi_mix_ZeroBias.clear();
	fourMuFit_mu3E_mix_ZeroBias.clear();
	fourMuFit_mu4Pt_mix_ZeroBias.clear();
	fourMuFit_mu4Eta_mix_ZeroBias.clear();
	fourMuFit_mu4Phi_mix_ZeroBias.clear();
	fourMuFit_mu4E_mix_ZeroBias.clear();
	mu3_Pt_mix_ZeroBias.clear();
	mu3_Eta_mix_ZeroBias.clear();
	mu3_Phi_mix_ZeroBias.clear();
	mu3_E_mix_ZeroBias.clear();
	mu4_Pt_mix_ZeroBias.clear();
	mu4_Eta_mix_ZeroBias.clear();
	mu4_Phi_mix_ZeroBias.clear();
	mu4_E_mix_ZeroBias.clear();
	mu3_d0_mix_ZeroBias.clear();
	mu3_d0err_mix_ZeroBias.clear();
	mu4_d0_mix_ZeroBias.clear();
	mu4_d0err_mix_ZeroBias.clear();
	mu3_dz_mix_ZeroBias.clear();
	mu3_dzerr_mix_ZeroBias.clear();
	mu4_dz_mix_ZeroBias.clear();
	mu4_dzerr_mix_ZeroBias.clear();
	mu3Charge_mix_ZeroBias.clear();
	mu4Charge_mix_ZeroBias.clear();
	mu1_Tight_mix_ZeroBias.clear();
	mu2_Tight_mix_ZeroBias.clear();
	mu3_Tight_mix_ZeroBias.clear();
	mu4_Tight_mix_ZeroBias.clear();
	mu1_Medium_mix_ZeroBias.clear();
	mu2_Medium_mix_ZeroBias.clear();
	mu3_Medium_mix_ZeroBias.clear();
	mu4_Medium_mix_ZeroBias.clear();

	fourMuFit_Mass_mix3evts_ZeroBias.clear();
	fourMuFit_Pt_mix3evts_ZeroBias.clear();
	fourMuFit_Eta_mix3evts_ZeroBias.clear();
	fourMuFit_Phi_mix3evts_ZeroBias.clear();
	fourMuFit_VtxProb_mix3evts_ZeroBias.clear();

	// Pruned particles are the one containing "important" stuff
	edm::Handle<reco::GenParticleCollection> pruned;
	iEvent.getByToken(genCands_, pruned);

	// Packed particles are all the status 1, so usable to remake jets
	// The navigation from status 1 to pruned is possible (the other direction should be made by hand)
	edm::Handle<pat::PackedGenParticleCollection> packed;
	iEvent.getByToken(packCands_,  packed);

	if ((isMC_ || OnlyGen_) && packed.isValid() && pruned.isValid()) {
		dimuon_pdgId  = 0;
		gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
		int foundit   = 0;

		for (size_t i=0; i<pruned->size(); i++) {
			int p_id = abs((*pruned)[i].pdgId());
			const reco::Candidate *aonia = &(*pruned)[i];
			if (( p_id == pdgid_ ) && (aonia->status() == 2)) {
				dimuon_pdgId = p_id;
				foundit++;
				for (size_t j=0; j<packed->size(); j++) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
					const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
					const reco::Candidate * d = &(*packed)[j];
					if ( motherInPrunedCollection != nullptr && (d->pdgId() == 13 ) && isAncestor(aonia , motherInPrunedCollection) ){
						gen_mu2_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
						foundit++;
					} 
					if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aonia , motherInPrunedCollection) ) {
						gen_mu1_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
						foundit++;
					}
					if ( foundit == 3 ) break;               
				}
				if ( foundit == 3 ) {
					gen_dimuon_p4 = gen_mu2_p4 + gen_mu1_p4;   // this should take into account FSR
					mother_pdgId  = GetAncestor(aonia)->pdgId();
					break;
				} else {
					foundit = 0;
					dimuon_pdgId = 0;
					mother_pdgId = 0;
					gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
				}            
			}  // if ( p_id
		} // for (size

		// sanity check
		//if ( ! dimuon_pdgId ) std::cout << "MixingRootupler: does not found the given decay " << run << "," << event << std::endl;
	}  // end if isMC

	float OniaMassMax_ = OniaMassCuts_[1];
	float OniaMassMin_ = OniaMassCuts_[0];

	// Kinematic fit

	edm::ESHandle<TransientTrackBuilder> theB; 
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 
	edm::ESHandle<MagneticField> bFieldHandle;
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

	int nGoodUpsilonCand = 0;
	float bestYMass = 1000;
	pat::CompositeCandidate DimuonCand_bestYMass;
	if ( ! OnlyGen_ 
			&& dimuons.isValid() && dimuons->size() > 0
			&& muons.isValid()) {
		for(pat::CompositeCandidateCollection::const_iterator dimuonCand=dimuons->begin();dimuonCand!= dimuons->end(); ++dimuonCand)
		{
			if (dimuonCand->mass() < OniaMassMin_ || dimuonCand->mass() > OniaMassMax_) continue;
			if (dimuonCand->daughter("muon1")->charge() == dimuonCand->daughter("muon2")->charge() ) continue;
			if (dimuonCand->daughter("muon1")->pt()<2.0 || dimuonCand->daughter("muon2")->pt()<2.0 ) continue;
			if (dimuonCand->daughter("muon1")->eta()>2.4|| dimuonCand->daughter("muon2")->eta()>2.4) continue;

			//dimuon refit. 
			//Here we use the KinematicParticleVertexFitter with muon mass. But in the Onia2MuMu skim, it was just KalmanVertexFitter. 
			//Fitted Vertex from both methods are the same (dimuonCand->userData<reco::Vertex>("commonVertex")->z() == mumu_vFit_vertex_noMC->position().z()), 
			//but fitted p4 and mass are different. 
			reco::TrackRef JpsiTk[2]={  //this is from Chib code
				( dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon1") ) )->innerTrack(),
				( dynamic_cast<const pat::Muon*>(dimuonCand->daughter("muon2") ) )->innerTrack()
			};

			//std::vector<reco::TransientTrack> MuMuTT;
			//MuMuTT.push_back((*theB).build(&JpsiTk[0]));
			//MuMuTT.push_back((*theB).build(&JpsiTk[1]));	
			reco::TransientTrack muon1TT(JpsiTk[0], &(*bFieldHandle));
			reco::TransientTrack muon2TT(JpsiTk[1], &(*bFieldHandle));
			//std::cout<<muon1TT.isValid()<<" "<<muon1TT.track().pt()<<std::endl;
			reco::TrackTransientTrack muon1TTT(JpsiTk[0], &(*bFieldHandle));
			reco::TrackTransientTrack muon2TTT(JpsiTk[1], &(*bFieldHandle));

			KinematicParticleFactoryFromTransientTrack pmumuFactory;
			std::vector<RefCountedKinematicParticle> mumuParticles;
			mumuParticles.push_back(pmumuFactory.particle(muon1TT,muonMass,float(0),float(0),muonSigma));
			mumuParticles.push_back(pmumuFactory.particle(muon2TT,muonMass,float(0),float(0),muonSigma));

			KinematicParticleVertexFitter mumufitter;
			RefCountedKinematicTree mumuVertexFitTree;
			//try {mumuVertexFitTree = mumufitter.fit(mumuParticles);}
			//catch (...) {
			//std::cout<<"mumu fit: PerigeeKinematicState::kinematic state passed is not valid!"<<std::endl;
			//continue;
			//}   
			mumuVertexFitTree = mumufitter.fit(mumuParticles);
			RefCountedKinematicParticle mumu_vFit_noMC;
			RefCountedKinematicVertex mumu_vFit_vertex_noMC;

			if (!(mumuVertexFitTree->isValid())) continue;
			mumuVertexFitTree->movePointerToTheTop();
			mumu_vFit_noMC = mumuVertexFitTree->currentParticle();
			mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex();

			//if (mumu_vFit_noMC->currentState().mass() < 8 || mumu_vFit_noMC->currentState().mass() > 12) continue;
			if (fabs(mumu_vFit_noMC->currentState().mass()-upsilon_mass_) > (3*1.105*sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) ))) continue;
			if (ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) < 0.005) continue;
			//apply trigger cut: 36 for upsilon, 73 for Jpsi
			//if ((trigger&triggerCuts_)==0) continue;
			//if (dimuonCand->pt() < 7) ccontinue; //another method: using mumufit_ instead of dimuon			
			nGoodUpsilonCand++;
			pat::CompositeCandidate thisDimuonCand = *dimuonCand;

			if (nGoodUpsilonCand==1) {
				fillUpsilonBestVertex(mumuVertexFitTree,thisDimuonCand,bFieldHandle,bs);
				//main result: 2 events (MuOnia, ZeroBias) mixing
				fourMuonMixFit_ZeroBias(thisDimuonCand, muons, zeroBiasMuons, bFieldHandle, bs, thePrimaryV);
				//cross check: 3 events (MuOnia, ZeroBias, ZeroBias) mixing
				if (zeroBiasMuons_previousEvent.size()>2) fourMuonMixFit_ZeroBiasTwoEvts(thisDimuonCand, zeroBiasMuons_previousEvent.at(zeroBiasMuons_previousEvent.size()-3), zeroBiasMuons_previousEvent.at(zeroBiasMuons_previousEvent.size()-2), zeroBiasMuons, bFieldHandle, bs, thePrimaryV);
				//No mixing. 4 muon fit in the original MuOnia event
				fourMuonFit(thisDimuonCand, muons, bFieldHandle, bs, thePrimaryV);
			}

			if (fabs(mumu_vFit_noMC->currentState().mass()-9.46)<bestYMass) {
				bestYMass=fabs(mumu_vFit_noMC->currentState().mass()-9.46);
				bestVertex_and_bestYMass = 0;
				if (nGoodUpsilonCand==1) bestVertex_and_bestYMass = 1;
				DimuonCand_bestYMass = *dimuonCand;
				fillUpsilonBestMass(mumuVertexFitTree,DimuonCand_bestYMass,bFieldHandle,bs);
			}
		} //end of Upsilon loop

	}

	if (nGoodUpsilonCand>0) onia_tree->Fill();

}


void  MixingRootupler::fillUpsilonBestVertex(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs) {
	mumuVertexFitTree->movePointerToTheTop();     
	RefCountedKinematicParticle mumu_vFit_noMC = mumuVertexFitTree->currentParticle();    
	RefCountedKinematicVertex mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex(); //fitted vertex is same as the commonVertex in the Onia2MuMu skim
	//KinematicParameters mymumupara=  mumu_vFit_noMC->currentState().kinematicParameters();
	//cout<<"mymumupara px="<<mymumupara.momentum().x()<<",py="<<mymumupara.momentum().y()<<", m="<<mumu_vFit_noMC->currentState().mass()<<endl;      
	//float mymumuonlyctau=GetcTau(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);
	//float mymumuonlyctauerr=GetcTauErr(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);     
	//cout<<"mymumuonlyctau="<<mymumuonlyctau<<endl;
	//mumufit_Ctau->push_back( mymumuonlyctau );
	//mumufit_Ctauerr->push_back( mymumuonlyctauerr );
	mumufit_Mass = mumu_vFit_noMC->currentState().mass();
	mumufit_MassErr = sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) )  ;
	mumufit_VtxCL = ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) ;
	mumufit_VtxCL2 = mumu_vFit_vertex_noMC->chiSquared() ;
	mumufit_DecayVtxX = mumu_vFit_vertex_noMC->position().x() ;
	mumufit_DecayVtxY = mumu_vFit_vertex_noMC->position().y() ;
	mumufit_DecayVtxZ = mumu_vFit_vertex_noMC->position().z() ;
	mumufit_DecayVtxXE = mumu_vFit_vertex_noMC->error().cxx() ;
	mumufit_DecayVtxYE = mumu_vFit_vertex_noMC->error().cyy() ;
	mumufit_DecayVtxZE = mumu_vFit_vertex_noMC->error().czz() ;
	mumufit_p4.SetXYZM( mumu_vFit_noMC->currentState().globalMomentum().x(), mumu_vFit_noMC->currentState().globalMomentum().y(), mumu_vFit_noMC->currentState().globalMomentum().z(), mumufit_Mass ); 
	//mumufit_mupIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon1));
	//mumufit_mumIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon2));
	//std::cout<<"mass0="<<dimuonCand->mass()<<", mass1="<<mumufit_Mass<<std::endl;
	//std::cout<<"vProb0="<<vProb<<", vProb1="<<mumufit_VtxCL<<std::endl;

	//raw dimuon and muon
	dimuon_p4.SetPtEtaPhiM(dimuonCand.pt(), dimuonCand.eta(), dimuonCand.phi(), dimuonCand.mass());
	reco::Candidate::LorentzVector vP = dimuonCand.daughter("muon1")->p4();
	reco::Candidate::LorentzVector vM = dimuonCand.daughter("muon2")->p4();
	//std::cout<<"muon charge"<<dimuonCand->daughter("muon1")->charge()<<" "<<dimuonCand->daughter("muon2")->charge()<<std::endl;
	//if ( dimuonCand->daughter("muon1")->charge() < 0) {
	// vP = dimuonCand->daughter("muon2")->p4();
	// vM = dimuonCand->daughter("muon1")->p4();
	//}
	mu1_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
	mu2_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
	mu1Charge = dimuonCand.daughter("muon1")->charge();
	mu2Charge = dimuonCand.daughter("muon2")->charge();

	reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
	reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
	reco::TrackTransientTrack muon1TTT(muTrack1_ref, &(*bFieldHandle));
	reco::TrackTransientTrack muon2TTT(muTrack2_ref, &(*bFieldHandle));
	mu1_d0 = -muon1TTT.dxy(bs);
	mu1_d0err = muon1TTT.d0Error();
	mu1_dz = muon1TTT.dz();
	mu1_dzerr = muon1TTT.dzError();
	mu2_d0 = -muon2TTT.dxy(bs);
	mu2_d0err = muon2TTT.d0Error();
	mu2_dz = muon2TTT.dz();
	mu2_dzerr = muon2TTT.dzError();
}


void  MixingRootupler::fillUpsilonBestMass(RefCountedKinematicTree mumuVertexFitTree, pat::CompositeCandidate dimuonCand, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs) {
	mumuVertexFitTree->movePointerToTheTop();
	RefCountedKinematicParticle mumu_vFit_noMC = mumuVertexFitTree->currentParticle();
	RefCountedKinematicVertex mumu_vFit_vertex_noMC = mumuVertexFitTree->currentDecayVertex(); //fitted vertex is same as the commonVertex in the Onia2MuMu skim
	//KinematicParameters mymumupara=  mumu_vFit_noMC->currentState().kinematicParameters();
	//cout<<"mymumupara px="<<mymumupara.momentum().x()<<",py="<<mymumupara.momentum().y()<<", m="<<mumu_vFit_noMC->currentState().mass()<<endl;      
	//float mymumuonlyctau=GetcTau(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);
	//float mymumuonlyctauerr=GetcTauErr(mumu_vFit_vertex_noMC,mumu_vFit_noMC,thePrimaryV,beamSpot);     
	//cout<<"mymumuonlyctau="<<mymumuonlyctau<<endl;
	//mumufit_Ctau->push_back( mymumuonlyctau );
	//mumufit_Ctauerr->push_back( mymumuonlyctauerr );
	mumufit_Mass_bestYMass = mumu_vFit_noMC->currentState().mass();
	mumufit_MassErr_bestYMass = sqrt( mumu_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6) )  ;
	mumufit_VtxCL_bestYMass = ChiSquaredProbability((double)(mumu_vFit_vertex_noMC->chiSquared()),(double)(mumu_vFit_vertex_noMC->degreesOfFreedom())) ;
	mumufit_VtxCL2_bestYMass = mumu_vFit_vertex_noMC->chiSquared() ;
	mumufit_DecayVtxX_bestYMass = mumu_vFit_vertex_noMC->position().x() ;
	mumufit_DecayVtxY_bestYMass = mumu_vFit_vertex_noMC->position().y() ;
	mumufit_DecayVtxZ_bestYMass = mumu_vFit_vertex_noMC->position().z() ;
	mumufit_DecayVtxXE_bestYMass = mumu_vFit_vertex_noMC->error().cxx() ;
	mumufit_DecayVtxYE_bestYMass = mumu_vFit_vertex_noMC->error().cyy() ;
	mumufit_DecayVtxZE_bestYMass = mumu_vFit_vertex_noMC->error().czz() ;
	mumufit_p4_bestYMass.SetXYZM( mumu_vFit_noMC->currentState().globalMomentum().x(), mumu_vFit_noMC->currentState().globalMomentum().y(), mumu_vFit_noMC->currentState().globalMomentum().z(), mumufit_Mass_bestYMass );
	//mumufit_mupIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon1));
	//mumufit_mumIdx->push_back( std::distance(thePATMuonHandle->begin(), iMuon2));
	//std::cout<<"mass0="<<dimuonCand->mass()<<", mass1="<<mumufit_Mass<<std::endl;
	//std::cout<<"vProb0="<<vProb<<", vProb1="<<mumufit_VtxCL<<std::endl;

	//raw dimuon and muon
	dimuon_p4_bestYMass.SetPtEtaPhiM(dimuonCand.pt(), dimuonCand.eta(), dimuonCand.phi(), dimuonCand.mass());
	reco::Candidate::LorentzVector vP = dimuonCand.daughter("muon1")->p4();
	reco::Candidate::LorentzVector vM = dimuonCand.daughter("muon2")->p4();
	//std::cout<<"muon charge"<<dimuonCand->daughter("muon1")->charge()<<" "<<dimuonCand->daughter("muon2")->charge()<<std::endl;
	//if ( dimuonCand->daughter("muon1")->charge() < 0) {
	// vP = dimuonCand->daughter("muon2")->p4();
	// vM = dimuonCand->daughter("muon1")->p4();
	//}
	mu1_p4_bestYMass.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
	mu2_p4_bestYMass.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
	mu1Charge_bestYMass = dimuonCand.daughter("muon1")->charge();
	mu2Charge_bestYMass = dimuonCand.daughter("muon2")->charge();

	reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
	reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
	reco::TrackTransientTrack muon1TTT(muTrack1_ref, &(*bFieldHandle));
	reco::TrackTransientTrack muon2TTT(muTrack2_ref, &(*bFieldHandle));
	mu1_d0_bestYMass = -muon1TTT.dxy(bs);
	mu1_d0err_bestYMass = muon1TTT.d0Error();
	mu1_dz_bestYMass = muon1TTT.dz();
	mu1_dzerr_bestYMass = muon1TTT.dzError();
	mu2_d0_bestYMass = -muon2TTT.dxy(bs);
	mu2_d0err_bestYMass = muon2TTT.d0Error();
	mu2_dz_bestYMass = muon2TTT.dz();
	mu2_dzerr_bestYMass = muon2TTT.dzError();
}


// ------------ method called once each job just before starting event loop  ------------
void MixingRootupler::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void MixingRootupler::endJob(const edm::Event & iEvent) {
}

// ------------ method called when starting to processes a run  ------------
void MixingRootupler::beginRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a run  ------------
void MixingRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
void MixingRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void MixingRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MixingRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

int MixingRootupler::mediumMuon(edm::View<pat::Muon>::const_iterator rmu) {
	int goodMediumMuon=0;

	bool goodGlob = rmu->isGlobalMuon() && 
		rmu->globalTrack()->normalizedChi2() < 3 && 
		rmu->combinedQuality().chi2LocalPosition < 12 && 
		rmu->combinedQuality().trkKink < 20; 
	if(  muon::isLooseMuon(*rmu) && 
			rmu->innerTrack()->validFraction() > 0.8 && 
			muon::segmentCompatibility(*rmu) > (goodGlob ? 0.303 : 0.451) 
	  ) goodMediumMuon=1;

	return goodMediumMuon;
}

int MixingRootupler::mediumMuon(std::vector<pat::Muon>::iterator rmu) {
	int goodMediumMuon=0;

	bool goodGlob = rmu->isGlobalMuon() &&   
		rmu->globalTrack()->normalizedChi2() < 3 &&
		rmu->combinedQuality().chi2LocalPosition < 12 &&
		rmu->combinedQuality().trkKink < 20;
	if(  muon::isLooseMuon(*rmu) &&              
			rmu->innerTrack()->validFraction() > 0.8 && 
			muon::segmentCompatibility(*rmu) > (goodGlob ? 0.303 : 0.451)  
	  ) goodMediumMuon=1;

	return goodMediumMuon;
}

int MixingRootupler::tightMuon(edm::View<pat::Muon>::const_iterator rmu, reco::Vertex vertex) {
	int goodTightMuon=0;

	if( rmu->isGlobalMuon()
			&& rmu->isPFMuon()
			&& rmu->globalTrack()->normalizedChi2()<10.0
			&& rmu->globalTrack()->hitPattern().numberOfValidMuonHits()>0
			&& rmu->numberOfMatchedStations()>1
			&& fabs(rmu->muonBestTrack()->dxy( vertex.position() ))<0.2
			&& fabs(rmu->muonBestTrack()->dz( vertex.position() )) < 0.5
			&& rmu->innerTrack()->hitPattern().numberOfValidPixelHits()>0
			&& rmu->track()->hitPattern().trackerLayersWithMeasurement()>5
	  )  goodTightMuon = 1;

	return goodTightMuon;
}

int MixingRootupler::tightMuon(std::vector<pat::Muon>::iterator rmu, reco::Vertex vertex) {
	int goodTightMuon=0;

	if( rmu->isGlobalMuon()
			&& rmu->isPFMuon()
			&& rmu->globalTrack()->normalizedChi2()<10.0
			&& rmu->globalTrack()->hitPattern().numberOfValidMuonHits()>0
			&& rmu->numberOfMatchedStations()>1
			&& fabs(rmu->muonBestTrack()->dxy( vertex.position() ))<0.2
			&& fabs(rmu->muonBestTrack()->dz( vertex.position() )) < 0.5
			&& rmu->innerTrack()->hitPattern().numberOfValidPixelHits()>0
			&& rmu->track()->hitPattern().trackerLayersWithMeasurement()>5
	  )  goodTightMuon = 1; 

	return goodTightMuon;
}



int MixingRootupler::fourMuonMixFit_ZeroBiasTwoEvts(pat::CompositeCandidate dimuonCand, std::vector<pat::Muon> muons_previousA_ZeroBias, std::vector<pat::Muon> muons_previousB_ZeroBias, const std::vector<pat::Muon>* muons_ZeroBias, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	int nGoodFourMuonMix=0;
	//fourMuFit_VtxProb_mix3evts_ZeroBias = -1;

	std::vector<pat::Muon> combinedMuons;
	for (std::vector<pat::Muon>::iterator muA = muons_previousA_ZeroBias.begin(); muA != muons_previousA_ZeroBias.end(); ++muA){
		combinedMuons.push_back(*muA);
	}
	for (std::vector<pat::Muon>::iterator muB = muons_previousB_ZeroBias.begin(); muB !=  muons_previousB_ZeroBias.end(); ++muB){
		combinedMuons.push_back(*muB);
	}


	//std::cout<<"muon size="<<muons_ZeroBias->size()<<", muon previous size="<<combinedMuons.size()<<std::endl;
	for (std::vector<pat::Muon>::const_iterator mu3 = muons_ZeroBias->begin(), mu3end = muons_ZeroBias->end(); mu3 != mu3end; ++mu3){
		//std::cout<<"outer loop: "<<mu3- muons_ZeroBias->begin()<<" "<<mu3end - mu3<<std::endl;
		for (std::vector<pat::Muon>::iterator mu4 = combinedMuons.begin(), mu4end=combinedMuons.end() ; mu4-combinedMuons.begin() < mu3end-1-mu3 && mu4!=mu4end; ++mu4) {
			//std::cout<<"inner loop: "<<mu4-combinedMuons.begin()<<std::endl;
			if (mu3->charge() == mu4->charge()) continue;
			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();
				if (fitFourMu->currentState().isValid() 
						//	&&  ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb_mix3evts_ZeroBias
					)
				{ //Get chib         
					fourMuFit_Mass_mix3evts_ZeroBias.push_back(fitFourMu->currentState().mass());
					TLorentzVector p4;
					p4.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fitFourMu->currentState().mass());
					fourMuFit_Pt_mix3evts_ZeroBias.push_back(p4.Pt());
					fourMuFit_Eta_mix3evts_ZeroBias.push_back(p4.Eta());
					fourMuFit_Phi_mix3evts_ZeroBias.push_back(p4.Phi());
					fourMuFit_VtxProb_mix3evts_ZeroBias.push_back(ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())));
				}
			}
		}
	}
	return nGoodFourMuonMix;
}

void MixingRootupler::fourMuonFit(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	std::vector<pat::Muon> theRestMuons;
	//fourMuFit_VtxProb = -1;
	//std::cout<<"mu1Index="<<mu1Index<<", mu2Index="<<mu2Index<<std::endl;	
	for (edm::View<pat::Muon>::const_iterator mu3 = muons->begin(), muend = muons->end(); mu3 != muend; ++mu3){
		if (mu3->pt()<2.0 || fabs(mu3->eta())>2.4)  continue;
		if (mu3-muons->begin() == dimuonCand.userInt("mu1Index"))  continue; 
		if (mu3-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		reco::GenParticleRef genMu3;
		if (isMC_) genMu3 = mu3->genParticleRef();
		if (!isMC_ || (isMC_ && !genMu3.isNonnull())) theRestMuons.push_back(*mu3);
		//std::cout<<"fill vector: "<<mu3->pt()<<std::endl;

		for(edm::View<pat::Muon>::const_iterator mu4 = mu3+1 ; mu4 != muend; ++mu4){ 
			if (mu4->pt()<2.0 || fabs(mu4->eta())>2.4)  continue;
			if (mu4-muons->begin() == dimuonCand.userInt("mu1Index")) continue; 
			if (mu4-muons->begin() == dimuonCand.userInt("mu2Index")) continue;
			reco::GenParticleRef genMu4;
			if (isMC_) genMu4 = mu4->genParticleRef();

			if (mu3->charge() == mu4->charge()) continue;
			/*if ( (tightMuon(muons->begin()+mu1Index)+
			  tightMuon(muons->begin()+mu2Index)+
			  tightMuon(mu3)+
			  tightMuon(mu4)) < 2 
			  ) continue;*/

			if (!isMC_ || (isMC_ && !genMu3.isNonnull() && !genMu4.isNonnull())) {
				TLorentzVector mu1p4,mu2p4,mu3p4,mu4p4,fourMup4;
				reco::Candidate::LorentzVector v1 = dimuonCand.daughter("muon1")->p4();
				reco::Candidate::LorentzVector v2 = dimuonCand.daughter("muon2")->p4();
				mu1p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
				mu2p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
				mu3p4.SetPtEtaPhiM(mu3->pt(), mu3->eta(), mu3->phi(), mu3->mass());
				mu4p4.SetPtEtaPhiM(mu4->pt(), mu4->eta(), mu4->phi(), mu4->mass());
				fourMup4 = mu1p4 + mu2p4 + mu3p4 + mu4p4;
				fourMuFit_Mass_allComb.push_back(fourMup4.M());
			}

			//std::cout<<"found good mu3mu4: "<<mu3->pt()<<" "<<mu4->pt()<<", mu1: "<<muon1TT.track().pt()<<", mu2: "<<muon2TT.track().pt()<<std::endl;
			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));

			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{    
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();

				if (fitFourMu->currentState().isValid()
						// &&	ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb)
					)
					{ //Get chib         
						fourMuFit_Mass.push_back(fitFourMu->currentState().mass());
						fourMuFit_MassErr.push_back(sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6)));
						TLorentzVector p4;
						p4.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fitFourMu->currentState().mass());
						fourMuFit_Pt.push_back(p4.Pt());
						fourMuFit_Eta.push_back(p4.Eta());
						fourMuFit_Phi.push_back(p4.Phi());
						fourMuFit_VtxX.push_back(FourMuDecayVertex->position().x());
						fourMuFit_VtxY.push_back(FourMuDecayVertex->position().y());
						fourMuFit_VtxZ.push_back(FourMuDecayVertex->position().z());
						fourMuFit_VtxProb.push_back(ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())));
						fourMuFit_Chi2.push_back(FourMuDecayVertex->chiSquared());
						fourMuFit_ndof.push_back(FourMuDecayVertex->degreesOfFreedom());

						//get first muon
						bool child = fourMuTree->movePointerToTheFirstChild();
						RefCountedKinematicParticle fitMu1 = fourMuTree->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu1->currentState().kinematicParameters().momentum().x(), fitMu1->currentState().kinematicParameters().momentum().y(), fitMu1->currentState().kinematicParameters().momentum().z(), fitMu1->currentState().mass() );
						fourMuFit_mu1Pt.push_back(p4.Pt());
						fourMuFit_mu1Eta.push_back(p4.Eta());
						fourMuFit_mu1Phi.push_back(p4.Phi());
						fourMuFit_mu1E.push_back(p4.E());

						//get second muon
						child = fourMuTree->movePointerToTheNextChild();
						RefCountedKinematicParticle fitMu2 = fourMuTree->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu2->currentState().kinematicParameters().momentum().x(), fitMu2->currentState().kinematicParameters().momentum().y(), fitMu2->currentState().kinematicParameters().momentum().z(), fitMu2->currentState().mass() );
						fourMuFit_mu2Pt.push_back(p4.Pt());
						fourMuFit_mu2Eta.push_back(p4.Eta());
						fourMuFit_mu2Phi.push_back(p4.Phi());
						fourMuFit_mu2E.push_back(p4.E());

						//get third muon
						child = fourMuTree->movePointerToTheNextChild();
						RefCountedKinematicParticle fitMu3 = fourMuTree->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu3->currentState().kinematicParameters().momentum().x(), fitMu3->currentState().kinematicParameters().momentum().y(), fitMu3->currentState().kinematicParameters().momentum().z(), fitMu3->currentState().mass() );
						fourMuFit_mu3Pt.push_back(p4.Pt());
						fourMuFit_mu3Eta.push_back(p4.Eta());
						fourMuFit_mu3Phi.push_back(p4.Phi());
						fourMuFit_mu3E.push_back(p4.E());

						//get fourth muon
						child = fourMuTree->movePointerToTheNextChild();
						RefCountedKinematicParticle fitMu4 = fourMuTree->currentParticle();
						if(!child) break;
						p4.SetXYZM( fitMu4->currentState().kinematicParameters().momentum().x(), fitMu4->currentState().kinematicParameters().momentum().y(), fitMu4->currentState().kinematicParameters().momentum().z(), fitMu4->currentState().mass() );
						fourMuFit_mu4Pt.push_back(p4.Pt());
						fourMuFit_mu4Eta.push_back(p4.Eta());
						fourMuFit_mu4Phi.push_back(p4.Phi());
						fourMuFit_mu4E.push_back(p4.E());

						//std::cout<<fourMuFit_Mass<<" "<<(fourMuFit_mu1p4+fourMuFit_mu2p4+fourMuFit_mu3p4+fourMuFit_mu4p4).M()<<std::endl;

						mu3_Pt.push_back(mu3->pt());
						mu3_Eta.push_back(mu3->eta());
						mu3_Phi.push_back(mu3->phi());
						mu3_E.push_back(mu3->energy());
						mu4_Pt.push_back(mu4->pt());
						mu4_Eta.push_back(mu4->eta());
						mu4_Phi.push_back(mu4->phi());
						mu4_E.push_back(mu4->energy());
						reco::TrackTransientTrack muon3TTT(muTrack3_ref, &(*bFieldHandle));
						reco::TrackTransientTrack muon4TTT(muTrack4_ref, &(*bFieldHandle));
						mu3_d0.push_back(-muon3TTT.dxy(bs));
						mu3_d0err.push_back(muon3TTT.d0Error());
						mu3_dz.push_back(muon3TTT.dz());
						mu3_dzerr.push_back(muon3TTT.dzError());
						mu4_d0.push_back(-muon4TTT.dxy(bs));
						mu4_d0err.push_back(muon4TTT.d0Error());
						mu4_dz.push_back(muon4TTT.dz());
						mu4_dzerr.push_back(muon4TTT.dzError());
						mu3Charge.push_back(mu3->charge());
						mu4Charge.push_back(mu4->charge());
						mu1_Tight.push_back(tightMuon(muons->begin()+dimuonCand.userInt("mu1Index"),thePrimaryV));
						mu2_Tight.push_back(tightMuon(muons->begin()+dimuonCand.userInt("mu2Index"),thePrimaryV));
						mu3_Tight.push_back(tightMuon(mu3,thePrimaryV));
						mu4_Tight.push_back(tightMuon(mu4,thePrimaryV));
						mu1_Medium.push_back(mediumMuon(muons->begin()+dimuonCand.userInt("mu1Index")));
						mu2_Medium.push_back(mediumMuon(muons->begin()+dimuonCand.userInt("mu2Index")));
						mu3_Medium.push_back(mediumMuon(mu3));
						mu4_Medium.push_back(mediumMuon(mu4));

						if (isMC_) { 
							reco::GenParticleRef genMu1 = (muons->begin()+dimuonCand.userInt("mu1Index"))->genParticleRef();
							reco::GenParticleRef genMu2 = (muons->begin()+dimuonCand.userInt("mu2Index"))->genParticleRef();
							//if (genMu1->motherRef()==genMu2->motherRef()) std::cout<<"genMu1->motherRef()->pdgId()="<<genMu1->motherRef()->pdgId()<<std::endl;
							//else std::cout<<"genMu1->motherRef()->pdgId()="<<genMu1->motherRef()->pdgId()<<", genMu2->motherRef()->pdgId()="<<genMu2->motherRef()->pdgId()
							//	<<"genMu1->motherRef()->motherRef()->pdgId()="<<genMu1->motherRef()->motherRef()->pdgId()<<", genMu2->motherRef()->motherRef()->pdgId()="<<genMu2->motherRef()->motherRef()->pdgId()<<std::endl;

							//reco::GenParticleRef genMu3 = mu3->genParticleRef();
							//reco::GenParticleRef genMu4 = mu4->genParticleRef();
							if (genMu3.isNonnull() ){
								mu3_pdgID.push_back(genMu3->pdgId());
								/*
									if (genMu3->numberOfMothers()>0){ 
									reco::GenParticleRef mom3 = genMu3->motherRef();
									if (mom3.isNonnull()) { 
									std::cout<<""<<"mom pdgID= "<<mom3->pdgId()<<std::endl;
									if (mom3==genMu1->motherRef()) std::cout<<"same mother"<<std::endl;
									}    
									else std::cout<<"mom non"<<std::endl;
									}    
									else std::cout<<"# mom = 0"<<std::endl;
									*/
							}
							else mu3_pdgID.push_back(0);
							if (genMu4.isNonnull() ) mu4_pdgID.push_back(genMu4->pdgId());
							else mu4_pdgID.push_back(0);
						}    
					}
			}
		}
	}
	muons_previousEvent.push_back(theRestMuons);
}


int MixingRootupler::fourMuonMixFit_ZeroBias(pat::CompositeCandidate dimuonCand, edm::Handle< edm::View<pat::Muon> > muons, const std::vector<pat::Muon>* muons_ZeroBias, edm::ESHandle<MagneticField> bFieldHandle, reco::BeamSpot bs, reco::Vertex thePrimaryV){
	int nGoodFourMuonMix=0;

	// prepare a vectpr which contains: 1. all muons from a Zerobias event; 2. muons from a MuONia event but excluded the 2 muons from Upsilon
	std::vector<pat::Muon> combinedMuons;
	for (std::vector<pat::Muon>::const_iterator mu3 = muons_ZeroBias->begin(), muend = muons_ZeroBias->end(); mu3 != muend; ++mu3){
		combinedMuons.push_back(*mu3);
		/*for((edm::View<pat::Muon>::const_iterator  mu4 = muons->begin() ; mu4 != muons->end(); ++mu4){
		  if (mu4->pt()<2.0 || fabs(mu4->eta())>2.4) continue;
		  if (mu4-muons->begin() == dimuonCand.userInt("mu1Index"))  continue;
		  if (mu4-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		  if (mu3->charge() == mu4->charge()) continue;

		  TLorentzVector mu1p4,mu2p4,mu3p4,mu4p4,fourMup4;
		  reco::Candidate::LorentzVector v1 = dimuonCand.daughter("muon1")->p4();
		  reco::Candidate::LorentzVector v2 = dimuonCand.daughter("muon2")->p4();
		  mu1p4.SetPtEtaPhiM(v1.pt(),v1.eta(),v1.phi(),v1.mass());
		  mu2p4.SetPtEtaPhiM(v2.pt(),v2.eta(),v2.phi(),v2.mass());
		  mu3p4.SetPtEtaPhiM(mu3->pt(), mu3->eta(), mu3->phi(), mu3->mass());
		  mu4p4.SetPtEtaPhiM(mu4->pt(), mu4->eta(), mu4->phi(), mu4->mass());
		  fourMup4 = mu1p4 + mu2p4 + mu3p4 + mu4p4;
		  fourMuFit_Mass_allComb_mix_ZeroBias.push_back(fourMup4.M());
		  }*/
	}
	//std::cout<<"ZeroBias muons Size="<<muons_ZeroBias->size()<<", MuOnia muons Size="<<muons->size()<<std::endl;
	int muons_ZeroBias_originalSize=muons_ZeroBias->size();
	for (edm::View<pat::Muon>::const_iterator mu = muons->begin(); mu !=  muons->end(); ++mu){
		if (mu->pt()<2.0 || fabs(mu->eta())>2.4) continue;
		if (mu-muons->begin() == dimuonCand.userInt("mu1Index"))  continue;
		if (mu-muons->begin() == dimuonCand.userInt("mu2Index"))  continue;
		reco::GenParticleRef genMu;
		if (isMC_) genMu = mu->genParticleRef();
		if (!isMC_ || (isMC_ && !genMu.isNonnull())) {
			combinedMuons.push_back(*mu);
		}
	}

	//fourMuFit_VtxProb_mix_ZeroBias = -1;
	//std::cout<<"combinedSize="<<combinedMuons.size()<<std::endl;
	for (std::vector<pat::Muon>::iterator mu3 = combinedMuons.begin(), muend = combinedMuons.end(); mu3-combinedMuons.begin() < muons_ZeroBias_originalSize; ++mu3){//loop muons from ZeroBias
		//std::cout<<"mu3-combinedMuons.begin() = "<<mu3-combinedMuons.begin()<<", muons_ZeroBias_originalSize="<<muons_ZeroBias_originalSize<<std::endl;

		// Option 1: do 3+1 only (pick 3 muons from MuOnia and 1 muon from ZeroBias). This is the most common scenario
		for (std::vector<pat::Muon>::iterator mu4 = combinedMuons.begin()+muons_ZeroBias_originalSize; mu4 != muend; ++mu4) //loop muons from MuOnia

			// Option 2: Do 3+1 and 2+2
			//for (std::vector<pat::Muon>::iterator mu4 = mu3+1 ; mu4 != muend; ++mu4)  //loop muons from both ZeroBias and MuOnia

			// Option 3: Do 2+2 only 
			//for (std::vector<pat::Muon>::iterator mu4 = mu3+1 ; mu4-combinedMuons.begin() < muons_ZeroBias_originalSize; ++mu4) //loop muons from ZeroBias

		{
			//std::cout<<"   small loop:"<< mu4-combinedMuons.begin()<<std::endl;
			if (mu3->charge() == mu4->charge()) continue;

			reco::TrackRef muTrack1_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon1") ) )->innerTrack();
			reco::TrackRef muTrack2_ref = ( dynamic_cast<const pat::Muon*>(dimuonCand.daughter("muon2") ) )->innerTrack();
			reco::TrackRef muTrack3_ref = mu3->track();
			reco::TrackRef muTrack4_ref = mu4->track();
			reco::TransientTrack muon1TT(muTrack1_ref, &(*bFieldHandle));
			reco::TransientTrack muon2TT(muTrack2_ref, &(*bFieldHandle));
			reco::TransientTrack muon3TT(muTrack3_ref, &(*bFieldHandle));
			reco::TransientTrack muon4TT(muTrack4_ref, &(*bFieldHandle));
			KinematicParticleFactoryFromTransientTrack fourMuFactory;
			std::vector<RefCountedKinematicParticle> fourMuParticles;
			fourMuParticles.push_back(fourMuFactory.particle (muon1TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon2TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon3TT, muonMass, float(0), float(0), muonSigma));
			fourMuParticles.push_back(fourMuFactory.particle (muon4TT, muonMass, float(0), float(0), muonSigma));
			KinematicConstrainedVertexFitter constVertexFitter;
			//fit w/ mass constraint
			//MultiTrackKinematicConstraint *upsilon_mtc = new  TwoTrackMassKinematicConstraint(upsilon_mass_);
			//RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles,upsilon_mtc);

			//fit w/o any mass constraint
			RefCountedKinematicTree fourMuTree = constVertexFitter.fit(fourMuParticles);

			if(!fourMuTree->isEmpty())
			{
				fourMuTree->movePointerToTheTop();
				RefCountedKinematicParticle fitFourMu = fourMuTree->currentParticle();
				RefCountedKinematicVertex FourMuDecayVertex = fourMuTree->currentDecayVertex();
				if (fitFourMu->currentState().isValid()
						//   &&  ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())) > fourMuFit_VtxProb_mix_ZeroBias
					)
				{ //Get chib         
					nGoodFourMuonMix = 1;
					fourMuFit_Mass_mix_ZeroBias.push_back(fitFourMu->currentState().mass());
					fourMuFit_MassErr_mix_ZeroBias.push_back(sqrt(fitFourMu->currentState().kinematicParametersError().matrix()(6,6)));
					TLorentzVector p4;
					p4.SetXYZM(fitFourMu->currentState().kinematicParameters().momentum().x(),fitFourMu->currentState().kinematicParameters().momentum().y(),fitFourMu->currentState().kinematicParameters().momentum().z(),fitFourMu->currentState().mass());
					fourMuFit_Pt_mix_ZeroBias.push_back(p4.Pt());
					fourMuFit_Eta_mix_ZeroBias.push_back(p4.Eta());
					fourMuFit_Phi_mix_ZeroBias.push_back(p4.Phi());
					fourMuFit_VtxX_mix_ZeroBias.push_back(FourMuDecayVertex->position().x());
					fourMuFit_VtxY_mix_ZeroBias.push_back(FourMuDecayVertex->position().y());
					fourMuFit_VtxZ_mix_ZeroBias.push_back(FourMuDecayVertex->position().z());
					fourMuFit_VtxProb_mix_ZeroBias.push_back(ChiSquaredProbability((double)(FourMuDecayVertex->chiSquared()),(double)(FourMuDecayVertex->degreesOfFreedom())));
					fourMuFit_Chi2_mix_ZeroBias.push_back(FourMuDecayVertex->chiSquared());
					fourMuFit_ndof_mix_ZeroBias.push_back(FourMuDecayVertex->degreesOfFreedom());
					if (mu4-combinedMuons.begin() >= muons_ZeroBias_originalSize) fourMuFit_3plus1_mix_ZeroBias.push_back(3);
					else fourMuFit_3plus1_mix_ZeroBias.push_back(2);

					//get first muon
					bool child = fourMuTree->movePointerToTheFirstChild();
					RefCountedKinematicParticle fitMu1 = fourMuTree->currentParticle();
					if(!child) break;
					p4.SetXYZM( fitMu1->currentState().kinematicParameters().momentum().x(), fitMu1->currentState().kinematicParameters().momentum().y(), fitMu1->currentState().kinematicParameters().momentum().z(), fitMu1->currentState().mass() );
					fourMuFit_mu1Pt_mix_ZeroBias.push_back(p4.Pt());
					fourMuFit_mu1Eta_mix_ZeroBias.push_back(p4.Eta());
					fourMuFit_mu1Phi_mix_ZeroBias.push_back(p4.Phi());
					fourMuFit_mu1E_mix_ZeroBias.push_back(p4.E());

					//get second muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu2 = fourMuTree->currentParticle();
					if(!child) break;
					p4.SetXYZM( fitMu2->currentState().kinematicParameters().momentum().x(), fitMu2->currentState().kinematicParameters().momentum().y(), fitMu2->currentState().kinematicParameters().momentum().z(), fitMu2->currentState().mass() );
					fourMuFit_mu2Pt_mix_ZeroBias.push_back(p4.Pt());
					fourMuFit_mu2Eta_mix_ZeroBias.push_back(p4.Eta());
					fourMuFit_mu2Phi_mix_ZeroBias.push_back(p4.Phi());
					fourMuFit_mu2E_mix_ZeroBias.push_back(p4.E());

					//get third muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu3 = fourMuTree->currentParticle();
					if(!child) break;
					p4.SetXYZM( fitMu3->currentState().kinematicParameters().momentum().x(), fitMu3->currentState().kinematicParameters().momentum().y(), fitMu3->currentState().kinematicParameters().momentum().z(), fitMu3->currentState().mass() );
					fourMuFit_mu3Pt_mix_ZeroBias.push_back(p4.Pt());
					fourMuFit_mu3Eta_mix_ZeroBias.push_back(p4.Eta());
					fourMuFit_mu3Phi_mix_ZeroBias.push_back(p4.Phi());
					fourMuFit_mu3E_mix_ZeroBias.push_back(p4.E());

					//get fourth muon
					child = fourMuTree->movePointerToTheNextChild();
					RefCountedKinematicParticle fitMu4 = fourMuTree->currentParticle();
					if(!child) break;
					p4.SetXYZM( fitMu4->currentState().kinematicParameters().momentum().x(), fitMu4->currentState().kinematicParameters().momentum().y(), fitMu4->currentState().kinematicParameters().momentum().z(), fitMu4->currentState().mass() );
					fourMuFit_mu4Pt_mix_ZeroBias.push_back(p4.Pt());
					fourMuFit_mu4Eta_mix_ZeroBias.push_back(p4.Eta());
					fourMuFit_mu4Phi_mix_ZeroBias.push_back(p4.Phi());
					fourMuFit_mu4E_mix_ZeroBias.push_back(p4.E());

					//save original information of the Zerobias muons
					mu3_Pt_mix_ZeroBias.push_back(mu3->pt());
					mu3_Eta_mix_ZeroBias.push_back(mu3->eta());
					mu3_Phi_mix_ZeroBias.push_back(mu3->phi());
					mu3_E_mix_ZeroBias.push_back(mu3->energy());
					mu4_Pt_mix_ZeroBias.push_back(mu4->pt());
					mu4_Eta_mix_ZeroBias.push_back(mu4->eta());
					mu4_Phi_mix_ZeroBias.push_back(mu4->phi());
					mu4_E_mix_ZeroBias.push_back(mu4->energy());
					reco::TrackTransientTrack muon3TTT(muTrack3_ref, &(*bFieldHandle));
					reco::TrackTransientTrack muon4TTT(muTrack4_ref, &(*bFieldHandle));
					mu3_d0_mix_ZeroBias.push_back(-muon3TTT.dxy(bs));
					mu3_d0err_mix_ZeroBias.push_back(muon3TTT.d0Error());
					mu3_dz_mix_ZeroBias.push_back(muon3TTT.dz());
					mu3_dzerr_mix_ZeroBias.push_back(muon3TTT.dzError());
					mu4_d0_mix_ZeroBias.push_back(-muon4TTT.dxy(bs));
					mu4_d0err_mix_ZeroBias.push_back(muon4TTT.d0Error());
					mu4_dz_mix_ZeroBias.push_back(muon4TTT.dz());
					mu4_dzerr_mix_ZeroBias.push_back(muon4TTT.dzError());
					mu3Charge_mix_ZeroBias.push_back(mu3->charge());
					mu4Charge_mix_ZeroBias.push_back(mu4->charge());
					mu1_Tight_mix_ZeroBias.push_back(tightMuon(muons->begin()+dimuonCand.userInt("mu1Index"),thePrimaryV));
					mu2_Tight_mix_ZeroBias.push_back(tightMuon(muons->begin()+dimuonCand.userInt("mu2Index"),thePrimaryV));
					mu3_Tight_mix_ZeroBias.push_back(tightMuon(mu3,thePrimaryV));
					mu4_Tight_mix_ZeroBias.push_back(tightMuon(mu4,thePrimaryV));
					mu1_Medium_mix_ZeroBias.push_back(mediumMuon(muons->begin()+dimuonCand.userInt("mu1Index")));
					mu2_Medium_mix_ZeroBias.push_back(mediumMuon(muons->begin()+dimuonCand.userInt("mu2Index")));
					mu3_Medium_mix_ZeroBias.push_back(mediumMuon(mu3));
					mu4_Medium_mix_ZeroBias.push_back(mediumMuon(mu4));

				}
			}
		}
	}
	return nGoodFourMuonMix;
}


//define this as a plug-in
DEFINE_FWK_MODULE(MixingRootupler);
