#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParametersError.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

// FIXME! this is horrible, but it allows to handle different kinds of C++ objects at the same time
class KinematicVertexFitter {

  public:
    KinematicVertexFitter() {};
    virtual ~KinematicVertexFitter() {};

    reco::TransientTrack getTransientTrack(const reco::Track& track) {    
      reco::TransientTrack transientTrack(track, paramField);
      return transientTrack;
    }

    reco::TransientTrack getTransientTrack(const reco::GsfTrack& track) {    
      reco::GsfTransientTrack gsfTransientTrack(track, paramField);
      reco::TransientTrack transientTrack(gsfTransientTrack.track(), paramField);
      return transientTrack;
    }

    bool checkIsLooseMuon(const reco::Muon & recoMu) {
      return muon::isLooseMuon(recoMu);
    }
  
    bool checkIsMediumMuon(const reco::Muon & recoMu) {
      return muon::isMediumMuon(recoMu);
    }

    bool checkIsTightMuon(const reco::Muon & recoMu, const reco::Vertex & vtx) {
      return muon::isTightMuon(recoMu, vtx);
    }

    double getLxy(const reco::Vertex & vtx1, const reco::Vertex & vtx2) {
      VertexDistanceXY vertTool;
      return vertTool.distance(vtx1, vtx2).value();
    }

    double getLxyError(const reco::Vertex & vtx1, const reco::Vertex & vtx2) {
      VertexDistanceXY vertTool;
      return vertTool.distance(vtx1, vtx2).error();
    }

    // charged candidates and packed candidates
    RefCountedKinematicTree Fit(const std::vector<reco::RecoChargedCandidate> & charged, const std::vector<pat::PackedCandidate> & packed){

      KinematicParticleFactoryFromTransientTrack pFactory;  
      std::vector<RefCountedKinematicParticle> XParticles;

      // loop over the RecoChargedCandidates
      for (std::vector<reco::RecoChargedCandidate>::const_iterator icc = charged.begin(); icc != charged.end(); ++icc){
        float pmass  = icc->mass();
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *(icc->track()) ), pmass, chi, ndf, pmasse));
      }

      // loop over the PackedCandidates, notice the different way to access the track
      for (std::vector<pat::PackedCandidate>::const_iterator ipc = packed.begin(); ipc != packed.end(); ++ipc){
        float pmass  = ipc->mass();
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *(ipc->bestTrack()) ), pmass, chi, ndf, pmasse));
      }

      KinematicConstrainedVertexFitter kvFitter;
      RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles); 
      
      return KinVtx;
        
    }

    // charged candidates and tracks
    RefCountedKinematicTree Fit(const std::vector<reco::RecoChargedCandidate> & charged, const std::vector<reco::Track> & packed, float trackmass){

      KinematicParticleFactoryFromTransientTrack pFactory;  
      std::vector<RefCountedKinematicParticle> XParticles;

      // loop over the RecoChargedCandidates
      for (std::vector<reco::RecoChargedCandidate>::const_iterator icc = charged.begin(); icc != charged.end(); ++icc){
        float pmass  = icc->mass();
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *(icc->track()) ), pmass, chi, ndf, pmasse));
      }

      // loop over the Tracks, notice the different way to access the track
      for (std::vector<reco::Track>::const_iterator ipc = packed.begin(); ipc != packed.end(); ++ipc){
        float pmass  = trackmass;
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *(ipc) ), pmass, chi, ndf, pmasse));
      }

      KinematicConstrainedVertexFitter kvFitter;
      RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles); 
      
      return KinVtx;
        
    }

    // tracks and tracks
    RefCountedKinematicTree Fit(const std::vector<reco::Track> & charged, const std::vector<reco::Track> & packed, float lepmass, float trackmass1, float trackmass2){

      KinematicParticleFactoryFromTransientTrack pFactory;  
      std::vector<RefCountedKinematicParticle> XParticles;

      // loop over the RecoChargedCandidates
      for (std::vector<reco::Track>::const_iterator icc = charged.begin(); icc != charged.end(); ++icc){
        float pmass  = lepmass;
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *(icc) ), pmass, chi, ndf, pmasse));
      }

      // loop over the Tracks, notice the different way to access the track
      for (std::vector<reco::Track>::const_iterator ipc = packed.begin(); ipc != packed.end(); ++ipc){
        float pmass  = ((ipc - packed.begin()) == 0 ) ? trackmass1:trackmass2;
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *(ipc) ), pmass, chi, ndf, pmasse));
      }

      KinematicConstrainedVertexFitter kvFitter;
      RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles); 
      
      return KinVtx;
        
    }

    // gsftracks and tracks
    RefCountedKinematicTree Fit(const std::vector<reco::GsfTrack> & charged, const std::vector<reco::Track> & packed, float trackmass){

      KinematicParticleFactoryFromTransientTrack pFactory;  
      std::vector<RefCountedKinematicParticle> XParticles;

      // loop over the RecoChargedCandidates
      for (std::vector<reco::GsfTrack>::const_iterator icc = charged.begin(); icc != charged.end(); ++icc){
        float pmass  = 0.0005109989461;
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *(icc) ), pmass, chi, ndf, pmasse));
      }

      // loop over the Tracks, notice the different way to access the track
      for (std::vector<reco::Track>::const_iterator ipc = packed.begin(); ipc != packed.end(); ++ipc){
        float pmass  = trackmass;
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *(ipc) ), pmass, chi, ndf, pmasse));
      }

      KinematicConstrainedVertexFitter kvFitter;
      RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles); 
      
      return KinVtx;
        
    }


    // pat::Electrons and packed candidates
    RefCountedKinematicTree Fit(std::vector<reco::Track> & eletracks, const std::vector<pat::PackedCandidate> & packed){

      KinematicParticleFactoryFromTransientTrack pFactory;  
      std::vector<RefCountedKinematicParticle> XParticles;

      // loop over the RecoChargedCandidates
      for (std::vector<reco::Track>::const_iterator iepair = eletracks.begin(); iepair != eletracks.end(); ++iepair){
        float pmass  = 0.0005109989461;
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *iepair ), pmass, chi, ndf, pmasse));
      }

      // loop over the PackedCandidates, notice the different way to access the track
      for (std::vector<pat::PackedCandidate>::const_iterator ipc = packed.begin(); ipc != packed.end(); ++ipc){
        float pmass  = ipc->mass();
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *(ipc->bestTrack()) ), pmass, chi, ndf, pmasse));
      }

      KinematicConstrainedVertexFitter kvFitter;
      RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles); 
      
      return KinVtx;
        
    }


    // packed candidates
    RefCountedKinematicTree Fit(const std::vector<pat::PackedCandidate> & packed){

      KinematicParticleFactoryFromTransientTrack pFactory;  
      std::vector<RefCountedKinematicParticle> XParticles;

      // loop over the PackedCandidates, notice the different way to access the track
      for (std::vector<pat::PackedCandidate>::const_iterator ipc = packed.begin(); ipc != packed.end(); ++ipc){
        float pmass  = ipc->mass();
        float pmasse = 1.e-6 * pmass;
        XParticles.push_back(pFactory.particle(getTransientTrack( *(ipc->bestTrack()) ), pmass, chi, ndf, pmasse));
      }

      KinematicConstrainedVertexFitter kvFitter;
      RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles); 
      
      return KinVtx;
        
    }


//*/    

  private:
    OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");
    // Insignificant mass sigma to avoid singularities in the covariance matrix.
    // initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered 
    float chi        = 0.;
    float ndf        = 0.;

};

