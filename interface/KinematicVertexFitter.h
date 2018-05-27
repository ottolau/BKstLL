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

// FIXME! this is horrible, but it allows to handle different kinds of C++ objects at the same time
class KinematicVertexFitter {

  public:
    KinematicVertexFitter() {};
    virtual ~KinematicVertexFitter() {};

    reco::TransientTrack getTransientTrack(const reco::Track& track) {    
      reco::TransientTrack transientTrack(track, paramField);
      return transientTrack;
    }

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

  private:
    OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");
    // Insignificant mass sigma to avoid singularities in the covariance matrix.
    // initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered 
    float chi        = 0.;
    float ndf        = 0.;

};

// 
// class KinematicVertexFitterGSF {
// 
//   public:
//     KinematicVertexFitterGSF() {};
//     virtual ~KinematicVertexFitterGSF() {};
// 
// 
// 
// 	edm::ESHandle<TransientTrackBuilder> ttrack_builder;
// 	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttrack_builder);
// 
// 
// 
// 			std::vector<reco::TransientTrack> trks{
// 				ttrack_builder->build(trk), ttrack_builder->build(ele->gsfTrack())};
// 
// 
// 
//     reco::TransientTrack getGsfTransientTrack(const reco::GsfTrack& track) {    
//       reco::GsfTransientTrack gsfTransientTrack(track, paramField);
//       reco::TransientTrack    transientTrack(gsfTransientTrack.candidatePtr(), paramField)
//       return transientTrack;
//     }
// 
//     RefCountedKinematicTree Fit(const std::vector<pat::Electron> & electrons, const std::vector<pat::PackedCandidate> & packed){
// 
//       KinematicParticleFactoryFromTransientTrack pFactory;  
//       std::vector<RefCountedKinematicParticle> XParticles;
// 
//       // loop over the Electrons
//       for (std::vector<pat::Electron>::const_iterator iele = electrons.begin(); iele != electrons.end(); ++iele){
//         float pmass  = iele->mass();
//         float pmasse = 1.e-6 * pmass;
//         XParticles.push_back(pFactory.particle(getGsfTransientTrack( *(iele->gsfTrack()) ), pmass, chi, ndf, pmasse));
//       }
// 
//       // loop over the PackedCandidates, notice the different way to access the track
//       for (std::vector<pat::PackedCandidate>::const_iterator ipc = packed.begin(); ipc != packed.end(); ++ipc){
//         float pmass  = ipc->mass();
//         float pmasse = 1.e-6 * pmass;
//         XParticles.push_back(pFactory.particle(getTransientTrack( *(ipc->bestTrack()) ), pmass, chi, ndf, pmasse));
//       }
// 
//       KinematicConstrainedVertexFitter kvFitter;
//       RefCountedKinematicTree KinVtx = kvFitter.fit(XParticles); 
//       
//       return KinVtx;
//         
//     }
// 
//   private:
//     OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");
//     // Insignificant mass sigma to avoid singularities in the covariance matrix.
//     // initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered 
//     float chi        = 0.;
//     float ndf        = 0.;
// 
// };
// 


