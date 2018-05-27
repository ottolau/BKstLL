#define G__DICTIONARY

#include "DataFormats/Common/interface/Wrapper.h"
#include "CMGTools/BKstLL/interface/KinematicVertexFitter.h"
#include "CMGTools/BKstLL/plugins/AddElectronTransientTrack.h"

namespace {
  struct CMG_BKstLL {
    KinematicVertexFitter KinVtx_;

//     std::pair<edm::Ptr<pat::Electron>,reco::TransientTrack> pppettk;
//     edm::Wrapper<std::pair<edm::Ptr<pat::Electron>,reco::TransientTrack> > wpppettk;
// 
//     std::vector<std::pair<edm::Ptr<pat::Electron>,reco::TransientTrack> > vpppettk;
//     edm::Wrapper<std::vector<std::pair<edm::Ptr<pat::Electron>,reco::TransientTrack> > > wvpppettk;

    std::pair<edm::Ptr<pat::Electron>,reco::Track> pppettk;
    edm::Wrapper<std::pair<edm::Ptr<pat::Electron>,reco::Track> > wpppettk;

    std::vector<std::pair<edm::Ptr<pat::Electron>,reco::Track> > vpppettk;
    edm::Wrapper<std::vector<std::pair<edm::Ptr<pat::Electron>,reco::Track> > > wvpppettk;

  };
}
