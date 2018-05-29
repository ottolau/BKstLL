import ROOT
from itertools import product, combinations
import math
from copy import deepcopy as dc

from PhysicsTools.Heppy.analyzers.core.Analyzer      import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle    import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar             import deltaR, deltaR2, bestMatch
from PhysicsTools.Heppy.physicsobjects.Muon          import Muon
from PhysicsTools.Heppy.physicsobjects.Electron      import Electron
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject

from CMGTools.BKstLL.analyzers.utils import isAncestor, displacement2D, displacement3D, makeRecoVertex # utility functions
from CMGTools.BKstLL.physicsobjects.BKLL import BKLL


from pdb import set_trace

##########################################################################################
# load custom library to ROOT. This contains the kinematic vertex fitter class
ROOT.gSystem.Load('libCMGToolsBKstLL')
from ROOT import KinematicVertexFitter as VertexFitter

class BKJPsiEEAnalyzer(Analyzer):
    '''
    '''
    def __init__(self, *args, **kwargs):
        super(BKJPsiEEAnalyzer, self).__init__(*args, **kwargs)
        # stuff I need to instantiate only once
        self.vtxfit = VertexFitter()
        # create a std::vector<reco::Track> to be passed to the fitter 
        self.tofit_el = ROOT.std.vector('reco::Track')()
        # create a std::vector<pat::PackedCandidate> to be passed to the fitter 
        self.tofit_pc = ROOT.std.vector('pat::PackedCandidate')()

    def declareHandles(self):
        super(BKJPsiEEAnalyzer, self).declareHandles()

        # miniAOD collections
        self.handles['electrons' ] = AutoHandle('slimmedElectrons'             , 'std::vector<pat::Electron>'       )
        self.handles['muons'     ] = AutoHandle('slimmedMuons'                 , 'std::vector<pat::Muon>'           )
        self.handles['losttracks'] = AutoHandle('lostTracks'                   , 'std::vector<pat::PackedCandidate>')
        self.handles['pfcands'   ] = AutoHandle('packedPFCandidates'           , 'std::vector<pat::PackedCandidate>')
        self.handles['pvs'       ] = AutoHandle('offlineSlimmedPrimaryVertices', 'std::vector<reco::Vertex>'        )
        
        # enriched miniAOD collections
        self.handles['eletrackMap'] = AutoHandle(('ttk', 'eleTtkMap', 'TTK'), 'std::vector<pair<edm::Ptr<pat::Electron>,reco::Track> >')
        
        # AOD collections
        self.handles['beamspot'  ] = AutoHandle('offlineBeamSpot', 'reco::BeamSpot')

    def beginLoop(self, setup):
        super(BKJPsiEEAnalyzer, self).beginLoop(setup)
        self.counters.addCounter('BKJPsiEEAnalyzer')
        count = self.counters.counter('BKJPsiEEAnalyzer')
        count.register('all events')
        count.register('>= 2 electrons')
        count.register('diele mass < 6')
        count.register('diele dz < 1')
        count.register('good diele vtx')
        count.register('B(KLL) mass < 6')

    def process(self, event):
        self.readCollections(event.input)

        self.counters.counter('BKJPsiEEAnalyzer').inc('all events')

        # vertex stuff
        event.pvs         = self.handles['pvs'     ].product()
        event.beamspot    = self.handles['beamspot'].product()

        # get the tracks
        allpf      = map(PhysicsObject, self.handles['pfcands'   ].product())
        losttracks = map(PhysicsObject, self.handles['losttracks'].product())

        # merge the track collections
        event.alltracks = sorted([tt for tt in allpf + losttracks if tt.charge() != 0 and abs(tt.pdgId()) not in (11,13)], key = lambda x : x.pt(), reverse = True)
        
        # get the offline electrons and muons
        event.electrons = map(Electron, self.handles['electrons'].product())
        event.muons     = map(Muon    , self.handles['muons'    ].product())

        # get the electron-track map BEFORE selections
        eletrks = self.handles['eletrackMap'].product()

        # check consistency between electron-track map and electrons
        if len(event.electrons) != len(eletrks):
            print 'run %d \tlumi %d \tevent %d' %(event.run, event.lumi, event.eventId)
            print 'different number of electrons and ele tracks, %d and %d respectively' %(len(event.electrons), len(eletrks))
            for ii in event.electrons:
                print ii, ii.gsfTrack().pt()
            for ii in eletrks:
                print ii.first.pt(), ii.first.eta(), ii.first.phi(), ii.second.pt()
            set_trace()
        
        for jj in range(len(event.electrons)):
            event.electrons[jj].trackFromGsfTrack = eletrks[jj].second
            event.electrons[jj].ptr = eletrks[jj].first
            if event.electrons[jj].gsfTrack().pt() != eletrks[jj].second.pt():
                set_trace()

        # preselect electrons
        event.electrons = [ele for ele in event.electrons if self.testEle(ele)]

        # at least two electrons  
        if len(event.electrons)<2: return False
        self.counters.counter('BKJPsiEEAnalyzer').inc('>= 2 electrons')

        # build all di-ele pairs
        dieles = [(ele1, ele2) for ele1, ele2 in combinations(event.electrons, 2)]

        # opposite sign di-muon
        dieles = [(ele1, ele2) for ele1, ele2 in dieles if ele1.charge() != ele2.charge()]

        # SAVE THEM ALL, ALSO NON RESONANT!
        # invariant mass http://cmslxr.fnal.gov/source/DataFormats/EgammaCandidates/interface/GsfElectron.h#0795
        # dieles = [(ele1, ele2) for ele1, ele2 in dieles if (ele1.p4(1) + ele2.p4(1)).mass() < 6.]
        dieles = [(ele1, ele2) for ele1, ele2 in dieles if (ele1.p4() + ele2.p4()).mass() < 6.]
        if not len(dieles): return False
        self.counters.counter('BKJPsiEEAnalyzer').inc('diele mass < 6')

        # select by dz
        dieles = [(ele1, ele2) for ele1, ele2 in dieles if abs(ele1.vz() - ele2.vz()) < 1.]
        if not len(dieles): return False
        self.counters.counter('BKJPsiEEAnalyzer').inc('diele dz < 1')
        
        selDieles = []
        
        # dieles with a good vertex
        for diele in dieles:
            # clear the vectors
            self.tofit_el.clear()
            self.tofit_pc.clear()
            # create the vector of bla bla FIXME!
            for il in diele:
                # set_trace()
                self.tofit_el.push_back(il.trackFromGsfTrack)
            # fit it!
            # set_trace()
            svtree = self.vtxfit.Fit(self.tofit_el, self.tofit_pc) # actual vertex fitting
            # check that the vertex is good
            if not svtree.get().isEmpty() and svtree.get().isValid():
                svtree.movePointerToTheTop()
                sv = svtree.currentDecayVertex().get()
                recoSv = makeRecoVertex(sv, kinVtxTrkSize=2) # need to do some gymastics
                selDieles.append(diele)
        
        if len(selDieles)==0: return False
        self.counters.counter('BKJPsiEEAnalyzer').inc('good diele vtx')

        cands = []
        
        # set_trace()

        # add a track
        for diele, tk in product(selDieles, event.alltracks):
            # try to fit a vertex only if the track's close to the di-ele cand
            if max([abs(ee.vz() - tk.vz()) for ee in diele]) > 1.5:
                continue
            # remove double countings 
            if deltaR(tk, diele[0])<0.01 or deltaR(tk, diele[1])<0.01:
                continue
            # pt and sanity checks on the track
            # set_trace()
            if not tk.bestTrack():
                continue
            if tk.pt()<0.8:
                continue

            # kaon mass
            tk.setMass(0.493677)
            # clear the vectors
            self.tofit_el.clear()
            self.tofit_pc.clear()
            # create a RecoChargedCandidate for each reconstructed lepton and flush it into the vector
            for il in diele:
                self.tofit_el.push_back(il.trackFromGsfTrack)
            
            # push the track into the vector
            self.tofit_pc.push_back(tk.physObj)

            # fit it!
            svtree = self.vtxfit.Fit(self.tofit_el, self.tofit_pc) # actual vertex fitting
            # check that the vertex is good
            if not svtree.get().isEmpty() and svtree.get().isValid():
                svtree.movePointerToTheTop()
                sv = svtree.currentDecayVertex().get()
                recoSv = makeRecoVertex(sv, kinVtxTrkSize=3) # need to do some gymastics
                cands.append(BKLL(diele[0], diele[1], tk, recoSv, event.beamspot))

        if not len(cands): return False
        self.counters.counter('BKJPsiEEAnalyzer').inc('good diele vtx')

        cands = [cand for cand in cands if cand.b().mass()<6]
        if not len(cands): return False
        self.counters.counter('BKJPsiEEAnalyzer').inc('B(KLL) mass < 6')

#         set_trace()
#         for jj in cands: print jj.b().mass(), jj.b().pt(), jj.b().eta(), jj.b().phi(), jj.ll().mass(), jj.ls2d(), jj.vtxprob(), jj.llcone(), jj.bcone()

        event.myB = sorted(cands, key = lambda x : x.vtxprob(), reverse=True)[0]
        
        return True

    
    def testEle(self, ele):
        return ele.pt()>2.        and \
               abs(ele.eta())<2.5 
    
    
    