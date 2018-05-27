import ROOT
from itertools import product, combinations
import math

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
    def declareHandles(self):
        super(BKJPsiEEAnalyzer, self).declareHandles()

        # miniAOD collections
#         self.handles['electrons' ] = AutoHandle('slimmedElectrons'                           , 'std::vector<pat::Electron>'                    )
#         self.handles['muons'     ] = AutoHandle('slimmedMuons'                               , 'std::vector<pat::Muon>'                        )
#         self.handles['losttracks'] = AutoHandle('lostTracks'                                 , 'std::vector<pat::PackedCandidate>'             )
#         self.handles['pfcands'   ] = AutoHandle('packedPFCandidates'                         , 'std::vector<pat::PackedCandidate>'             )
#         self.handles['pvs'       ] = AutoHandle(('offlineSlimmedPrimaryVertices', '', 'PAT' ), 'std::vector<reco::Vertex>'                     )

        # AOD collections
        self.handles['muons'     ] = AutoHandle(('muons'                 , '', 'RECO'), 'vector<reco::Muon>'       )
        self.handles['electrons' ] = AutoHandle(('gedGsfElectrons'       , '', 'RECO'), 'vector<reco::GsfElectron>')
        self.handles['tracks'    ] = AutoHandle(('generalTracks'         , '', 'RECO'), 'vector<reco::Track>'      )
        self.handles['pvs'       ] = AutoHandle(('offlinePrimaryVertices', '', 'RECO'), 'vector<reco::Vertex>'     )

        self.handles['beamspot'  ] = AutoHandle(('offlineBeamSpot'       , '', 'RECO'), 'reco::BeamSpot'           )


    def beginLoop(self, setup):
        super(BKJPsiEEAnalyzer, self).beginLoop(setup)
        self.counters.addCounter('BKJPsiEEAnalyzer')
        count = self.counters.counter('BKJPsiEEAnalyzer')
        count.register('all events')
        count.register('>= 2 electrons')
        count.register('diele mass < 6')
        count.register('good diele vtx')
        count.register('B(KLL) mass < 6')

        # stuff I need to instantiate only once
        self.vtxfit = VertexFitter()
        # create a std::vector<reco::RecoChargedCandidate> to be passed to the fitter 
        self.tofit_cc = ROOT.std.vector('reco::RecoChargedCandidate')()
        # create a std::vector<pat::PackedCandidate> to be passed to the fitter 
        self.tofit_pc = ROOT.std.vector('pat::PackedCandidate')()

    def process(self, event):
        self.readCollections(event.input)

        self.counters.counter('BKJPsiEEAnalyzer').inc('all events')

        # vertex stuff
        event.pvs         = self.handles['pvs'     ].product()
        event.beamspot    = self.handles['beamspot'].product()

        # get the tracks
#         allpf      = map(PhysicsObject, self.handles['pfcands'   ].product())
#         losttracks = map(PhysicsObject, self.handles['losttracks'].product())

        # merge the track collections
#         event.alltracks = sorted([tt for tt in allpf + losttracks if tt.charge() != 0 and abs(tt.pdgId()) not in (11,13)], key = lambda x : x.pt(), reverse = True)
        
        # get the offline electrons and muons
#         event.electrons = map(Electron, self.handles['electrons'].product())
#         event.muons     = map(Muon    , self.handles['muons'    ].product())
        
        # preselect electrons
#         event.electrons = [ele for ele in event.electrons if self.testEle(ele)]
        
  
        # AOD tracks
        tracks    = self.handles['tracks'].product()
        seltracks = [tk for tk in tracks if tk.pt()>0.8 and abs(tk.eta())<2.5]

        # AOD electrons
        electrons    = self.handles['electrons'].product()
        selelectrons = [el for el in electrons if self.testEle(el)]
        
        if len(selelectrons)<2: return False
        self.counters.counter('BKJPsiEEAnalyzer').inc('>= 2 electrons')

        # build all di-ele pairs
        dieles = [(ele1, ele2) for ele1, ele2 in combinations(selelectrons, 2)]

        # opposite sign di-muon
        dieles = [(ele1, ele2) for ele1, ele2 in dieles if ele1.charge() != ele2.charge()]

        # SAVE THEM ALL, ALSO NON RESONANT!
        # invariant mass http://cmslxr.fnal.gov/source/DataFormats/EgammaCandidates/interface/GsfElectron.h#0795
        set_trace()
        dieles = [(ele1, ele2) for ele1, ele2 in dieles if (ele1.p4(1) + ele2.p4(1)).mass() < 6.]
        if not len(dieles): return False
        self.counters.counter('BKJPsiEEAnalyzer').inc('diele mass < 6')
        
        selDieles = []
        
        # dieles with a good vertex
        for diele in dieles:
            # clear the vectors
            self.tofit_cc.clear()
            self.tofit_pc.clear()
            # create a RecoChargedCandidate for each reconstructed lepton and flush it into the vector
            for il in diele:
                # if the reco particle is a displaced thing, it does not have the p4() method, so let's build it 
                myp4 = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzE4D<double> >')(il.px(), il.py(), il.pz(), math.sqrt(il.mass()**2 + il.px()**2 + il.py()**2 + il.pz()**2))
                ic = ROOT.reco.RecoChargedCandidate() # instantiate a dummy RecoChargedCandidate
                ic.setCharge(il.charge())             # assign the correct charge
                ic.setP4(myp4)                        # assign the correct p4
                set_trace()
                try:
                    ic.setTrack(il.bestTrackRef())        # set the correct TrackRef
                except:
                    set_trace()
                
                if ic.track().isNonnull():            # check that the track is valid
                    self.tofit_cc.push_back(ic)
            
#             set_trace()
            # further sanity check: two *distinct* tracks
            if self.tofit_cc.size()==2 and self.tofit_cc[0].bestTrackRef() != self.tofit_cc[1].bestTrackRef():
                # fit it!
                svtree = self.vtxfit.Fit(self.tofit_cc, self.tofit_pc) # actual vertex fitting
                # check that the vertex is good
                if not svtree.get().isEmpty() and svtree.get().isValid():
                    svtree.movePointerToTheTop()
                    sv = svtree.currentDecayVertex().get()
                    recoSv = makeRecoVertex(sv, kinVtxTrkSize=2) # need to do some gymastics
                    selDieles.append(diele)
        
        if len(selDieles)==0: return False
        self.counters.counter('BKJPsiEEAnalyzer').inc('good diele vtx')

        cands = []
        
        # add a track
        for diele, tk in product(selDieles, seltracks):
            # remove double countings 
            if deltaR(tk, diele[0])<0.01 or deltaR(tk, diele[1])<0.01:
                continue
            # pt and sanity checks on the track
#             set_trace()
            if not tk.bestTrack():
                continue
            if tk.pt()<0.8:
                continue

            # kaon mass
            tk.setMass(0.493677)
            # clear the vectors
            self.tofit_cc.clear()
            self.tofit_pc.clear()
            # create a RecoChargedCandidate for each reconstructed lepton and flush it into the vector
            for il in diele:
                # if the reco particle is a displaced thing, it does not have the p4() method, so let's build it 
                myp4 = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzE4D<double> >')(il.px(), il.py(), il.pz(), math.sqrt(il.mass()**2 + il.px()**2 + il.py()**2 + il.pz()**2))
                ic = ROOT.reco.RecoChargedCandidate() # instantiate a dummy RecoChargedCandidate
                ic.setCharge(il.charge())             # assign the correct charge
                ic.setP4(myp4)                        # assign the correct p4                
                ic.setTrack(il.bestTrackRef())        # set the correct TrackRef
                if ic.track().isNonnull():            # check that the track is valid
                    self.tofit_cc.push_back(ic)
            
            # push the track into the vector
            self.tofit_pc.push_back(tk.physObj)

            # fit it!
            svtree = self.vtxfit.Fit(self.tofit_cc, self.tofit_pc) # actual vertex fitting
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
    
    
    