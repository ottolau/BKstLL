import ROOT
from itertools import product, combinations
import math

from PhysicsTools.Heppy.analyzers.core.Analyzer      import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle    import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar             import deltaR, deltaR2, bestMatch
from PhysicsTools.Heppy.physicsobjects.Muon          import Muon
from PhysicsTools.Heppy.physicsobjects.Electron      import Electron
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject
from CMGTools.BKstLL.analyzers.utils import displacement2D, displacement3D, makeRecoVertex # utility functions
from CMGTools.BKstLL.physicsobjects.BsJPsiPhi import BsJPsiPhi

from pdb import set_trace

##########################################################################################
# load custom library to ROOT. This contains the kinematic vertex fitter class
ROOT.gSystem.Load('libCMGToolsBKstLL')
from ROOT import KinematicVertexFitter as VertexFitter

class BsJPsiPhiGenAnalyzer_Skim_AOD(Analyzer):
    '''
    '''
    def declareHandles(self):
        super(BsJPsiPhiGenAnalyzer_Skim_AOD, self).declareHandles()

        self.handles['L1muons'] = AutoHandle(
            ('gmtStage2Digis', 'Muon'),
            'BXVector<l1t::Muon>'
        )

        self.handles['electrons'] = AutoHandle(
            #'slimmedElectrons',
            #'selectedPatElectrons',
            #'std::vector<pat::Electron>'
            'gedGsfElectrons',
            'std::vector<reco::GsfElectron>'
        )

        self.handles['muons'] = AutoHandle(
            #'slimmedMuons',
            #'selectedPatMuons',
            #'patMuons',
            #'std::vector<pat::Muon>'
            #'std::vector<reco::Muon>'
            'muons',
            'std::vector<reco::Muon>'
        )

        self.handles['tracks'] = AutoHandle(
            'generalTracks',
            #'packedPFCandidates',
            'std::vector<reco::Track>'
            #'std::vector<reco::TrackCollection>'
            #vector<reco::Track>                   "generalTracks"
        )
 
        self.handles['gsftracks'] = AutoHandle(
            'electronGsfTracks',
            #'packedPFCandidates',
            'std::vector<reco::GsfTrack>'
            #'std::vector<reco::TrackCollection>'
            #vector<reco::Track>                   "generalTracks"
        )
       
        #self.handles['tracks'] = Handle("std::vector")

        self.mchandles['GenParticles'] = AutoHandle(
            'genParticles',
            #('GenParticles', '', 'HLT'),
            'std::vector<reco::GenParticle>'
        )

        self.mchandles['packedGenParticles'] = AutoHandle(
            'packedGenParticles',
            'std::vector<pat::PackedGenParticle>'
        )

        self.handles['beamspot'  ] = AutoHandle(
            ('offlineBeamSpot'              , '', 'RECO'),
            'reco::BeamSpot'
        )


    def beginLoop(self, setup):
        super(BsJPsiPhiGenAnalyzer_Skim_AOD, self).beginLoop(setup)
        self.counters.addCounter('BsJPsiPhiGenAnalyzer_Skim_AOD')
        count = self.counters.counter('BsJPsiPhiGenAnalyzer_Skim_AOD')
        count.register('all events')
        count.register('has a good gen Bs->JPsiPhi')

        # stuff I need to instantiate only once
        self.vtxfit = VertexFitter()
        # create a std::vector<reco::RecoChargedCandidate> to be passed to the fitter 
        self.tofit_cc = ROOT.std.vector('reco::RecoChargedCandidate')()
        # create a std::vector<pat::PackedCandidate> to be passed to the fitter 
        self.tofit_pc = ROOT.std.vector('reco::Track')()

    def process(self, event):
        self.readCollections(event.input)

        self.counters.counter('BsJPsiPhiGenAnalyzer_Skim_AOD').inc('all events')

        # vertex stuff
        event.beamspot    = self.handles['beamspot'].product()

        # get the tracks
        #tracks = map(PhysicsObject, self.handles['tracks'].product())
        tracks = self.handles['tracks'].product()

        # merge the track collections
        event.alltracks = sorted([tt for tt in tracks if tt.charge() != 0], key = lambda x : x.pt(), reverse = True)

        gsftracks = self.handles['gsftracks'].product()
        event.gsftracks = sorted([tt for tt in gsftracks], key = lambda x: x.pt(), reverse = True)


        event.muons     = self.handles['muons'    ].product()
        event.electrons = self.handles['electrons'].product()

        pruned_gen_particles = self.mchandles['GenParticles'].product()
        #packed_gen_particles = self.mchandles['packedGenParticles'].product()
        packed_gen_particles = [ip for ip in pruned_gen_particles if ip.status()==1]

        all_gen_particles = [ip for ip in pruned_gen_particles] #+ [ip for ip in packed_gen_particles]
         
        # HOOK RM
        event.pruned_gen_particles = pruned_gen_particles
        event.packed_gen_particles = packed_gen_particles
        event.all_gen_particles    = all_gen_particles

        # match gen mu to offline mu
        genmus = [ii for ii in all_gen_particles if abs(ii.pdgId())==13 and ii.status()==1]
        for imu in genmus:
            bm, dr = bestMatch(imu, event.muons)
            if dr<0.3:
                imu.reco = bm
        if len(genmus)>0:
            event.thetagmu = sorted(genmus, key = lambda x : x.pt(), reverse = True)[0]

        # match gen ele to offline ele
        geneles = [ii for ii in all_gen_particles if abs(ii.pdgId())==11 and ii.status()==1]
        for iele in geneles:
            bm, dr = bestMatch(iele, event.alltracks)
            if dr<0.3:
                iele.reco = bm

        event.gen_bmesons  = [pp for pp in pruned_gen_particles if abs(pp.pdgId()) > 500 and abs(pp.pdgId()) < 600 and pp.isPromptDecayed()]

        event.gen_bsmesons = [pp for pp in pruned_gen_particles if abs(pp.pdgId())==531]
        
        # walk down the lineage of the B mesons and find the final state muons and charged particles
        for ip in event.gen_bsmesons + event.gen_bmesons:
            if getattr(self.cfg_ana, 'verbose', False):
                print 'PdgId : %s   pt : %s  eta : %s   phi : %s' %(ip.pdgId(), ip.pt(), ip.eta(), ip.phi())    
                print '     daughters'
            finaldaughters = []
            finalcharged   = []
            finalmuons     = []
            for ipp in packed_gen_particles:
                mother = ipp.mother(0)
                if mother and self.isAncestor(ip, mother):
                    if abs(ipp.pdgId())==13:
                        finalmuons.append(ipp)
                    if abs(ipp.charge())==1:
                        finalcharged.append(ipp)
                    if getattr(self.cfg_ana, 'verbose', False):
                        print '     PdgId : %s   pt : %s  eta : %s   phi : %s' %(ipp.pdgId(), ipp.pt(), ipp.eta(), ipp.phi())
                    finaldaughters.append(ipp)
                    if getattr(self.cfg_ana, 'verbose', False):
                        print '     PdgId : %s   pt : %s  eta : %s   phi : %s' %(ipp.pdgId(), ipp.pt(), ipp.eta(), ipp.phi())
            ip.finalmuondaughters    = sorted(finalmuons    , key = lambda x : x.pt(), reverse = True)
            ip.finalchargeddaughters = sorted(finalcharged  , key = lambda x : x.pt(), reverse = True)
            ip.finaldaughters        = sorted(finaldaughters, key = lambda x : x.pt(), reverse = True)

        flav = self.cfg_ana.flavour
        for ibs in event.gen_bsmesons:

            # positive-charged leptons from Jpsi
            lps = [ip for ip in ibs.finaldaughters if ip.pdgId()==-abs(flav)]
            # negative-charged leptons from Jpsi
            lms = [ip for ip in ibs.finaldaughters if ip.pdgId()== abs(flav)]

            for ilep in lps+lms:
                # avoid matching the same particle twice, use the charge to distinguish (charge flip not accounted...)
                tomatch = [jj for jj in event.alltracks if jj.charge()==ilep.charge()]
                bm, dr = bestMatch(ilep, tomatch)
                if dr<0.3:
                    ilep.reco = bm
                tomatch = [jj for jj in event.gsftracks if jj.charge()==ilep.charge()]
                bm, dr = bestMatch(ilep, tomatch)
                if dr<0.3:
                    ilep.gsf = bm
            event.lp = lps[0]
            event.lm = lpm[0]
            break # yeah, only one at a time, mate!
            

        self.counters.counter('BsJPsiPhiGenAnalyzer_Skim_AOD').inc('has a good gen Bs->JPsiPhi')
        
        return True

    @staticmethod
    def isAncestor(a, p):
        if a == p :
            return True
        for i in xrange(0,p.numberOfMothers()):
            if BsJPsiPhiGenAnalyzer_Skim_AOD.isAncestor(a,p.mother(i)):
                return True
        return False

    
    
    
    
    
    
    
