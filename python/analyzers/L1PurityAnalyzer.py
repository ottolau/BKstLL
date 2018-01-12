import ROOT
from itertools import product, combinations
import math

from PhysicsTools.Heppy.analyzers.core.Analyzer   import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar          import deltaR, deltaR2

from pdb import set_trace

class L1PurityAnalyzer(Analyzer):
    '''
    '''
    def declareHandles(self):
        super(L1PurityAnalyzer, self).declareHandles()

        self.handles['taus'] = AutoHandle(
            'slimmedTaus',
            'std::vector<pat::Tau>'
        )

        self.handles['electrons'] = AutoHandle(
            'slimmedElectrons',
            'std::vector<pat::Electron>'
        )

        self.handles['muons'] = AutoHandle(
            'slimmedMuons',
            'std::vector<pat::Muon>'
        )

        self.handles['jets'] = AutoHandle(
            'slimmedJets',
            'std::vector<pat::Jet>'
        )

        self.handles['L1muons'] = AutoHandle(
            ('gmtStage2Digis', 'Muon'),
            'BXVector<l1t::Muon>'
        )

        self.mchandles['prunedGenParticles'] = AutoHandle(
            'prunedGenParticles',
            'std::vector<reco::GenParticle>'
        )

        self.mchandles['packedGenParticles'] = AutoHandle(
            'packedGenParticles',
            'std::vector<pat::PackedGenParticle>'
        )

        self.mchandles['genInfo'] = AutoHandle(
            'generator',
            'GenEventInfoProduct'
        )

    def beginLoop(self, setup):
        super(L1PurityAnalyzer, self).beginLoop(setup)
        self.counters.addCounter('L1Purity')
        count = self.counters.counter('L1Purity')
        count.register('all events')

    def process(self, event):
        self.readCollections(event.input)

        # attach qscale (pt hat) to the event
        event.qscale = self.mchandles['genInfo'].product().qScale()

        # find the gen B mesons from the hard scattering
        pruned_gen_particles = self.mchandles['prunedGenParticles'].product()
        packed_gen_particles = self.mchandles['packedGenParticles'].product()
        
        event.gen_bmesons = [pp for pp in pruned_gen_particles if abs(pp.pdgId()) > 500 and abs(pp.pdgId()) < 600 and pp.isPromptDecayed()]

        # walk down the lineage of the B mesons and find the final state muons and charged particles
        for ip in event.gen_bmesons :
            if getattr(self.cfg_ana, 'verbose', False):
                print 'PdgId : %s   pt : %s  eta : %s   phi : %s' %(ip.pdgId(), ip.pt(), ip.eta(), ip.phi())    
                print '     daughters'
            finalmuons   = []
            finalcharged = []
            for ipp in packed_gen_particles:
                mother = ipp.mother(0)
                if mother and self.isAncestor(ip, mother):
                    if abs(ipp.pdgId())==13:
                        finalmuons.append(ipp)
                    if abs(ipp.charge())==1:
                        finalcharged.append(ipp)
                    if getattr(self.cfg_ana, 'verbose', False):
                        print '     PdgId : %s   pt : %s  eta : %s   phi : %s' %(ipp.pdgId(), ipp.pt(), ipp.eta(), ipp.phi())
            ip.finalmuondaughters    = finalmuons
            ip.finalchargeddaughters = finalcharged

        # now find the L1 muons from BX = 0
        L1muons_allbx = self.handles['L1muons'].product()

        L1_muons = []
    
        for jj in range(L1muons_allbx.size(0)):
            L1_muons.append(L1muons_allbx.at(0,jj))
        
        event.L1_muons = L1_muons

        return True

    @staticmethod
    def isAncestor(a, p):
        if a == p :
            return True
        for i in xrange(0,p.numberOfMothers()):
            if isAncestor(a,p.mother(i)):
                return True
        return False

    
    
    
    
    
    
    