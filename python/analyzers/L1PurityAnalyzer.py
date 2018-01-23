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

        event.bqs = [pp for pp in [jj for jj in packed_gen_particles] + [kk for kk in pruned_gen_particles] if abs(pp.pdgId())==5]

        event.gen_bmesons.sort(key = lambda x : x.pt(), reverse = True)

#         if len(event.gen_bmesons): 
#             for i in event.bqs: 
#                 print 'pdgid %d  \t\tstatus %d\tpt %.2f\teta %.2f\tphi %.2f\t# of moms %d\tmom(0) pdgid %d' %(i.pdgId(), i.status(), i.pt(), i.eta(), i.phi(), i.numberOfMothers(), i.mother(0).pdgId())
            #import pdb ; pdb.set_trace()

        # bb pairs
        event.gluonsplit = False
        event.gluon = None
        for b1, b2 in combinations(event.bqs, 2):
            mothers1 = [b1.mother(ii) for ii in range(b1.numberOfMothers())]
            mothers2 = [b2.mother(ii) for ii in range(b2.numberOfMothers())]
            common = [mm for mm in mothers2 if mm in mothers1]
            if len(common)==1 and common[0].pdgId()==21:
                event.gluonsplit = True
                event.gluon = common[0]
                break
#             if len(common): import pdb ; pdb.set_trace()
        
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
            ip.finalmuondaughters    = sorted(finalmuons, key = lambda x : x.pt(), reverse = True)
            ip.finalchargeddaughters = sorted(finalcharged, key = lambda x : x.pt(), reverse = True)

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
            if L1PurityAnalyzer.isAncestor(a,p.mother(i)):
                return True
        return False

    
    
    
    
    
    
    