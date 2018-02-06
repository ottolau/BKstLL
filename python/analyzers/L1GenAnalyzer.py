import ROOT
from itertools import product, combinations
import math

from PhysicsTools.Heppy.analyzers.core.Analyzer       import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle     import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar              import deltaR, deltaR2
from PhysicsTools.Heppy.physicsobjects.PhysicsObjects import Jet

from pdb import set_trace

class L1GenAnalyzer(Analyzer):
    '''
    '''
    def declareHandles(self):
        super(L1GenAnalyzer, self).declareHandles()

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

        self.handles['jets'] = AutoHandle(
            'slimmedJets',
            'std::vector<pat::Jet>'
        )

    def beginLoop(self, setup):
        super(L1GenAnalyzer, self).beginLoop(setup)
        self.counters.addCounter('L1Gen')
        count = self.counters.counter('L1Gen')
        count.register('all events')
        count.register('no PU vertex has pt hat higher than qsquare')

    def process(self, event):
        self.readCollections(event.input)

        self.counters.counter('L1Gen').inc('all events')

        # all jets
        miniaodjets = [Jet(jet) for jet in self.handles['jets'].product()]

        # attach qscale (pt hat) to the event
        event.qscale = self.mchandles['genInfo'].product().qScale()
        
        # filter events where the pt hat of a PU vertex is higher than that of the hard scattering
        if getattr(self.cfg_ana, 'filterPU', True):
            if len(event.pileUpVertex_ptHat)>0:
                if max(event.pileUpVertex_ptHat) > event.qscale and 'QCD' in self.cfg_comp.name:
                    return False
            self.counters.counter('L1Gen').inc('no PU vertex has pt hat higher than qsquare')

        # find the gen B mesons from the hard scattering
        pruned_gen_particles = self.mchandles['prunedGenParticles'].product()
        packed_gen_particles = self.mchandles['packedGenParticles'].product()
 
        all_gen_particles = [ip for ip in pruned_gen_particles] + [ip for ip in packed_gen_particles]
         
        event.gen_bmesons         = [pp for pp in pruned_gen_particles if abs(pp.pdgId()) > 500 and abs(pp.pdgId()) < 600 and pp.isPromptDecayed()]
        event.gen_dmesons         = [pp for pp in pruned_gen_particles if abs(pp.pdgId()) > 400 and abs(pp.pdgId()) < 500 and pp.isPromptDecayed()]
        event.gen_prompt_jpsis    = [pp for pp in pruned_gen_particles if abs(pp.pdgId()) in [443, 100443, 30443, 9000443, 9010443, 9020443] and pp.isPromptDecayed()]
        event.gen_prompt_upsilons = [pp for pp in pruned_gen_particles if abs(pp.pdgId()) in [553, 100553, 200553, 300553] and pp.isPromptDecayed()]
        event.gen_vbosons         = event.genVBosons
        event.gen_topquarks       = event.gentopquarks

#         if len(event.gen_vbosons) and len(event.gen_dmesons):
#             import pdb ; pdb.set_trace()

        # is this double parton scattering, WOW!
        # https://indico.cern.ch/event/629232/contributions/2807337/attachments/1579753/2495883/Epiphany2018-Maciula.pdf
        # https://indico.cern.ch/event/75247/contributions/1246824/attachments/1048811/1495114/jackson.pdf
#         if len(event.gen_bmesons) and len(event.gen_dmesons):
#             import pdb ; pdb.set_trace()
            
        # walk down the lineage of the B mesons and find the final state muons and charged particles
        for ip in event.gen_bmesons + event.gen_dmesons + event.gen_prompt_jpsis + event.gen_prompt_upsilons + event.gen_vbosons + event.gen_topquarks:
            if getattr(self.cfg_ana, 'verbose', False):
                print 'PdgId : %s   pt : %s  eta : %s   phi : %s' %(ip.pdgId(), ip.pt(), ip.eta(), ip.phi())    
                print '     daughters'
            finaldaus    = []
            finalmuons   = []
            finalcharged = []
            for ipp in packed_gen_particles:
                mother = ipp.mother(0)
                if mother and self.isAncestor(ip, mother):
                    finaldaus.append(ipp)
                    if abs(ipp.pdgId())==13:
                        finalmuons.append(ipp)
                    if abs(ipp.charge())==1:
                        finalcharged.append(ipp)
                    if getattr(self.cfg_ana, 'verbose', False):
                        print '     PdgId : %s   pt : %s  eta : %s   phi : %s' %(ipp.pdgId(), ipp.pt(), ipp.eta(), ipp.phi())
            ip.finaldaughters        = sorted(finaldaus   , key = lambda x : x.pt(), reverse = True)
            ip.finalmuondaughters    = sorted(finalmuons  , key = lambda x : x.pt(), reverse = True)
            ip.finalchargeddaughters = sorted(finalcharged, key = lambda x : x.pt(), reverse = True)

        # now find the L1 muons from BX = 0
        L1muons_allbx = self.handles['L1muons'].product()

        L1_muons = []
    
        for jj in range(L1muons_allbx.size(0)):
            L1_muons.append(L1muons_allbx.at(0,jj))
        
        event.L1_muons = L1_muons

        # append a reco jet (if it exists) to each L1 muon
        for mu in event.L1_muons:
            jets = sorted([jj for jj in event.jets if deltaR2(jj.eta(), jj.phi(), mu.etaAtVtx(), mu.phiAtVtx()) < 0.16], key = lambda x : deltaR2(x, mu))
            if len(jets):
                mu.jet = jets[0]

        return True

    @staticmethod
    def isAncestor(a, p):
        if a == p :
            return True
        for i in xrange(0,p.numberOfMothers()):
            if L1GenAnalyzer.isAncestor(a,p.mother(i)):
                return True
        return False

    
    
    
    
    
    
    