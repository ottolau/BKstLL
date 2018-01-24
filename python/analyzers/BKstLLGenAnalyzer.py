import ROOT
from itertools import product, combinations
import math

from PhysicsTools.Heppy.analyzers.core.Analyzer   import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar          import deltaR, deltaR2

from pdb import set_trace

class BKstLLGenAnalyzer(Analyzer):
    '''
    '''
    def declareHandles(self):
        super(BKstLLGenAnalyzer, self).declareHandles()

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
        super(BKstLLGenAnalyzer, self).beginLoop(setup)
        self.counters.addCounter('BKstLLGenAnalyzer')
        count = self.counters.counter('BKstLLGenAnalyzer')
        count.register('all events')
        count.register('has a good gen B0->K*LL')

    def process(self, event):
        self.readCollections(event.input)

        self.counters.counter('BKstLLGenAnalyzer').inc('all events')

        # attach qscale (pt hat) to the event
        event.qscale = self.mchandles['genInfo'].product().qScale()

        # find the gen B mesons from the hard scattering
        pruned_gen_particles = self.mchandles['prunedGenParticles'].product()
        packed_gen_particles = self.mchandles['packedGenParticles'].product()
        
        event.gen_bmesons  = [pp for pp in pruned_gen_particles if abs(pp.pdgId()) > 500 and abs(pp.pdgId()) < 600 and pp.isPromptDecayed()]

        event.gen_b0mesons = [pp for pp in pruned_gen_particles if abs(pp.pdgId())==511]
        
        # walk down the lineage of the B mesons and find the final state muons and charged particles
        for ip in event.gen_b0mesons + event.gen_bmesons:
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

        for ib0 in event.gen_b0mesons:
            isit, lp, lm, pi, k = self.isKstLL(ib0, abs(self.cfg_ana.flavour)) 
            if isit:
                event.kstll    = ib0
                event.kstll.lp = lp
                event.kstll.lm = lm
                event.kstll.pi = pi
                event.kstll.k  = k 
                break # yeah, only one at a time, mate!
        
        if not hasattr(event, 'kstll'):
            return False
        self.counters.counter('BKstLLGenAnalyzer').inc('has a good gen B0->K*LL')
        
        toclean = None
        
        if event.kstll.isPromptDecayed():
            toclean = event.kstll
        
        else:
            for ip in event.gen_bmesons:
                if self.isAncestor(ip, event.kstll):
                    toclean = ip
                    break
        
#         if toclean is None:
#             print 'nothing to clean lumi %d, ev %d' %(event.lumi, event.eventId)
#             return False
#             import pdb ; pdb.set_trace()
# 
#         self.counters.counter('BKstLLGenAnalyzer').inc('no fuck ups')
        
        
#         elif abs(event.kstll.mother(0).pdgId())>500 and abs(event.kstll.mother(0).pdgId())<600:
#             if toclean = event.kstll.mother(0)
#         elif abs(event.kstll.mother(0).mother(0).pdgId())>500 and abs(event.kstll.mother(0).mother(0).pdgId())<600:
#         else:
#             print 'nothing to clean lumi %d, ev %d' %(event.lumi, event.eventId)
#             import pdb ; pdb.set_trace()
            
        if toclean is not None:
            event.clean_gen_bmesons = [ib for ib in event.gen_bmesons if ib!=toclean]
        else:
            event.clean_gen_bmesons = event.gen_bmesons 
            
#         import pdb ; pdb.set_trace()
                
        # now find the L1 muons from BX = 0
        L1muons_allbx = self.handles['L1muons'].product()

        L1_muons = []
    
        for jj in range(L1muons_allbx.size(0)):
            L1_muons.append(L1muons_allbx.at(0,jj))
        
        event.L1_muons = L1_muons

        return True


    @staticmethod
    def isKstLL(b0meson, flav=13):
#         isB0 = (b0meson.pdgId()==511)
#         isAntiB0 = (not isB0)
        # positive-charged leptons
        lps = [ip for ip in b0meson.finaldaughters if ip.pdgId()==-abs(flav)]
        # negative-charged leptons
        lms = [ip for ip in b0meson.finaldaughters if ip.pdgId()== abs(flav)]
        # pions from K*
        pis = [ip for ip in b0meson.finaldaughters if abs(ip.pdgId())== 211]
        # kaons from K*
        ks  = [ip for ip in b0meson.finaldaughters if abs(ip.pdgId())== 321]
        
        if len(lps) == len(lms) == len(pis) == len(ks) == 1:
            return True, lps[0], lms[0], pis[0], ks[0]
        else:
            return False, None, None, None, None        
        
    @staticmethod
    def isAncestor(a, p):
        if a == p :
            return True
        for i in xrange(0,p.numberOfMothers()):
            if BKstLLGenAnalyzer.isAncestor(a,p.mother(i)):
                return True
        return False

    
    
    
    
    
    
    