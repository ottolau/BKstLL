import ROOT
from itertools import product, combinations
import math

from PhysicsTools.Heppy.analyzers.core.Analyzer       import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle     import AutoHandle
from PhysicsTools.Heppy.physicsobjects.PhysicsObjects import Jet, Muon
from PhysicsTools.HeppyCore.utils.deltar              import deltaR, deltaR2

from pdb import set_trace

class L1RateAnalyzer(Analyzer):
    '''
    '''
    def declareHandles(self):
        super(L1RateAnalyzer, self).declareHandles()

        self.handles['L1muons'] = AutoHandle(
            ('gmtStage2Digis', 'Muon'),
            'BXVector<l1t::Muon>'
        )

        self.handles['ugt'] = AutoHandle(
            ('gtStage2Digis'),
            'BXVector<GlobalAlgBlk>'
        )

        self.handles['jets'] = AutoHandle(
            'slimmedJets',
            'std::vector<pat::Jet>'
        )

        self.handles['muons'] = AutoHandle(
            'slimmedMuons',
            'std::vector<pat::Muon>'
        )

    def beginLoop(self, setup):
        super(L1RateAnalyzer, self).beginLoop(setup)
        self.counters.addCounter('L1Rate')
        count = self.counters.counter('L1Rate')
        count.register('all events')

    def process(self, event):
        self.readCollections(event.input)

        # all jets
        miniaodjets = [Jet(jet) for jet in self.handles['jets'].product()]

        # ugt decision
        event.ugt = self.handles['ugt'].product().at(0,0)
        
        # now find the L1 muons from BX = 0
        L1muons_allbx = self.handles['L1muons'].product()

        if getattr(self.cfg_ana, 'onlyBX0', True):
            L1_muons = []
            for jj in range(L1muons_allbx.size(0)):
                L1_muons.append(L1muons_allbx.at(0,jj))
        
            event.L1_muons = L1_muons

        else:
            event.L1_muons = [mu for mu in L1muons_allbx]
            if len(event.L1_muons): import pdb ; pdb.set_trace()
        
        # append a reco jet (if it exists) to each L1 muon
        for mu in event.L1_muons:
            jets = sorted([jj for jj in event.jets if deltaR2(jj.eta(), jj.phi(), mu.etaAtVtx(), mu.phiAtVtx()) < 0.16], key = lambda x : deltaR2(x, mu))
            if len(jets):
                mu.jet = jets[0]
        
        event.L1_muons_eta_0p5_q8  = sorted([mu for mu in event.L1_muons if abs(mu.eta())<0.5       and mu.hwQual()>=8 ], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_0p5_q12 = sorted([mu for mu in event.L1_muons if abs(mu.eta())<0.5       and mu.hwQual()>=12], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_0p8_q8  = sorted([mu for mu in event.L1_muons if abs(mu.eta())<0.8       and mu.hwQual()>=8 ], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_0p8_q12 = sorted([mu for mu in event.L1_muons if abs(mu.eta())<0.8       and mu.hwQual()>=12], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_1p0_q8  = sorted([mu for mu in event.L1_muons if abs(mu.eta())<1.0       and mu.hwQual()>=8 ], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_1p0_q12 = sorted([mu for mu in event.L1_muons if abs(mu.eta())<1.0       and mu.hwQual()>=12], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_1p5_q8  = sorted([mu for mu in event.L1_muons if abs(mu.eta())<1.5       and mu.hwQual()>=8 ], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_1p5_q12 = sorted([mu for mu in event.L1_muons if abs(mu.eta())<1.5       and mu.hwQual()>=12], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_2p1_q8  = sorted([mu for mu in event.L1_muons if abs(mu.eta())<2.1043125 and mu.hwQual()>=8 ], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_2p1_q12 = sorted([mu for mu in event.L1_muons if abs(mu.eta())<2.1043125 and mu.hwQual()>=12], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_2p5_q8  = sorted([mu for mu in event.L1_muons if abs(mu.eta())<2.5       and mu.hwQual()>=8 ], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_2p5_q12 = sorted([mu for mu in event.L1_muons if abs(mu.eta())<2.5       and mu.hwQual()>=12], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_inf_q8  = sorted([mu for mu in event.L1_muons if                             mu.hwQual()>=8 ], key = lambda mu : mu.pt(), reverse = True)
        event.L1_muons_eta_inf_q12 = sorted([mu for mu in event.L1_muons if                             mu.hwQual()>=12], key = lambda mu : mu.pt(), reverse = True)

#         # 29  L1_SingleMu22er2p1
#         if ugt.getAlgoDecisionFinal(29) and len([mu for mu in event.L1_muons_eta_2p1_q12 if mu.pt()>=22.])==0:
#             import pdb ; pdb.set_trace()

#         if any([mm.eta() == 2.1 for mm in event.L1_muons]):
#             print 'cazz...'
#             import pdb ; pdb.set_trace()
        
#         if len(event.L1_muons_eta_2p5_q12)>1:
#             for mm in event.L1_muons_eta_2p5_q12: print mm.pt(), mm.eta(), mm.phi(), mm.hwQual()
#             import pdb ; pdb.set_trace()

#         if event.run==305081: 

        event.fakeL1mu = True
        event.deltapt  = -9999.
        event.isiso    = -99.
        event.lowptjet = False

        if len(event.L1_muons_eta_2p1_q12)>0 and event.L1_muons_eta_2p1_q12[0].pt()>=22.:
        
            print 'Jets, unfiltered'
            for jet in miniaodjets: 
                print jet, '\tpfJetId', jet.jetID("POG_PFID_Loose") , '\tpuJetId', jet.puJetId(), '\tbtag', jet.btag('pfCombinedInclusiveSecondaryVertexV2BJetTags')
                if deltaR2(jet.eta(), jet.phi(), event.L1_muons_eta_2p1_q12[0].etaAtVtx(), event.L1_muons_eta_2p1_q12[0].phiAtVtx())<0.25:
                    if jet.pt()<20.:
                        event.lowptjet = True

            print 'L1 muons'
            for mu in event.L1_muons_eta_2p1_q12:
                print 'pt %.2f\teta %.2f\tphi %.2f\tetaAtVtx %.2f\tphiAtVtx %.2f' %(mu.pt(), mu.eta(), mu.phi(), mu.etaAtVtx(), mu.phiAtVtx()) 
                        
            print 'offline muons'
            for mm in [Muon(ii) for ii in self.handles['muons'].product()]:
                print mm
                
                if deltaR2(mm.eta(), mm.phi(), event.L1_muons_eta_2p1_q12[0].etaAtVtx(), event.L1_muons_eta_2p1_q12[0].phiAtVtx())<0.25:
                    event.fakeL1mu = False
                    event.isiso    = int(mm.relIso()<1.)
                    event.deltapt  = event.L1_muons_eta_2p1_q12[0].pt() - mm.pt()
                     
            # jj = miniaodjets[1]
            # deltaR(jj.eta(), jj.phi(), mu.etaAtVtx(), mu.phiAtVtx())
            
#             import pdb ; pdb.set_trace()

        return True

