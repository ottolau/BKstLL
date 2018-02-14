import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from CMGTools.BKstLL.analyzers.L1Seeds import single_muon, di_muon, tri_muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from itertools import combinations

class BKstLLGenTreeProducer(L1PurityTreeProducerBase):

    '''
    '''

    def declareVariables(self, setup):
        '''
        '''
        self.bookEvent(self.tree)

        # gen pt hat
        self.var(self.tree, 'qscale')
        
        # gen tagging muon
        self.bookGenParticle(self.tree, 'tag_mu')
        self.bookParticle   (self.tree, 'tag_mu_reco')
                
        # gen level B mesons from the hard interaction, sorted by pt
        self.var(self.tree, 'nbmesons')
        
        self.bookGenParticle(self.tree, 'bd')

        self.bookGenParticle(self.tree, 'b1')
        self.bookGenParticle(self.tree, 'b2')
        self.bookGenParticle(self.tree, 'b3')

        self.bookJet(self.tree, 'b1_jet', fill_extra=True)
        self.bookJet(self.tree, 'b2_jet', fill_extra=True)
        self.bookJet(self.tree, 'b3_jet', fill_extra=True)

        # bbbar topology
        self.var(self.tree, 'dr_btag_bprobe'   )
        self.var(self.tree, 'deta_btag_bprobe' )
        self.var(self.tree, 'dphi_btag_bprobe' )
        self.var(self.tree, 'dr_mutag_bprobe'  )
        self.var(self.tree, 'deta_mutag_bprobe')
        self.var(self.tree, 'dphi_mutag_bprobe')

        # first two gen level muons from B mesons, sorted by pt
        self.bookGenParticle(self.tree, 'bd_lp') ; self.bookParticle(self.tree, 'bd_lp_reco')
        self.bookGenParticle(self.tree, 'bd_lm') ; self.bookParticle(self.tree, 'bd_lm_reco')
        self.bookGenParticle(self.tree, 'bd_pi') ; self.bookParticle(self.tree, 'bd_pi_reco')
        self.bookGenParticle(self.tree, 'bd_k' ) ; self.bookParticle(self.tree, 'bd_k_reco' )
        self.bookJet(self.tree, 'bd_jet', fill_extra=True)
        self.var(self.tree, 'bd_dr' )

        # L1 seeds. Bool, either fired or not
        self.var(self.tree, 'L1_SingleMu_22_eta2p1_Q12'          )
        self.var(self.tree, 'L1_SingleMu_25_Q12'                 )
        self.var(self.tree, 'L1_SingleMu_25_eta1p0_Q12'          )
  
        self.var(self.tree, 'L1_SingleMu_7_eta1p0_Q8'            )
        self.var(self.tree, 'L1_SingleMu_7_eta1p0_Q12'           )
        self.var(self.tree, 'L1_SingleMu_7_eta1p5_Q8'            )
        self.var(self.tree, 'L1_SingleMu_7_eta1p5_Q12'           )
        self.var(self.tree, 'L1_SingleMu_7_eta2p1_Q8'            )
        self.var(self.tree, 'L1_SingleMu_7_eta2p1_Q12'           )
        self.var(self.tree, 'L1_SingleMu_7_eta2p5_Q8'            )
        self.var(self.tree, 'L1_SingleMu_7_eta2p5_Q12'           )
  
        self.var(self.tree, 'L1_SingleMu_10_eta1p0_Q8'           )
        self.var(self.tree, 'L1_SingleMu_10_eta1p0_Q12'          )
        self.var(self.tree, 'L1_SingleMu_10_eta1p5_Q8'           )
        self.var(self.tree, 'L1_SingleMu_10_eta1p5_Q12'          )
        self.var(self.tree, 'L1_SingleMu_10_eta2p1_Q8'           )
        self.var(self.tree, 'L1_SingleMu_10_eta2p1_Q12'          )
        self.var(self.tree, 'L1_SingleMu_10_eta2p5_Q8'           )
        self.var(self.tree, 'L1_SingleMu_10_eta2p5_Q12'          )
  
        self.var(self.tree, 'L1_SingleMu_15_eta1p0_Q8'           )
        self.var(self.tree, 'L1_SingleMu_15_eta1p0_Q12'          )
        self.var(self.tree, 'L1_SingleMu_15_eta1p5_Q8'           )
        self.var(self.tree, 'L1_SingleMu_15_eta1p5_Q12'          )
        self.var(self.tree, 'L1_SingleMu_15_eta2p1_Q8'           )
        self.var(self.tree, 'L1_SingleMu_15_eta2p1_Q12'          )
        self.var(self.tree, 'L1_SingleMu_15_eta2p5_Q8'           )
        self.var(self.tree, 'L1_SingleMu_15_eta2p5_Q12'          )
  
        self.var(self.tree, 'L1_DoubleMu_15_7_Q8'                )
        self.var(self.tree, 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4'  )
        self.var(self.tree, 'L1_DoubleMu4_SQ_OS_dR_Max1p2'       )
        self.var(self.tree, 'L1_DoubleMu4p5_SQ_OS_dR_Max1p2'     )
        self.var(self.tree, 'L1_DoubleMu5_3er2p1_SQ_OS_dR_Max1p4')

        # is the given L1 seed fired by one B meson?
        self.var(self.tree, 'matched_L1_SingleMu_22_eta2p1_Q12'          )
        self.var(self.tree, 'matched_L1_SingleMu_25_Q12'                 )
        self.var(self.tree, 'matched_L1_SingleMu_25_eta1p0_Q12'          )
  
        self.var(self.tree, 'matched_L1_SingleMu_7_eta1p0_Q8'            )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta1p0_Q12'           )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta1p5_Q8'            )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta1p5_Q12'           )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta2p1_Q8'            )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta2p1_Q12'           )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta2p5_Q8'            )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta2p5_Q12'           )
  
        self.var(self.tree, 'matched_L1_SingleMu_10_eta1p0_Q8'           )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta1p0_Q12'          )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta1p5_Q8'           )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta1p5_Q12'          )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta2p1_Q8'           )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta2p1_Q12'          )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta2p5_Q8'           )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta2p5_Q12'          )
  
        self.var(self.tree, 'matched_L1_SingleMu_15_eta1p0_Q8'           )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta1p0_Q12'          )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta1p5_Q8'           )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta1p5_Q12'          )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta2p1_Q8'           )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta2p1_Q12'          )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta2p5_Q8'           )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta2p5_Q12'          )
  
        self.var(self.tree, 'matched_L1_DoubleMu_15_7_Q8'                )
        self.var(self.tree, 'matched_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4'  )
        self.var(self.tree, 'matched_L1_DoubleMu4_SQ_OS_dR_Max1p2'       )
        self.var(self.tree, 'matched_L1_DoubleMu4p5_SQ_OS_dR_Max1p2'     )
        self.var(self.tree, 'matched_L1_DoubleMu5_3er2p1_SQ_OS_dR_Max1p4')

    def process(self, event):
        '''
        '''
        self.readCollections(event.input)
        self.tree.reset()

        if not eval(self.skimFunction):
            return False

        self.fillEvent(self.tree, event)
        self.fill(self.tree, 'qscale', event.qscale)
        
        self.fill(self.tree, 'nbmesons', len(event.gen_bmesons))

        # gen tagging muon
        if hasattr(event, 'thetagmu'):
            self.fillGenParticle(self.tree, 'tag_mu', event.thetagmu)
            if hasattr(event.thetagmu, 'reco'):
                self.fillParticle(self.tree, 'tag_mu_reco', event.thetagmu.reco)

        self.fillGenParticle(self.tree, 'bd', event.kstll)
        self.fill(self.tree, 'bd_dr', event.kstll.dr)

        for jet in event.jets:
            if deltaR(event.kstll, jet) < 0.4:
                self.fillJet(self.tree, 'bd_jet', jet, fill_extra=True)

        self.fillGenParticle(self.tree, 'bd_lp', event.kstll.lp)
                
        if hasattr(event.kstll.lp, 'reco'):
            self.fillParticle(self.tree, 'bd_lp_reco', event.kstll.lp.reco)       
        self.fillGenParticle(self.tree, 'bd_lm', event.kstll.lm)
        if hasattr(event.kstll.lm, 'reco'):
            self.fillParticle(self.tree, 'bd_lm_reco', event.kstll.lm.reco)
        self.fillGenParticle(self.tree, 'bd_pi', event.kstll.pi)
        if hasattr(event.kstll.pi, 'reco'):
            self.fillParticle(self.tree, 'bd_pi_reco', event.kstll.pi.reco)
        self.fillGenParticle(self.tree, 'bd_k' , event.kstll.k )
        if hasattr(event.kstll.k, 'reco'):
            self.fillParticle(self.tree, 'bd_k_reco', event.kstll.k.reco)

        for i, ib in enumerate(event.clean_gen_bmesons[:3]):
            self.fillGenParticle(self.tree, 'b%d' %(i+1), ib)
            for jet in event.jets:
                if deltaR(ib, jet) < 0.4:
                    self.fillJet(self.tree, 'b%d_jet' %(i+1), jet, fill_extra=True)

        if len(event.clean_gen_bmesons)>0:
            self.fill(self.tree, 'dr_btag_bprobe'   , deltaR  (event.clean_gen_bmesons[0], event.kstll)             )
            self.fill(self.tree, 'deta_btag_bprobe' , abs     (event.clean_gen_bmesons[0].eta() - event.kstll.eta()))
            self.fill(self.tree, 'dphi_btag_bprobe' , deltaPhi(event.clean_gen_bmesons[0].phi(), event.kstll.phi()) )
        if hasattr(event, 'thetagmu'):
            self.fill(self.tree, 'dr_mutag_bprobe'  , deltaR  (event.thetagmu, event.kstll)                         )
            self.fill(self.tree, 'deta_mutag_bprobe', abs     (event.thetagmu.eta() - event.kstll.eta())            )
            self.fill(self.tree, 'dphi_mutag_bprobe', deltaPhi(event.thetagmu.phi(), event.kstll.phi())             )
                        
        fired, matched, index = single_muon(event.L1_muons, 22, 2.1, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_22_eta2p1_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_22_eta2p1_Q12', int(matched) * (index+1))
        
        fired, matched, index = single_muon(event.L1_muons, 25, 2.5, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_25_Q12'       , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_25_Q12'       , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 25, 1.0, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_25_eta1p0_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_25_eta1p0_Q12', int(matched) * (index+1))

        fired, matched, index = single_muon(event.L1_muons,  7, 1.0,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta1p0_Q8'  , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta1p0_Q8'  , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons,  7, 1.0, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta1p0_Q12' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta1p0_Q12' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons,  7, 1.5,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta1p5_Q8'  , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta1p5_Q8'  , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons,  7, 1.5, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta1p5_Q12' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta1p5_Q12' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons,  7, 2.1,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta2p1_Q8'  , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta2p1_Q8'  , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons,  7, 2.1, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta2p1_Q12' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta2p1_Q12' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons,  7, 2.5,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta2p5_Q8'  , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta2p5_Q8'  , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons,  7, 2.5, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta2p5_Q12' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta2p5_Q12' , int(matched) * (index+1))

        fired, matched, index = single_muon(event.L1_muons, 10, 1.0,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta1p0_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta1p0_Q8' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 10, 1.0, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta1p0_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta1p0_Q12', int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 10, 1.5,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta1p5_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta1p5_Q8' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 10, 1.5, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta1p5_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta1p5_Q12', int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 10, 2.1,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta2p1_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta2p1_Q8' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 10, 2.1, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta2p1_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta2p1_Q12', int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 10, 2.5,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta2p5_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta2p5_Q8' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 10, 2.5, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta2p5_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta2p5_Q12', int(matched) * (index+1))

        fired, matched, index = single_muon(event.L1_muons, 15, 1.0,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta1p0_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta1p0_Q8' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 15, 1.0, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta1p0_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta1p0_Q12', int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 15, 1.5,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta1p5_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta1p5_Q8' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 15, 1.5, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta1p5_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta1p5_Q12', int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 15, 2.1,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta2p1_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta2p1_Q8' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 15, 2.1, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta2p1_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta2p1_Q12', int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 15, 2.5,  8, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta2p5_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta2p5_Q8' , int(matched) * (index+1))
        fired, matched, index = single_muon(event.L1_muons, 15, 2.5, 12, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta2p5_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta2p5_Q12', int(matched) * (index+1))

        fired, matched, index = di_muon    (event.L1_muons, 15  , 7  , qual1= 8, qual2= 8, matches=event.clean_gen_bmesons);                                        self.fill(self.tree, 'L1_DoubleMu_15_7_Q8'                , fired) ; self.fill(self.tree, 'matched_L1_DoubleMu_15_7_Q8'                , int(matched) * (index+1))
        fired, matched, index = di_muon    (event.L1_muons,  0  , 0  , eta1=1.5, eta2=1.5, qual1=12, qual2=12, maxDr=1.4, sign=0, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4'  , fired) ; self.fill(self.tree, 'matched_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4'  , int(matched) * (index+1))
        fired, matched, index = di_muon    (event.L1_muons,  4  , 4  , qual1=12, qual2=12, maxDr=1.2, sign=0, matches=event.clean_gen_bmesons);                     self.fill(self.tree, 'L1_DoubleMu4_SQ_OS_dR_Max1p2'       , fired) ; self.fill(self.tree, 'matched_L1_DoubleMu4_SQ_OS_dR_Max1p2'       , int(matched) * (index+1))
        fired, matched, index = di_muon    (event.L1_muons,  4.5, 4.5, qual1=12, qual2=12, maxDr=1.2, sign=0, matches=event.clean_gen_bmesons);                     self.fill(self.tree, 'L1_DoubleMu4p5_SQ_OS_dR_Max1p2'     , fired) ; self.fill(self.tree, 'matched_L1_DoubleMu4p5_SQ_OS_dR_Max1p2'     , int(matched) * (index+1))
        fired, matched, index = di_muon    (event.L1_muons,  5  , 3  , eta1=2.1, eta2=2.1, qual1=12, qual2=12, maxDr=1.4, sign=0, matches=event.clean_gen_bmesons); self.fill(self.tree, 'L1_DoubleMu5_3er2p1_SQ_OS_dR_Max1p4', fired) ; self.fill(self.tree, 'matched_L1_DoubleMu5_3er2p1_SQ_OS_dR_Max1p4', int(matched) * (index+1))
        
        self.fillTree(event)

