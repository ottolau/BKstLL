import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from CMGTools.BKstLL.analyzers.L1Seeds import single_muon, di_muon, tri_muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR
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
                
        # gen level B mesons from the hard interaction, sorted by pt
        self.var(self.tree, 'nbmesons')
        
        self.bookGenParticle(self.tree, 'bd')

        self.bookGenParticle(self.tree, 'b1')
        self.bookGenParticle(self.tree, 'b2')
        self.bookGenParticle(self.tree, 'b3')

        # first two gen level muons from B mesons, sorted by pt
        self.bookGenParticle(self.tree, 'bd_lp')
        self.bookGenParticle(self.tree, 'bd_lm')
        self.bookGenParticle(self.tree, 'bd_pi')
        self.bookGenParticle(self.tree, 'bd_k' )


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
        
        self.fillGenParticle(self.tree, 'bd', event.kstll)

        self.fillGenParticle(self.tree, 'bd_lp', event.kstll.lp)
        self.fillGenParticle(self.tree, 'bd_lm', event.kstll.lm)
        self.fillGenParticle(self.tree, 'bd_pi', event.kstll.pi)
        self.fillGenParticle(self.tree, 'bd_k' , event.kstll.k )

        for i, ib in enumerate(event.clean_gen_bmesons[:3]):
            self.fillGenParticle(self.tree, 'b%d' %(i+1), ib)

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

