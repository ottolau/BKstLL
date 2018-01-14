import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from CMGTools.BKstLL.analyzers.L1Seeds import single_muon, di_muon, tri_muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR
from itertools import combinations

class L1PurityTreeProducer(L1PurityTreeProducerBase):

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
        self.bookGenParticle(self.tree, 'b1')
        self.bookGenParticle(self.tree, 'b2')
        self.bookGenParticle(self.tree, 'b3')
        self.bookGenParticle(self.tree, 'b4')

        # first two gen level muons from B mesons, sorted by pt
        self.bookParticle(self.tree, 'b1_m1') ; self.bookParticle(self.tree, 'b1_m2')
        self.bookParticle(self.tree, 'b2_m1') ; self.bookParticle(self.tree, 'b2_m2')
        self.bookParticle(self.tree, 'b3_m1') ; self.bookParticle(self.tree, 'b3_m2')
        self.bookParticle(self.tree, 'b4_m1') ; self.bookParticle(self.tree, 'b4_m2')

        # number of muons / charged particles from each B meson in the final state 
        self.var(self.tree, 'b1_nmuons'  )
        self.var(self.tree, 'b1_ncharged')
        self.var(self.tree, 'b2_nmuons'  )
        self.var(self.tree, 'b2_ncharged')
        self.var(self.tree, 'b3_nmuons'  )
        self.var(self.tree, 'b3_ncharged')
        self.var(self.tree, 'b4_nmuons'  )
        self.var(self.tree, 'b4_ncharged')

        # dR between B meson pairs
        self.var(self.tree, 'dr_12')
        self.var(self.tree, 'dr_13')
        self.var(self.tree, 'dr_14')
        self.var(self.tree, 'dr_23')
        self.var(self.tree, 'dr_24')
        self.var(self.tree, 'dr_34')

        # L1 seeds. Bool, either fired or not
        self.var(self.tree, 'L1_SingleMu_22_eta2p1_Q12'        )
        self.var(self.tree, 'L1_SingleMu_25_Q12'               )
        self.var(self.tree, 'L1_SingleMu_25_eta1p0_Q12'        )

        self.var(self.tree, 'L1_SingleMu_7_eta1p0_Q8'          )
        self.var(self.tree, 'L1_SingleMu_7_eta1p0_Q12'         )
        self.var(self.tree, 'L1_SingleMu_7_eta1p5_Q8'          )
        self.var(self.tree, 'L1_SingleMu_7_eta1p5_Q12'         )
        self.var(self.tree, 'L1_SingleMu_7_eta2p1_Q8'          )
        self.var(self.tree, 'L1_SingleMu_7_eta2p1_Q12'         )
        self.var(self.tree, 'L1_SingleMu_7_eta2p5_Q8'          )
        self.var(self.tree, 'L1_SingleMu_7_eta2p5_Q12'         )

        self.var(self.tree, 'L1_SingleMu_10_eta1p0_Q8'         )
        self.var(self.tree, 'L1_SingleMu_10_eta1p0_Q12'        )
        self.var(self.tree, 'L1_SingleMu_10_eta1p5_Q8'         )
        self.var(self.tree, 'L1_SingleMu_10_eta1p5_Q12'        )
        self.var(self.tree, 'L1_SingleMu_10_eta2p1_Q8'         )
        self.var(self.tree, 'L1_SingleMu_10_eta2p1_Q12'        )
        self.var(self.tree, 'L1_SingleMu_10_eta2p5_Q8'         )
        self.var(self.tree, 'L1_SingleMu_10_eta2p5_Q12'        )

        self.var(self.tree, 'L1_SingleMu_15_eta1p0_Q8'         )
        self.var(self.tree, 'L1_SingleMu_15_eta1p0_Q12'        )
        self.var(self.tree, 'L1_SingleMu_15_eta1p5_Q8'         )
        self.var(self.tree, 'L1_SingleMu_15_eta1p5_Q12'        )
        self.var(self.tree, 'L1_SingleMu_15_eta2p1_Q8'         )
        self.var(self.tree, 'L1_SingleMu_15_eta2p1_Q12'        )
        self.var(self.tree, 'L1_SingleMu_15_eta2p5_Q8'         )
        self.var(self.tree, 'L1_SingleMu_15_eta2p5_Q12'        )

        self.var(self.tree, 'L1_DoubleMu_15_7_Q8'              )
        self.var(self.tree, 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4')
        self.var(self.tree, 'L1_DoubleMu4_SQ_OS_dR_Max1p2'     )
        self.var(self.tree, 'L1_DoubleMu4p5_SQ_OS_dR_Max1p2'   )

        # is the given L1 seed fired by one B meson?
        self.var(self.tree, 'matched_L1_SingleMu_22_eta2p1_Q12'        )
        self.var(self.tree, 'matched_L1_SingleMu_25_Q12'               )
        self.var(self.tree, 'matched_L1_SingleMu_25_eta1p0_Q12'        )

        self.var(self.tree, 'matched_L1_SingleMu_7_eta1p0_Q8'          )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta1p0_Q12'         )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta1p5_Q8'          )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta1p5_Q12'         )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta2p1_Q8'          )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta2p1_Q12'         )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta2p5_Q8'          )
        self.var(self.tree, 'matched_L1_SingleMu_7_eta2p5_Q12'         )

        self.var(self.tree, 'matched_L1_SingleMu_10_eta1p0_Q8'         )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta1p0_Q12'        )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta1p5_Q8'         )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta1p5_Q12'        )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta2p1_Q8'         )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta2p1_Q12'        )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta2p5_Q8'         )
        self.var(self.tree, 'matched_L1_SingleMu_10_eta2p5_Q12'        )

        self.var(self.tree, 'matched_L1_SingleMu_15_eta1p0_Q8'         )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta1p0_Q12'        )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta1p5_Q8'         )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta1p5_Q12'        )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta2p1_Q8'         )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta2p1_Q12'        )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta2p5_Q8'         )
        self.var(self.tree, 'matched_L1_SingleMu_15_eta2p5_Q12'        )

        self.var(self.tree, 'matched_L1_DoubleMu_15_7_Q8'              )
        self.var(self.tree, 'matched_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4')
        self.var(self.tree, 'matched_L1_DoubleMu4_SQ_OS_dR_Max1p2'     )
        self.var(self.tree, 'matched_L1_DoubleMu4p5_SQ_OS_dR_Max1p2'   )

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
        
        for i, ib in enumerate(event.gen_bmesons[:4]):
            self.fillGenParticle(self.tree, 'b%d' %(i+1), ib)
            self.fill(self.tree, 'b%d_nmuons'   %(i+1), len(ib.finalmuondaughters   ))
            self.fill(self.tree, 'b%d_ncharged' %(i+1), len(ib.finalchargeddaughters))
            for j, im in enumerate(ib.finalmuondaughters[:2]):
                self.fillParticle(self.tree, 'b%d_m%d' %(i+1,j+1), im)
                
        if len(event.gen_bmesons)>=2:
            for a,b in combinations(range(min(4, len(event.gen_bmesons))), 2): 
                self.fill(self.tree, 'dr_%d%d' %(a+1, b+1), deltaR(event.gen_bmesons[a], event.gen_bmesons[b]))

        fired, matched = single_muon(event.L1_muons, 22, 2.1, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_22_eta2p1_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_22_eta2p1_Q12', matched)
        fired, matched = single_muon(event.L1_muons, 25, 2.5, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_25_Q12'       , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_25_Q12'       , matched)
        fired, matched = single_muon(event.L1_muons, 25, 1.0, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_25_eta1p0_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_25_eta1p0_Q12', matched)

        fired, matched = single_muon(event.L1_muons,  7, 1.0,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta1p0_Q8'  , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta1p0_Q8'  , matched)
        fired, matched = single_muon(event.L1_muons,  7, 1.0, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta1p0_Q12' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta1p0_Q12' , matched)
        fired, matched = single_muon(event.L1_muons,  7, 1.5,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta1p5_Q8'  , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta1p5_Q8'  , matched)
        fired, matched = single_muon(event.L1_muons,  7, 1.5, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta1p5_Q12' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta1p5_Q12' , matched)
        fired, matched = single_muon(event.L1_muons,  7, 2.1,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta2p1_Q8'  , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta2p1_Q8'  , matched)
        fired, matched = single_muon(event.L1_muons,  7, 2.1, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta2p1_Q12' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta2p1_Q12' , matched)
        fired, matched = single_muon(event.L1_muons,  7, 2.5,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta2p5_Q8'  , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta2p5_Q8'  , matched)
        fired, matched = single_muon(event.L1_muons,  7, 2.5, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_7_eta2p5_Q12' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_7_eta2p5_Q12' , matched)

        fired, matched = single_muon(event.L1_muons, 10, 1.0,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta1p0_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta1p0_Q8' , matched)
        fired, matched = single_muon(event.L1_muons, 10, 1.0, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta1p0_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta1p0_Q12', matched)
        fired, matched = single_muon(event.L1_muons, 10, 1.5,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta1p5_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta1p5_Q8' , matched)
        fired, matched = single_muon(event.L1_muons, 10, 1.5, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta1p5_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta1p5_Q12', matched)
        fired, matched = single_muon(event.L1_muons, 10, 2.1,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta2p1_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta2p1_Q8' , matched)
        fired, matched = single_muon(event.L1_muons, 10, 2.1, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta2p1_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta2p1_Q12', matched)
        fired, matched = single_muon(event.L1_muons, 10, 2.5,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta2p5_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta2p5_Q8' , matched)
        fired, matched = single_muon(event.L1_muons, 10, 2.5, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_10_eta2p5_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_10_eta2p5_Q12', matched)

        fired, matched = single_muon(event.L1_muons, 15, 1.0,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta1p0_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta1p0_Q8' , matched)
        fired, matched = single_muon(event.L1_muons, 15, 1.0, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta1p0_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta1p0_Q12', matched)
        fired, matched = single_muon(event.L1_muons, 15, 1.5,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta1p5_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta1p5_Q8' , matched)
        fired, matched = single_muon(event.L1_muons, 15, 1.5, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta1p5_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta1p5_Q12', matched)
        fired, matched = single_muon(event.L1_muons, 15, 2.1,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta2p1_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta2p1_Q8' , matched)
        fired, matched = single_muon(event.L1_muons, 15, 2.1, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta2p1_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta2p1_Q12', matched)
        fired, matched = single_muon(event.L1_muons, 15, 2.5,  8, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta2p5_Q8' , fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta2p5_Q8' , matched)
        fired, matched = single_muon(event.L1_muons, 15, 2.5, 12, matches=event.gen_bmesons); self.fill(self.tree, 'L1_SingleMu_15_eta2p5_Q12', fired) ; self.fill(self.tree, 'matched_L1_SingleMu_15_eta2p5_Q12', matched)

        fired, matched = di_muon    (event.L1_muons, 15  , 7  , qual1= 8, qual2= 8, matches=event.gen_bmesons);                                        self.fill(self.tree, 'L1_DoubleMu_15_7_Q8'              , fired) ; self.fill(self.tree, 'matched_L1_DoubleMu_15_7_Q8'              , matched)
        fired, matched = di_muon    (event.L1_muons,  0  , 0  , eta1=1.5, eta2=1.5, qual1=12, qual2=12, maxDr=1.4, sign=0, matches=event.gen_bmesons); self.fill(self.tree, 'L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4', fired) ; self.fill(self.tree, 'matched_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4', matched)
        fired, matched = di_muon    (event.L1_muons,  4  , 4  , qual1=12, qual2=12, maxDr=1.2, sign=0, matches=event.gen_bmesons);                     self.fill(self.tree, 'L1_DoubleMu4_SQ_OS_dR_Max1p2'     , fired) ; self.fill(self.tree, 'matched_L1_DoubleMu4_SQ_OS_dR_Max1p2'     , matched)
        fired, matched = di_muon    (event.L1_muons,  4.5, 4.5, qual1=12, qual2=12, maxDr=1.2, sign=0, matches=event.gen_bmesons);                     self.fill(self.tree, 'L1_DoubleMu4p5_SQ_OS_dR_Max1p2'   , fired) ; self.fill(self.tree, 'matched_L1_DoubleMu4p5_SQ_OS_dR_Max1p2'   , matched)

        self.fillTree(event)

