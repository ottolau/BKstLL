import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from itertools import combinations

class L1RateTreeProducer(L1PurityTreeProducerBase):

    '''
    '''

    def declareVariables(self, setup):
        '''
        '''
        self.bookEvent(self.tree)

        # level-1 leading muon
        self.bookParticle(self.tree, 'mu_eta_0p5_q8' )
        self.bookParticle(self.tree, 'mu_eta_0p5_q12')
        self.bookParticle(self.tree, 'mu_eta_0p8_q8' )
        self.bookParticle(self.tree, 'mu_eta_0p8_q12')
        self.bookParticle(self.tree, 'mu_eta_1p0_q8' )
        self.bookParticle(self.tree, 'mu_eta_1p0_q12')
        self.bookParticle(self.tree, 'mu_eta_1p5_q8' )
        self.bookParticle(self.tree, 'mu_eta_1p5_q12')
        self.bookParticle(self.tree, 'mu_eta_2p1_q8' )
        self.bookParticle(self.tree, 'mu_eta_2p1_q12')
        self.bookParticle(self.tree, 'mu_eta_2p5_q8' )
        self.bookParticle(self.tree, 'mu_eta_2p5_q12')
        self.bookParticle(self.tree, 'mu_eta_inf_q8' )
        self.bookParticle(self.tree, 'mu_eta_inf_q12')
        
        # book control seed
        self.var(self.tree, 'L1_SingleMu22er2p1')

    def process(self, event):
        '''
        '''
        self.readCollections(event.input)
        self.tree.reset()

        if not eval(self.skimFunction):
            return False

        self.fillEvent(self.tree, event)
        if len(event.L1_muons_eta_0p5_q8 ): self.fillParticle(self.tree, 'mu_eta_0p5_q8' , event.L1_muons_eta_0p5_q8 [0])
        if len(event.L1_muons_eta_0p5_q12): self.fillParticle(self.tree, 'mu_eta_0p5_q12', event.L1_muons_eta_0p5_q12[0])
        if len(event.L1_muons_eta_0p8_q8 ): self.fillParticle(self.tree, 'mu_eta_0p8_q8' , event.L1_muons_eta_0p8_q8 [0])
        if len(event.L1_muons_eta_0p8_q12): self.fillParticle(self.tree, 'mu_eta_0p8_q12', event.L1_muons_eta_0p8_q12[0])
        if len(event.L1_muons_eta_1p0_q8 ): self.fillParticle(self.tree, 'mu_eta_1p0_q8' , event.L1_muons_eta_1p0_q8 [0])
        if len(event.L1_muons_eta_1p0_q12): self.fillParticle(self.tree, 'mu_eta_1p0_q12', event.L1_muons_eta_1p0_q12[0])
        if len(event.L1_muons_eta_1p5_q8 ): self.fillParticle(self.tree, 'mu_eta_1p5_q8' , event.L1_muons_eta_1p5_q8 [0])
        if len(event.L1_muons_eta_1p5_q12): self.fillParticle(self.tree, 'mu_eta_1p5_q12', event.L1_muons_eta_1p5_q12[0])
        if len(event.L1_muons_eta_2p1_q8 ): self.fillParticle(self.tree, 'mu_eta_2p1_q8' , event.L1_muons_eta_2p1_q8 [0])
        if len(event.L1_muons_eta_2p1_q12): self.fillParticle(self.tree, 'mu_eta_2p1_q12', event.L1_muons_eta_2p1_q12[0])
        if len(event.L1_muons_eta_2p5_q8 ): self.fillParticle(self.tree, 'mu_eta_2p5_q8' , event.L1_muons_eta_2p5_q8 [0])
        if len(event.L1_muons_eta_2p5_q12): self.fillParticle(self.tree, 'mu_eta_2p5_q12', event.L1_muons_eta_2p5_q12[0])
        if len(event.L1_muons_eta_inf_q8 ): self.fillParticle(self.tree, 'mu_eta_inf_q8' , event.L1_muons_eta_inf_q8 [0])
        if len(event.L1_muons_eta_inf_q12): self.fillParticle(self.tree, 'mu_eta_inf_q12', event.L1_muons_eta_inf_q12[0])
        
        self.fill(self.tree, 'L1_SingleMu22er2p1', event.ugt.getAlgoDecisionFinal(29))
        
        self.fillTree(event)

