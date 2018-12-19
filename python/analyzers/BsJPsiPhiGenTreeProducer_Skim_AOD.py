import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from CMGTools.BKstLL.analyzers.L1Seeds import single_muon, di_muon, tri_muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from itertools import combinations
from ROOT import KinematicVertexFitter as VertexFitter

class BsJPsiPhiGenTreeProducer_Skim_AOD(L1PurityTreeProducerBase):

    '''
    '''

    def declareVariables(self, setup):
        '''
        '''

        
        # first two gen level muons from B mesons, sorted by pt
        self.bookGenParticle(self.tree, 'bs_lp') ; self.bookParticle(self.tree, 'bs_lp_reco')
        self.bookGenParticle(self.tree, 'bs_lm') ; self.bookParticle(self.tree, 'bs_lm_reco')
        self.bookParticle(self.tree, 'bs_lp_gsf')
        self.bookParticle(self.tree, 'bs_lm_gsf')



    def process(self, event):
        '''
        '''
        self.readCollections(event.input)
        self.tree.reset()

        if not eval(self.skimFunction):
            return False

        #self.fillEvent(self.tree, event)

        #self.fill(self.tree, 'qscale', event.qscale)
        self.fillGenParticle(self.tree, 'bs_lp', event.lp)
        if hasattr(event.lp, 'reco'):
            self.fillParticle(self.tree, 'bs_lp_reco', event.lp.reco)       
        self.fillGenParticle(self.tree, 'bs_lm', event.lm)
        if hasattr(event.lm, 'reco'):
            self.fillParticle(self.tree, 'bs_lm_reco', event.lm.reco)
        if hasattr(event.lp, 'gsf'):
            self.fillParticle(self.tree, 'bs_lp_gsf', event.lp.gsf)       
        if hasattr(event.lm, 'gsf'):
            self.fillParticle(self.tree, 'bs_lm_gsf', event.lm.gsf)       
 

        self.fillTree(event)


