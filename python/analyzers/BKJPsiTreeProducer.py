import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from CMGTools.BKstLL.analyzers.L1Seeds import single_muon, di_muon, tri_muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from itertools import combinations

class BKJPsiTreeProducer(L1PurityTreeProducerBase):

    '''
    '''

    def declareVariables(self, setup):
        '''
        '''
        self.bookEvent(self.tree)
        
        self.bookParticle(self.tree, 'l_1')
        self.bookParticle(self.tree, 'l_2')
        self.bookParticle(self.tree, 'k')
        self.bookParticle(self.tree, 'b')
        self.bookParticle(self.tree, 'll')
        
        self.var(self.tree, 'vtx_ls')
        self.var(self.tree, 'vtx_prob')

    def process(self, event):
        '''
        '''
        self.readCollections(event.input)
        self.tree.reset()

        if not eval(self.skimFunction):
            return False

        self.fillEvent(self.tree, event)

        self.fillParticle(self.tree, 'l_1', event.myB.l0())
        self.fillParticle(self.tree, 'l_2', event.myB.l1())
        self.fillParticle(self.tree, 'k'  , event.myB.k ())
        self.fillParticle(self.tree, 'b'  , event.myB.b() )
        self.fillParticle(self.tree, 'll' , event.myB.ll())

        self.fill(self.tree, 'vtx_ls'  , event.myB.ls2d()   )
        self.fill(self.tree, 'vtx_prob', event.myB.vtxprob())
        
        self.fillTree(event)


