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
        
        self.bookParticle(self.tree, 'l1')
        self.bookParticle(self.tree, 'l2')
        self.bookParticle(self.tree, 'k')
        self.bookParticle(self.tree, 'b')
        self.bookParticle(self.tree, 'll')

        self.var(self.tree, 'llcone')
        self.var(self.tree, 'bcone')
        
        self.var(self.tree, 'vtx_ls')
        self.var(self.tree, 'vtx_prob')
        self.var(self.tree, 'vtx_cos')

        if getattr(self.cfg_ana, 'addTagMu', False):
            self.bookParticle(self.tree, 'tag_mu')

    def process(self, event):
        '''
        '''
        self.readCollections(event.input)
        self.tree.reset()

        if not eval(self.skimFunction):
            return False

        self.fillEvent(self.tree, event)

        self.fillParticle(self.tree, 'l1', event.myB.l0())
        self.fillParticle(self.tree, 'l2', event.myB.l1())
        self.fillParticle(self.tree, 'k'  , event.myB.k ())
        self.fillParticle(self.tree, 'b'  , event.myB.b() )
        self.fillParticle(self.tree, 'll' , event.myB.ll())

        self.fill(self.tree, 'llcone', event.myB.llcone())
        self.fill(self.tree, 'bcone' , event.myB.bcone())

        self.fill(self.tree, 'vtx_ls'  , event.myB.ls2d()   )
        self.fill(self.tree, 'vtx_prob', event.myB.vtxprob())
        self.fill(self.tree, 'vtx_cos' , event.myB.vtxcos() )

        if getattr(self.cfg_ana, 'addTagMu', False):
            if len(event.muons):
                muons = sorted(event.muons, key = lambda x : x.pt(), reverse=True)
                self.fillParticle(self.tree, 'tag_mu', muons[0])
                    
        self.fillTree(event)


