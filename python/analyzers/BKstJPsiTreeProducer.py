import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from CMGTools.BKstLL.analyzers.L1Seeds import single_muon, di_muon, tri_muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from itertools import combinations

class BKstJPsiTreeProducer(L1PurityTreeProducerBase):

    '''
    '''

    def declareVariables(self, setup):
        '''
        '''
        self.bookEvent(self.tree)
        
        self.bookParticle(self.tree, 'l1')
        self.bookParticle(self.tree, 'l2')
        self.bookParticle(self.tree, 'pi')
        self.bookParticle(self.tree, 'kst')
        self.bookParticle(self.tree, 'k')
        self.bookParticle(self.tree, 'b')
        self.bookParticle(self.tree, 'll')

        self.var(self.tree, 'llcone')
        self.var(self.tree, 'bcone')
        self.var(self.tree, 'kstcone')
        
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

#         if abs(event.myB.l0().mass() - 0.000511)> 0.0001:
#             print 'problem with electron mass'
#             import pdb ; pdb.set_trace()
# 
#         if abs(event.myB.pi().mass() - 0.13957061) > 0.1:
#             print 'problem with pi mass'
#             import pdb ; pdb.set_trace()
#         
#         if abs(event.myB.k().mass() - 0.493677) > 0.1:
#             print 'problem with k mass'
#             import pdb ; pdb.set_trace()

        self.fillParticle(self.tree, 'l1' , event.myB.l0 ())
        self.fillParticle(self.tree, 'l2' , event.myB.l1 ())
        self.fillParticle(self.tree, 'pi' , event.myB.pi ())
        self.fillParticle(self.tree, 'kst', event.myB.kst())
        self.fillParticle(self.tree, 'k'  , event.myB.k  ())
        self.fillParticle(self.tree, 'b'  , event.myB.b  ())
        self.fillParticle(self.tree, 'll' , event.myB.ll ())

        self.fill(self.tree, 'llcone' , event.myB.llcone ())
        self.fill(self.tree, 'bcone'  , event.myB.bcone  ())
        self.fill(self.tree, 'kstcone', event.myB.kstcone())

        self.fill(self.tree, 'vtx_ls'  , event.myB.ls2d()   )
        self.fill(self.tree, 'vtx_prob', event.myB.vtxprob())
        self.fill(self.tree, 'vtx_cos' , event.myB.vtxcos() )

        if getattr(self.cfg_ana, 'addTagMu', False):
            if len(event.muons):
                muons = sorted(event.muons, key = lambda x : x.pt(), reverse=True)
                self.fillParticle(self.tree, 'tag_mu', muons[0])
                    
        self.fillTree(event)


