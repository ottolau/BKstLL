import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from CMGTools.BKstLL.analyzers.L1Seeds import single_muon, di_muon, tri_muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from itertools import combinations

class BKstJPsiEEGenTreeProducer_AOD(L1PurityTreeProducerBase):

    '''
    '''

    def declareVariables(self, setup):
        '''
        '''
        #self.bookEvent(self.tree)

        # gen pt hat
        #self.var(self.tree, 'qscale')
        
        # gen tagging muon
        self.bookGenParticle(self.tree, 'tag_mu')
        self.bookParticle   (self.tree, 'tag_mu_reco')
                
        # gen level B mesons from the hard interaction, sorted by pt
        self.var(self.tree, 'nbmesons')
        
        self.bookGenParticle(self.tree, 'b0')

        self.bookGenParticle(self.tree, 'b1')
        self.bookGenParticle(self.tree, 'b2')
        self.bookGenParticle(self.tree, 'b3')

        # bbbar topology
        self.var(self.tree, 'dr_btag_bprobe'   )
        self.var(self.tree, 'deta_btag_bprobe' )
        self.var(self.tree, 'dphi_btag_bprobe' )
        self.var(self.tree, 'dr_mutag_bprobe'  )
        self.var(self.tree, 'deta_mutag_bprobe')
        self.var(self.tree, 'dphi_mutag_bprobe')

        # first two gen level muons from B mesons, sorted by pt
        self.bookGenParticle(self.tree, 'b0_lp') ; self.bookParticle(self.tree, 'b0_lp_reco')
        self.bookGenParticle(self.tree, 'b0_lm') ; self.bookParticle(self.tree, 'b0_lm_reco')
        self.bookGenParticle(self.tree, 'b0_pi') ; self.bookParticle(self.tree, 'b0_pi_reco')
        self.bookGenParticle(self.tree, 'b0_k') ; self.bookParticle(self.tree, 'b0_k_reco' )
        self.var(self.tree, 'b0_dr' )

        self.var(self.tree, 'tag_mu_reco_dxy')
        self.var(self.tree, 'tag_mu_reco_dz')
        self.var(self.tree, 'tag_mu_reco_dxyError')
        self.var(self.tree, 'tag_mu_reco_dzError')
        self.var(self.tree, 'tag_mu_reco_vz')

        self.var(self.tree, 'foundB0')
        self.var(self.tree, 'b0_ll_reco_dR')
        self.var(self.tree, 'b0_lp_reco_dxy')
        self.var(self.tree, 'b0_lp_reco_dz')
        self.var(self.tree, 'b0_lp_reco_dxyError')
        self.var(self.tree, 'b0_lp_reco_dzError')
        self.var(self.tree, 'b0_lp_reco_chi2')
        self.var(self.tree, 'b0_lp_reco_ndof')
        self.var(self.tree, 'b0_lp_reco_vz')

        self.var(self.tree, 'b0_lm_reco_dxy')
        self.var(self.tree, 'b0_lm_reco_dz')
        self.var(self.tree, 'b0_lm_reco_dxyError')
        self.var(self.tree, 'b0_lm_reco_dzError')
        self.var(self.tree, 'b0_lm_reco_chi2')
        self.var(self.tree, 'b0_lm_reco_ndof')
        self.var(self.tree, 'b0_lm_reco_vz')

        self.var(self.tree, 'b0_kk_reco_dR')
        self.var(self.tree, 'b0_pi_reco_dxy')
        self.var(self.tree, 'b0_pi_reco_dz')
        self.var(self.tree, 'b0_pi_reco_dxyError')
        self.var(self.tree, 'b0_pi_reco_dzError')
        self.var(self.tree, 'b0_pi_reco_chi2')
        self.var(self.tree, 'b0_pi_reco_ndof')
        self.var(self.tree, 'b0_pi_reco_vz')

        self.var(self.tree, 'b0_k_reco_dxy')
        self.var(self.tree, 'b0_k_reco_dz')
        self.var(self.tree, 'b0_k_reco_dxyError')
        self.var(self.tree, 'b0_k_reco_dzError')
        self.var(self.tree, 'b0_k_reco_chi2')
        self.var(self.tree, 'b0_k_reco_ndof')
        self.var(self.tree, 'b0_k_reco_vz')

        self.var(self.tree, 'b0_reco_chi2')
        self.var(self.tree, 'b0_reco_ndof')
        self.var(self.tree, 'b0_reco_cosAngle')
        self.var(self.tree, 'b0_reco_x')
        self.var(self.tree, 'b0_reco_y')
        self.var(self.tree, 'b0_reco_z')
        self.var(self.tree, 'b0_reco_xError')
        self.var(self.tree, 'b0_reco_yError')
        self.var(self.tree, 'b0_reco_zError')
        self.var(self.tree, 'b0_reco_vtx')
        self.var(self.tree, 'b0_reco_vty')
        self.var(self.tree, 'b0_reco_vtz')


    def process(self, event):
        '''
        '''
        self.readCollections(event.input)
        self.tree.reset()

        if not eval(self.skimFunction):
            return False

        #self.fillEvent(self.tree, event)
        #self.fill(self.tree, 'qscale', event.qscale)
        
        self.fill(self.tree, 'nbmesons', len(event.gen_bmesons))

        # gen tagging muon
        if hasattr(event, 'thetagmu'):
            self.fillGenParticle(self.tree, 'tag_mu', event.thetagmu)
            if hasattr(event.thetagmu, 'reco'):
                self.fillParticle(self.tree, 'tag_mu_reco', event.thetagmu.reco)

        self.fillGenParticle(self.tree, 'b0', event.kstll)
        self.fill(self.tree, 'b0_dr', event.kstll.dr)

        self.fillGenParticle(self.tree, 'b0_lp', event.kstll.lp)
                
        if hasattr(event.kstll.lp, 'reco'):
            self.fillParticle(self.tree, 'b0_lp_reco', event.kstll.lp.reco)       
        self.fillGenParticle(self.tree, 'b0_lm', event.kstll.lm)
        if hasattr(event.kstll.lm, 'reco'):
            self.fillParticle(self.tree, 'b0_lm_reco', event.kstll.lm.reco)
        self.fillGenParticle(self.tree, 'b0_pi', event.kstll.pi)
        if hasattr(event.kstll.pi, 'reco'):
            self.fillParticle(self.tree, 'b0_pi_reco', event.kstll.pi.reco)
        self.fillGenParticle(self.tree, 'b0_k' , event.kstll.k )
        if hasattr(event.kstll.k, 'reco'):
            self.fillParticle(self.tree, 'b0_k_reco', event.kstll.k.reco)

        if hasattr(event.thetagmu, 'reco'):
            self.fill(self.tree, 'tag_mu_reco_dxy', event.thetagmu.reco.bestTrack().dxy(event.pv.position()))
            self.fill(self.tree, 'tag_mu_reco_dz', event.thetagmu.reco.bestTrack().dz(event.pv.position()))
            self.fill(self.tree, 'tag_mu_reco_dxyError', event.thetagmu.reco.bestTrack().dxyError())
            self.fill(self.tree, 'tag_mu_reco_dzError', event.thetagmu.reco.bestTrack().dzError())
            self.fill(self.tree, 'tag_mu_reco_vz', event.thetagmu.reco.bestTrack().vz())
        
        if hasattr(event.kstll.lp, 'reco') and hasattr(event.kstll.lm, 'reco') and hasattr(event.kstll.pi, 'reco') and hasattr(event.kstll.k, 'reco'):
            if hasattr(event, 'myB'):
                self.fill(self.tree, 'foundB0', 1)
            else:
                self.fill(self.tree, 'foundB0', 0.5)
        else:
            self.fill(self.tree, 'foundB0', 0)

        if hasattr(event, 'myB'):
            pv = ROOT.math.XYZPoint(event.myB.bs().x(),event.myB.bs().y(),event.myB.bs().z())

            #self.fill(self.tree, 'foundB0', 1)
            self.fill(self.tree, 'b0_ll_reco_dR', event.myB.llcone())
            self.fill(self.tree, 'b0_lp_reco_dxy', event.myB.l0().dxy(pv))
            self.fill(self.tree, 'b0_lp_reco_dz', event.myB.l0().dz(pv))
            self.fill(self.tree, 'b0_lp_reco_dxyError', event.myB.l0().dxyError())
            self.fill(self.tree, 'b0_lp_reco_dzError', event.myB.l0().dzError())
            self.fill(self.tree, 'b0_lp_reco_chi2', event.myB.l0().chi2())
            self.fill(self.tree, 'b0_lp_reco_ndof', event.myB.l0().ndof())
            self.fill(self.tree, 'b0_lp_reco_vz', event.myB.l0().vz())

            self.fill(self.tree, 'b0_lm_reco_dxy', event.myB.l1().dxy(pv))
            self.fill(self.tree, 'b0_lm_reco_dz', event.myB.l1().dz(pv))
            self.fill(self.tree, 'b0_lm_reco_dxyError', event.myB.l1().dxyError())
            self.fill(self.tree, 'b0_lm_reco_dzError', event.myB.l1().dzError())
            self.fill(self.tree, 'b0_lm_reco_chi2', event.myB.l1().chi2())
            self.fill(self.tree, 'b0_lm_reco_ndof', event.myB.l1().ndof())
            self.fill(self.tree, 'b0_lm_reco_vz', event.myB.l1().vz())

            self.fill(self.tree, 'b0_kk_reco_dR', event.myB.kstcone())
            self.fill(self.tree, 'b0_pi_reco_dxy', event.myB.pi().dxy(pv))
            self.fill(self.tree, 'b0_pi_reco_dz', event.myB.pi().dz(pv))
            self.fill(self.tree, 'b0_pi_reco_dxyError', event.myB.pi().dxyError())
            self.fill(self.tree, 'b0_pi_reco_dzError', event.myB.pi().dzError())
            self.fill(self.tree, 'b0_pi_reco_chi2', event.myB.pi().chi2())
            self.fill(self.tree, 'b0_pi_reco_ndof', event.myB.pi().ndof())
            self.fill(self.tree, 'b0_pi_reco_vz', event.myB.pi().vz())

            self.fill(self.tree, 'b0_k_reco_dxy', event.myB.k().dxy(pv))
            self.fill(self.tree, 'b0_k_reco_dz', event.myB.k().dz(pv))
            self.fill(self.tree, 'b0_k_reco_dxyError', event.myB.k().dxyError())
            self.fill(self.tree, 'b0_k_reco_dzError', event.myB.k().dzError())
            self.fill(self.tree, 'b0_k_reco_chi2', event.myB.k().chi2())
            self.fill(self.tree, 'b0_k_reco_ndof', event.myB.k().ndof())
            self.fill(self.tree, 'b0_k_reco_vz', event.myB.k().vz())

            self.fill(self.tree, 'b0_reco_chi2', event.myB.vtx().chi2())
            self.fill(self.tree, 'b0_reco_ndof', event.myB.vtx().ndof())
            self.fill(self.tree, 'b0_reco_cosAngle', event.myB.vtxcos())
            self.fill(self.tree, 'b0_reco_x', event.myB.vtx().x())
            self.fill(self.tree, 'b0_reco_y', event.myB.vtx().y())
            self.fill(self.tree, 'b0_reco_z', event.myB.vtx().z())
            self.fill(self.tree, 'b0_reco_xError', event.myB.vtx().xError())
            self.fill(self.tree, 'b0_reco_yError', event.myB.vtx().yError())
            self.fill(self.tree, 'b0_reco_zError', event.myB.vtx().zError())
            self.fill(self.tree, 'b0_reco_vtx', event.myB.bs().x())
            self.fill(self.tree, 'b0_reco_vty', event.myB.bs().y())
            self.fill(self.tree, 'b0_reco_vtz', event.myB.bs().z())


        #else: self.fill(self.tree, 'foundB0', 0)
        


        for i, ib in enumerate(event.clean_gen_bmesons[:3]):
            self.fillGenParticle(self.tree, 'b%d' %(i+1), ib)

        if len(event.clean_gen_bmesons)>0:
            self.fill(self.tree, 'dr_btag_bprobe'   , deltaR  (event.clean_gen_bmesons[0], event.kstll)             )
            self.fill(self.tree, 'deta_btag_bprobe' , abs     (event.clean_gen_bmesons[0].eta() - event.kstll.eta()))
            self.fill(self.tree, 'dphi_btag_bprobe' , deltaPhi(event.clean_gen_bmesons[0].phi(), event.kstll.phi()) )
        if hasattr(event, 'thetagmu'):
            self.fill(self.tree, 'dr_mutag_bprobe'  , deltaR  (event.thetagmu, event.kstll)                         )
            self.fill(self.tree, 'deta_mutag_bprobe', abs     (event.thetagmu.eta() - event.kstll.eta())            )
            self.fill(self.tree, 'dphi_mutag_bprobe', deltaPhi(event.thetagmu.phi(), event.kstll.phi())             )


        self.fillTree(event)


