import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from CMGTools.BKstLL.analyzers.L1Seeds import single_muon, di_muon, tri_muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from itertools import combinations

class BKstJPsiMMGenTreeProducer(L1PurityTreeProducerBase):

    '''
    '''

    def declareVariables(self, setup):
        '''
        '''
        self.bookEvent(self.tree)

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

        self.var(self.tree, 'foundB0')
        self.var(self.tree, 'b0_ll_reco_dR')
        self.var(self.tree, 'b0_lp_reco_dxy')
        self.var(self.tree, 'b0_lp_reco_dz')
        self.var(self.tree, 'b0_lp_reco_dB')
        self.var(self.tree, 'b0_lp_reco_SIP')
        self.var(self.tree, 'b0_lp_reco_chi2')
        self.var(self.tree, 'b0_lp_reco_ndof')
        self.var(self.tree, 'b0_lp_reco_idSelectorBit')
        self.var(self.tree, 'b0_lp_reco_idPatMVA')
        self.var(self.tree, 'b0_lp_reco_idPatSoftMVA')

        self.var(self.tree, 'b0_lm_reco_dxy')
        self.var(self.tree, 'b0_lm_reco_dz')
        self.var(self.tree, 'b0_lm_reco_dB')
        self.var(self.tree, 'b0_lm_reco_SIP')
        self.var(self.tree, 'b0_lm_reco_chi2')
        self.var(self.tree, 'b0_lm_reco_ndof')
        self.var(self.tree, 'b0_lm_reco_idSelectorBit')
        self.var(self.tree, 'b0_lm_reco_idPatMVA')
        self.var(self.tree, 'b0_lm_reco_idPatSoftMVA')

        self.var(self.tree, 'b0_kk_reco_dR')
        self.var(self.tree, 'b0_pi_reco_dxy')
        self.var(self.tree, 'b0_pi_reco_dz')
        self.var(self.tree, 'b0_pi_reco_dxyError')
        self.var(self.tree, 'b0_pi_reco_dzError')
        self.var(self.tree, 'b0_pi_reco_chi2')
        self.var(self.tree, 'b0_pi_reco_ndof')

        self.var(self.tree, 'b0_k_reco_dxy')
        self.var(self.tree, 'b0_k_reco_dz')
        self.var(self.tree, 'b0_k_reco_dxyError')
        self.var(self.tree, 'b0_k_reco_dzError')
        self.var(self.tree, 'b0_k_reco_chi2')
        self.var(self.tree, 'b0_k_reco_ndof')

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

        self.fillEvent(self.tree, event)
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

        
        if hasattr(event.kstll.lp, 'reco') and hasattr(event.kstll.lm, 'reco') and hasattr(event.kstll.pi, 'reco') and hasattr(event.kstll.k, 'reco'):
            if hasattr(event, 'myB'):
                self.fill(self.tree, 'foundB0', 1)
            else:
                self.fill(self.tree, 'foundB0', 0.5)
        else:
            self.fill(self.tree, 'foundB0', 0)

        if hasattr(event, 'myB'):
            pv = ROOT.math.XYZPoint(event.myB.bs().x0(),event.myB.bs().y0(),event.myB.bs().z0())

            #self.fill(self.tree, 'foundB0', 1)
            self.fill(self.tree, 'b0_ll_reco_dR', event.myB.llcone())
            self.fill(self.tree, 'b0_lp_reco_dxy', event.myB.l0().muonBestTrack().dxy(event.myB.bs()))
            self.fill(self.tree, 'b0_lp_reco_dz', event.myB.l0().muonBestTrack().dz(pv))
            self.fill(self.tree, 'b0_lp_reco_SIP', abs(event.myB.l0().dB(ROOT.pat.Muon.PV3D))/event.myB.l0().edB(ROOT.pat.Muon.PV3D))
            self.fill(self.tree, 'b0_lp_reco_chi2', event.myB.l0().muonBestTrack().chi2())
            self.fill(self.tree, 'b0_lp_reco_ndof', event.myB.l0().muonBestTrack().ndof())
            self.fill(self.tree, 'b0_lp_reco_idSelectorBit', event.myB.l0().selectors())
            self.fill(self.tree, 'b0_lp_reco_idPatMVA', event.myB.l0().mvaValue())
            self.fill(self.tree, 'b0_lp_reco_idPatSoftMVA', event.myB.l0().softMvaValue())

            self.fill(self.tree, 'b0_lm_reco_dxy', event.myB.l1().muonBestTrack().dxy(event.myB.bs()))
            self.fill(self.tree, 'b0_lm_reco_dz', event.myB.l1().muonBestTrack().dz(pv))
            self.fill(self.tree, 'b0_lm_reco_SIP', abs(event.myB.l1().dB(ROOT.pat.Muon.PV3D))/event.myB.l1().edB(ROOT.pat.Muon.PV3D))
            self.fill(self.tree, 'b0_lm_reco_chi2', event.myB.l1().muonBestTrack().chi2())
            self.fill(self.tree, 'b0_lm_reco_ndof', event.myB.l1().muonBestTrack().ndof())
            self.fill(self.tree, 'b0_lm_reco_idSelectorBit', event.myB.l1().selectors())
            self.fill(self.tree, 'b0_lm_reco_idPatMVA', event.myB.l1().mvaValue())
            self.fill(self.tree, 'b0_lm_reco_idPatSoftMVA', event.myB.l1().softMvaValue())

            self.fill(self.tree, 'b0_kk_reco_dR', event.myB.kstcone())
            self.fill(self.tree, 'b0_pi_reco_dxy', event.myB.pi().bestTrack().dxy(event.myB.bs()))
            self.fill(self.tree, 'b0_pi_reco_dz', event.myB.pi().bestTrack().dz(pv))
            self.fill(self.tree, 'b0_pi_reco_dxyError', event.myB.pi().bestTrack().dxyError())
            self.fill(self.tree, 'b0_pi_reco_dzError', event.myB.pi().bestTrack().dzError())
            self.fill(self.tree, 'b0_pi_reco_chi2', event.myB.pi().bestTrack().chi2())
            self.fill(self.tree, 'b0_pi_reco_ndof', event.myB.pi().bestTrack().ndof())
            self.fill(self.tree, 'b0_k_reco_dxy', event.myB.k().bestTrack().dxy(event.myB.bs()))
            self.fill(self.tree, 'b0_k_reco_dz', event.myB.k().bestTrack().dz(pv))
            self.fill(self.tree, 'b0_k_reco_dxyError', event.myB.k().bestTrack().dxyError())
            self.fill(self.tree, 'b0_k_reco_dzError', event.myB.k().bestTrack().dzError())
            self.fill(self.tree, 'b0_k_reco_chi2', event.myB.k().bestTrack().chi2())
            self.fill(self.tree, 'b0_k_reco_ndof', event.myB.k().bestTrack().ndof())

            self.fill(self.tree, 'b0_reco_chi2', event.myB.vtx().chi2())
            self.fill(self.tree, 'b0_reco_ndof', event.myB.vtx().ndof())
            self.fill(self.tree, 'b0_reco_cosAngle', event.myB.vtxcos())
            self.fill(self.tree, 'b0_reco_x', event.myB.vtx().x())
            self.fill(self.tree, 'b0_reco_y', event.myB.vtx().y())
            self.fill(self.tree, 'b0_reco_z', event.myB.vtx().z())
            self.fill(self.tree, 'b0_reco_xError', event.myB.vtx().xError())
            self.fill(self.tree, 'b0_reco_yError', event.myB.vtx().yError())
            self.fill(self.tree, 'b0_reco_zError', event.myB.vtx().zError())
            self.fill(self.tree, 'b0_reco_vtx', event.myB.bs().x0())
            self.fill(self.tree, 'b0_reco_vty', event.myB.bs().y0())
            self.fill(self.tree, 'b0_reco_vtz', event.myB.bs().z0())


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


