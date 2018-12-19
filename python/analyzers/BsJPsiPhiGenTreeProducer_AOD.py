import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from CMGTools.BKstLL.analyzers.L1Seeds import single_muon, di_muon, tri_muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from itertools import combinations
from ROOT import KinematicVertexFitter as VertexFitter

class BsJPsiPhiGenTreeProducer_AOD(L1PurityTreeProducerBase):

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
        
        self.bookGenParticle(self.tree, 'bs')

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
        self.bookGenParticle(self.tree, 'bs_lp') ; self.bookParticle(self.tree, 'bs_lp_reco')
        self.bookGenParticle(self.tree, 'bs_lm') ; self.bookParticle(self.tree, 'bs_lm_reco')
        self.bookGenParticle(self.tree, 'bs_kp') ; self.bookParticle(self.tree, 'bs_kp_reco')
        self.bookGenParticle(self.tree, 'bs_km') ; self.bookParticle(self.tree, 'bs_km_reco' )
        self.var(self.tree, 'bs_dr' )

        self.var(self.tree, 'foundBs')
        self.var(self.tree, 'bs_ll_reco_dR')
        self.var(self.tree, 'bs_lp_reco_dxy')
        self.var(self.tree, 'bs_lp_reco_dz')
        self.var(self.tree, 'bs_lp_reco_dxyError')
        self.var(self.tree, 'bs_lp_reco_dzError')
        self.var(self.tree, 'bs_lp_reco_chi2')
        self.var(self.tree, 'bs_lp_reco_ndof')
        #self.var(self.tree, 'bs_lp_reco_dB')
        #self.var(self.tree, 'bs_lp_reco_SIP')
        self.var(self.tree, 'bs_lp_reco_isLooseMuon')
        self.var(self.tree, 'bs_lp_reco_isMediumMuon')
        self.var(self.tree, 'bs_lp_reco_isTightMuon')

        self.var(self.tree, 'bs_lm_reco_dxy')
        self.var(self.tree, 'bs_lm_reco_dz')
        self.var(self.tree, 'bs_lm_reco_dxyError')
        self.var(self.tree, 'bs_lm_reco_dzError')
        #self.var(self.tree, 'bs_lm_reco_dB')
        #self.var(self.tree, 'bs_lm_reco_SIP')
        self.var(self.tree, 'bs_lm_reco_chi2')
        self.var(self.tree, 'bs_lm_reco_ndof')
        self.var(self.tree, 'bs_lm_reco_isLooseMuon')
        self.var(self.tree, 'bs_lm_reco_isMediumMuon')
        self.var(self.tree, 'bs_lm_reco_isTightMuon')

        self.var(self.tree, 'bs_kk_reco_dR')
        self.var(self.tree, 'bs_kp_reco_dxy')
        self.var(self.tree, 'bs_kp_reco_dz')
        self.var(self.tree, 'bs_kp_reco_dxyError')
        self.var(self.tree, 'bs_kp_reco_dzError')
        self.var(self.tree, 'bs_kp_reco_chi2')
        self.var(self.tree, 'bs_kp_reco_ndof')

        self.var(self.tree, 'bs_km_reco_dxy')
        self.var(self.tree, 'bs_km_reco_dz')
        self.var(self.tree, 'bs_km_reco_dxyError')
        self.var(self.tree, 'bs_km_reco_dzError')
        self.var(self.tree, 'bs_km_reco_chi2')
        self.var(self.tree, 'bs_km_reco_ndof')

        self.var(self.tree, 'bs_reco_chi2')
        self.var(self.tree, 'bs_reco_ndof')
        self.var(self.tree, 'bs_reco_cosAngle')
        self.var(self.tree, 'bs_reco_x')
        self.var(self.tree, 'bs_reco_y')
        self.var(self.tree, 'bs_reco_z')
        self.var(self.tree, 'bs_reco_xError')
        self.var(self.tree, 'bs_reco_yError')
        self.var(self.tree, 'bs_reco_zError')
        self.var(self.tree, 'bs_reco_vtx')
        self.var(self.tree, 'bs_reco_vty')
        self.var(self.tree, 'bs_reco_vtz')


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

        self.fillGenParticle(self.tree, 'bs', event.jpsiphi)
        self.fill(self.tree, 'bs_dr', event.jpsiphi.dr)

        self.fillGenParticle(self.tree, 'bs_lp', event.jpsiphi.lp)
                
        if hasattr(event.jpsiphi.lp, 'reco'):
            self.fillParticle(self.tree, 'bs_lp_reco', event.jpsiphi.lp.reco)       
        self.fillGenParticle(self.tree, 'bs_lm', event.jpsiphi.lm)
        if hasattr(event.jpsiphi.lm, 'reco'):
            self.fillParticle(self.tree, 'bs_lm_reco', event.jpsiphi.lm.reco)
        self.fillGenParticle(self.tree, 'bs_kp', event.jpsiphi.kp)
        if hasattr(event.jpsiphi.kp, 'reco'):
            self.fillParticle(self.tree, 'bs_kp_reco', event.jpsiphi.kp.reco)
        self.fillGenParticle(self.tree, 'bs_km' , event.jpsiphi.km )
        if hasattr(event.jpsiphi.km, 'reco'):
            self.fillParticle(self.tree, 'bs_km_reco', event.jpsiphi.km.reco)
    
        if hasattr(event, 'myB'):
            self.vtxfit = VertexFitter()

            pv = ROOT.math.XYZPoint(event.myB.bs().x0(),event.myB.bs().y0(),event.myB.bs().z0())
            vtxpoint = ROOT.reco.Vertex.Point(event.myB.bs().x0(),event.myB.bs().y0(),event.myB.bs().z0())
            vtxerror = event.myB.vtx().error()
            recoVtx = ROOT.reco.Vertex(vtxpoint, vtxerror)

            self.fill(self.tree, 'foundBs', 1)
            self.fill(self.tree, 'bs_ll_reco_dR', event.myB.llcone())
            self.fill(self.tree, 'bs_lp_reco_dxy', event.myB.l0().muonBestTrack().dxy(event.myB.bs()))
            self.fill(self.tree, 'bs_lp_reco_dz', event.myB.l0().muonBestTrack().dz(pv))
            self.fill(self.tree, 'bs_lp_reco_dxyError', event.myB.l0().muonBestTrack().dxyError())
            self.fill(self.tree, 'bs_lp_reco_dzError', event.myB.l0().muonBestTrack().dzError())
            #self.fill(self.tree, 'bs_lp_reco_SIP', abs(event.myB.l0().dB(ROOT.pat.Muon.PV3D))/event.myB.l0().edB(ROOT.pat.Muon.PV3D))
            self.fill(self.tree, 'bs_lp_reco_chi2', event.myB.l0().muonBestTrack().chi2())
            self.fill(self.tree, 'bs_lp_reco_ndof', event.myB.l0().muonBestTrack().ndof())

            # Only works for >CMSSW_9_4_X
            #self.fill(self.tree, 'bs_lp_reco_isLooseMuon', 1 if event.myB.l0().passed(0) else 0)
            #self.fill(self.tree, 'bs_lp_reco_isMediumMuon', 1 if event.myB.l0().passed(1) else 0)
            #self.fill(self.tree, 'bs_lp_reco_isTightMuon', 1 if event.myB.l0().passed(3) else 0)

            self.fill(self.tree, 'bs_lp_reco_isLooseMuon', 1 if self.vtxfit.checkIsLooseMuon(event.myB.l0()) else 0)
            self.fill(self.tree, 'bs_lp_reco_isMediumMuon', 1 if self.vtxfit.checkIsMediumMuon(event.myB.l0()) else 0)
            self.fill(self.tree, 'bs_lp_reco_isTightMuon', 1 if self.vtxfit.checkIsTightMuon(event.myB.l0(), recoVtx) else 0)

            self.fill(self.tree, 'bs_lm_reco_dxy', event.myB.l1().muonBestTrack().dxy(event.myB.bs()))
            self.fill(self.tree, 'bs_lm_reco_dz', event.myB.l1().muonBestTrack().dz(pv))
            self.fill(self.tree, 'bs_lm_reco_dxyError', event.myB.l1().muonBestTrack().dxyError())
            self.fill(self.tree, 'bs_lm_reco_dzError', event.myB.l1().muonBestTrack().dzError())
            #self.fill(self.tree, 'bs_lm_reco_SIP', abs(event.myB.l1().dB(ROOT.pat.Muon.PV3D))/event.myB.l1().edB(ROOT.pat.Muon.PV3D))
            self.fill(self.tree, 'bs_lm_reco_chi2', event.myB.l1().muonBestTrack().chi2())
            self.fill(self.tree, 'bs_lm_reco_ndof', event.myB.l1().muonBestTrack().ndof())

            # Only works for >CMSSW_9_4_X
            #self.fill(self.tree, 'bs_lm_reco_isLooseMuon', 1 if event.myB.l1().passed(0) else 0)
            #self.fill(self.tree, 'bs_lm_reco_isMediumMuon', 1 if event.myB.l1().passed(1) else 0)
            #self.fill(self.tree, 'bs_lm_reco_isTightMuon', 1 if event.myB.l1().passed(3) else 0)

            self.fill(self.tree, 'bs_lm_reco_isLooseMuon', 1 if self.vtxfit.checkIsLooseMuon(event.myB.l1()) else 0)
            self.fill(self.tree, 'bs_lm_reco_isMediumMuon', 1 if self.vtxfit.checkIsLooseMuon(event.myB.l1()) else 0)
            self.fill(self.tree, 'bs_lm_reco_isTightMuon', 1 if self.vtxfit.checkIsTightMuon(event.myB.l1(), recoVtx) else 0)

            self.fill(self.tree, 'bs_kk_reco_dR', event.myB.phicone())
            self.fill(self.tree, 'bs_kp_reco_dxy', event.myB.kp().dxy(event.myB.bs()))
            self.fill(self.tree, 'bs_kp_reco_dz', event.myB.kp().dz(pv))
            self.fill(self.tree, 'bs_kp_reco_dxyError', event.myB.kp().dxyError())
            self.fill(self.tree, 'bs_kp_reco_dzError', event.myB.kp().dzError())
            self.fill(self.tree, 'bs_kp_reco_chi2', event.myB.kp().chi2())
            self.fill(self.tree, 'bs_kp_reco_ndof', event.myB.kp().ndof())
            self.fill(self.tree, 'bs_km_reco_dxy', event.myB.km().dxy(event.myB.bs()))
            self.fill(self.tree, 'bs_km_reco_dz', event.myB.km().dz(pv))
            self.fill(self.tree, 'bs_km_reco_dxyError', event.myB.km().dxyError())
            self.fill(self.tree, 'bs_km_reco_dzError', event.myB.km().dzError())
            self.fill(self.tree, 'bs_km_reco_chi2', event.myB.km().chi2())
            self.fill(self.tree, 'bs_km_reco_ndof', event.myB.km().ndof())

            self.fill(self.tree, 'bs_reco_chi2', event.myB.vtx().chi2())
            self.fill(self.tree, 'bs_reco_ndof', event.myB.vtx().ndof())
            self.fill(self.tree, 'bs_reco_cosAngle', event.myB.vtxcos())
            self.fill(self.tree, 'bs_reco_x', event.myB.vtx().x())
            self.fill(self.tree, 'bs_reco_y', event.myB.vtx().y())
            self.fill(self.tree, 'bs_reco_z', event.myB.vtx().z())
            self.fill(self.tree, 'bs_reco_xError', event.myB.vtx().xError())
            self.fill(self.tree, 'bs_reco_yError', event.myB.vtx().yError())
            self.fill(self.tree, 'bs_reco_zError', event.myB.vtx().zError())
            self.fill(self.tree, 'bs_reco_vtx', event.myB.bs().x0())
            self.fill(self.tree, 'bs_reco_vty', event.myB.bs().y0())
            self.fill(self.tree, 'bs_reco_vtz', event.myB.bs().z0())


        else: self.fill(self.tree, 'foundBs', 0)
        


        for i, ib in enumerate(event.clean_gen_bmesons[:3]):
            self.fillGenParticle(self.tree, 'b%d' %(i+1), ib)

        if len(event.clean_gen_bmesons)>0:
            self.fill(self.tree, 'dr_btag_bprobe'   , deltaR  (event.clean_gen_bmesons[0], event.jpsiphi)             )
            self.fill(self.tree, 'deta_btag_bprobe' , abs     (event.clean_gen_bmesons[0].eta() - event.jpsiphi.eta()))

        self.fillTree(event)


