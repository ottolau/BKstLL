import ROOT
from CMGTools.BKstLL.analyzers.L1PurityTreeProducerBase import L1PurityTreeProducerBase
from CMGTools.BKstLL.analyzers.L1Seeds import single_muon, di_muon, tri_muon
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi, bestMatch
from itertools import combinations

class L1GenTreeProducer(L1PurityTreeProducerBase):

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
        self.var(self.tree, 'ndmesons')
        self.var(self.tree, 'njpsis'  )
        self.var(self.tree, 'nupsilon')
        self.var(self.tree, 'nvboson' )
        self.var(self.tree, 'ntops'   )

        # L1 seeds. Bool, either fired or not
        self.var(self.tree, 'L1_SingleMu_22_eta2p1_Q12')

        # passing L1 muons
        self.bookParticle(self.tree, 'mu_eta_2p1_q12_pt22_mu1')
        self.bookParticle(self.tree, 'mu_eta_2p1_q12_pt22_mu2')
        self.bookParticle(self.tree, 'mu_eta_2p1_q12_pt22_mu3')
        self.bookParticle(self.tree, 'mu_eta_2p1_q12_pt22_mu4')

        # the particle that fired the seed
        self.bookGenParticle(self.tree, 'mu_eta_2p1_q12_pt22_particle1')
        self.bookGenParticle(self.tree, 'mu_eta_2p1_q12_pt22_particle2')
        self.bookGenParticle(self.tree, 'mu_eta_2p1_q12_pt22_particle3')
        self.bookGenParticle(self.tree, 'mu_eta_2p1_q12_pt22_particle4')

        # jets matching to L1 muons
        self.bookJet(self.tree, 'mu_eta_2p1_q12_pt22_jet1')
        self.bookJet(self.tree, 'mu_eta_2p1_q12_pt22_jet2')
        self.bookJet(self.tree, 'mu_eta_2p1_q12_pt22_jet3')
        self.bookJet(self.tree, 'mu_eta_2p1_q12_pt22_jet4')

    def process(self, event):
        '''
        '''
        self.readCollections(event.input)
        self.tree.reset()

        if not eval(self.skimFunction):
            return False

        self.fillEvent(self.tree, event)
        self.fill(self.tree, 'qscale', event.qscale)
        
        self.fill(self.tree, 'nbmesons', len(event.gen_bmesons        ))
        self.fill(self.tree, 'ndmesons', len(event.gen_dmesons        ))
        self.fill(self.tree, 'njpsis'  , len(event.gen_prompt_jpsis   ))
        self.fill(self.tree, 'nupsilon', len(event.gen_prompt_upsilons))
        self.fill(self.tree, 'nvboson' , len(event.gen_vbosons        ))
        self.fill(self.tree, 'ntops'   , len(event.gen_topquarks      ))
           
        tomatch = event.gen_bmesons + event.gen_dmesons + event.gen_prompt_jpsis + event.gen_prompt_upsilons + event.gen_vbosons + event.gen_topquarks
        tomatch.sort(key = lambda x : x.pt(), reverse = True)

        fired, matched, index, matched_particle = single_muon(event.L1_muons, 22, 2.1043125, 12, matches=tomatch) 

        if not fired                                    : self.fill(self.tree, 'L1_SingleMu_22_eta2p1_Q12', -1)
        elif index<0                                    : self.fill(self.tree, 'L1_SingleMu_22_eta2p1_Q12',  0)
        elif tomatch[index] in event.gen_bmesons        : self.fill(self.tree, 'L1_SingleMu_22_eta2p1_Q12',  1)
        elif tomatch[index] in event.gen_dmesons        : self.fill(self.tree, 'L1_SingleMu_22_eta2p1_Q12',  2)
        elif tomatch[index] in event.gen_prompt_jpsis   : self.fill(self.tree, 'L1_SingleMu_22_eta2p1_Q12',  3)
        elif tomatch[index] in event.gen_prompt_upsilons: self.fill(self.tree, 'L1_SingleMu_22_eta2p1_Q12',  4)
        elif tomatch[index] in event.gen_vbosons        : self.fill(self.tree, 'L1_SingleMu_22_eta2p1_Q12',  5)
        elif tomatch[index] in event.gen_topquarks      : self.fill(self.tree, 'L1_SingleMu_22_eta2p1_Q12',  6)
        else                                            : pass

        # save the information of the jet matching the L1 muon
        L1_muons_eta_2p1_q12_pt22 = sorted([mu for mu in event.L1_muons if abs(mu.eta())<2.1043125 and mu.hwQual()>=12 and mu.pt()>= 22.], key = lambda mu : mu.pt(), reverse = True)
        for ii, mu in enumerate(L1_muons_eta_2p1_q12_pt22[:4]):
            self.fillParticle(self.tree, 'mu_eta_2p1_q12_pt22_mu%d' %(ii+1), mu)
            if hasattr(mu, 'jet'):
                self.fillJet(self.tree, 'mu_eta_2p1_q12_pt22_jet%d' %(ii+1), mu.jet)
            if index>=0 and hasattr(tomatch[index], 'finalchargeddaughters'):
                if tomatch[index].finalchargeddaughters:
                    self.fillGenParticle(self.tree, 'mu_eta_2p1_q12_pt22_particle%d' %(ii+1), bestMatch(mu, tomatch[index].finalchargeddaughters))
                    
        self.fillTree(event)

