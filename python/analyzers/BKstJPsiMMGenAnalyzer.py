import ROOT
from itertools import product, combinations
import math

from PhysicsTools.Heppy.analyzers.core.Analyzer      import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle    import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar             import deltaR, deltaR2, bestMatch
from PhysicsTools.Heppy.physicsobjects.Muon          import Muon
from PhysicsTools.Heppy.physicsobjects.Electron      import Electron
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject
from CMGTools.BKstLL.analyzers.utils import displacement2D, displacement3D, makeRecoVertex # utility functions
from CMGTools.BKstLL.physicsobjects.BKstLL import BKstLL

from pdb import set_trace

##########################################################################################
# load custom library to ROOT. This contains the kinematic vertex fitter class
ROOT.gSystem.Load('libCMGToolsBKstLL')
from ROOT import KinematicVertexFitter as VertexFitter

class BKstJPsiMMGenAnalyzer(Analyzer):
    '''
    '''
    def declareHandles(self):
        super(BKstJPsiMMGenAnalyzer, self).declareHandles()

        self.handles['L1muons'] = AutoHandle(
            ('gmtStage2Digis', 'Muon'),
            'BXVector<l1t::Muon>'
        )

        self.handles['electrons'] = AutoHandle(
            'slimmedElectrons',
            'std::vector<pat::Electron>'
        )

        self.handles['muons'] = AutoHandle(
            'slimmedMuons',
            'std::vector<pat::Muon>'
        )

        self.handles['losttracks'] = AutoHandle(
            'lostTracks',
            'std::vector<pat::PackedCandidate>'
        )

        self.handles['pfcands'] = AutoHandle(
            'packedPFCandidates',
            'std::vector<pat::PackedCandidate>'
        )

        self.mchandles['prunedGenParticles'] = AutoHandle(
            'prunedGenParticles',
            'std::vector<reco::GenParticle>'
        )

        self.mchandles['packedGenParticles'] = AutoHandle(
            'packedGenParticles',
            'std::vector<pat::PackedGenParticle>'
        )

        self.handles['beamspot'  ] = AutoHandle(
            ('offlineBeamSpot'              , '', 'RECO'),
            'reco::BeamSpot'
        )


    def beginLoop(self, setup):
        super(BKstJPsiMMGenAnalyzer, self).beginLoop(setup)
        self.counters.addCounter('BKstJPsiMMGenAnalyzer')
        count = self.counters.counter('BKstJPsiMMGenAnalyzer')
        count.register('all events')
        count.register('has a good gen B0->KstJPsiMM')

        # stuff I need to instantiate only once
        self.vtxfit = VertexFitter()
        # create a std::vector<reco::RecoChargedCandidate> to be passed to the fitter 
        self.tofit_cc = ROOT.std.vector('reco::RecoChargedCandidate')()
        # create a std::vector<pat::PackedCandidate> to be passed to the fitter 
        self.tofit_pc = ROOT.std.vector('pat::PackedCandidate')()

    def process(self, event):
        self.readCollections(event.input)

        self.counters.counter('BKstJPsiMMGenAnalyzer').inc('all events')

        # vertex stuff
        event.beamspot    = self.handles['beamspot'].product()

        # get the tracks
        allpf      = map(PhysicsObject, self.handles['pfcands'   ].product())
        losttracks = map(PhysicsObject, self.handles['losttracks'].product())

        # merge the track collections
        event.alltracks = sorted([tt for tt in allpf + losttracks if tt.charge() != 0 and abs(tt.pdgId()) not in (11,13)], key = lambda x : x.pt(), reverse = True)

        # get the offline electrons and muons
        event.electrons = map(Electron, self.handles['electrons'].product())
        event.muons     = map(Muon    , self.handles['muons'    ].product())

        pruned_gen_particles = self.mchandles['prunedGenParticles'].product()
        packed_gen_particles = self.mchandles['packedGenParticles'].product()
 
        all_gen_particles = [ip for ip in pruned_gen_particles] + [ip for ip in packed_gen_particles]
         
        # HOOK RM
        event.pruned_gen_particles = pruned_gen_particles
        event.packed_gen_particles = packed_gen_particles
        event.all_gen_particles    = all_gen_particles

        # match gen mu to offline mu
        genmus = [ii for ii in all_gen_particles if abs(ii.pdgId())==13 and ii.status()==1]
        for imu in genmus:
            bm, dr = bestMatch(imu, event.muons)
            if dr<0.3:
                imu.reco = bm
        if len(genmus)>0:
            event.thetagmu = sorted(genmus, key = lambda x : x.pt(), reverse = True)[0]

        # match gen ele to offline ele
        geneles = [ii for ii in all_gen_particles if abs(ii.pdgId())==11 and ii.status()==1]
        for iele in geneles:
            bm, dr = bestMatch(iele, event.electrons)
            if dr<0.3:
                iele.reco = bm
        
        event.gen_bmesons  = [pp for pp in pruned_gen_particles if abs(pp.pdgId()) > 500 and abs(pp.pdgId()) < 600 and pp.isPromptDecayed()]

        event.gen_b0mesons = [pp for pp in pruned_gen_particles if abs(pp.pdgId())==511]
        
        # walk down the lineage of the B mesons and find the final state muons and charged particles
        for ip in event.gen_b0mesons + event.gen_bmesons:
            if getattr(self.cfg_ana, 'verbose', False):
                print 'PdgId : %s   pt : %s  eta : %s   phi : %s' %(ip.pdgId(), ip.pt(), ip.eta(), ip.phi())    
                print '     daughters'
            finaldaughters = []
            finalcharged   = []
            finalmuons     = []
            for ipp in packed_gen_particles:
                mother = ipp.mother(0)
                if mother and self.isAncestor(ip, mother):
                    if abs(ipp.pdgId())==13:
                        finalmuons.append(ipp)
                    if abs(ipp.charge())==1:
                        finalcharged.append(ipp)
                    if getattr(self.cfg_ana, 'verbose', False):
                        print '     PdgId : %s   pt : %s  eta : %s   phi : %s' %(ipp.pdgId(), ipp.pt(), ipp.eta(), ipp.phi())
                    finaldaughters.append(ipp)
                    if getattr(self.cfg_ana, 'verbose', False):
                        print '     PdgId : %s   pt : %s  eta : %s   phi : %s' %(ipp.pdgId(), ipp.pt(), ipp.eta(), ipp.phi())
            ip.finalmuondaughters    = sorted(finalmuons    , key = lambda x : x.pt(), reverse = True)
            ip.finalchargeddaughters = sorted(finalcharged  , key = lambda x : x.pt(), reverse = True)
            ip.finaldaughters        = sorted(finaldaughters, key = lambda x : x.pt(), reverse = True)

        for ib0 in event.gen_b0mesons:
            if   abs(self.cfg_ana.flavour)==11 : togenmatchleptons = event.electrons
            elif abs(self.cfg_ana.flavour)==13 : togenmatchleptons = event.muons
            else                               : print 'you can only pick either pdgId 11 or 13' ; raise
            isit, lp, lm, pi, k = self.isKstLL(ib0, togenmatchleptons, event.alltracks, abs(self.cfg_ana.flavour)) 
            if isit:
                event.kstll    = ib0
                event.kstll.lp = lp
                event.kstll.lm = lm
                event.kstll.pi = pi
                event.kstll.k = k
                event.kstll.dr = min([deltaR(event.kstll, jj) for jj in [lp, lm, pi, k]])


                # if all tracks are reconstructed
                if hasattr(event.kstll.lp, 'reco') and hasattr(event.kstll.lm, 'reco') and hasattr(event.kstll.pi, 'reco') and hasattr(event.kstll.k, 'reco') and event.kstll.pi.reco.bestTrack() and event.kstll.k.reco.bestTrack():

                    # vertex fit
                    # clear the vectors
                    self.tofit_cc.clear()
                    self.tofit_pc.clear()
                    dimu = (event.kstll.lp.reco, event.kstll.lm.reco)
                    # create a RecoChargedCandidate for each reconstructed lepton and flush it into the vector
                    for il in dimu:
                        # if the reco particle is a displaced thing, it does not have the p4() method, so let's build it 
                        myp4 = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzE4D<double> >')(il.px(), il.py(), il.pz(), math.sqrt(il.mass()**2 + il.px()**2 + il.py()**2 + il.pz()**2))
                        ic = ROOT.reco.RecoChargedCandidate() # instantiate a dummy RecoChargedCandidate
                        ic.setCharge(il.charge())             # assign the correct charge
                        ic.setP4(myp4)                        # assign the correct p4                
                        ic.setTrack(il.track())               # set the correct TrackRef
                        if ic.track().isNonnull():            # check that the track is valid
                            self.tofit_cc.push_back(ic)
                    #set_trace()
                    m_k = 0.493677
                    # push the track into the vector
                    self.tofit_pc.push_back(event.kstll.pi.reco.physObj)
                    self.tofit_pc.push_back(event.kstll.k.reco.physObj)

                    # fit it!
                    try:
                        svtree = self.vtxfit.Fit(self.tofit_cc, self.tofit_pc) # actual vertex fitting
                    except:
                        set_trace()
                    # check that the vertex is good
                    if not svtree.get().isEmpty() and svtree.get().isValid():
                        svtree.movePointerToTheTop()
                        sv = svtree.currentDecayVertex().get()
                        recoSv = makeRecoVertex(sv, kinVtxTrkSize=4) # need to do some gymastics
                        event.kstll.sv = recoSv
                        event.myB = BKstLL(dimu[0], dimu[1], event.kstll.pi.reco, event.kstll.k.reco, recoSv, event.beamspot)


                break # yeah, only one at a time, mate!
        
#         if hasattr(event.kstll.lp, 'reco') and hasattr(event.kstll.lm, 'reco'):
#             if event.kstll.lp.reco.pt() ==  event.kstll.lm.reco.pt():
#                 import pdb ; pdb.set_trace()
            
        if not hasattr(event, 'kstll'):
            return False
        self.counters.counter('BKstJPsiMMGenAnalyzer').inc('has a good gen B0->KstJPsiMM')
        
        toclean = None
        
        if event.kstll.isPromptDecayed():
            toclean = event.kstll
        
        else:
            for ip in event.gen_bmesons:
                if self.isAncestor(ip, event.kstll):
                    toclean = ip
                    break
        
#         if toclean is None:
#             print 'nothing to clean lumi %d, ev %d' %(event.lumi, event.eventId)
#             return False
#             import pdb ; pdb.set_trace()
# 
#         self.counters.counter('BKstJPsiMMGenAnalyzer').inc('no fuck ups')
        
        
#         elif abs(event.kstll.mother(0).pdgId())>500 and abs(event.kstll.mother(0).pdgId())<600:
#             if toclean = event.kstll.mother(0)
#         elif abs(event.kstll.mother(0).mother(0).pdgId())>500 and abs(event.kstll.mother(0).mother(0).pdgId())<600:
#         else:
#             print 'nothing to clean lumi %d, ev %d' %(event.lumi, event.eventId)
#             import pdb ; pdb.set_trace()
            
        if toclean is not None:
            event.clean_gen_bmesons = [ib for ib in event.gen_bmesons if ib!=toclean]
        else:
            event.clean_gen_bmesons = event.gen_bmesons 
            
#         import pdb ; pdb.set_trace()
                
        # now find the L1 muons from BX = 0
        L1muons_allbx = self.handles['L1muons'].product()

        L1_muons = []
    
        for jj in range(L1muons_allbx.size(0)):
            L1_muons.append(L1muons_allbx.at(0,jj))
        
        event.L1_muons = L1_muons

        event.selectedLeptons = [] # don't really care about isolated leptons, right?
        
        return True

    @staticmethod
    def isKstLL(b0meson, togenmatchleptons, togenmatchhadrons, flav=13):
#         isB0 = (b0meson.pdgId()==511)
#         isAntiB0 = (not isB0)
        # positive-charged leptons
        lps = [ip for ip in b0meson.finaldaughters if ip.pdgId()==-abs(flav)]
        # negative-charged leptons
        lms = [ip for ip in b0meson.finaldaughters if ip.pdgId()== abs(flav)]
        # pions from K*
        pis = [ip for ip in b0meson.finaldaughters if abs(ip.pdgId())== 211]
        # kaons from K*
        ks  = [ip for ip in b0meson.finaldaughters if abs(ip.pdgId())== 321]

        for ilep in lps+lms:
            # avoid matching the same particle twice, use the charge to distinguish (charge flip not accounted...)
            tomatch = [jj for jj in togenmatchleptons if jj.charge()==ilep.charge()]
            bm, dr = bestMatch(ilep, tomatch)
            if dr<0.3:
                ilep.reco = bm

        for ihad in pis+ks:
            # avoid matching the same particle twice, use the charge to distinguish (charge flip not accounted...)
            tomatch = [jj for jj in togenmatchhadrons if jj.charge()==ihad.charge()]
            bm, dr = bestMatch(ihad, tomatch)
            if dr<0.3:
                ihad.reco = bm
        
        if len(lps) == len(lms) == len(pis) == len(ks) == 1:
            return True, lps[0], lms[0], pis[0], ks[0]
        else:
            return False, None, None, None, None  



    @staticmethod
    def isAncestor(a, p):
        if a == p :
            return True
        for i in xrange(0,p.numberOfMothers()):
            if BKstJPsiMMGenAnalyzer.isAncestor(a,p.mother(i)):
                return True
        return False

    
    
    
    
    
    
    
