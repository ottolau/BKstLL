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

class BKstJPsiEEGenAnalyzer_AOD(Analyzer):
    '''
    '''
    def declareHandles(self):
        super(BKstJPsiEEGenAnalyzer_AOD, self).declareHandles()

        self.handles['L1muons'] = AutoHandle(
            ('gmtStage2Digis', 'Muon'),
            'BXVector<l1t::Muon>'
        )

        self.handles['electrons'] = AutoHandle(
            #'selectedPatElectrons',
            #'std::vector<pat::Electron>'
            'gedGsfElectrons',
            'std::vector<reco::GsfElectron>'

        )

        self.handles['muons'] = AutoHandle(
            #'selectedPatMuons',
            #'std::vector<pat::Muon>'
            'muons',
            'std::vector<reco::Muon>'

        )

        self.handles['tracks'] = AutoHandle(
            'generalTracks',
            'std::vector<reco::Track>'
        )
        
        self.mchandles['GenParticles'] = AutoHandle(
            'genParticles',
            #('GenParticles', '', 'HLT'),
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

        self.handles['pvs'] = AutoHandle(
            'offlinePrimaryVertices',
            'std::vector<reco::Vertex>'
        )

        #self.handles['eletrackMap'] = AutoHandle(('ttk', 'eleTtkMap', 'TTK'), 'std::vector<pair<edm::Ptr<pat::Electron>,reco::Track> >')

    def beginLoop(self, setup):
        super(BKstJPsiEEGenAnalyzer_AOD, self).beginLoop(setup)
        self.counters.addCounter('BKstJPsiEEGenAnalyzer_AOD')
        count = self.counters.counter('BKstJPsiEEGenAnalyzer_AOD')
        count.register('all events')
        count.register('has a good gen B0->KstJPsiEE')

        # stuff I need to instantiate only once
        self.vtxfit = VertexFitter()
        # create a std::vector<reco::RecoChargedCandidate> to be passed to the fitter 
        #self.tofit_cc = ROOT.std.vector('reco::RecoChargedCandidate')()
        self.tofit_cc = ROOT.std.vector('reco::Track')()
        # create a std::vector<pat::PackedCandidate> to be passed to the fitter 
        self.tofit_pc = ROOT.std.vector('reco::Track')()

    def process(self, event):
        self.readCollections(event.input)

        self.counters.counter('BKstJPsiEEGenAnalyzer_AOD').inc('all events')
        
        # vertex stuff
        event.beamspot    = self.handles['beamspot'].product()
        event.pvs = self.handles['pvs'].product()

        # get the tracks
        #tracks = map(PhysicsObject, self.handles['tracks'].product())
        tracks = self.handles['tracks'].product()

        # merge the track collections
        event.alltracks = sorted([tt for tt in tracks if tt.charge() != 0], key = lambda x : x.pt(), reverse = True)

        #print len(event.alltracks)

        # get the offline electrons and muons
        #event.electrons = map(Electron, self.handles['electrons'].product())
        #event.muons     = map(Muon    , self.handles['muons'    ].product())

        event.muons     = self.handles['muons'    ].product()
        event.electrons = self.handles['electrons'].product()

        pruned_gen_particles = self.mchandles['GenParticles'].product()
        #packed_gen_particles = self.mchandles['packedGenParticles'].product()
        packed_gen_particles = [ip for ip in pruned_gen_particles if ip.status()==1]

        all_gen_particles = [ip for ip in pruned_gen_particles] #+ [ip for ip in packed_gen_particles]

        """
        # get the electron-track map BEFORE selections
        eletrks = self.handles['eletrackMap'].product()

        # check consistency between electron-track map and electrons
        if len(event.electrons) != len(eletrks):
            print 'run %d \tlumi %d \tevent %d' %(event.run, event.lumi, event.eventId)
            print 'different number of electrons and ele tracks, %d and %d respectively' %(len(event.electrons), len(eletrks))
            for ii in event.electrons:
                print ii, ii.gsfTrack().pt()
            for ii in eletrks:
                print ii.first.pt(), ii.first.eta(), ii.first.phi(), ii.second.pt()
            set_trace()
        
        for jj in range(len(event.electrons)):
            event.electrons[jj].trackFromGsfTrack = eletrks[jj].second
            event.electrons[jj].ptr = eletrks[jj].first
            if event.electrons[jj].gsfTrack().pt() != eletrks[jj].second.pt():
                set_trace()
        """

         
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
        
        #event.pv = ROOT.math.XYZPoint(0.,0.,0.)
        thetagmuPvDz = 1.e+10
        if hasattr(event.thetagmu, 'reco'):
            for pv in event.pvs:
                if pv.isFake(): continue # work for only AOD
                if abs(pv.z() - event.thetagmu.reco.vz()) < thetagmuPvDz:
                    thetagmuPvDz = abs(pv.z() - event.thetagmu.reco.vz())
                    #event.pv.SetXYZ(pv.x(), pv.y(), pv.z()) 
                    event.pv = pv
        else:
            event.pv = event.pvs[0]
        #event.pv = event.pvs[0]
        if len(genmus)>0 and hasattr(event.thetagmu, 'reco'):
            if (event.thetagmu.reco.track() is not None) and (not event.thetagmu.reco.track().isNull()):
                event.thetagmu.dxy = event.thetagmu.reco.track().dxy(event.pv.position())
                event.thetagmu.dz = event.thetagmu.reco.track().dz(event.pv.position())
                event.thetagmu.dxyError = event.thetagmu.reco.track().dxyError()
                event.thetagmu.dzError = event.thetagmu.reco.track().dzError()
                event.thetagmu.vz = event.thetagmu.reco.track().vz()
            else:
                event.thetagmu.dxy = -99
                event.thetagmu.dz = -99
                event.thetagmu.dxyError = -99
                event.thetagmu.dzError = -99
                event.thetagmu.vz = -99

        # match gen ele to offline ele
        geneles = [ii for ii in all_gen_particles if abs(ii.pdgId())==11 and ii.status()==1]
        for iele in geneles:
            bm, dr = bestMatch(iele, event.electrons)
            if dr<0.3:
                iele.reco = bm
        
        event.gen_bmesons  = [pp for pp in pruned_gen_particles if abs(pp.pdgId()) > 500 and abs(pp.pdgId()) < 600 and pp.isPromptDecayed()]

        event.gen_b0mesons = [pp for pp in pruned_gen_particles if abs(pp.pdgId())==511]
        #print len(event.gen_bmesons)
        #print len(event.gen_b0mesons)

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
            #isit, lp, lm, pi, k = self.isKstLL(ib0, event.alltracks, event.alltracks, abs(self.cfg_ana.flavour)) 
            if isit:
                event.kstll    = ib0
                event.kstll.lp = lp
                event.kstll.lm = lm
                event.kstll.pi = pi
                event.kstll.k = k
                event.kstll.dr = min([deltaR(event.kstll, jj) for jj in [lp, lm, pi, k]])


                # if all tracks are reconstructed
                if hasattr(event.kstll.lp, 'reco') and hasattr(event.kstll.lm, 'reco') and hasattr(event.kstll.pi, 'reco') and hasattr(event.kstll.k, 'reco'):

                    # vertex fit
                    # clear the vectors
                    self.tofit_cc.clear()
                    self.tofit_pc.clear()
                    dimu = (event.kstll.lp.reco, event.kstll.lm.reco)
                    # create a RecoChargedCandidate for each reconstructed lepton and flush it into the vector
                    #for il in dimu:
                        # if the reco particle is a displaced thing, it does not have the p4() method, so let's build it 
                        #myp4 = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzE4D<double> >')(il.px(), il.py(), il.pz(), math.sqrt(il.mass()**2 + il.px()**2 + il.py()**2 + il.pz()**2))
                        #ic = ROOT.reco.RecoChargedCandidate() # instantiate a dummy RecoChargedCandidate
                        #ic.setCharge(il.charge())             # assign the correct charge
                        #ic.setP4(myp4)                        # assign the correct p4                
                        #ic.setTrack(il.track())               # set the correct TrackRef
                        #if ic.track().isNonnull():            # check that the track is valid
                            #self.tofit_cc.push_back(ic)
                    #set_trace()
                    self.tofit_cc.push_back(event.kstll.lp.reco)
                    self.tofit_cc.push_back(event.kstll.lm.reco)
                    m_k = 0.493677
                    m_pi = 0.13957018
                    m_ele = 0.0005109989461
                    # push the track into the vector
                    self.tofit_pc.push_back(event.kstll.pi.reco)
                    self.tofit_pc.push_back(event.kstll.k.reco)

                    # fit it!
                    try:
                        svtree = self.vtxfit.Fit(self.tofit_cc, self.tofit_pc, m_ele, m_pi, m_k) # actual vertex fitting
                    except:
                        set_trace()
                    # check that the vertex is good
                    if not svtree.get().isEmpty() and svtree.get().isValid():
                        svtree.movePointerToTheTop()
                        sv = svtree.currentDecayVertex().get()
                        recoSv = makeRecoVertex(sv, kinVtxTrkSize=4) # need to do some gymastics
                        event.kstll.sv = recoSv
                        event.myB = BKstLL(dimu[0], dimu[1], event.kstll.pi.reco, event.kstll.k.reco, recoSv, event.pv)

                        #Dispvtx = ROOT.GlobalPoint(-1.0*(event.pv.x() - recoSv.x()), -1.0*(event.pv.y() - recoSv.y()), 0.0)
                        #Disperr = sv.error()
                        #event.Lxy = Dispvtx.perp()
                        #event.LxyError = ROOT.TMath.Sqrt(Disperr.rerr(Dispvtx))

                        event.Lxy = self.vtxfit.getLxy(recoSv, event.pv)
                        event.LxyError = self.vtxfit.getLxyError(recoSv, event.pv)

                        #print 'foundB0'

                break # yeah, only one at a time, mate!
        
#         if hasattr(event.kstll.lp, 'reco') and hasattr(event.kstll.lm, 'reco'):
#             if event.kstll.lp.reco.pt() ==  event.kstll.lm.reco.pt():
#                 import pdb ; pdb.set_trace()
            
        if not hasattr(event, 'kstll'):
            print 'no kstll'
            return False
        self.counters.counter('BKstJPsiEEGenAnalyzer_AOD').inc('has a good gen B0->KstJPsiEE')
        
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
#         self.counters.counter('BKstJPsiEEGenAnalyzer_AOD').inc('no fuck ups')
        
        
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
        #print flav
        # positive-charged leptons
        lps = [ip for ip in b0meson.finaldaughters if ip.pdgId()==-abs(flav)]
        # negative-charged leptons
        lms = [ip for ip in b0meson.finaldaughters if ip.pdgId()== abs(flav)]
        # pions from K*
        pis = [ip for ip in b0meson.finaldaughters if abs(ip.pdgId())== 211]
        # kaons from K*
        ks  = [ip for ip in b0meson.finaldaughters if abs(ip.pdgId())== 321]
        #for ii in b0meson.finaldaughters:
            #print ii.pdgId()
        #print len(lps), len(lms), len(pis), len(ks)
        for ilep in lps+lms:
            # avoid matching the same particle twice, use the charge to distinguish (charge flip not accounted...)
            tomatch = [jj for jj in togenmatchhadrons if jj.charge()==ilep.charge()]
            bm, dr = bestMatch(ilep, tomatch)
            #print ilep.pdgId(), dr
            #if dr<0.3:
            ilep.reco = bm

        # match with gsfEelctron
        for ilep in lps+lms:
            # avoid matching the same particle twice, use the charge to distinguish (charge flip not accounted...)
            tomatch = [jj for jj in togenmatchleptons if jj.charge()==ilep.charge()]
            bm, dr = bestMatch(ilep, tomatch)
            #print ilep.pdgId(), dr
            if dr<0.3:
                ilep.PF = bm

        for ihad in pis+ks:
            # avoid matching the same particle twice, use the charge to distinguish (charge flip not accounted...)
            tomatch = [jj for jj in togenmatchhadrons if jj.charge()==ihad.charge()]
            bm, dr = bestMatch(ihad, tomatch)
            #print ihad.pdgId(), dr
            #if dr<0.3:
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
            if BKstJPsiEEGenAnalyzer_AOD.isAncestor(a,p.mother(i)):
                return True
        return False

    
    
    
    
    
    
    
