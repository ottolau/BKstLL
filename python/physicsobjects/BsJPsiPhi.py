import ROOT
import math

from itertools import combinations, product
from copy import deepcopy as dc

from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from PhysicsTools.Heppy.physicsobjects.PhysicsObject import PhysicsObject
from CMGTools.BKstLL.analyzers.utils import makeRecoVertex # utility functions

from ROOT import TVector3, Math

from ROOT import KinematicVertexFitter as VertexFitter

# PDG masses
global m_pi_ch, m_k_ch, m_k_st_zero
m_pi_ch     = 0.13957061
m_k_ch      = 0.493677
m_k_st_zero = 0.89176

class Phi(object):
    '''
    Builds a K*(892) candidate out of two tracks.
    Assigns \pi and K mass hypothesis, find the assignment that returns the
    invariant mass closer to the nominal K*(892) mass.
    As an option, it can be required that the two tracks make a vertex.
    '''
    
    def __init__(self, tk1, tk2, requireVtx=False):
        self.tk1_ = tk1 
        self.tk2_ = tk2 

        self.arbitrate()

        self.vtx     = None
        self.vtxprob = -99.

        if requireVtx:
            # stuff I need to instantiate only once
            self.vtxfit = VertexFitter()
            # create a std::vector<reco::Track> to be passed to the fitter 
            self.tofit = ROOT.std.vector('reco::Track')()
            # fill the vector
            self.tofit.push_back(self.tk1_.physObj)
            self.tofit.push_back(self.tk2_.physObj)
            
            # fit the vertex and save the information
            self.makeVtx_()
                
    def makeVtx_(self):
        # fit it!
        svtree = self.vtxfit.Fit(self.tofit) # actual vertex fitting
        # check that the vertex is good
        if not svtree.get().isEmpty() and svtree.get().isValid():
            svtree.movePointerToTheTop()
            sv = svtree.currentDecayVertex().get()
            recoSv = makeRecoVertex(sv, kinVtxTrkSize=2) # need to do some gymastics
        self.vtx = recoSv
        self.vtxprob = ROOT.TMath.Prob(self.vtx.chi2(), int(self.vtx().ndof()))
   
    def arbitrate(self):
        
        # assign the mass hypotheses
        if self.tk1_.charge() > 0:
            self.kp_ = PhysicsObject(ROOT.pat.PackedCandidate(self.tk1_.physObj))
            self.km_ = PhysicsObject(ROOT.pat.PackedCandidate(self.tk2_.physObj))
        else:
            self.kp_ = PhysicsObject(ROOT.pat.PackedCandidate(self.tk2_.physObj))
            self.km_ = PhysicsObject(ROOT.pat.PackedCandidate(self.tk1_.physObj))
        
        #self.kp_.setMass(m_k_ch)
        #self.km_.setMass(m_k_ch)

        #self.kp_.setPdgId(self.kp_.charge() * 321)
        #self.km_.setPdgId(self.km_ .charge() * 321)

    def charge(self):
        return self.tk1_.charge() + self.tk2_.charge()
    
    def kp(self):
        return self.kp_

    def km(self):
        return self.km_
    
    def p4(self):
        return self.km().p4() + self.kp().p4()
    
    def mass(self):
        return self.p4().mass()

    def eta(self):
        return self.p4().eta()

    def phi(self):
        return self.p4().phi()

    def e(self):
        return self.p4().e()

    def pt(self):
        return self.p4().pt()

    def p(self):
        return math.sqrt(self.p4().px()**2 + self.p4().py()**2 + self.p4().pz()**2)

    def px(self):
        return self.p4().px()

    def py(self):
        return self.p4().py()

    def pz(self):
        return self.p4().pz()
    
    def pdgId():
        # FIXME! get the correct pdgId based on the kaon an pion charges
        return 333
    


class BsJPsiPhi(object):

    ''' 
    B->Kst(Kpi)LL object.
    Add as many methods as convenient, to pre-compute relevant 
    observables, e.g. angular variables
    '''

    def __init__(self, l0, l1, kp, km, vtx=None, beamspot=None):
        self.l0_    = l0   # leading pt lepton
        self.l1_    = l1   # trailing pt displaced lepton
        self.kp_    = kp   # kaon
        self.km_    = km   # pion
        self.vtx_   = vtx  # vertex
        self.bs_    = beamspot 

        m_mu = 0.1056583745
        self.l0p4_ = ROOT.TLorentzVector()
        self.l1p4_ = ROOT.TLorentzVector()
        self.l0p4_.SetPtEtaPhiM(l0.pt(), l0.eta(), l0.phi(), m_mu)
        self.l1p4_.SetPtEtaPhiM(l1.pt(), l1.eta(), l1.phi(), m_mu)

        m_k = 0.493677
        self.kpp4_ = ROOT.TLorentzVector()
        self.kmp4_ = ROOT.TLorentzVector()
        self.kpp4_.SetPtEtaPhiM(kp.pt(), kp.eta(), kp.phi(), m_k)
        self.kmp4_.SetPtEtaPhiM(km.pt(), km.eta(), km.phi(), m_k)

        #self.kp_.setMass(m_k_ch)
        #self.km_.setMass(m_k_ch)

        point = ROOT.reco.Vertex.Point(
            self.bs_.position().x(),
            self.bs_.position().y(),
            self.bs_.position().z(),
        )
        error = self.bs_.covariance3D()
        chi2 = 0.
        ndof = 0.
        bsvtx = ROOT.reco.Vertex(point, error, chi2, ndof, 3) # size? say 3? does it matter?

        self.bsvtx_ = bsvtx 

    # objects
    def l0(self):
        return self.l0_

    def l1(self):
        return self.l1_

    def kp(self):
        return self.kp_

    def km(self):
        return self.km_

    def p4(self):
        return self.kpp4_ + self.kmp4_ + self.l0p4_ + self.l1p4_ 

    def mass(self):
        return self.p4().mass()

    def eta(self):
        return self.p4().eta()

    def phi(self):
        return self.p4().phi()

    def e(self):
        return self.p4().e()

    def pt(self):
        return self.p4().pt()

    def p(self):
        return math.sqrt(self.p4().px()**2 + self.p4().py()**2 + self.p4().pz()**2)

    def px(self):
        return self.p4().px()

    def py(self):
        return self.p4().py()

    def pz(self):
        return self.p4().pz()
    
    def pdgId():
        # FIXME! get the correct pdgId based on the kaon an pion charges
        return 531

    def Phi(self):
        return self.kpp4_ + self.kmp4_

    def charge(self):
        return self.kp().charge() + self.km().charge() + self.l0().charge() + self.l1.charge()

    def bs(self):
        return self.bs_
            
    def vtx(self):
        return self.vtx_

    def b(self):
        return self.l0p4_ + self.l1p4_ + self.kpp4_ + self.kmp4_
   
    def ll(self):
        return self.l0p4_ + self.l1p4_

    def bcone(self):
        return max([deltaR(ii, self.b()) for ii in [self.l0(), self.l1(), self.kp(), self.km()]])

    def llcone(self):
        return deltaR(self.l0(), self.l1())

    def phicone(self):
        return deltaR(self.kp(), self.km())

    def disp2d(self):
        ''' return 2D displacement from beamspot '''
        return ROOT.VertexDistanceXY().distance(self.vtx(), self.bsvtx_)
        
    def ls2d(self):
        ''' return 2D L/sigma '''
        return self.disp2d().significance()

    def vtxprob(self):
        return ROOT.TMath.Prob(self.vtx().chi2(), int(self.vtx().ndof()))

    def vtxchi2(self):
        return self.vtx().chi2()

    def vtxndof(self):
        return int(self.vtx().ndof()) 

    def vtxcos(self):
        perp = ROOT.math.XYZVector(self.b().Px(),
                                   self.b().Py(),
                                   0.)
        
        dxybs = ROOT.GlobalPoint(-1*((self.bs().x0() - self.vtx().x()) + (self.vtx().z() - self.bs().z0()) * self.bs().dxdz()), 
                                 -1*((self.bs().y0() - self.vtx().y()) + (self.vtx().z() - self.bs().z0()) * self.bs().dydz()),
                                  0)
        
        vperp = ROOT.math.XYZVector(dxybs.x(), dxybs.y(), 0.)
        
        cos = vperp.Dot(perp)/(vperp.R()*perp.R())
        
        return cos

