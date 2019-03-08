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

class Kst(object):
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
        pi1 = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzM4D<double> >')(self.tk1_.px(), self.tk1_.py(), self.tk1_.pz(), m_pi_ch)
        pi2 = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzM4D<double> >')(self.tk2_.px(), self.tk2_.py(), self.tk2_.pz(), m_pi_ch)
        k1  = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzM4D<double> >')(self.tk1_.px(), self.tk1_.py(), self.tk1_.pz(), m_k_ch)
        k2  = ROOT.Math.LorentzVector('<ROOT::Math::PxPyPzM4D<double> >')(self.tk2_.px(), self.tk2_.py(), self.tk2_.pz(), m_k_ch)
        
        # assign the mass hypotheses
        if abs((pi1 + k2).mass() - m_k_st_zero) < abs((pi2 + k1).mass() - m_k_st_zero):
            self.pi_ = PhysicsObject(ROOT.pat.PackedCandidate(self.tk1_.physObj))
            self.k_  = PhysicsObject(ROOT.pat.PackedCandidate(self.tk2_.physObj))
        else:
            self.pi_ = PhysicsObject(ROOT.pat.PackedCandidate(self.tk2_.physObj))
            self.k_  = PhysicsObject(ROOT.pat.PackedCandidate(self.tk1_.physObj))
        
        self.pi_.setMass(m_pi_ch)
        self.k_ .setMass(m_k_ch)

        self.pi_.setPdgId(self.pi_.charge() * 211)
        self.k_ .setPdgId(self.k_ .charge() * 321)

    def charge(self):
        return self.tk1_.charge() + self.tk2_.charge()
    
    def pi(self):
        return self.pi_

    def k(self):
        return self.k_
    
    def p4(self):
        return self.pi().p4() + self.k().p4()
    
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
        return 313 if self.k().charge()>0 else -313
    


class BKstLL(object):

    ''' 
    B->Kst(Kpi)LL object.
    Add as many methods as convenient, to pre-compute relevant 
    observables, e.g. angular variables
    '''

    def __init__(self, l0, l1, k, pi, vtx=None, beamspot=None):
        self.l0_    = l0   # leading pt lepton
        self.l1_    = l1   # trailing pt displaced lepton
        self.k_     = k    # kaon
        self.pi_    = pi   # pion
        self.vtx_   = vtx  # vertex
        self.bs_    = beamspot 

        m_ele = 0.0005109989461
        self.l0p4_ = ROOT.TLorentzVector()
        self.l1p4_ = ROOT.TLorentzVector()
        self.l0p4_.SetPtEtaPhiM(l0.pt(), l0.eta(), l0.phi(), m_ele)
        self.l1p4_.SetPtEtaPhiM(l1.pt(), l1.eta(), l1.phi(), m_ele)

        m_k = 0.493677
        m_pi = 0.13957018
        self.kp4_ = ROOT.TLorentzVector()
        self.pip4_ = ROOT.TLorentzVector()
        self.kp4_.SetPtEtaPhiM(k.pt(), k.eta(), k.phi(), m_k)
        self.pip4_.SetPtEtaPhiM(pi.pt(), pi.eta(), pi.phi(), m_pi)

        #self.k_ .setMass(m_k_ch)
        #self.pi_.setMass(m_pi_ch)

        """
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
        """
        self.bsvtx_ = self.bs_

    # objects
    def l0(self):
        return self.l0_

    def l1(self):
        return self.l1_

    def k(self):
        return self.k_

    def pi(self):
        return self.pi_

    def p4(self):
        return self.pip4_ + self.kp4_ + self.l0p4_ + self.l1p4_ 

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
        return 511 if self.k().charge()>0 else -511

    def kst(self):
        return self.pip4_ + self.kp4_

    def charge(self):
        return self.k().charge() + self.pi().charge() + self.l0().charge() + self.l1.charge()

    def bs(self):
        return self.bs_
            
    def vtx(self):
        return self.vtx_

    def b(self):
        return self.l0p4_ + self.l1p4_ + self.kp4_ + self.pip4_ 
   
    def ll(self):
        return self.l0p4_ + self.l1p4_

    def bcone(self):
        return max([deltaR(ii, self.b()) for ii in [self.l0(), self.l1(), self.k(), self.pi()]])

    def llcone(self):
        return deltaR(self.l0(), self.l1())

    def kstcone(self):
        return deltaR(self.pi(), self.k())

    def disp2d(self):
        ''' return 2D displacement from beamspot '''
        return ROOT.VertexDistanceXY().distance(self.vtx(), self.bsvtx_)
        
    def ls2d(self):
        ''' return 2D L/sigma '''
        return self.disp2d().significance()

    def vtxprob(self):
        return ROOT.TMath.Prob(self.vtx().chi2(), int(self.vtx().ndof()))

    
    def vtxcos(self):
        perp = ROOT.math.XYZVector(self.b().Px(),
                                   self.b().Py(),
                                   0.)
        
        #dxybs = ROOT.GlobalPoint(-1*((self.bs().x0() - self.vtx().x()) + (self.vtx().z() - self.bs().z0()) * self.bs().dxdz()), 
        #                         -1*((self.bs().y0() - self.vtx().y()) + (self.vtx().z() - self.bs().z0()) * self.bs().dydz()),
        #                          0)
 
        dxybs = ROOT.GlobalPoint(-1*(self.bs().x() - self.vtx().x()), 
                                 -1*(self.bs().y() - self.vtx().y()),
                                  0)
       
        vperp = ROOT.math.XYZVector(dxybs.x(), dxybs.y(), 0.)
        
        cos = vperp.Dot(perp)/(vperp.R()*perp.R())
        
        return cos

class BsPhiLL(object):

    ''' 
    B->Kst(Kpi)LL object.
    Add as many methods as convenient, to pre-compute relevant 
    observables, e.g. angular variables
    '''

    def __init__(self, l0, l1, k, pi, vtx=None, beamspot=None):
        self.l0_    = l0   # leading pt lepton
        self.l1_    = l1   # trailing pt displaced lepton
        self.k_     = k    # kaon
        self.pi_    = pi   # pion
        self.vtx_   = vtx  # vertex
        self.bs_    = beamspot 

        m_ele = 0.0005109989461
        self.l0p4_ = ROOT.TLorentzVector()
        self.l1p4_ = ROOT.TLorentzVector()
        self.l0p4_.SetPtEtaPhiM(l0.pt(), l0.eta(), l0.phi(), m_ele)
        self.l1p4_.SetPtEtaPhiM(l1.pt(), l1.eta(), l1.phi(), m_ele)

        m_k = 0.493677
        m_pi = 0.13957018
        self.kp4_ = ROOT.TLorentzVector()
        self.pip4_ = ROOT.TLorentzVector()
        self.kp4_.SetPtEtaPhiM(k.pt(), k.eta(), k.phi(), m_k)
        self.pip4_.SetPtEtaPhiM(pi.pt(), pi.eta(), pi.phi(), m_k)

        #self.k_ .setMass(m_k_ch)
        #self.pi_.setMass(m_pi_ch)

        """
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
        """
        self.bsvtx_ = self.bs_


    # objects
    def l0(self):
        return self.l0_

    def l1(self):
        return self.l1_

    def k(self):
        return self.k_

    def pi(self):
        return self.pi_

    def p4(self):
        return self.pip4_ + self.kp4_ + self.l0p4_ + self.l1p4_ 

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
        return 531 if self.k().charge()>0 else -531

    def kst(self):
        return self.pip4_ + self.kp4_

    def charge(self):
        return self.k().charge() + self.pi().charge() + self.l0().charge() + self.l1.charge()

    def bs(self):
        return self.bs_
            
    def vtx(self):
        return self.vtx_

    def b(self):
        return self.l0p4_ + self.l1p4_ + self.kp4_ + self.pip4_ 
   
    def ll(self):
        return self.l0p4_ + self.l1p4_

    def bcone(self):
        return max([deltaR(ii, self.b()) for ii in [self.l0(), self.l1(), self.k(), self.pi()]])

    def llcone(self):
        return deltaR(self.l0(), self.l1())

    def kstcone(self):
        return deltaR(self.pi(), self.k())

    def disp2d(self):
        ''' return 2D displacement from beamspot '''
        return ROOT.VertexDistanceXY().distance(self.vtx(), self.bsvtx_)
        
    def ls2d(self):
        ''' return 2D L/sigma '''
        return self.disp2d().significance()

    def vtxprob(self):
        return ROOT.TMath.Prob(self.vtx().chi2(), int(self.vtx().ndof()))

    
    def vtxcos(self):
        perp = ROOT.math.XYZVector(self.b().Px(),
                                   self.b().Py(),
                                   0.)
        
        #dxybs = ROOT.GlobalPoint(-1*((self.bs().x0() - self.vtx().x()) + (self.vtx().z() - self.bs().z0()) * self.bs().dxdz()), 
        #                         -1*((self.bs().y0() - self.vtx().y()) + (self.vtx().z() - self.bs().z0()) * self.bs().dydz()),
        #                          0)
 
        dxybs = ROOT.GlobalPoint(-1*(self.bs().x() - self.vtx().x()), 
                                 -1*(self.bs().y() - self.vtx().y()),
                                  0)
       
        vperp = ROOT.math.XYZVector(dxybs.x(), dxybs.y(), 0.)
        
        cos = vperp.Dot(perp)/(vperp.R()*perp.R())
        
        return cos

