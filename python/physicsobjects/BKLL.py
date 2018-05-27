import ROOT
import math

from itertools import combinations, product
# from copy import deepcopy as dc

from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from ROOT import TVector3, Math

from ROOT import KinematicVertexFitter as VertexFitter

class BKLL(object):

    ''' 
    B->KLL object.
    Add as many methods as convenient, to pre-compute relevant 
    observables, e.g. angular variables
    '''

    def __init__(self, l0, l1, k, vtx=None, beamspot=None):
        self.l0_    = l0   # leading pt lepton
        self.l1_    = l1   # trailing pt displaced lepton
        self.k_     = k    # kaon
        self.vtx_   = vtx  # vertex
        self.bs_    = beamspot 

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

    def k(self):
        return self.k_

    def charge(self):
        return self.k().charge() + self.l0().charge() + self.l1.charge()
    
    def vtx(self):
        return self.vtx_

    def b(self):
        return self.l0().p4() + self.l1().p4() + self.k().p4() 
   
    def ll(self):
        return self.l0().p4() + self.l1().p4()

    def bcone(self):
        return max([deltaR(ii, self.b()) for ii in [self.l0(), self.l1(), self.k()]])

    def llcone(self):
        return deltaR(self.l0(), self.l1())

    def disp2d(self):
        ''' return 2D displacement from beamspot '''
        return ROOT.VertexDistanceXY().distance(self.vtx(), self.bsvtx_)
        
    def ls2d(self):
        ''' return 2D L/sigma '''
        return self.disp2d().significance()

    def vtxprob(self):
        return ROOT.TMath.Prob(self.vtx().chi2(), int(self.vtx().ndof()))


'''
        # calculate decay length significance w.r.t. the beamspot
        tauperp = ROOT.math.XYZVector(event.tau3mu.refitTau.px(),
                                      event.tau3mu.refitTau.py(),
                                      0.)
        
        displacementFromBeamspotTau = ROOT.GlobalPoint(-1*((event.bs.x0() - tauvtx.x()) + (tauvtx.z() - event.bs.z0()) * event.bs.dxdz()), 
                                                       -1*((event.bs.y0() - tauvtx.y()) + (tauvtx.z() - event.bs.z0()) * event.bs.dydz()),
                                                        0)
        
        vperptau = ROOT.math.XYZVector(displacementFromBeamspotTau.x(), displacementFromBeamspotTau.y(), 0.)
        
        event.tau3muRefit.refittedVertex.cos = vperptau.Dot(tauperp)/(vperptau.R()*tauperp.R())
'''