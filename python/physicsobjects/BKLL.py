import math

from itertools import combinations, product
# from copy import deepcopy as dc

from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from ROOT import TVector3, Math


class BKLL(object):

    ''' 
    B->KLL object.
    Add as many methods as convenient, to pre-compute relevant 
    observables, e.g. angular variables
    '''

    def __init__(self, l0, l1, k, vtx=None):
        self.l0_  = l0   # leading pt lepton
        self.l1_  = l1   # trailing pt displaced lepton
        self.k_   = k    # kaon
        self.vtx_ = vtx  # vertex

    # objects
    def l0(self):
        return self.l0_

    def l1(self):
        return self.l1_

    def k(self):
        return self.k_

    def vtx(self):
        return self.vtx_

    def b(self):
        return l0.p4() + l1.p4() + k.p4() 
   
    def ll(self):
        return l0.p4() + l1.p4()
