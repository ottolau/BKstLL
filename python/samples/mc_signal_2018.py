import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator

creator = ComponentCreator()

BdKstEE_ttk = creator.makeMCComponentFromLocal(
    'BdKstEE_ttk', 
    'XXX', 
    path = '/eos/uscms/store/user/klau/BKstPsiEEMuFilter/', 
    #pattern='.*root', 
    pattern='.*ttk_debug.*root',
    xSec=1.0, 
)



