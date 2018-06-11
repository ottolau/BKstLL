#! /usr/bin/env python

import ROOT
ROOT.gSystem.Load('libCMGToolsBKstLL')

import sys
from DataFormats.FWLite import Events, Handle

# Make VarParsing object
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideAboutPythonConfigFile#VarParsing_Example
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()

# Events takes either
# - single file name
# - list of file names
# - VarParsing options

# use Varparsing object
events = Events ('output_v2.root')

# create handle outside of loop
handle_ele  = Handle ("std::vector<pat::Electron>")
handle_map  = Handle ("vector<pair<edm::Ptr<pat::Electron>,reco::Track>>")

# for now, label is just a tuple of strings that is initialized just
# like and edm::InputTag
label_ele = ("slimmedElectrons")
label_map = ("ttk", "eleTtkMap", "TTK")

# Create histograms, etc.
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background

# loop over events
for ii, event in enumerate(events):
    if event.eventAuxiliary().run() != 316060:
         continue
    if event.eventAuxiliary().luminosityBlock() != 306:
         continue
    if event.eventAuxiliary().event() != 291325724:
         continue
    print ii
    # use getByLabel, just like in cmsRun
    event.getByLabel (label_ele, handle_ele)
    event.getByLabel (label_map, handle_map)
    # get the product
    eles = handle_ele.product()
    maps = handle_map.product()
    
    if len(eles)>0:
        import pdb ; pdb.set_trace()
    










