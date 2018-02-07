import ROOT
import numpy as np
from array import array
from copy import deepcopy as dc
from glob import glob

ROOT.TH1.SetDefaultSumw2()

# QCD cross sections
cross_sections = [
    5887580000          , # QCD_Pt10to15    
    1837410000          , # QCD_Pt15to30    
     140932000          , # QCD_Pt30to50    
      19204300          , # QCD_Pt50to80    
       2762530          , # QCD_Pt80to120   
        471100          , # QCD_Pt120to170  
        117276          , # QCD_Pt170to300  
          7823          , # QCD_Pt300to470  
           648.2        , # QCD_Pt470to600  
           186.9        , # QCD_Pt600to800  
            32.293      , # QCD_Pt800to1000 
             9.4183     , # QCD_Pt1000to1400
             0.84265    , # QCD_Pt1400to1800
             0.114943   , # QCD_Pt1800to2400
             0.00682981 , # QCD_Pt2400to3200
             0.000165445, # QCD_Pt3200toInf 
]

bins = [
       10.,
       15., 
       30., 
       50., 
       80., 
      120., 
      170., 
      300., 
      470., 
      600., 
      800., 
     1000., 
     1400., 
     1800., 
     2400., 
     3200., 
    14000.,
]


chain = ROOT.TChain('tree')

# samples = glob('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/*Chunk*/L1PurityTreeProducer/tree.root')

# samples = [
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt15to30/QCD_Pt_15to30/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt30to50/QCD_Pt_30to50/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt50to80/QCD_Pt_50to80/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt80to120/QCD_Pt_80to120/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt120to170/QCD_Pt_120to170/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt170to300/QCD_Pt_170to300/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt300to470/QCD_Pt_300to470/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt470to600/QCD_Pt_470to600/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt800to1000/QCD_Pt_800to1000/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt1000to1400/QCD_Pt_1000to1400/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt1400to1800/QCD_Pt_1400to1800/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt1800to2400/QCD_Pt_1800to2400/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt2400to3200/QCD_Pt_2400to3200/L1GenTreeProducer/tree.root',
#     '/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/pt3200toInf/QCD_Pt_3200toInf/L1GenTreeProducer/tree.root',
# ]

samples = glob('/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v2/*Chunk*/L1GenTreeProducer/tree.root')

print 'loading samples ...'
for sample in samples:
    chain.Add(sample)
'... done'

# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_15to30/L1PurityTreeProducer/tree.root'    )
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_30to50/L1PurityTreeProducer/tree.root'    )
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_50to80/L1PurityTreeProducer/tree.root'    )
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_80to120/L1PurityTreeProducer/tree.root'   )
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_120to170/L1PurityTreeProducer/tree.root'  )
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_170to300/L1PurityTreeProducer/tree.root'  )
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_300to470/L1PurityTreeProducer/tree.root'  )
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_470to600/L1PurityTreeProducer/tree.root'  )
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_600to800/L1PurityTreeProducer/tree.root'  )
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_800to1000/L1PurityTreeProducer/tree.root' )
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_1000to1400/L1PurityTreeProducer/tree.root')
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_1400to1800/L1PurityTreeProducer/tree.root')
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_1800to2400/L1PurityTreeProducer/tree.root')
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_2400to3200/L1PurityTreeProducer/tree.root')
# chain.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/QCD_Pt_3200toInf/L1PurityTreeProducer/tree.root' )

# f1 = ROOT.TFile.Open('purity_tuple_v5.root', 'read')
# f1.cd()
# t1 = f1.Get('tree')

bins   = array('d', bins[1:])
values = array('d', cross_sections[1:])

flat_pthat   = ROOT.TH1F('flat_pthat'  , 'flat_pthat'  , len(bins)-1, bins)
real_pthat   = ROOT.TH1F('real_pthat'  , 'real_pthat'  , len(bins)-1, bins)
weight_pthat = ROOT.TH1F('weight_pthat', 'weight_pthat', len(bins)-1, bins)

chain.Draw('qscale >> flat_pthat')

for ii in range(len(values)):
    real_pthat.SetBinContent(ii+1, values[ii])

flat_pthat.Scale(1./flat_pthat.Integral())
real_pthat.Scale(1./real_pthat.Integral())

weights = []

for ii in range(len(values)):
    weight = real_pthat.GetBinContent(ii+1) / flat_pthat.GetBinContent(ii+1)
    weight_pthat.SetBinContent(ii+1, weight)
    weights.append(weight)

weight_pthat.Draw()
ROOT.gPad.Update()


# reweigh = dc(flat_pthat)
# reweigh.SetTitle('weights')

# reweight.Divide(real_pthat)
# 
# reweight.Draw()



# In [1]: weights
# Out[1]:
# [10.919193711641835,
#  0.81795030487855,
#  0.1102147655818603,
#  0.025108659340147128,
#  0.003993077196987641,
#  0.0009911563312897424,
#  3.553529868585685e-05,
#  3.83303535686388e-06,
#  7.893770288183407e-07,
#  3.715447520938835e-07,
#  1.0828089252003986e-07,
#  9.739178737950565e-09,
#  1.3409608412733624e-09,
#  8.162270429079056e-11,
#  1.9712296083927894e-12]