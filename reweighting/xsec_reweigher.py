import ROOT
import numpy as np
from array import array
from copy import deepcopy as dc
from glob import glob
from collections import OrderedDict
import pandas, root_numpy
import root_pandas

ROOT.TH1.SetDefaultSumw2()

cross_sections = OrderedDict()
# QCD cross sections
# cross_sections['QCD_Pt_10to15'    ] = 5887580000.         
cross_sections['QCD_Pt_15to30'    ] = 1837410000.         
cross_sections['QCD_Pt_30to50'    ] =  140932000.         
cross_sections['QCD_Pt_50to80'    ] =   19204300.         
cross_sections['QCD_Pt_80to120'   ] =    2762530.         
cross_sections['QCD_Pt_120to170'  ] =     471100.         
cross_sections['QCD_Pt_170to300'  ] =     117276.         
cross_sections['QCD_Pt_300to470'  ] =       7823.         
cross_sections['QCD_Pt_470to600'  ] =        648.2        
cross_sections['QCD_Pt_600to800'  ] =        186.9        
cross_sections['QCD_Pt_800to1000' ] =         32.293      
cross_sections['QCD_Pt_1000to1400'] =          9.4183     
cross_sections['QCD_Pt_1400to1800'] =          0.84265    
cross_sections['QCD_Pt_1800to2400'] =          0.114943   
cross_sections['QCD_Pt_2400to3200'] =          0.00682981 
cross_sections['QCD_Pt_3200toInf' ] =          0.000165445
# EWK cross section
cross_sections['DYJetsM10to50'    ] =      18610.  
cross_sections['DYJetsM50'        ] = 3. *  2008.  
cross_sections['WJetsToLNu'       ] = 3. * 20508.9 

# find all samples in the production round
allsamples = glob('/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v5/*Chunk*/L1GenTreeProducer/tree.root')

lowpts = cross_sections.keys()[:3]

for k, v in cross_sections.iteritems():
    print '\nworking on', k

    # restrict to those that are relevant
    nowsamples = [ii for ii in allsamples if k in ii]

    if len(nowsamples) == 0:
        print '... no root files! Skipping!'
        continue

    if k in lowpts:
        maxchunks = len(nowsamples)     
    else:
        maxchunks = 15
    
    chain = ROOT.TChain('tree')

    print 'loading %d samples ...' %len(nowsamples)
    for i, sample in enumerate(nowsamples[:maxchunks]):
        print '\t\t', i
        chain.Add(sample)
    '... done'

    nevents = chain.GetEntries()
    
    equivalent_lumi = float(nevents) / float(v)
        
    print '\tnsamples %d \tnevents %d\txsec %f' %(len(nowsamples), nevents, v)
    print '\tequivalent lumi = %.6f [pb^-1]' %equivalent_lumi
    
    # normalise to 1 pb-1
    weight = 1./equivalent_lumi
    print '\tweight to get 1pb-1 equivalent yield = %.6f' %weight
    
    print 'loading dataset %s...' %k
    dataset = pandas.DataFrame( root_numpy.root2array(nowsamples[:maxchunks], 'tree') )
    print '\t...done'

    print 'computing weights and adding new column to the dataset...'
    dataset['weight'] = weight
    print '\t...done'

    print 'staging out the enriched dataset...'
    dataset.to_root('/eos/cms/store/group/phys_tau/BKstLL/L1Gen_v5/merged/%s/tree.root' %k, key='tree', store_index=False)
    print '\t...done'

    
    
