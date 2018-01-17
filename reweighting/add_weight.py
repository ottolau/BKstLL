import ROOT
import numpy as np
import pandas, root_numpy
import root_pandas
from glob import glob

@np.vectorize
def ptHatWeight(qscale):
    if   qscale > 3200.: return 1.8956843890843723e-12
    elif qscale > 2400.: return 8.160483583860695e-11
    elif qscale > 1800.: return 1.3557068341766136e-09
    elif qscale > 1400.: return 9.737048034862061e-09
    elif qscale > 1000.: return 1.0825719822589222e-07
    elif qscale >  800.: return 3.7146342104002605e-07
    elif qscale >  600.: return 7.89204251750761e-07
    elif qscale >  470.: return 3.832196524686683e-06
    elif qscale >  300.: return 3.5527522981988425e-05
    elif qscale >  170.: return 0.0009909394626551465
    elif qscale >  120.: return 0.003992203477653263
    elif qscale >   80.: return 0.025103166027266692
    elif qscale >   50.: return 0.11019064968756105
    elif qscale >   30.: return 0.831341988631457
    elif qscale >   15.: return 10.916804834723132
    else               : return None

samples = glob('/eos/cms/store/group/phys_tau/BKstLL/L1Purity_v3/*Chunk*/L1PurityTreeProducer/tree.root')

for ii in samples:

    print 'loading dataset %s...' %ii
    dataset = pandas.DataFrame( root_numpy.root2array(ii, 'tree') )
    print '\t...done'

    print 'computing weights and adding new column to the dataset...'
    dataset['qscale_weight'] = ptHatWeight(dataset.qscale)
    print '\t...done'

    print 'staging out the enriched dataset...'
    dataset.to_root(ii.replace('tree.root', 'tree_enriched.root'), key='tree', store_index=False)
    print '\t...done'

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
