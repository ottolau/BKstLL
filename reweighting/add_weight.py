import ROOT
import numpy as np
import pandas, root_numpy
import root_pandas
from glob import glob

@np.vectorize
def ptHatWeight(qscale):
    if   qscale > 3200.: return 1.9712296083927894e-12
    elif qscale > 2400.: return 8.162270429079056e-11
    elif qscale > 1800.: return 1.3409608412733624e-09
    elif qscale > 1400.: return 9.739178737950565e-09
    elif qscale > 1000.: return 1.0828089252003986e-07
    elif qscale >  800.: return 3.715447520938835e-07
    elif qscale >  600.: return 7.893770288183407e-07
    elif qscale >  470.: return 3.83303535686388e-06
    elif qscale >  300.: return 3.553529868585685e-05
    elif qscale >  170.: return 0.0009911563312897424
    elif qscale >  120.: return 0.003993077196987641
    elif qscale >   80.: return 0.025108659340147128
    elif qscale >   50.: return 0.1102147655818603
    elif qscale >   30.: return 0.81795030487855
    elif qscale >   15.: return 10.919193711641835
    else               : return None

samples = glob('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/*Chunk*/L1PurityTreeProducer/tree.root')

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

