import ROOT

ROOT.TH1.SetDefaultSumw2()

# ROOT.gROOT.SetBatch(True)

print 'loading chain...'
t1 = ROOT.TChain('tree')
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_15to30/L1PurityTreeProducer/tree_enriched.root'    )
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_30to50/L1PurityTreeProducer/tree_enriched.root'    )
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_50to80/L1PurityTreeProducer/tree_enriched.root'    )
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_80to120/L1PurityTreeProducer/tree_enriched.root'   ) 
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_120to170/L1PurityTreeProducer/tree_enriched.root'  ) 
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_170to300/L1PurityTreeProducer/tree_enriched.root'  )
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_300to470/L1PurityTreeProducer/tree_enriched.root'  )
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_470to600/L1PurityTreeProducer/tree_enriched.root'  )
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_600to800/L1PurityTreeProducer/tree_enriched.root'  )
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_800to1000/L1PurityTreeProducer/tree_enriched.root' )
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_1000to1400/L1PurityTreeProducer/tree_enriched.root')
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_1400to1800/L1PurityTreeProducer/tree_enriched.root')
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_1800to2400/L1PurityTreeProducer/tree_enriched.root')
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_2400to3200/L1PurityTreeProducer/tree_enriched.root')
t1.Add('/eos/cms/store/group/phys_tau/BKstLL/L1Purity/QCD_Pt_3200toInf/L1PurityTreeProducer/tree_enriched.root' )
print '\t ...done'

h1 = ROOT.TH1F('h1', '', 2, 0, 2)
h1.SetMinimum(0.)

print '\n\n'
print 'B meson purity out of events firing the trigger'

t1.Draw('(matched_L1_SingleMu_22_eta2p1_Q12         & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_22_eta2p1_Q12         & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_22_eta2p1_Q12          %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_22_eta2p1_q12.pdf'        )        
t1.Draw('(matched_L1_SingleMu_25_Q12                & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_25_Q12                & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_25_Q12                 %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_25_q12.pdf'               )               
t1.Draw('(matched_L1_SingleMu_25_eta1p0_Q12         & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_25_eta1p0_Q12         & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_25_eta1p0_Q12          %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_25_eta1p0_q12.pdf'        )        
t1.Draw('(matched_L1_SingleMu_10_eta1p0_Q8          & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_10_eta1p0_Q8          & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_10_eta1p0_Q8           %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_10_eta1p0_q8.pdf'         )         
t1.Draw('(matched_L1_SingleMu_10_eta1p0_Q12         & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_10_eta1p0_Q12         & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_10_eta1p0_Q12          %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_10_eta1p0_q12.pdf'        )        
t1.Draw('(matched_L1_SingleMu_10_eta1p5_Q8          & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_10_eta1p5_Q8          & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_10_eta1p5_Q8           %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_10_eta1p5_q8.pdf'         )         
t1.Draw('(matched_L1_SingleMu_10_eta1p5_Q12         & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_10_eta1p5_Q12         & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_10_eta1p5_Q12          %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_10_eta1p5_q12.pdf'        )        
t1.Draw('(matched_L1_SingleMu_15_eta1p0_Q8          & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_15_eta1p0_Q8          & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_15_eta1p0_Q8           %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_15_eta1p0_q8.pdf'         )         
t1.Draw('(matched_L1_SingleMu_15_eta1p0_Q12         & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_15_eta1p0_Q12         & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_15_eta1p0_Q12          %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_15_eta1p0_q12.pdf'        )        
t1.Draw('(matched_L1_SingleMu_15_eta1p5_Q8          & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_15_eta1p5_Q8          & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_15_eta1p5_Q8           %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_15_eta1p5_q8.pdf'         )         
t1.Draw('(matched_L1_SingleMu_15_eta1p5_Q12         & nbmesons>=2) >> h1', '(qscale_weight) * (L1_SingleMu_15_eta1p5_Q12         & qscale > 20)', 'norm hist e') ; print 'Purity of L1_SingleMu_15_eta1p5_Q12          %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_singlemu_15_eta1p5_q12.pdf'        )        
t1.Draw('(matched_L1_DoubleMu_15_7_Q8               & nbmesons>=2) >> h1', '(qscale_weight) * (L1_DoubleMu_15_7_Q8               & qscale > 20)', 'norm hist e') ; print 'Purity of L1_DoubleMu_15_7_Q8                %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_doublemu_15_7_q8.pdf'              )              
t1.Draw('(matched_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 & nbmesons>=2) >> h1', '(qscale_weight) * (L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 & qscale > 20)', 'norm hist e') ; print 'Purity of L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4  %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_doublemu0er1p5_sq_os_dr_max1p4.pdf')
t1.Draw('(matched_L1_DoubleMu4_SQ_OS_dR_Max1p2      & nbmesons>=2) >> h1', '(qscale_weight) * (L1_DoubleMu4_SQ_OS_dR_Max1p2      & qscale > 20)', 'norm hist e') ; print 'Purity of L1_DoubleMu4_SQ_OS_dR_Max1p2       %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_doublemu4_sq_os_dr_max1p2.pdf'     )     
t1.Draw('(matched_L1_DoubleMu4p5_SQ_OS_dR_Max1p2    & nbmesons>=2) >> h1', '(qscale_weight) * (L1_DoubleMu4p5_SQ_OS_dR_Max1p2    & qscale > 20)', 'norm hist e') ; print 'Purity of L1_DoubleMu4p5_SQ_OS_dR_Max1p2     %.2f +/- %.2f' %(100.*h1.GetMean(), 100.*h1.GetMeanError()) ; ROOT.gPad.SaveAs('purity_l1_doublemu4p5_sq_os_dr_max1p2.pdf'   )   




h2 = ROOT.TH1F('dr12_', 'dr12_', 25, 0, 6)
t1.Draw('dr_12 >> dr12_', '(qscale_weight) * (nbmesons==2 & qscale>20)'                                            , 'norm hist e') ; ROOT.gPad.SaveAs('dr_12_coarsebinning.pdf')
t1.Draw('dr_12 >> dr12_', '(qscale_weight) * (nbmesons==2 & qscale>20 & matched_L1_SingleMu_25_Q12)'               , 'norm hist e') ; ROOT.gPad.SaveAs('dr_12_coarsebinning_l1_singlemu_25_q12.pdf')
t1.Draw('dr_12 >> dr12_', '(qscale_weight) * (nbmesons==2 & qscale>20 & matched_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4)', 'norm hist e') ; ROOT.gPad.SaveAs('dr_12_coarsebinning_l1_doublemu0er1p5_sq_os_dr_max1p4.pdf')




