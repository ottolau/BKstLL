import ROOT

ROOT.TH1.SetDefaultSumw2()

ROOT.gStyle.SetOptStat(0)

ROOT.gROOT.SetBatch(True)

f1 = ROOT.TFile.Open('gluon/QCD_Pt_15to30/L1PurityTreeProducer/tree.root', 'read')
f1.cd()
t1 = f1.Get('tree')

f2 = ROOT.TFile.Open('gluon/QCD_Pt_80to120/L1PurityTreeProducer/tree.root', 'read')
f2.cd()
t2 = f2.Get('tree')

f3 = ROOT.TFile.Open('gluon/QCD_Pt_3200toInf/L1PurityTreeProducer/tree.root', 'read')
f3.cd()
t3 = f3.Get('tree')

h1_1 = ROOT.TH1F('dr_inclusive_lowpt', '15 < #hat{pt} < 30 GeV' , 50, 0, 7); h1_1.SetLineColor(ROOT.kBlack)
h1_2 = ROOT.TH1F('dr_gs_lowpt'       , '15 < #hat{pt} < 30 GeV' , 50, 0, 7); h1_2.SetLineColor(ROOT.kRed  )
h1_3 = ROOT.TH1F('dr_other_lowpt'    , '15 < #hat{pt} < 30 GeV' , 50, 0, 7); h1_3.SetLineColor(ROOT.kBlue )

h2_1 = ROOT.TH1F('dr_inclusive_medpt', '80 < #hat{pt} < 120 GeV', 50, 0, 7); h2_1.SetLineColor(ROOT.kBlack)
h2_2 = ROOT.TH1F('dr_gs_medpt'       , '80 < #hat{pt} < 120 GeV', 50, 0, 7); h2_2.SetLineColor(ROOT.kRed  )
h2_3 = ROOT.TH1F('dr_other_medpt'    , '80 < #hat{pt} < 120 GeV', 50, 0, 7); h2_3.SetLineColor(ROOT.kBlue )

h3_1 = ROOT.TH1F('dr_inclusive_hipt' , '#hat{pt} > 3200 GeV'    , 50, 0, 7); h3_1.SetLineColor(ROOT.kBlack)
h3_2 = ROOT.TH1F('dr_gs_hipt'        , '#hat{pt} > 3200 GeV'    , 50, 0, 7); h3_2.SetLineColor(ROOT.kRed  )
h3_3 = ROOT.TH1F('dr_other_hipt'     , '#hat{pt} > 3200 GeV'    , 50, 0, 7); h3_3.SetLineColor(ROOT.kBlue )

for hh in [h1_1, h1_2, h1_3, h2_1, h2_2, h2_3, h3_1, h3_2, h3_3]:
    hh.GetXaxis().SetTitle('dR(B_{1},B_{2})')
    hh.GetYaxis().SetTitle('a.u.')


t1.Draw('dr_12 >> dr_inclusive_lowpt', 'nbmesons==2'         )
t1.Draw('dr_12 >> dr_gs_lowpt'       , 'nbmesons==2 & is_gs' )
t1.Draw('dr_12 >> dr_other_lowpt'    , 'nbmesons==2 & !is_gs')

t2.Draw('dr_12 >> dr_inclusive_medpt', 'nbmesons==2'         )
t2.Draw('dr_12 >> dr_gs_medpt'       , 'nbmesons==2 & is_gs' )
t2.Draw('dr_12 >> dr_other_medpt'    , 'nbmesons==2 & !is_gs')

t3.Draw('dr_12 >> dr_inclusive_hipt' , 'nbmesons==2'         )
t3.Draw('dr_12 >> dr_gs_hipt'        , 'nbmesons==2 & is_gs' )
t3.Draw('dr_12 >> dr_other_hipt'     , 'nbmesons==2 & !is_gs')


leg = ROOT.TLegend(.6,.7,.88,.88)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
# leg.SetTextFont(42)
# leg.SetTextSize(0.035)
leg.AddEntry(h1_1, 'inclusive production', 'L')
leg.AddEntry(h1_2, 'gluon splitting'     , 'L')
leg.AddEntry(h1_3, 'other modes'         , 'L')


h1_1.Draw('hist')
h1_2.Draw('hist same')
h1_3.Draw('hist same')
leg.Draw('same')

ROOT.gPad.SaveAs('dr_bb_lowpt.pdf')

h2_1.Draw('hist')
h2_2.Draw('hist same')
h2_3.Draw('hist same')
leg.Draw('same')

ROOT.gPad.SaveAs('dr_bb_medpt.pdf')

h3_1.Draw('hist')
h3_2.Draw('hist same')
h3_3.Draw('hist same')
leg.Draw('same')

ROOT.gPad.SaveAs('dr_bb_hipt.pdf')






h_gpt_1  = ROOT.TH1F('gpt_lowpt' , '' , 100,  0, 3500); h_gpt_1 .SetLineColor(ROOT.kBlack)
h_gpt_2  = ROOT.TH1F('gpt_medpt' , '' , 100,  0, 3500); h_gpt_2 .SetLineColor(ROOT.kRed  )
h_gpt_3  = ROOT.TH1F('gpt_hipt'  , '' , 100,  0, 3500); h_gpt_3 .SetLineColor(ROOT.kBlue )

h_geta_1 = ROOT.TH1F('geta_lowpt', '' ,  25, -6,    6); h_geta_1.SetLineColor(ROOT.kBlack)
h_geta_2 = ROOT.TH1F('geta_medpt', '' ,  25, -6,    6); h_geta_2.SetLineColor(ROOT.kRed  )
h_geta_3 = ROOT.TH1F('geta_hipt' , '' ,  25, -6,    6); h_geta_3.SetLineColor(ROOT.kBlue )

for hh in [h_gpt_1, h_gpt_2, h_gpt_3]:
    hh.GetXaxis().SetTitle('gluon p_{T} [GeV]')
    hh.GetYaxis().SetTitle('dN/dp_{T}')

for hh in [h_geta_1, h_geta_2, h_geta_3]:
    hh.GetXaxis().SetTitle('gluon #eta')
    hh.GetYaxis().SetTitle('dN/d#eta')


t1.Draw('g_pt  >> gpt_lowpt' , 'nbmesons==2 & is_gs')
t1.Draw('g_eta >> geta_lowpt', 'nbmesons==2 & is_gs')

t2.Draw('g_pt  >> gpt_medpt' , 'nbmesons==2 & is_gs')
t2.Draw('g_eta >> geta_medpt', 'nbmesons==2 & is_gs')

t3.Draw('g_pt  >> gpt_hipt'  , 'nbmesons==2 & is_gs')
t3.Draw('g_eta >> geta_hipt' , 'nbmesons==2 & is_gs')

leg = ROOT.TLegend(.6,.7,.88,.88)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
# leg.SetTextFont(42)
# leg.SetTextSize(0.035)
leg.AddEntry(h_gpt_1, '15 < #hat{pt} < 30 GeV'  , 'L')
leg.AddEntry(h_gpt_2, '80 < #hat{pt} < 120 GeV' , 'L')
leg.AddEntry(h_gpt_3, '#hat{pt} > 3200 GeV (x5)', 'L')


for hh in [h_gpt_1, h_gpt_2, h_gpt_3, h_geta_1, h_geta_2, h_geta_3]:
    hh.Scale(1./hh.Integral())

h_gpt_1.Draw('hist')
h_gpt_2.Draw('hist same')
h_gpt_3.Scale(5.)
h_gpt_3.Draw('hist same')
leg.Draw('same')

# ROOT.gPad.SetLogx(True)
ROOT.gPad.Update()
ROOT.gPad.SaveAs('gluon_pt.pdf')


leg = ROOT.TLegend(.6,.7,.88,.88)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
# leg.SetTextFont(42)
# leg.SetTextSize(0.035)
leg.AddEntry(h_gpt_1, '15 < #hat{pt} < 30 GeV' , 'L')
leg.AddEntry(h_gpt_2, '80 < #hat{pt} < 120 GeV', 'L')
leg.AddEntry(h_gpt_3, '#hat{pt} > 3200 GeV'    , 'L')

h_geta_3.Draw('hist')
h_geta_1.Draw('hist same')
h_geta_2.Draw('hist same')
leg.Draw('same')

ROOT.gPad.SetLogx(False)
ROOT.gPad.Update()
ROOT.gPad.SaveAs('gluon_eta.pdf')




