import ROOT
from itertools import combinations
from math import cos, cosh, sqrt
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi, bestMatch

def single_muon(muons, pt, eta=2.5, qual=8, matches=[], cone=0.3, onlymuons=False, atvtx=True):
    
    cone2 = cone * cone
    
    results = []
    
    for mu in muons:
        matched = False
        dRmin = 999.
        if mu.hwQual()   < qual: continue
        if mu.pt()       < pt  : continue
        if abs(mu.eta()) > eta : continue
        
        mu_for_matching = mu.clone()
        if atvtx:
            mu_for_matching.eta = mu.etaAtVtx
            mu_for_matching.phi = mu.phiAtVtx
        
        match_index = -99
                
        for ii, imatch in enumerate(matches):
            imatch.matched = False
            matched = False
            if onlymuons:
                best_match, dRmin = bestMatch(mu_for_matching, imatch.finalmuondaughters)            
            else:
                best_match, dRmin = bestMatch(mu_for_matching, imatch.finalchargeddaughters)
            matched = dRmin<cone2
            if matched: 
                match_index = ii
                imatch.matched = True
                break
        
        results.append((1, int(matched), match_index))
    
    if len(results)>0:
        results.sort(key = lambda x : (x[0], x[1]), reverse=True)
        return int(results[0][0]), int(results[0][1]), int(results[0][2])

    return 0, 0, -1
    

def di_muon(muons, pt1, pt2, eta1=2.5, eta2=2.5, 
            qual1=8, qual2=8, minMass=-1., 
            maxMass=1.E10, minDr=-1., maxDr=1.E10, sign=-1, 
            matches=[], cone=0.3, onlymuons=False, atvtx=True):
    
    cone2 = cone*cone
    
    if len(muons)<2:
        return 0, 0, -1
    
    results = []
    
    muons.sort(key = lambda x : x.pt(), reverse = True)
    
    for mu1, mu2 in combinations(muons, 2):
        matched  = False
        if mu1.hwQual() < qual1: continue
        if mu2.hwQual() < qual2: continue

        if mu1.pt() < pt1: continue
        if mu2.pt() < pt2: continue
    
        if abs(mu1.eta()) > eta1: continue
        if abs(mu2.eta()) > eta2: continue
        
        if abs(mu1.charge() + mu2.charge())!=sign and sign>=0: continue

        mu1_p4_atVtx = ROOT.TLorentzVector()
        mu1_p4_atVtx.SetPtEtaPhiM(mu1.pt(), mu1.etaAtVtx(), mu1.phiAtVtx(), 0.105658)
    
        mu2_p4_atVtx = ROOT.TLorentzVector()
        mu2_p4_atVtx.SetPtEtaPhiM(mu2.pt(), mu2.etaAtVtx(), mu2.phiAtVtx(), 0.105658)
        
        mass        = (mu1_p4_atVtx + mu2_p4_atVtx).M()
        mass_approx = sqrt(2. * mu1.pt() * mu2.pt() * (cosh(mu1.etaAtVtx() - mu2.etaAtVtx()) - cos(mu1.phiAtVtx() - mu2.phiAtVtx())))
        dR          = mu1_p4_atVtx.DeltaR(mu2_p4_atVtx)
        
        if mass_approx < minMass: continue 
        if mass_approx > maxMass: continue 
        
        if dR < minDr: continue
        if dR > maxDr: continue

        mu1_for_matching = mu1.clone()
        mu2_for_matching = mu2.clone()
        if atvtx:
#             print 'pre  matching mu1: eta %.3f \t phi %.3f phi' %(mu1_for_matching.eta(), mu1_for_matching.phi())
#             print '              mu2: eta %.3f \t phi %.3f phi' %(mu2_for_matching.eta(), mu2_for_matching.phi())
            mu1_for_matching.eta = mu1.etaAtVtx
            mu1_for_matching.phi = mu1.phiAtVtx
            mu2_for_matching.eta = mu2.etaAtVtx
            mu2_for_matching.phi = mu2.phiAtVtx
#             print 'post matching mu1: eta %.3f \t phi %.3f phi' %(mu1_for_matching.eta(), mu1_for_matching.phi())
#             print '              mu2: eta %.3f \t phi %.3f phi' %(mu2_for_matching.eta(), mu2_for_matching.phi())
                
        match_index = -99

        for ii, imatch in enumerate(matches):
            imatch.matched = True
            goodmatches = []
            matched = False
            dRmin1 = 999.
            dRmin2 = 999.
            
            if onlymuons:
                best_match1, dRmin1 = bestMatch(mu1_for_matching, imatch.finalmuondaughters)  
                if dRmin1 < cone2:
                    goodmatches.append(best_match1)       
                best_match2, dRmin2 = bestMatch(mu2_for_matching, imatch.finalmuondaughters)  
                if dRmin2 < cone2 and best_match2 != best_match1:
                    goodmatches.append(best_match2)       
            else:
                best_match1, dRmin1 = bestMatch(mu1_for_matching, imatch.finalchargeddaughters)  
                if dRmin1 < cone2:
                    goodmatches.append(best_match1)       
                best_match2, dRmin2 = bestMatch(mu2_for_matching, imatch.finalchargeddaughters)  
                if dRmin2 < cone2 and best_match2 != best_match1:
                    goodmatches.append(best_match2)       

            matched = len(goodmatches)>1
            if matched: 
                match_index = ii
                imatch.matched = True
                break
        
#         if matched:
#             for i in goodmatches: print i.pdgId(), i.pt(), i.eta(), i.phi()
#             print ''
#             for j in muons: print j.pt(), j.eta(), j.phi()
#             import pdb ; pdb.set_trace()

        results.append((1, int(matched), match_index))
        
    if len(results)>0:
        results.sort(key = lambda x : (x[0], x[1]), reverse=True)
#         if len(results)>1 and any(ii==(True, False) for ii in results) and any(ii==(True, True) for ii in results): 
#             import pdb ; pdb.set_trace()
        return int(results[0][0]), int(results[0][1]), int(results[0][2])
    
    return 0, 0, -1


def tri_muon(muons, pt1, pt2, pt3, eta1=2.5, eta2=2.5, eta3=2.5, 
             qual1=8, qual2=8, qual3=8, minMass=-1., 
             maxMass=1.E10, minDr=-1., maxDr=1.E10, sign=-1,
             matches=[], cone=0.3, onlymuons=False):
    
    cone2 = cone * cone
    
    passed = False

    results = []

    if len(muons) < 3:
        return 0, 0, -1

    muons.sort(key = lambda x : x.pt(), reverse = True)    
    triplets = combinations(muons, 3)
    
    for mu1, mu2, mu3 in triplets:
        matched  = False
        
        if mu1.hwQual() < qual1: continue
        if mu2.hwQual() < qual2: continue
        if mu3.hwQual() < qual3: continue

        if mu1.pt() < pt1: continue
        if mu2.pt() < pt2: continue
        if mu3.pt() < pt3: continue
    
        if abs(mu1.eta()) > eta1: continue
        if abs(mu2.eta()) > eta2: continue
        if abs(mu3.eta()) > eta2: continue
    
        passed = di_muon([mu1, mu2, mu3], pt1, pt2, eta1, eta2, 
                         qual1, qual2, minMass, 
                         maxMass, minDr, maxDr, sign)


        if not passed: continue

        mu1_for_matching = mu1.clone()
        mu2_for_matching = mu2.clone()
        mu3_for_matching = mu3.clone()
        if atvtx:
            mu1_for_matching.eta = mu1.etaAtVtx
            mu1_for_matching.phi = mu1.phiAtVtx
            mu2_for_matching.eta = mu2.etaAtVtx
            mu2_for_matching.phi = mu2.phiAtVtx
            mu3_for_matching.eta = mu3.etaAtVtx
            mu3_for_matching.phi = mu3.phiAtVtx

        match_index = -99

        for imatch in matches:
            imatch.matched = False
            goodmatches = []
            matched = False
            dRmin1 = 999.
            dRmin2 = 999.
            dRmin3 = 999.
            if onlymuons:
                best_match1, dRmin1 = bestMatch(mu1_for_matching, imatch.finalmuondaughters)  
                if dRmin1 < cone2:
                    goodmatches.append(best_match1)       
                best_match2, dRmin2 = bestMatch(mu2_for_matching, imatch.finalmuondaughters)  
                if dRmin2 < cone2 and best_match2 != best_match1:
                    goodmatches.append(best_match2)       
                best_match3, dRmin3 = bestMatch(mu3_for_matching, imatch.finalmuondaughters)  
                if dRmin3 < cone2 and best_match3 != best_match2 and best_match3 != best_match1:
                    goodmatches.append(best_match3)       
            else:
                best_match1, dRmin1 = bestMatch(mu1_for_matching, imatch.finalchargeddaughters)  
                if dRmin1 < cone2:
                    goodmatches.append(best_match1)       
                best_match2, dRmin2 = bestMatch(mu2_for_matching, imatch.finalchargeddaughters)  
                if dRmin2 < cone2 and best_match2 != best_match1:
                    goodmatches.append(best_match2)       
                best_match3, dRmin3 = bestMatch(mu3_for_matching, imatch.finalchargeddaughters)  
                if dRmin3 < cone2 and best_match3 != best_match2 and best_match3 != best_match1:
                    goodmatches.append(best_match3)       

            matched = len(goodmatches)>2
            if matched: 
                match_index = ii
                imatch.matched = True
                break
        
        results.append((1, int(matched), match_index))
        
    if len(results)>0:
        results.sort(key = lambda x : (x[0], x[1]), reverse=True)
        return int(results[0][0]), int(results[0][1]), int(results[0][2])
    
    return 0, 0
