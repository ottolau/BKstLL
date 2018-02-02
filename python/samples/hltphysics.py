import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator

json = '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_9_2_2_minimal_recipe/src/CMGTools/BKstLL/python/samples/ephemeral_hlt_physics_run305636_pu55.txt'

creator = ComponentCreator()

EphemeralHLTPhysics_2017F     = creator.makeDataComponent(
    "EphemeralHLTPhysics_2017F", 
    "/EphemeralHLTPhysics1/Run2017F-PromptReco-v1/MINIAOD", 
    "CMS", 
    ".*root", 
    json,
    useAAA=True
)

EphemeralHLTPhysics_2017F.files = [
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/027FB71B-7CBE-E711-ADB8-02163E01A6C2.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/02BE7273-75BE-E711-9D90-02163E019B9E.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/0451082D-87BE-E711-A36A-02163E013793.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/08385575-7ABE-E711-ACA2-02163E01A7A7.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/18508822-77BE-E711-A4BC-02163E019DF5.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/1E3D3A50-9ABE-E711-8B14-02163E01455D.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/308DD555-8DBE-E711-B91C-02163E0123AE.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/3C57078D-85BE-E711-8BB4-02163E011D4D.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/42C8086E-90BE-E711-B298-02163E019CF2.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/5012AE5C-95BE-E711-AF25-02163E013825.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/50F1E726-8ABE-E711-9E51-02163E01A25F.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/52CA4328-84BE-E711-BA79-02163E019DBA.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/62AE896B-75BE-E711-93DB-02163E01A5D1.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/6CA80765-81BE-E711-B63A-02163E019C78.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/6E213335-74BE-E711-A121-02163E0134AA.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/70E464D4-88BE-E711-9C9A-02163E01A240.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/78E58D60-7FBE-E711-9E1A-02163E01A4E6.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/7AB34339-94BE-E711-856B-02163E011EE5.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/8012CEC1-80BE-E711-92EE-02163E01A4C7.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/82A344A5-9DBE-E711-BC10-02163E019D9E.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/8E24DBAC-17BF-E711-BCC4-02163E01A65B.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/90D28640-92BE-E711-8E2A-02163E01A41F.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/92BDEDB4-B0BE-E711-A40A-02163E0144CA.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/9C920966-7DBE-E711-8BF6-02163E01A25F.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/A0E48D38-8FBE-E711-A8EF-02163E01426A.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/BE887035-A2BE-E711-8C09-02163E014310.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/C81D6925-7ABE-E711-A514-02163E0145DC.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/CE08A8E7-9BBE-E711-8DB0-02163E01A5C6.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/CEA3AF99-97BE-E711-80D7-02163E011ACC.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/ECD1ED74-82BE-E711-9F0E-02163E011A04.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics1/MINIAOD/PromptReco-v1/000/305/636/00000/EED2C652-82BE-E711-9D93-02163E012317.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/06025449-2BBE-E711-A29C-02163E019CC9.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/060EFB7F-37BE-E711-AB2E-02163E0146E8.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/0815CDE2-17BE-E711-B71A-02163E01A738.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/1089A03B-1FBE-E711-A76A-02163E019D66.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/16B19FA2-21BE-E711-89F4-02163E0143E7.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/22E9F280-0ABF-E711-B56B-02163E01A681.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/26B4138B-20BE-E711-9938-02163E01A7A7.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/3018D11F-2EBE-E711-A16B-02163E014537.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/3A45D932-22BE-E711-A5FE-02163E019DB0.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/3C0A39E1-24BE-E711-BD43-02163E01A282.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/3C83AD34-B5BF-E711-B76A-02163E0145ED.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/52F0E90C-33BE-E711-8AA5-02163E011B4C.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/601D22A8-15BE-E711-9A79-02163E014141.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/62A3979E-1CBE-E711-9643-02163E01427A.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/6CBC17DF-24BE-E711-AD58-02163E019CA2.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/7CDCCE79-16BE-E711-9BCC-02163E01A4BC.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/7CEBECFD-1EBE-E711-8D6D-02163E0140EB.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/8232982C-1BBE-E711-A422-02163E0137EE.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/8234AABF-19BE-E711-890F-02163E019DB2.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/86820664-26BE-E711-A5C4-02163E01A2B5.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/86A3A4A3-2CBE-E711-A444-02163E019DDA.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/92F6316D-26BE-E711-BFBB-02163E01A3B1.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/96F2088A-23BE-E711-A665-02163E0123AE.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/A44D7EB0-19BE-E711-A4E6-02163E013462.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/A4A640F9-0FBE-E711-B6D5-02163E019CF6.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/AC7BFC5B-28BE-E711-8657-02163E013825.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/D232D1B7-29BE-E711-B90D-02163E01A5EA.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/DCFB687D-3FBE-E711-BCB9-02163E01435C.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/E2CB5B6A-16BE-E711-A57C-02163E019C97.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics2/MINIAOD/PromptReco-v1/000/305/636/00000/E6B4C793-11BE-E711-A01F-02163E011ECF.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/02A12792-78BE-E711-8B7F-02163E014141.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/04B34850-85BE-E711-A0EC-02163E011E17.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/085F79EB-03BF-E711-9316-02163E01A3D7.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/0A846B01-71BE-E711-903D-02163E01A305.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/1633475F-72BE-E711-BDA7-02163E0141B1.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/1A452B78-7DBE-E711-B86C-02163E011C66.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/2A5AB970-6DBE-E711-9D56-02163E019C01.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/3805B751-85BE-E711-A490-02163E01A597.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/3E4C0554-6FBE-E711-BBF0-02163E01A22D.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/42198E62-75BE-E711-9FA0-02163E01A747.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/466136AC-73BE-E711-9913-02163E01A632.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/5267DD61-82BE-E711-B154-02163E019E26.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/5422A111-84BE-E711-868A-02163E0144DE.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/54CFF201-7CBE-E711-92B7-02163E019BA1.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/5A351172-72BE-E711-A551-02163E014280.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/6E414024-67BE-E711-8039-02163E011E00.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/7629584B-6FBE-E711-8729-02163E013940.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/7AD6A361-7ABE-E711-BCA1-02163E01A5CC.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/8855DEC2-6ABE-E711-83B0-02163E011A69.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/909CBC9A-65BE-E711-98EC-02163E0136CA.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/927667C3-88BE-E711-BD71-02163E01A4B6.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/AAA46A2D-8ABE-E711-AF13-02163E012557.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/AEE7CD02-74BE-E711-A6A2-02163E019B80.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/AEFB2C2D-77BE-E711-B484-02163E01A78F.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/CCA38D2D-87BE-E711-9AAF-02163E0133CC.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/D477D377-7FBE-E711-8F4D-02163E01A22D.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/D8A61822-ADBE-E711-AD49-02163E01A214.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/DCDE9B54-75BE-E711-854D-02163E0123C6.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/DE948324-8DBE-E711-96A3-02163E01A595.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/E40FC757-95BE-E711-ABB0-02163E01A2E5.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/E6ACF40E-6CBE-E711-982A-02163E013573.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/E6C99775-90BE-E711-961E-02163E0133F1.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/EE5C78AF-80BE-E711-95F4-02163E01A70B.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics3/MINIAOD/PromptReco-v1/000/305/636/00000/F80BA302-7CBE-E711-AE40-02163E01A6E1.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/04025C65-FABD-E711-BDD5-02163E019D20.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/1489E7B3-06BE-E711-94A9-02163E019CA9.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/1C20F15E-03BE-E711-A401-02163E019DD0.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/24E20FAF-29BE-E711-BD97-02163E01A59E.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/3E2F470B-20BE-E711-A5BF-02163E014529.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/4867B74E-C4BF-E711-9AB0-02163E019B80.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/5032BB4A-11BE-E711-BCFD-02163E019B63.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/52C5271E-0ABE-E711-AFFE-02163E01A552.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/546CDA6A-F8BD-E711-92FF-02163E01A5B8.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/5CDE1A04-18BE-E711-A7A5-02163E01A60D.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/60F0B7CC-19BE-E711-AA30-02163E0128FC.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/74557259-13BE-E711-A2F1-02163E01A36C.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/766EA7DE-FCBD-E711-A973-02163E014108.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/808E9DDA-FEBD-E711-B8BE-02163E01A685.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/94DC3E90-0EBE-E711-A6A7-02163E0142B0.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/98609150-05BE-E711-8F67-02163E0122EE.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/9E1D8150-F8BD-E711-8E30-02163E011BD0.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/ACE46889-02BE-E711-91A6-02163E019C83.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/BEF08645-FBBD-E711-A9EC-02163E012067.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/C0082D4C-F3BD-E711-964C-02163E013647.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/C8DAC6CD-16BE-E711-BF8C-02163E01A76A.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/CCBA066C-F5BD-E711-9DCF-02163E019BD8.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/D2052C2E-0DBE-E711-926D-02163E01411D.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/D6C166D5-FCBD-E711-A45A-02163E0133F1.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/D6D32962-F0BD-E711-AE5F-02163E01A2C4.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/DAF79BE1-CBBE-E711-A574-02163E013503.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/DECFCBE3-FFBD-E711-8806-02163E014467.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/FA6BD15F-F3BD-E711-8AC4-02163E01A4F3.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/FC615B82-1CBE-E711-BF3C-02163E01A76C.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/FE4E109F-15BE-E711-992B-02163E019CDF.root',
    'root://cms-xrd-global.cern.ch//store/data/Run2017F/EphemeralHLTPhysics4/MINIAOD/PromptReco-v1/000/305/636/00000/FE8A650E-10BE-E711-B788-02163E011C1F.root',
]

