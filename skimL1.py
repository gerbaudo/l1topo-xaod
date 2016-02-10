#!/usr/bin/env python
import sys
from math import pi,sqrt,cos
from ctypes import *
import array

import ROOT
ROOT.gROOT.Macro( '$ROOTCOREDIR/scripts/load_packages.C' )
#ROOT.gROOT.ProcessLine(".L L1MuonRoI_cxx.so")
ROOT.gROOT.LoadMacro( "L1MuonRoI.cxx+" )

# Initialize the xAOD infrastructure:
if(not ROOT.xAOD.Init().isSuccess()): print "Failed xAOD.Init()"

inputFiles = sys.argv[1].split(',')
print "inputFiles = ", inputFiles
fileName=inputFiles[0]
if len(fileName) == 0 :
    print "Please provide input"
    exit()

process_all_events = False
filter_events = False
verbose = False

######### Read input tree

treeName = "CollectionTree" # default when making transient tree anyway

ch = ROOT.TChain(treeName)
for i in range(1,len(sys.argv)) :
    ch.Add(sys.argv[i])
    print sys.argv[i]

t = ROOT.xAOD.MakeTransientTree( ch ) #f, treeName)

######### Initialize TrigDecision

ROOT.gROOT.ProcessLine(
"  TrigConf::xAODConfigTool configTool( \"xAODConfigTool\" );\
  ToolHandle<TrigConf::ITrigConfigTool> handle(&configTool);\
  handle->initialize().ignore();\
  Trig::TrigDecisionTool  trigDecTool(\"TrigDecTool\");\
  trigDecTool.setProperty(\"ConfigTool\",handle).ignore();\
  trigDecTool.setProperty(\"TrigDecisionKey\",\"xTrigDecision\").ignore();\
  trigDecTool.setProperty(\"OutputLevel\", MSG::INFO).ignore();\
  trigDecTool.initialize().ignore();" )




trigList = [
    "L1_MU4","L1_MU6",
    # "L1_MU10","L1_MU11","L1_MU15","L1_MU20",
    # "L1_2MU4","L1_2MU6","L1_2MU10","L1_MU10_2MU4","L1_MU10_2MU6","L1_LFV-MU",
    "L1_LFV-MU-topo",
    # "L1_3MU6","L1_3MU4","L1_MU6_3MU4",
    # "HLT_mu4",  "HLT_mu6", "HLT_mu10", "HLT_mu18","HLT_mu24",
    # "HLT_2mu14",
    # "HLT_2mu10",
    # "HLT_mu18_mu8noL1", "HLT_mu18_2mu4noL1",

    # "HLT_3mu4",
    # "HLT_mu4_2mu6",
    # "HLT_3mu6",
    # "HLT_3mu6_msonly",
    # "HLT_3mu4_bTau", "HLT_3mu6_bTau",

    # "HLT_mu20_msonly_mu6noL1_msonly_nscan05",
    # "HLT_mu20_mu6noL1_nscan03",

    # #"HLT_mu11_2mu4noL1_nscan03",
    # "HLT_mu11_2mu4noL1_nscan03_L1MU11_2MU6",
    # 'HLT_mu11_llns_2mu4noL1_nscan03_L1MU11_2MU6',
    # "HLT_mu11_L1MU10_2mu4noL1_nscan03_L1MU10_2MU6",

    # "HLT_noalg_eb_L1PhysicsLow_noPS", #L1_2MU6 , L1_3MU4, L1_MU20,
    # "HLT_noalg_eb_L1PhysicsHigh_noPS", # 3MU6
    # "HLT_eb_low_L1RD2_FILLED", # L1_MU6' L1_MU4_J12'
    # "HLT_eb_high_L1RD2_FILLED"  # L1_2MU4', L1_MU15 L1_MU6_2MU4', 'L1_MU6_J20'  'L1_EM7_MU10', 'L1_EM15_MU4',
]

from OTriggerHLT import OTriggerHLT
HLT = OTriggerHLT()


counters = {}
for trig in trigList :
    counters[trig] = 0

passTrig = {}

fOut = ROOT.TFile("tmptrig.root","RECREATE")

tOut = ROOT.TTree("trig","truth")
passTrig = {}
for trig in trigList :
    passTrig[trig] = array.array("i",(0 for i in range(0,2)))
    tOut.Branch(str(trig), passTrig[trig],str(trig)+"/I")

eventNumber = array.array("i",(0 for i in range(0,2)))
tOut.Branch("eventNumber", eventNumber,"eventNumber/I")
runNumber = array.array("i",(0 for i in range(0,2)))
tOut.Branch("runNumber", runNumber,"runNumber/I")

#l1muons = ROOT.vector('TLorentzVector')(10)
l1muons = ROOT.vector('L1MuonRoI')(10)
tOut.Branch('l1muons',  l1muons  )
l1muonsn= array.array("i",(0 for i in range(0,1)))
tOut.Branch('l1muons_n',  l1muonsn, "l1muons_n/I"  )

topomuons = ROOT.vector('TLorentzVector')(10)
tOut.Branch('topomuons',  topomuons  )
topomuonsn= array.array("i",(0 for i in range(0,1)))
tOut.Branch('topomuons_n',  topomuonsn, "topomuons_n/I"  )

recomuons = ROOT.vector('TLorentzVector')(10)
tOut.Branch('recomuons',  recomuons  )
recomuonsn= array.array("i",(0 for i in range(0,1)))
tOut.Branch('recomuons_n',  recomuonsn, "recomuons_n/I"  )
recotype = ROOT.vector('int')(10)
tOut.Branch('recotype',  recotype  )


################################################################

print( "Number of input events: %s" % t.GetEntries() )
numEntriesToProcess = max([10, t.GetEntries()]) if process_all_events else min([10, t.GetEntries()])
print("About to process %d entries" % numEntriesToProcess)
mev2gev=1.0E-3

def get_l1met_p4(l1met, verbose=False):
    mev2gev = 1.0E-3
    metx = l1met.exMiss()
    mety = l1met.eyMiss()
    met = l1met.energyT()
    metp4 = ROOT.TLorentzVector(0.0,0.0,0.0,0.0)
    if verbose:
        print "met from input: ",l1met.energyT()
        print "met computed: ",sqrt(metx*metx + mety*mety)
    metp4.SetPxPyPzE(metx*mev2gev, mety*mev2gev, 0.0, met*mev2gev)
    return metp4

def print_l1met(l1met, l1metp4):
    print("l1met: ex\tey\tet")
    print("       %d\t%d\t%4.f" % (int(mev2gev*l1met.exMiss()),
                                   int(mev2gev*l1met.eyMiss()),
                                   l1metp4.Et()))

def print_l1jets(l1jets):
    print("l1jet[%d]: et8x8\tet4x4\teta\tphi"%len(l1jets))
    for iJet, l1jet in enumerate(l1jets):
        # format from TopoInputEvent::dump()
        # Et1 Et2 eta phi etaDouble phiDouble
        # print("l1jet[%d]: et8x8 %.2f et4x4 %.2f eta %.2f phi %.2f"%
        #       (iJet, l1jet.et8x8(), l1jet.et4x4(), l1jet.eta(), l1jet.phi()))
        print("[%03d]: %d\t%d\t%d\t%d"%(iJet,
                                        int(mev2gev*l1jet.et8x8()), int(mev2gev*l1jet.et4x4()),
                                        int(10*l1jet.eta()), int(10*l1jet.phi())))
def print_l1emtaus(l1emtaus):
    print("jets[%d]: et emIso eta phi"%len(l1emtaus))
    for iEmtau, l1emtau in enumerate(l1emtaus): # todo check duplication
        # format from TopoInputEvent::dump()
        # Et isolation eta phi etaDouble phiDouble
        # no etaDouble, phiDouble
        # http://acode-browser.usatlas.bnl.gov/lxr/source/atlas/Event/xAOD/xAODTrigger/xAODTrigger/versions/EmTauRoI_v2.h
        # print("l1emtau[%d]: et %.2f emIso %.2f eta %.2f phi %.2f" %
        #       (iEmtau, l1emtau.eT(), l1emtau.emIsol(), l1emtau.eta(), l1emtau.phi()))
        print("[%03d]: %d\t%d\t%d\t%d" %
              (iEmtau,
               int(mev2gev*l1emtau.eT()), int(mev2gev*l1emtau.emIsol()),
               int(10*l1emtau.eta()), int(10*l1emtau.phi())))
        # Q for xAOD devs & Joerg:
        # - there is a factor of 2 between the RAW values and the xAOD ones (both
        # - the type of l1emtau.isol() is 'str' why?
        # - the RAW iso() seems to be 2*xAOD::emIsol()
        # Q for myself:
        # - try to figure out the int/float rounding issues (due to python? to xAOD conversion?)

def skip_run_event(r, e):
    return (r!=287924
            or
            e not in [178424911,178424911,178496695,178662528,178432898,178525775])

def l1muonroi2l1muon(l1muonroi):
    "make an l1mu with p4 to be stored as output"
    mu = l1muonroi
    #if mu.isVetoed() : continue # add this check in the list comprehension?
    l1 = ROOT.L1MuonRoI()
    l1.p4.SetPtEtaPhiM( mu.thrValue(), mu.eta(), mu.phi(), 105.65)
    l1.roiWord= mu.roiWord()
    l1.thrValue= mu.thrValue()
    l1.isVetoed = mu.isVetoed()
    return l1

def toporoi2mu(toporoi):
    """from a topo(mu?) roi, make a TLV to be stored as output

    Ask Olya:
    - is this a muon toporoi, or a generic toporoi?
    - why isn't the pt decoded in L1MuonRoI.h?
    """
    mu = toporoi
    pt1 = mu.pT()
    if   pt1 == 1 : pt1 = 4000
    elif pt1 == 2 : pt1 = 6000
    else          : pt1 = 10000
    topo = ROOT.TLorentzVector()
    topo.SetPtEtaPhiM( pt1, mu.eta(), mu.phi(), 105.65)
    return topo

def filter_muons(muons):
    "filter the reco muons from ntuples based on their pt+type, and add the algo_type attribute (to be stored)"
    def add_algo(mu):
        comb = ROOT.xAOD.Muon.Combined
        segm = ROOT.xAOD.Muon.SegmentTagged
        stan = ROOT.xAOD.Muon.MuonStandAlone
        mu_type = mu.muonType()
        mu.algo_type = (0 if mu_type==comb else
                        1 if mu_type==segm else
                        2 if mu_type==stan else
                        None)
        return mu
    muons = [add_algo(m) for m in muons]
    return [m for m in muons if m.pt()>=2.5 and m.algo_type is not None]

def muon2recomuon(mu):
    "convert a reco muon to a TLV for the output ntuple"
    mu_mass = 105.65
    recomuon = ROOT.TLorentzVector()
    recomuon.SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), mu_mass)
    return recomuon

for entry in xrange(numEntriesToProcess):
    t.GetEntry( entry )
    HLT.setEvent(t)
    emulated = HLT.emulateDecision("L1_LFV-MU")

    eventNumber[0] = t.EventInfo.eventNumber()
    runNumber[0] = t.EventInfo.runNumber()
    l1met = t.LVL1EnergySumRoI
    l1metp4 = get_l1met_p4(l1met, verbose)
    l1jets = t.LVL1JetRoIs
    eventN = t.EventInfo.eventNumber()
    if filter_events and skip_run_event(runNumber[0], eventN):
        continue
    print("%d / %d" % (runNumber[0], eventNumber[0]))
    # met->Ex() << "  " << met->Ey() << "  " << met->Et()
    # print("l1met: ex %.2f ey %2.f et %2.f" % (l1met.exMiss(), l1met.eyMiss(), l1metp4.Et()))
    if verbose:
        print_l1met(l1met, l1metp4)
        print_l1jets(l1jets)

    l1emtaus= t.LVL1EmTauRoIs
    if verbose:
        print_l1emtaus(l1emtaus)

    for trig in trigList :
        passTrig[trig][0] = 0
        if ROOT.trigDecTool.isPassed( trig ):
            passTrig[trig][0] = 1
    passTrig["L1_LFV-MU-topo"][0] = emulated

    l1muons.clear()
    for l1mu in t.LVL1MuonRoIs:
        l1muons.push_back(l1muonroi2l1muon(l1mu))
    l1muonsn[0] = len(l1muons)

    topomuons.clear()
    for toporoi in HLT.topoRoIs:
        topomuons.push_back(toporoi2mu(toporoi))
    topomuonsn[0] = len(topomuons)

    recomuons.clear()
    recotype.clear()
    ntup_recomuons = sorted(filter_muons(t.Muons), key=lambda m : m.pt(), reverse=True)
    for mu in ntup_recomuons:
        recomuons.push_back(muon2recomuon(mu))
        recotype.push_back(mu.algo_type)
    recomuonsn[0] = len(recomuons)

    tOut.Fill()

print "Write out trig ntuple ", fOut.GetName() , " with ", tOut.GetEntries(), " events"
fOut.Write()
fOut.Close()
exit(0)
