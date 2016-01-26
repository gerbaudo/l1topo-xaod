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
    "L1_MU4","L1_MU6","L1_MU10","L1_MU11","L1_MU15","L1_MU20",
    "L1_2MU4","L1_2MU6","L1_2MU10","L1_MU10_2MU4","L1_MU10_2MU6","L1_LFV-MU",
    "L1_LFV-MU-topo",
    "L1_3MU6","L1_3MU4","L1_MU6_3MU4",
    "HLT_mu4",  "HLT_mu6", "HLT_mu10", "HLT_mu18","HLT_mu24",
    "HLT_2mu14",
    "HLT_2mu10",
    "HLT_mu18_mu8noL1", "HLT_mu18_2mu4noL1",
    
    "HLT_3mu4",
    "HLT_mu4_2mu6",
    "HLT_3mu6",
    "HLT_3mu6_msonly",
    "HLT_3mu4_bTau", "HLT_3mu6_bTau", 

    "HLT_mu20_msonly_mu6noL1_msonly_nscan05",
    "HLT_mu20_mu6noL1_nscan03",

    #"HLT_mu11_2mu4noL1_nscan03",
    "HLT_mu11_2mu4noL1_nscan03_L1MU11_2MU6",
    'HLT_mu11_llns_2mu4noL1_nscan03_L1MU11_2MU6',
    "HLT_mu11_L1MU10_2mu4noL1_nscan03_L1MU10_2MU6",

    "HLT_noalg_eb_L1PhysicsLow_noPS", #L1_2MU6 , L1_3MU4, L1_MU20,
    "HLT_noalg_eb_L1PhysicsHigh_noPS", # 3MU6
    "HLT_eb_low_L1RD2_FILLED", # L1_MU6' L1_MU4_J12'
    "HLT_eb_high_L1RD2_FILLED"  # L1_2MU4', L1_MU15 L1_MU6_2MU4', 'L1_MU6_J20'  'L1_EM7_MU10', 'L1_EM15_MU4',
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
for entry in xrange( t.GetEntries() ):
    t.GetEntry( entry )
    HLT.setEvent(t)
    emulated = HLT.emulateDecision("L1_LFV-MU")
    
    eventNumber[0] = t.EventInfo.eventNumber()
    runNumber[0] = t.EventInfo.runNumber()
    for trig in trigList :
        passTrig[trig][0] = 0
        if  ROOT.trigDecTool.isPassed( trig ) : passTrig[trig][0] = 1 
    passTrig["L1_LFV-MU-topo"][0] = emulated
        
    l1muons.clear()
    for imu in xrange(len(t.LVL1MuonRoIs)) :
        mu = t.LVL1MuonRoIs[imu]
	#if mu.isVetoed() : continue
        l1 = ROOT.L1MuonRoI()
        l1.p4.SetPtEtaPhiM( mu.thrValue(), mu.eta(), mu.phi(), 105.65)
        l1.roiWord= mu.roiWord()
        l1.thrValue= mu.thrValue()
        l1.isVetoed = mu.isVetoed()
        l1muons.push_back(l1)
    l1muonsn[0] = len(l1muons)

    topomuons.clear()
    for imu in xrange(len(HLT.topoRoIs)) :
        mu = HLT.topoRoIs[imu]
        pt1 = mu.pT()
        if   pt1 == 1 : pt1 = 4000
        elif pt1 == 2 : pt1 = 6000
        else          : pt1 = 10000
        topo = ROOT.TLorentzVector()
        topo.SetPtEtaPhiM( pt1, mu.eta(), mu.phi(), 105.65)
        topomuons.push_back(topo)
    topomuonsn[0] = len(topomuons)

    pts = []
    for imu in xrange(len(t.Muons)) :
        mu = t.Muons[imu] 
        mutype = 0
        if mu.muonType() == ROOT.xAOD.Muon.Combined :
            mutype = 0
        elif mu.muonType() == ROOT.xAOD.Muon.SegmentTagged :
            mutype = 1
        elif mu.muonType() == ROOT.xAOD.Muon.MuonStandAlone :
            mutype = 2
        else :
            continue
        mupt = mu.pt()
        if mupt<2.5 : continue
        reco = ROOT.TLorentzVector()
        reco.SetPtEtaPhiM( mupt, mu.eta(), mu.phi(), 105.65)
        #print "r ", reco.Pt()
        pts.append([reco,mutype])

    #print pts
    pts.sort( key=lambda x: x[0].Pt(), reverse=True)
    recomuons.clear()
    for x in pts :
        recomuons.push_back(x[0])
        recotype.push_back(x[1])
        #print "x ", x.Pt()

    recomuonsn[0] = len(recomuons)

    tOut.Fill()
    
print "Write out trig ntuple ", fOut.GetName() , " with ", tOut.GetEntries(), " events"
fOut.Write()
fOut.Close()
exit(0)
