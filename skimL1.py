#!/usr/bin/env python
import sys
from math import pi,sqrt,cos
import array

import ROOT

mev2gev = 1.0E-3

def main():
    # init: load packages, setup xAOD infrastructure
    ROOT.gROOT.SetBatch(1)
    ROOT.gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')
    ROOT.gROOT.LoadMacro('L1MuonRoI.cxx+')
    if(not ROOT.xAOD.Init().isSuccess()): print "Failed xAOD.Init()"

    # Read input tree
    inputFiles = sys.argv[1].split(',')
    print "inputFiles = ", inputFiles
    fileName=inputFiles[0]
    if len(fileName) == 0 :
        print "Please provide input"
        exit()

    process_all_events = True # False
    filter_events = False # True
    verbose = False # True

    treeName = "CollectionTree" # default when making transient tree anyway
    ch = ROOT.TChain(treeName)
    for input_file in sys.argv[1:]:
        ch.Add(input_file)
        print input_file
    t = ROOT.xAOD.MakeTransientTree( ch ) #f, treeName)

    initialize_trigger_decision_tool()
    trigList = get_muon_triggers_to_process()

    from OTriggerHLT import OTriggerHLT
    HLT = OTriggerHLT()
    counters = dict((trig, 0) for trig in trigList)
    passTrig = {}

    # create the output tree with its branches
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

    histos = book_histos()
    counters = book_counters()

    # Start processing events
    print( "Number of input events: %s" % t.GetEntries() )
    numEntriesToProcess = max([10, t.GetEntries()]) if process_all_events else min([10, t.GetEntries()])
    print("About to process %d entries" % numEntriesToProcess)

    for entry in xrange(numEntriesToProcess):
        t.GetEntry( entry )
        HLT.setEvent(t)
        emulated = HLT.emulateDecision("L1_LFV-MU")

        eventNumber[0] = t.EventInfo.eventNumber()
        runNumber[0] = t.EventInfo.runNumber()
        l1met = t.LVL1EnergySumRoI
        l1metp4 = get_l1met_p4(l1met)
        l1jets = t.LVL1JetRoIs
        eventN = t.EventInfo.eventNumber()
        if filter_events and skip_run_event(runNumber[0], eventN):
            continue
        if verbose: # True
            print(10*'-')
            print("%d / %d" % (runNumber[0], eventNumber[0]))
            # continue
            # print_l1met(l1met, l1metp4)
            # print_l1jets(l1jets)

        l1emtaus= t.LVL1EmTauRoIs
        if verbose:
            print_l1emtaus(l1emtaus)
        check_L1_LAR_EM      (l1emtaus=l1emtaus, tdt=ROOT.trigDecTool, histos=histos, counters=counters)
        check_L1_JPSI_1M5_EM7(l1emtaus=l1emtaus, tdt=ROOT.trigDecTool, histos=histos, counters=counters)

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

        check_L1_BTAG_MU4J15(l1mus=topomuons, l1jets=l1jets, tdt=ROOT.trigDecTool, histograms=histos, counters=counters)
        continue

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
    print(10*'-')
    print('Summary')
    print('\t'.join(['trig', 'pass-tdt', 'pass-emul']))
    for trig, counts in counters.iteritems():
        print("%10s %6d   %6d"%(trig, counts['pass'], counts['emul']))

def get_l1met_p4(l1met, verbose=False):
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
    print("l1emtau[%d]: et emIso eta phi"%len(l1emtaus))
    twice = 2.0 # Murrough says that these clusters are in 0.5MeV unit?
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

def get_muon_triggers_to_process():
    return [
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
def initialize_trigger_decision_tool():
    cmd = """TrigConf::xAODConfigTool configTool( "xAODConfigTool" );
ToolHandle<TrigConf::ITrigConfigTool> handle(&configTool);
handle->initialize().ignore();
Trig::TrigDecisionTool  trigDecTool("TrigDecTool");
trigDecTool.setProperty("ConfigTool",handle).ignore();
trigDecTool.setProperty("TrigDecisionKey","xTrigDecision").ignore();
trigDecTool.setProperty("OutputLevel", MSG::INFO).ignore();
trigDecTool.initialize().ignore();
"""
    ROOT.gROOT.ProcessLine(cmd)

def book_histos():
    histos = dict()
    # L1_LAR-EM
    meT, MeT, net = 0, 100, 50
    meta, Meta, neta = -4, +4, 40
    mphi, Mphi, nphi = -pi, +pi, 40
    histos['h_any_l1em_et'      ] = ROOT.TH1D('h_any_l1em_et'      ,';ET [GeV]', net, meT, MeT)
    histos['h_pass_l1em_et'     ] = ROOT.TH1D('h_pass_l1em_et'     ,';ET [GeV]', net, meT, MeT)
    histos['h_emul_l1em_et'     ] = ROOT.TH1D('h_emul_l1em_et'     ,';ET [GeV]', net, meT, MeT)
    histos['h_any_l1em_eta_phi' ] = ROOT.TH2D('h_any_l1em_eta_phi' ,';eta;phi', neta, meta, Meta, nphi, mphi, Mphi)
    histos['h_pass_l1em_eta_phi'] = ROOT.TH2D('h_pass_l1em_eta_phi',';eta;phi', neta, meta, Meta, nphi, mphi, Mphi)
    histos['h_emul_l1em_eta_phi'] = ROOT.TH2D('h_emul_l1em_eta_phi',';eta;phi', neta, meta, Meta, nphi, mphi, Mphi)
    # L1_JPSI-1M5-EM7
    nem, mem, Mem = 11, -0.5, 10.5 # N em clusters
    nm, mm, Mm = 40, 0.0, 40.0 # inv mass
    histos['h_any_l1jpsi_mult_m' ] = ROOT.TH2D('h_any_l1jpsi_mult_m'  ,';mult;m[GeV]', nem, mem, Mem, nm, mm, Mm)
    histos['h_pass_l1jpsi_mult_m'] = ROOT.TH2D('h_pass_l1jpsi_mult_m' ,';mult;m[GeV]', nem, mem, Mem, nm, mm, Mm)
    histos['h_emul_l1jpsi_mult_m'] = ROOT.TH2D('h_emul_l1jpsi_mult_m' ,';mult;m[GeV]', nem, mem, Mem, nm, mm, Mm)
    binning = (11, -0.5, 10.5, 4, 0.0, 4.0)
    histos['h_l1jpsi_emall_mult_match'] = ROOT.TH2D('h_l1jpsi_emall_mult_match', ';emall multiplicity; hw/emul', *binning)
    set_hw_emul_bin_labels(histos['h_l1jpsi_emall_mult_match'])

    # L1_BTAG_MU4J15
    npairs, mpairs, Mpairs = 11, -0.5, 10.5 # N jet-mu pairs
    ndr, mdr, Mdr = 60, 0.0, 6.0 # deltaR(j,mu)
    other_args  = (';mult;#Delta R (jet, #mu)', 11, -0.5, 10.5, 60, 0.0, 6.0)
    histos['h_any_l1btag_mult_dr' ] = ROOT.TH2D('h_any_l1btag_mult_dr'  , *other_args)
    histos['h_pass_l1btag_mult_dr'] = ROOT.TH2D('h_pass_l1btag_mult_dr' , *other_args)
    histos['h_emul_l1btag_mult_dr'] = ROOT.TH2D('h_emul_l1btag_mult_dr' , *other_args)
    binning = (11, -0.5, 10.5, 4, 0.0, 4.0)
    histos['h_l1btag_mu4ab_mult_match' ] = ROOT.TH2D('h_l1btag_mu4ab_mult_match' , ';mu4ab multiplicity; hw/emul', *binning)
    histos['h_l1btag_cj15ab_mult_match'] = ROOT.TH2D('h_l1btag_cj15ab_mult_match', ';cj15ab multiplicity; hw/emul', *binning)
    set_hw_emul_bin_labels(histos['h_l1btag_mu4ab_mult_match' ])
    set_hw_emul_bin_labels(histos['h_l1btag_cj15ab_mult_match'])
    return histos

def hw_emul_binc(pass_hw, pass_em):
    "bin center assuming that the binning is 4 bins in [0.0, 4.0]"
    return (0.5 if pass_hw and pass_em else
            1.5 if pass_hw else
            2.5 if pass_em else
            3.5)

def hw_emul_label(pass_hw, pass_em):
    return ('1/1' if pass_hw and pass_em else
            '1/0' if pass_hw else
            '0/1' if pass_em else
            '0/0')

def set_hw_emul_bin_labels(h):
    yax = h.GetYaxis()
    for pass_hw in [True, False]:
        for pass_em in [True, False]:
            yax.SetBinLabel(yax.FindBin(hw_emul_binc(pass_hw, pass_em)),
                            hw_emul_label(pass_hw, pass_em))

def book_counters():
    counters = dict((t, {'any':0,
                         'pass':0,
                         'emul':0})
                    for t in ['L1_LAR-EM', 'L1_JPSI-1M5-EM7', 'L1_BTAG_MU4J15'])
    return counters

def check_L1_LAR_EM(l1emtaus=[], tdt=None, histos={}, counters={}):
    """check the L1_LAR-EM trigger in two ways:
    - emulation vs. tdt -> counters
    - pass/fail one of the two requirements (et, position) -> histograms
    """
    h_any_et       = histos['h_any_l1em_et'      ]
    h_pass_et      = histos['h_pass_l1em_et'     ]
    h_emul_et      = histos['h_emul_l1em_et'     ]
    h_any_eta_phi  = histos['h_any_l1em_eta_phi' ]
    h_pass_eta_phi = histos['h_pass_l1em_eta_phi']
    h_emul_eta_phi = histos['h_emul_l1em_eta_phi']
    twice = 1.0 # todo checkme 2.0 # Murrough says that these clusters are in 0.5MeV unit?
    def in_et(l): # todo: implement in int?
        return twice*mev2gev*l.eT() > 10
    def in_rectangle(l): # todo: implement this in int space? (rather than float)
        return (0.4<l.eta()<1.9) and (1.8<l.phi()<2.2)
    passed   = tdt.isPassed('L1_LAR-EM')
    emulated = any(l1em for l1em in l1emtaus if in_et(l1em) and in_rectangle(l1em))
    counters['L1_LAR-EM']['any' ] += 1
    counters['L1_LAR-EM']['pass'] += (1 if passed else 0)
    counters['L1_LAR-EM']['emul'] += (1 if emulated else 0)
    for iL1em, l1em in enumerate(l1emtaus):
        et = twice*mev2gev*l1em.eT()
        eta = l1em.eta()
        phi = l1em.phi()
        pass_et = in_et(l1em)
        pass_rc = in_rectangle(l1em)
        if pass_rc: h_any_et.Fill(et)
        if pass_et: h_any_eta_phi.Fill(eta, phi)
        if passed:
            if pass_rc: h_pass_et.Fill(et)
            if pass_et: h_pass_eta_phi.Fill(eta, phi)
        if emulated:
            if pass_rc: h_emul_et.Fill(et)
            if pass_et: h_emul_eta_phi.Fill(eta, phi)
        # print("  [%03d]: %.2f\t%.2f\t%.2f" % (iL1em, et, eta, phi))

def check_L1_JPSI_1M5_EM7(l1emtaus=[], tdt=None, histos={}, counters={}):
    em7s1 = sorted([em for em in l1emtaus if mev2gev*em.eT()>7.0], # todo check 'twice*'
                   key=lambda em : em.eT(), reverse=True
                   )[:1]
    emall = sorted(l1emtaus, key=lambda em : em.eT(), reverse=True)
    # in the real algo, is the pair (e1,e1) allowed?
    def addp4(o):
        o.p4 = ROOT.TLorentzVector(0.0,0.0,0.0,0.0)
        o.p4.SetPtEtaPhiM(mev2gev*o.eT(), o.eta(), o.phi(), 0.0)
        return o
    masses = [(addp4(e1).p4+addp4(e2).p4).M() for e1 in em7s1 for e2 in emall]
    # masses = [(j1+j2).Mag() for j1 in em7s1 for j2 in itertools.combinations(jets, 2)]
    passed   = tdt.isPassed('L1_JPSI-1M5-EM7')
    emulated = any(1.0<m and m<5.0 for m in masses)
    emulated_m = next(m for m in masses if 1.0<m and m<5.0) if emulated else None
    counters['L1_JPSI-1M5-EM7']['any' ] += 1
    counters['L1_JPSI-1M5-EM7']['pass'] += (1 if passed else 0)
    counters['L1_JPSI-1M5-EM7']['emul'] += (1 if emulated else 0)
    h_any  = histos['h_any_l1jpsi_mult_m' ]
    h_pass = histos['h_pass_l1jpsi_mult_m']
    h_emul = histos['h_emul_l1jpsi_mult_m']
    mult = len(masses)
    max_m = h_any.GetYaxis().GetXmax()
    for m in masses:
        m = min([m, max_m]) # overflow
        h_any .Fill(mult, m)
        if passed: h_pass.Fill(mult, m)
    if emulated: h_emul.Fill(mult, emulated_m)
    def shift_mult_overflow(v, h):
        max_emall_mult =h.GetXaxis().GetXmax()-1.0 # todo use binwidth?
        return min([v, max_emall_mult])
    h_match = histos['h_l1jpsi_emall_mult_match']
    h_match.Fill(shift_mult_overflow(len(emall), h_match), hw_emul_binc(passed, emulated))

def check_L1_BTAG_MU4J15(l1mus=[], l1jets=[], tdt=None, histograms={}, counters={}):
    def addp4(o):
        o.p4 = ROOT.TLorentzVector(0.0,0.0,0.0,0.0)
        o.p4.SetPtEtaPhiM(mev2gev*o.et8x8(), o.eta(), o.phi(), 0.0)
        return o
    mu4ab = [m for m in l1mus if m.Pt()*mev2gev>4.0] # todo: ask if the code has > or >= (here 1 bit: can diff)
    cj15ab = [addp4(j) for j in l1jets if j.et8x8()*mev2gev>15.0 and abs(j.eta())<2.6]
    drs = [j.p4.DeltaR(m) for j in cj15ab for m in mu4ab]
    trigger = 'L1_BTAG_MU4J15'
    passed   = tdt.isPassed(trigger)
    emulated = any(0.0<dr and dr<0.4 for dr in drs)
    emulated_dr = next(dr for dr in drs if 0.0<dr and dr<0.4) if emulated else None
    counts = counters[trigger]
    counts['any' ] += 1
    counts['pass'] += (1 if passed else 0)
    counts['emul'] += (1 if emulated else 0)
    h_any  = histograms['h_any_l1btag_mult_dr' ]
    h_pass = histograms['h_pass_l1btag_mult_dr']
    h_emul = histograms['h_emul_l1btag_mult_dr']
    max_mult = h_any.GetXaxis().GetXmax()-0.5
    max_dr   = h_any.GetYaxis().GetXmax()
    mult = min([max_mult, len(drs)]) # overflows
    for dr in drs:
        dr = min([dr, max_dr]) # overflow
        h_any .Fill(mult, dr)
        if passed: h_pass.Fill(mult, dr)
    if emulated:
        h_emul.Fill(mult, emulated_dr)
        print("L1_BTAG_MU4J15 pass emul w/ dr %.2f, mult %.2f"%(emulated_dr, mult))
    def shift_mult_overflow(v, h):
        max_mult =h.GetXaxis().GetXmax()-1.0 # todo use binwidth?
        return min([v, max_mult])
    h_match_mu4ab  = histograms['h_l1btag_mu4ab_mult_match' ]
    h_match_cj15ab = histograms['h_l1btag_cj15ab_mult_match']
    h_match_mu4ab .Fill(shift_mult_overflow(len(mu4ab) , h_match_mu4ab) , hw_emul_binc(passed, emulated))
    h_match_cj15ab.Fill(shift_mult_overflow(len(cj15ab), h_match_cj15ab), hw_emul_binc(passed, emulated))

if __name__=='__main__':
    main()
