from math import pi,sqrt,cos
import ctypes
import ROOT
from ROOT import *
class OTriggerHLT() :
    def __init__(self) :
        #mask = ROOT.UInt_t(4193281)
##         ROOT.gROOT.ProcessLine('#include "ElectronPhotonSelectorTools/AsgPhotonIsEMSelector.h"')
##         self.offlinePhotonEMSelector = ROOT.AsgPhotonIsEMSelector("offlinePhotonEMSelector" )
##         #print ROOT.egammaPID.PhotonLoose
##         self.offlinePhotonEMSelector.setIsemValue(ROOT.egammaPID.PhotonLoose )
##         self.offlinePhotonEMSelector.setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150518/PhotonIsEMLooseSelectorCutDefs.conf")
##         self.offlinePhotonEMSelector.initialize()
##         #self.offlinePhotonEMSelector.print()
        
##         self.triggerPhotonEMSelector = ROOT.AsgPhotonIsEMSelector("triggerPhotonEMSelector" )
##         self.triggerPhotonEMSelector.setIsemValue(ROOT.egammaPID.PhotonLoose )
##         self.triggerPhotonEMSelector.setProperty("ConfigFile","ElectronPhotonSelectorTools/trigger/mc15_20150429/PhotonIsEMLooseSelectorCutDefs.conf")
##         self.triggerPhotonEMSelector.initialize();


        #ROOT.gROOT.ProcessLine('#include "ElectronPhotonSelectorTools/AsgPhotonIsEMSelector.h"')
        ROOT.gROOT.ProcessLine('#include "MuCTPIPhase0Upgrade/MIOCTPhase0TopoDataDecoder.h"')
        ROOT.gROOT.ProcessLine('#include "MuCTPIPhase0Upgrade/MIOCTPhase0TopoDataEncoderxAOD.h"')
        self.MuCTPiToTopoEncoder = ROOT.MuCTPIPhase0Upgrade.MIOCTPhase0TopoDataEncoderxAOD(ROOT.MuCTPIPhase0Upgrade.miniroi2cands16bitFinal1 ) #, ROOT.MuCTPIPhase0Upgrade.DEBUG )        
        self.MuCTPiToTopoDecoder = ROOT.MuCTPIPhase0Upgrade.MIOCTPhase0TopoDataDecoder(ROOT.MuCTPIPhase0Upgrade.miniroi2cands16bitFinal1 ) #, ROOT.MuCTPIPhase0Upgrade.DEBUG)        

        self.LVL1EmTauRoIs = []
        self.offlinePhotons = []
        self.trigPhotons = []

        self.LVL1MuonRoIs = []    
        self.offlineMuons = []
        self.trigMuons = []
        self.topoRoIs = []

        self.L1EM_scale = 2.

        self.triggerDispatch = {
            #'HLT_mu4':  lambda x : self.passMuonPtCut(x,4000),

            'L1_MU4'  : lambda x : self.passL1Muon(3900),
            'L1_MU6'  : lambda x : self.passL1Muon(5900),
            'L1_MU10' : lambda x : self.passL1Muon(9900),
            'L1_MU11' : lambda x : self.passL1Muon(10900),
            'L1_MU15' : lambda x : self.passL1Muon(14900),
            'L1_MU20' : lambda x : self.passL1Muon(19900),

            'L1_2MU4'      : lambda x : self.passL1diMuon(3900,3900),
            'L1_2MU6'      : lambda x : self.passL1diMuon(5900,5900),
            'L1_MU6_2MU4' : lambda x : self.passL1diMuon(5900,3900),
            'L1_MU10_2MU6' : lambda x : self.passL1diMuon(9900,5900),
            'L1_LFV-MU-F' : lambda x : self.passL1diMuondR(9900,5900,1.0),
            'L1_LFV-MU' : lambda x : self.passL1MuCTPidiMuondR(9900,5900,1.0),
            'L1_3MU4'      : lambda x : self.passL1triMuon(3900,3900, 3900),
            'L1_3MU6'      : lambda x : self.passL1triMuon(5900,5900, 5900),
            'L1_MU6_3MU4'  : lambda x : self.passL1triMuon(5900,3900, 3900),
            'L1_2MU6_3MU4' : lambda x : self.passL1triMuon(5900,5900, 3900),
            'L1_3MU4_LFV-MU-F' : lambda x : self.passL1triMuon(3900,3900, 3900) and self.passL1diMuondR(3900,3900,0.4),
            'L1_MU6_3MU4_LFV-MU-F' : lambda x : self.passL1triMuon(5900,3900, 3900) and self.passL1diMuondR(5900,3900,0.4),
            'L1_3MU4_LFV-MU' : lambda x : self.passL1triMuon(3900,3900, 3900) and self.passL1MuCTPidiMuondR(3900,3900,0.4),
            'L1_MU6_3MU4_LFV-MU' : lambda x : self.passL1triMuon(5900,3900, 3900) and self.passL1MuCTPidiMuondR(5900,3900,0.4),

            'L1_EM7'  : lambda x : self.passL1EM( 7000),
            'L1_EM8'  : lambda x : self.passL1EM( 8000),
            'L1_EM8I' : lambda x : self.passL1EM( 8000, require_iso=True),
            'L1_EM10' : lambda x : self.passL1EM( 10000),
            'L1_EM12' : lambda x : self.passL1EM( 12000),
            'L1_EM15' : lambda x : self.passL1EM( 15000),

            'L1_EM15_MU4' : lambda x : (self.passL1Muon(3900) and self.passL1EM( 15000)),
            'L1_EM7_MU10' : lambda x : (self.passL1Muon(9900) and self.passL1EM( 7000)),
            'L1_EM8I_MU10': lambda x : (self.passL1Muon(9900) and self.passL1EM( 8000, require_iso=True)),
            'L1_LFV-EM8I' : lambda x : self.passL1_EMMU_LFV(8000,9900, 0.4, 0.3, require_iso=True),
            'L1_LFV-EM15I' : lambda x : self.passL1_EMMU_LFV(15000,3900, 0.4, 0.3, require_iso=True),
            
            'HLT_mu10'             :  lambda x : self.passMuonPtCut(10000) and self.passL1Muon(9900),
            'HLT_mu20'             :  lambda x :  self.passMuonPtCut(20000) and self.passL1Muon(19900),
            'HLT_g10_etcut'        :  lambda x : (self.passPhotonEtCut(10000)),
            'HLT_g10_etcut_L1EM7'  :  lambda x : (self.passPhotonEtCut(10000)  and self.emulateDecision('L1_EM7')),
            'HLT_g20_etcut_L1EM12' :  lambda x : (self.passPhotonEtCut(20000) and self.emulateDecision('L1_EM12')),
            'HLT_g10_loose'        :  lambda x : (self.passPhotonLooseEtCut(10000) and self.emulateDecision('L1_EM7')),
            'HLT_g10_looseNoHad'   :  lambda x : (self.passPhotonLooseNoHadEtCut(10000) and self.emulateDecision('L1_EM7')),

            'HLT_g10_etcut_mu10'   :  lambda x : self.passMuonPtCut(10000) and self.passPhotonEtCut(10000) and  self.emulateDecision('L1_EM7_MU10'),
            'HLT_g10_etcut_mu10_taumass'   : lambda x : (self.passTauMass(10000,10000) and self.emulateDecision('L1_LFV-EM8I')),        
            'HLT_g10_etcut_mu10_L1LFV-EM8I': lambda x : (self.passPhotonEtCut(10000)
                                                         and  self.passMuonPtCut(10000)
                                                         and self.emulateDecision('L1_LFV-EM8I')), 
            'HLT_g10_etcut_L1EM7_mu10_taumass': lambda x : (self.passTauMass(10000,10000)
                                                            and self.emulateDecision('L1_EM7_MU10')),
            'HLT_g10_looseNoHad_mu10_taumass' : lambda x : (self.passTauMass(10000,10000)
                                                            and self.passPhotonLooseNoHadEtCut(10000)
                                                            and self.emulateDecision('L1_EM7_MU10')),
            'HLT_g10_loose_mu10_taumass' : lambda x : (self.passTauMass(10000,10000)
                                                       and self.passPhotonLooseEtCut(10000)
                                                       and self.emulateDecision('L1_EM7_MU10')),
            
            }


##################################################################
    def deltaPhi( self, phi1, phi2):
        PHI = abs(phi1-phi2)
        while PHI >= pi : PHI = PHI - 2.*pi
        while PHI < -pi : PHI = PHI + 2.*pi
        return PHI

    def getPhotons(self, ch) :
        self.LVL1EmTauRoIs = []
        if hasattr(ch, "LVL1EmTauRoIs") :
            self.LVL1EmTauRoIs = ch.LVL1EmTauRoIs
        self.trigPhotons = []
        if hasattr(ch, "HLT_xAOD__PhotonContainer_egamma_Photons") :
            for iph in xrange(len(ch.HLT_xAOD__PhotonContainer_egamma_Photons )) :
                ph = ch.HLT_xAOD__PhotonContainer_egamma_Photons[iph]
                duplicate = False
                for tph in self.trigPhotons :
                    if ph.p4().DeltaR(tph.p4() ) < 0.03 :
                        duplicate =  True
                        break                
                if not duplicate :
                    self.trigPhotons.append(ph)
                
        self.offlinePhotons =  []
        if hasattr(ch, "Photons") :
            for iph in xrange(len(ch.Photons )) :
                ph = ch.Photons[iph]
                duplicate = False
                for tph in self.offlinePhotons :
                    if ph.p4().DeltaR(tph.p4() ) < 0.03 :
                        duplicate =  True
                        break                
                if not duplicate :
                    self.offlinePhotons.append(ph)


    def passPhotonEtCut(self, ETCUT ) :
        for iph in xrange(len(self.trigPhotons)) :
            ph = self.trigPhotons[iph]
            #print " Found Photon ", ph.pt(), ph.caloCluster().et()
            #if ph.caloCluster().et() >= ETCUT :
            if ph.pt() >= ETCUT :
                return True
        return False

    def passPhotonLooseEtCut(self, ETCUT ) :
        for iph in xrange(len(self.trigPhotons)) :
            ph = self.trigPhotons[iph]
            #if ph.caloCluster().et() >= ETCUT :
            if ph.pt() >= ETCUT :
                loose = self.triggerPhotonEMSelector.accept(ph);        
                if loose :
                    #print "found loose photon", ph.pt(), " isEM",self.triggerPhotonEMSelector.IsemValue()
                    return True
                #else :
                #    print " photon ",ph.pt()," is lost isEM = ", self.triggerPhotonEMSelector.IsemValue()
        return False
    
    def passPhotonLooseNoHadEtCut(self,  ETCUT ) :
        for iph in xrange(len(self.trigPhotons)) :
            ph = self.trigPhotons[iph]
            #if ph.caloCluster().et() >= ETCUT :
            if ph.pt() >= ETCUT :
                loose = self.triggerPhotonEMSelector.accept(ph);        
                ph_isEM = self.triggerPhotonEMSelector.IsemValue()
                if (ph_isEM &4192256)==0 :
                    return True
        return False

    def passL1EM(self, PTCUT, require_iso = False) :
        if len(self.LVL1EmTauRoIs) == 0 : return False
        for l1_em in  self.LVL1EmTauRoIs :
            if type(l1_em).__name__ == 'xAOD::EmTauRoI_v2' and l1_em.roiType() != 1 : continue
            #if l1_em.emClus() <= PTCUT : continue
            if float(l1_em.eT())/self.L1EM_scale <= PTCUT : continue
            if require_iso:
                if float(l1_em.eT())/self.L1EM_scale > 50000 or l1_em.emIsol()/1000. <=  2.:  #max(1.0,-1.8+l1_em.emClus()/8000.)  :
                    return True
            else :
                return True
        return False

                
##################################################################
    def getMuons(self, ch) :
        self.LVL1MuonRoIs = []
        #if hasattr(ch, "LVLMuonRoIs") :
        self.LVL1MuonRoIs = ch.LVL1MuonRoIs
        self.trigMuons = []
        if hasattr(ch, "HLT_xAOD__MuonContainer_MuonEFInfo") :
            for imu in xrange(len(ch.HLT_xAOD__MuonContainer_MuonEFInfo)):
                mu = ch.HLT_xAOD__MuonContainer_MuonEFInfo[imu]
                duplicate = False
                for tmu in self.trigMuons :
                    if mu.p4().DeltaR(tmu.p4() ) < 0.005 :
                        duplicate =  True
                        break                
                if not duplicate :
                    self.trigMuons.append(mu)
                    #print "add TrigMuon ", mu.pt()

        self.offlineMuons = []
        if hasattr(ch, "Muons") :
            for imu in xrange(len(ch.Muons )) :
                mu = ch.Muons[imu]
                duplicate = False
                for tmu in self.offlineMuons :
                    if mu.p4().DeltaR(tmu.p4() ) < 0.005 :
                        duplicate =  True
                        break                
                if not duplicate :
                    self.offlineMuons.append(mu)
                    #print "add offline muon ", mu.pt()

    def passMuonPtCut(self, PTCUT) :
        for imuon in xrange(len(self.trigMuons)) :
            muonFull=self.trigMuons[imuon]
            if not muonFull.primaryTrackParticle() : continue
            muon = muonFull.trackParticle(ROOT.xAOD.Muon.ExtrapolatedMuonSpectrometerTrackParticle)
            if not muon : continue        
            if abs(muon.pt()) > PTCUT :
                #print " Found HLT muon ", muon.pt()
                return True
        return False


    def passL1Muon(self, PTCUT) :
        if len(self.LVL1MuonRoIs) == 0 : return False
        for l1_mu in  self.LVL1MuonRoIs :
            if l1_mu.isVetoed() : continue
            if l1_mu.thrValue() > PTCUT : return True
        return False

    def passL1diMuon(self, PTCUT_high, PTCUT_low) :
        pts = sorted( [l1_mu.thrValue() for l1_mu in  self.LVL1MuonRoIs if not l1_mu.isVetoed() ], reverse =True )
        #pts = sorted( [l1_mu.thrValue() for l1_mu in  self.LVL1MuonRoIs  ], reverse =True )
        return (len(pts)>1 and pts[0] > PTCUT_high and pts[1]>PTCUT_low)

    def passL1triMuon(self, PTCUT_high, PTCUT_low, PTCUT_third) :
        pts = sorted( [l1_mu.thrValue() for l1_mu in  self.LVL1MuonRoIs if not l1_mu.isVetoed() ], reverse =True )
        #pts = sorted( [l1_mu.thrValue() for l1_mu in  self.LVL1MuonRoIs  ], reverse =True )
        return (len(pts)>2 and pts[0] > PTCUT_high and pts[1]>PTCUT_low and pts[2]>PTCUT_third)

    def passL1diMuondR(self, PTCUT_high, PTCUT_low, DR) :
        for imu1 in  xrange(0,len(self.LVL1MuonRoIs)) :
            if self.LVL1MuonRoIs[imu1].isVetoed() : continue
            pt1 = self.LVL1MuonRoIs[imu1].thrValue()
            if  pt1 <= PTCUT_low : continue            
            for imu2 in  xrange(imu1+1,len(self.LVL1MuonRoIs)) :
                if self.LVL1MuonRoIs[imu2].isVetoed() : continue
                pt2 = self.LVL1MuonRoIs[imu2].thrValue()
                if  pt2 <= PTCUT_low : continue            
                if pt1 <=  PTCUT_high and pt2 <=  PTCUT_high : continue

                dPhi = self.deltaPhi( self.LVL1MuonRoIs[imu1].phi(),  self.LVL1MuonRoIs[imu2].phi())
                dEta =  self.LVL1MuonRoIs[imu1].eta() - self.LVL1MuonRoIs[imu2].eta()
                if dPhi*dPhi + dEta*dEta  < DR*DR :
                    return True
        return False
                
    def passL1MuCTPidiMuondR(self, PTCUT_high, PTCUT_low, DR) :
        self.MuCTPiToTopoEncoder.clearRoIs()
        dataWords = self.MuCTPiToTopoEncoder.encode(self.LVL1MuonRoIs)
        self.topoRoIs = self.MuCTPiToTopoDecoder.decode(dataWords)

        for imu1 in  xrange(0,len(self.topoRoIs)) :
            mu1 = self.topoRoIs[imu1]
            
            pt1 = mu1.pT()
            if   pt1 == 1 : pt1 = 4000
            elif pt1 == 2 : pt1 = 6000
            else          : pt1 = 10000

            if  pt1 <= PTCUT_low : continue            
            for imu2 in  xrange(imu1+1,len(self.topoRoIs)) :
                mu2 = self.topoRoIs[imu2]
                pt2 = mu2.pT()
                if   pt2 == 1 : pt2 = 4000
                elif pt2 == 2 : pt2 = 6000
                else          : pt2 = 10000
                if  pt2 <= PTCUT_low : continue            
                if pt1 <=  PTCUT_high and pt2 <=  PTCUT_high : continue

                dPhi = self.deltaPhi( mu1.phi(),  mu2.phi())
                dEta =  mu1.eta() - mu2.eta()
                if dPhi*dPhi + dEta*dEta  < DR*DR :
                    return True
        dataWords.clear()
        self.topoRoIs.clear()
        return False
                
                
                

##################################################################
# LFV part
    def passL1_EMMU_LFV(self, EMPTCUT, MUPTCUT, DETA, DPHI, require_iso = False) :
        if len(self.LVL1EmTauRoIs) == 0 : return False
        if len(self.LVL1MuonRoIs) == 0 : return False
        l1EMRoIs = []
        for l1_em in  self.LVL1EmTauRoIs :
            if l1_em.roiType() != 1 : continue
            if float(l1_em.eT())/self.L1EM_scale <= EMPTCUT : continue
            #if l1_em.emClus() <= EMPTCUT : continue
            if require_iso:
                #if l1_em.emClus() > 50000 or l1_em.emIsol()/1000. <=  2.:  #max(1.0,-1.8+l1_em.emClus()/8000.)  :
                if float(l1_em.eT())/self.L1EM_scale > 50000 or l1_em.emIsol()/1000. <=  2.:  #max(1.0,-1.8+l1_em.emClus()/8000.)  :
                    l1EMRoIs.append(l1_em)
            else :
                l1EMRoIs.append(l1_em)

        l1MURoIs = []
        for l1_mu in  self.LVL1MuonRoIs :
            if l1_mu.isVetoed() : continue
            if l1_mu.thrValue() > MUPTCUT :
                l1MURoIs.append(l1_mu)

        # last step - find closeby pair
        for l1_em in l1EMRoIs :
            for l1_mu in l1MURoIs :
                dEta =  abs(l1_em.eta() - l1_mu.eta())
                if dEta > DETA : continue
                dPhi = abs(self.deltaPhi(l1_em.phi(), l1_mu.phi()))
                if dPhi > DPHI : continue
                #print " Accept ", l1_em.emClus(), "\t", l1_mu.thrValue(), "\t dEta=",dEta, "\t dPhi=",dPhi
                return True
        # nothing is found
        return False


    def passTauMass(self, EMPTCUT, MUPTCUT ) :
        gCand = []
        for iph in xrange(len(self.trigPhotons)) :
            ph = self.trigPhotons[iph]
            #if ph.caloCluster().et() >= EMPTCUT :
            if ph.pt() >= EMPTCUT :
                gCand.append(ph)
        muCand = []
        for imuon in xrange(len(self.trigMuons)) :
            muon=self.trigMuons[imuon]
            if muon.pt() > MUPTCUT : muCand.append(muon)

        for ph in gCand :
            for mu in muCand :
                if self.deltaPhi(ph.phi(), mu.phi()) > 1.5 : continue
                if ph.p4().DeltaR(mu.p4()) > 0.6 : continue
                tau = ph.p4() + mu.p4()
                if tau.M() > 2500 : continue
                return True
        # nothing was found
        return False


##################################################################
    def setEvent(self, ch ) :
        self.getPhotons(ch)
        self.getMuons(ch)

            
            
####################################################

        
    def emulateDecision(self, trig) :
        if trig in self.triggerDispatch.keys() :
            return self.triggerDispatch[trig](trig)
        else :
            print "Unknown trigger ", trig
            return False
            
            
