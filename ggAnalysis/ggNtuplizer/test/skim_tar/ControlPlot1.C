#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "math.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TInterpreter.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string> //remove space for successful compilation
#include <ostream> //remove space for successful compilation
#include <utility>
#include <vector>
using namespace std;
class ControlPlot {
public :
  ControlPlot(std::string& fname, bool muon);
  void Loop(std::vector<double>& puwt, std::vector<double> W_events,
	    std::vector<double>& DY_events);
  void SavePlot();
private:
  void                   InitTree();
  void                   plotFill(std::string name, float x, int nx, 
				  float nxmin, float nxmax, double weight=1);
  TH1F*                  nplot1(std::string name);
  float                  weightCalc(TH1F *Histo, std::string& outputName,
				    int njet,  std::vector<double>& W_events, 
				    std::vector<double>& DY_events);

  std::string                       fname_;
  TFile                            *f_Double_, *fout_;
  TH1F                             *HistoTot_, *h_mass_;
  TTree                            *Run_Tree_;
  bool                              Muon_;
  int                               njets;
  std::map<std::string, TH1F*>      myMap1_;

  // Declaration of leaf types
  Int_t                             run;
  Long64_t                          event;
  Int_t                             lumis;
  Bool_t                            isData;
  Int_t                             nVtx;
  Int_t                             nTrksPV;
  Bool_t                            isPVGood;
  Bool_t                            hasGoodVtx;
  Float_t                           vtx;
  Float_t                           vty;
  Float_t                           vtz;
  Float_t                           rho;
  Float_t                           rhoCentral;
  ULong64_t                         HLTEleMuX;
  ULong64_t                         HLTPho;
  ULong64_t                         HLTJet;
  ULong64_t                         HLTEleMuXIsPrescaled;
  ULong64_t                         HLTPhoIsPrescaled;
  ULong64_t                         HLTJetIsPrescaled;
  std::vector<float>               *pdf;
  Float_t                           pthat;
  Float_t                           processID;
  Float_t                           genWeight;
  Float_t                           num_gen_jets;
  Float_t                           genHT;
  TString                           *EventTag;
  Int_t                             nPUInfo;
  std::vector<int>                 *nPU;
  std::vector<int>                 *puBX;
  std::vector<float>               *puTrue;
  Int_t                             nMC;
  std::vector<int>                 *mcPID;
  std::vector<float>               *mcVtx;
  std::vector<float>               *mcVty;
  std::vector<float>               *mcVtz;
  std::vector<float>               *mcPt;
  std::vector<float>               *mcMass;
  std::vector<float>               *mcEta;
  std::vector<float>               *mcPhi;
  std::vector<float>               *mcE;
  std::vector<float>               *mcEt;
  std::vector<int>                 *mcGMomPID;
  std::vector<int>                 *mcMomPID;
  std::vector<float>               *mcMomPt;
  std::vector<float>               *mcMomMass;
  std::vector<float>               *mcMomEta;
  std::vector<float>               *mcMomPhi;
  std::vector<int>                 *mcIndex;
  std::vector<unsigned short>      *mcStatusFlag;
  std::vector<int>                 *mcParentage;
  std::vector<int>                 *mcStatus;
  std::vector<float>               *mcCalIsoDR03;
  std::vector<float>               *mcTrkIsoDR03;
  std::vector<float>               *mcCalIsoDR04;
  std::vector<float>               *mcTrkIsoDR04;
  std::vector<float>               *mcHadronPt;
  std::vector<float>               *mcHadronEta;
  std::vector<float>               *mcHadronPhi;
  std::vector<float>               *mcHadronMass;
  std::vector<float>               *mcHadronE;
  std::vector<float>               *mcHadronGMomPID;
  std::vector<float>               *mcHadronMomPID;
  std::vector<float>               *mcHadronMomPt;
  std::vector<float>               *mcHadronMomMass;
  std::vector<float>               *mcHadronMomEta;
  std::vector<float>               *mcHadronMomPhi;

  Int_t                             metFilters;
  Float_t                           genMET;
  Float_t                           genMETPhi;
  Float_t                           pfMET;
  Float_t                           pfMETPhi;
  Float_t                           pfMETsumEt;
  Float_t                           pfMETmEtSig;
  Float_t                           pfMETSig;
  Float_t                           pfMET_T1JERUp;
  Float_t                           pfMET_T1JERDo;
  Float_t                           pfMET_T1JESUp;
  Float_t                           pfMET_T1JESDo;
  Float_t                           pfMET_T1MESUp;
  Float_t                           pfMET_T1MESDo;
  Float_t                           pfMET_T1EESUp;
  Float_t                           pfMET_T1EESDo;
  Float_t                           pfMET_T1PESUp;
  Float_t                           pfMET_T1PESDo;
  Float_t                           pfMET_T1TESUp;
  Float_t                           pfMET_T1TESDo;
  Float_t                           pfMET_T1UESUp;
  Float_t                           pfMET_T1UESDo;
  Int_t                             nEle;
  std::vector<int>                 *eleCharge;
  std::vector<int>                 *eleChargeConsistent;
  std::vector<float>               *eleEn;
  std::vector<float>               *eleSCEn;
  std::vector<float>               *eleESEn;
  std::vector<float>               *eleESEnP1;
  std::vector<float>               *eleESEnP2;
  std::vector<float>               *eleD0;
  std::vector<float>               *eleDz;
  std::vector<float>               *elePt;
  std::vector<float>               *eleEta;
  std::vector<float>               *elePhi;
  std::vector<float>               *eleR9;
  std::vector<float>               *eleCalibPt;
  std::vector<float>               *eleCalibEn;
  std::vector<float>               *eleSCEta;
  std::vector<float>               *eleSCPhi;
  std::vector<float>               *eleSCRawEn;
  std::vector<float>               *eleSCEtaWidth;
  std::vector<float>               *eleSCPhiWidth;
  std::vector<float>               *eleHoverE;
  std::vector<float>               *eleEoverP;
  std::vector<float>               *eleEoverPout;
  std::vector<float>               *eleEoverPInv;
  std::vector<float>               *eleBrem;
  std::vector<float>               *eledEtaAtVtx;
  std::vector<float>               *eledPhiAtVtx;
  std::vector<float>               *eledEtaAtCalo;
  std::vector<float>               *eleSigmaIEtaIEta;
  std::vector<float>               *eleSigmaIEtaIPhi;
  std::vector<float>               *eleSigmaIPhiIPhi;
  std::vector<float>               *eleSigmaIEtaIEtaFull5x5;
  std::vector<float>               *eleSigmaIPhiIPhiFull5x5;
  std::vector<int>                 *eleConvVeto;
  std::vector<int>                 *eleMissHits;
  std::vector<float>               *eleESEffSigmaRR;
  std::vector<float>               *elePFChIso;
  std::vector<float>               *elePFPhoIso;
  std::vector<float>               *elePFNeuIso;
  std::vector<float>               *elePFPUIso;
  std::vector<float>               *elePFClusEcalIso;
  std::vector<float>               *elePFClusHcalIso;
  std::vector<float>               *elePFMiniIso;
  std::vector<float>               *eleIDMVANonTrg;
  std::vector<float>               *eleIDMVA;
  std::vector<float>               *eleIDMVATrg;
  std::vector<float>               *eledEtaseedAtVtx;
  std::vector<float>               *eleE1x5;
  std::vector<float>               *eleE2x5;
  std::vector<float>               *eleE5x5;
  std::vector<float>               *eleE1x5Full5x5;
  std::vector<float>               *eleE2x5Full5x5;
  std::vector<float>               *eleE5x5Full5x5;
  std::vector<float>               *eleR9Full5x5;
  std::vector<int>                 *eleEcalDrivenSeed;
  std::vector<float>               *eleDr03EcalRecHitSumEt;
  std::vector<float>               *eleDr03HcalDepth1TowerSumEt;
  std::vector<float>               *eleDr03HcalDepth2TowerSumEt;
  std::vector<float>               *eleDr03HcalTowerSumEt;
  std::vector<float>               *eleDr03TkSumPt;
  std::vector<float>               *elecaloEnergy;
  std::vector<float>               *eleTrkdxy;
  std::vector<float>               *eleKFHits;
  std::vector<float>               *eleKFChi2;
  std::vector<std::vector<float> > *eleGSFPt;
  std::vector<std::vector<float> > *eleGSFEta;
  std::vector<std::vector<float> > *eleGSFPhi;
  std::vector<std::vector<float> > *eleGSFCharge;
  std::vector<std::vector<int> > *eleGSFHits;
  std::vector<std::vector<int> > *eleGSFMissHits;
  std::vector<std::vector<int> > *eleGSFNHitsMax;
  std::vector<std::vector<float> > *eleGSFVtxProb;
  std::vector<std::vector<float> > *eleGSFlxyPV;
  std::vector<std::vector<float> > *eleGSFlxyBS;
  std::vector<std::vector<float> > *eleBCEn;
  std::vector<std::vector<float> > *eleBCEta;
  std::vector<std::vector<float> > *eleBCPhi;
  std::vector<std::vector<float> > *eleBCS25;
  std::vector<std::vector<float> > *eleBCS15;
  std::vector<std::vector<float> > *eleBCSieie;
  std::vector<std::vector<float> > *eleBCSieip;
  std::vector<std::vector<float> > *eleBCSipip;
  std::vector<int>                 *eleFiredTrgs;
  std::vector<unsigned short>      *eleIDbit;
  Int_t                             nMu;
  std::vector<float>               *muPt;
  std::vector<float>               *muEn;
  std::vector<float>               *muEta;
  std::vector<float>               *muPhi;
  std::vector<int>                 *muCharge;
  std::vector<int>                 *muType;
  std::vector<bool>                *muIsLooseID;
  std::vector<bool>                *muIsMediumID;
  std::vector<bool>                *muIsTightID;
  std::vector<bool>                *muIsSoftID;
  std::vector<bool>                *muIsHighPtID;
  std::vector<float>               *muD0;
  std::vector<float>               *muDz;
  std::vector<float>               *muChi2NDF;
  std::vector<float>               *muInnerD0;
  std::vector<float>               *muInnerDz;
  std::vector<int>                 *muTrkLayers;
  std::vector<int>                 *muPixelLayers;
  std::vector<int>                 *muPixelHits;
  std::vector<int>                 *muMuonHits;
  std::vector<int>                 *muStations;
  std::vector<int>                 *muMatches;
  std::vector<int>                 *muTrkQuality;
  std::vector<float>               *muIsoTrk;
  std::vector<float>               *muPFChIso;
  std::vector<float>               *muPFPhoIso;
  std::vector<float>               *muPFNeuIso;
  std::vector<float>               *muPFPUIso;
  std::vector<float>               *muPFMiniIso;
  std::vector<int>                 *muFiredTrgs;
  std::vector<float>               *muInnervalidFraction;
  std::vector<float>               *musegmentCompatibility;
  std::vector<float>               *muchi2LocalPosition;
  std::vector<float>               *mutrkKink;
  std::vector<float>               *muBestTrkPtError;
  std::vector<float>               *muBestTrkPt;
  std::vector<UShort_t>            *muIDbit;
  Int_t                             nTau;
  std::vector<bool>                *taupfTausDiscriminationByDecayModeFinding;
  std::vector<bool>                *taupfTausDiscriminationByDecayModeFindingNewDMs;
  std::vector<bool>                *tauByMVA6VLooseElectronRejection;
  std::vector<bool>                *tauByMVA6LooseElectronRejection;
  std::vector<bool>                *tauByMVA6MediumElectronRejection;
  std::vector<bool>                *tauByMVA6TightElectronRejection;
  std::vector<bool>                *tauByMVA6VTightElectronRejection;
  std::vector<bool>                *tauByLooseMuonRejection3;
  std::vector<bool>                *tauByTightMuonRejection3;
  std::vector<bool>                *tauByLooseCombinedIsolationDeltaBetaCorr3Hits;
  std::vector<bool>                *tauByMediumCombinedIsolationDeltaBetaCorr3Hits;
  std::vector<bool>                *tauByTightCombinedIsolationDeltaBetaCorr3Hits;
  std::vector<float>               *tauCombinedIsolationDeltaBetaCorrRaw3Hits;
  std::vector<float>               *tauByIsolationMVArun2v1DBnewDMwLTraw;
  std::vector<float>               *tauByIsolationMVArun2v1DBoldDMwLTraw;
  std::vector<float>               *tauByIsolationMVArun2v1PWnewDMwLTraw;
  std::vector<float>               *tauByIsolationMVArun2v1PWoldDMwLTraw;
  std::vector<bool>                *tauByVTightIsolationMVArun2v1DBnewDMwLT;
  std::vector<bool>                *tauByVTightIsolationMVArun2v1DBoldDMwLT;
  std::vector<bool>                *tauByVTightIsolationMVArun2v1PWnewDMwLT;
  std::vector<bool>                *tauByVTightIsolationMVArun2v1PWoldDMwLT;
  std::vector<bool>                *tauByTightIsolationMVArun2v1DBnewDMwLT;
  std::vector<bool>                *tauByTightIsolationMVArun2v1DBoldDMwLT;
  std::vector<bool>                *tauByTightIsolationMVArun2v1PWnewDMwLT;
  std::vector<bool>                *tauByTightIsolationMVArun2v1PWoldDMwLT;
  std::vector<bool>                *tauByMediumIsolationMVArun2v1DBnewDMwLT;
  std::vector<bool>                *tauByMediumIsolationMVArun2v1DBoldDMwLT;
  std::vector<bool>                *tauByMediumIsolationMVArun2v1PWnewDMwLT;
  std::vector<bool>                *tauByMediumIsolationMVArun2v1PWoldDMwLT;
  std::vector<bool>                *tauByLooseIsolationMVArun2v1DBnewDMwLT;
  std::vector<bool>                *tauByLooseIsolationMVArun2v1DBoldDMwLT;
  std::vector<bool>                *tauByLooseIsolationMVArun2v1PWnewDMwLT;
  std::vector<bool>                *tauByLooseIsolationMVArun2v1PWoldDMwLT;
  std::vector<bool>                *tauByVLooseIsolationMVArun2v1DBnewDMwLT;
  std::vector<bool>                *tauByVLooseIsolationMVArun2v1DBoldDMwLT;
  std::vector<bool>                *tauByVLooseIsolationMVArun2v1PWnewDMwLT;
  std::vector<bool>                *tauByVLooseIsolationMVArun2v1PWoldDMwLT;
  std::vector<float>               *tauEta;
  std::vector<float>               *tauPhi;
  std::vector<float>               *tauPt;
  std::vector<float>               *tauEt;
  std::vector<float>               *tauCharge;
  std::vector<float>               *tauP;
  std::vector<float>               *tauPx;
  std::vector<float>               *tauPy;
  std::vector<float>               *tauPz;
  std::vector<float>               *tauVz;
  std::vector<float>               *tauEnergy;
  std::vector<float>               *tauMass;
  std::vector<float>               *tauDxy;
  std::vector<float>               *tauZImpact;
  std::vector<int>                 *tauDecayMode;
  std::vector<bool>                *tauLeadChargedHadronExists;
  std::vector<float>               *tauLeadChargedHadronEta;
  std::vector<float>               *tauLeadChargedHadronPhi;
  std::vector<float>               *tauLeadChargedHadronPt;
  std::vector<float>               *tauChargedIsoPtSum;
  std::vector<float>               *tauNeutralIsoPtSum;
  std::vector<float>               *tauPuCorrPtSum;
  std::vector<int>                 *tauNumSignalPFChargedHadrCands;
  std::vector<int>                 *tauNumSignalPFNeutrHadrCands;
  std::vector<int>                 *tauNumSignalPFGammaCands;
  std::vector<int>                 *tauNumSignalPFCands;
  std::vector<int>                 *tauNumIsolationPFChargedHadrCands;
  std::vector<int>                 *tauNumIsolationPFNeutrHadrCands;
  std::vector<int>                 *tauNumIsolationPFGammaCands;
  std::vector<int>                 *tauNumIsolationPFCands;
  std::vector<float>               *taufootprintCorrection;
  std::vector<float>               *tauphotonPtSumOutsideSignalCone;
  std::vector<float>               *taudz;
  std::vector<float>               *taudxy;
  Int_t                             nJet;
  std::vector<float>               *jetPt;
  std::vector<float>               *jetEn;
  std::vector<float>               *jetEta;
  std::vector<float>               *jetPhi;
  std::vector<float>               *jetRawPt;
  std::vector<float>               *jetRawEn;
  std::vector<float>               *jetMt;
  std::vector<float>               *jetArea;
  std::vector<float>               *jetLeadTrackPt;
  std::vector<float>               *jetLeadTrackEta;
  std::vector<float>               *jetLeadTrackPhi;
  std::vector<int>                 *jetLepTrackPID;
  std::vector<float>               *jetLepTrackPt;
  std::vector<float>               *jetLepTrackEta;
  std::vector<float>               *jetLepTrackPhi;
  std::vector<float>               *jetpfCombinedInclusiveSecondaryVertexV2BJetTags;
  std::vector<float>               *jetJetProbabilityBJetTags;
  std::vector<float>               *jetpfCombinedMVAV2BJetTags;
  std::vector<int>                 *jetPartonID;
  std::vector<int>                 *jetHadFlvr;
  std::vector<int>                 *jetGenJetIndex;
  std::vector<float>               *jetGenJetEn;
  std::vector<float>               *jetGenJetPt;
  std::vector<float>               *jetGenJetEta;
  std::vector<float>               *jetGenJetPhi;
  std::vector<int>                 *jetGenPartonID;
  std::vector<float>               *jetGenEn;
  std::vector<float>               *jetGenPt;
  std::vector<float>               *jetGenEta;
  std::vector<float>               *jetGenPhi;
  std::vector<int>                 *jetGenPartonMomID;
  std::vector<bool>                *jetPFLooseId;
  std::vector<float>               *jetPUidFullDiscriminant;
  std::vector<float>               *jetJECUnc;
  std::vector<int>                 *jetFiredTrgs;
  std::vector<float>               *jetCHF;
  std::vector<float>               *jetCSV2BJetTags;
  std::vector<float>               *jetNHF;
  std::vector<float>               *jetCEF;
  std::vector<float>               *jetNEF;
  std::vector<int>                 *jetNCH;
  std::vector<float>               *jetVtxPt;
  std::vector<float>               *jetVtxMass;
  std::vector<float>               *jetVtxNtrks;
  std::vector<float>               *jetVtx3DVal;
  std::vector<float>               *jetVtx3DSig;
};

std::vector<double> W_EvenetMultiplicity() {
  std::vector<double> W_events;

  TFile * myFile_W0 = new TFile("WJetsToLNu.root");
  TH1F * Histo_W0 = (TH1F*) myFile_W0->Get("hcount");
  W_events.push_back(Histo_W0->GetBinContent(2));
  myFile_W0->Close();
  TFile * myFile_W1 = new TFile("W1JetsToLNu.root");
  TH1F * Histo_W1 = (TH1F*) myFile_W1->Get("hcount");
  W_events.push_back(Histo_W1->GetBinContent(2));
  myFile_W1->Close();
  TFile * myFile_W2 = new TFile("W2JetsToLNu.root");
  TH1F * Histo_W2 = (TH1F*) myFile_W2->Get("hcount");
  W_events.push_back(Histo_W2->GetBinContent(2));
  myFile_W2->Close();
  TFile * myFile_W3 = new TFile("W3JetsToLNu.root");
  TH1F * Histo_W3 = (TH1F*) myFile_W3->Get("hcount");
  W_events.push_back(Histo_W3->GetBinContent(2));
  myFile_W3->Close();
  TFile * myFile_W4 = new TFile("W4JetsToLNu.root");
  TH1F * Histo_W4 = (TH1F*) myFile_W4->Get("hcount");
  W_events.push_back(Histo_W4->GetBinContent(2));
  myFile_W4->Close();
  return W_events ;
}

std::vector<double> DY_EvenetMultiplicity() {
  std::vector<double> DY_events;

  TFile * myFile_DY0 = new TFile("DYJetsToLL.root");
  TH1F * Histo_DY0 = (TH1F*) myFile_DY0->Get("hcount");
  DY_events.push_back(Histo_DY0->GetBinContent(2));
  myFile_DY0->Close();
  TFile * myFile_DY1 = new TFile("DY1JetsToLL.root");
  TH1F * Histo_DY1 = (TH1F*) myFile_DY1->Get("hcount");
  DY_events.push_back(Histo_DY1->GetBinContent(2));
  myFile_DY1->Close();
  TFile * myFile_DY2 = new TFile("DY2JetsToLL.root");
  TH1F * Histo_DY2 = (TH1F*) myFile_DY2->Get("hcount");
  DY_events.push_back(Histo_DY2->GetBinContent(2));
  myFile_DY2->Close();  
  TFile * myFile_DY3 = new TFile("DY3JetsToLL.root");
  TH1F * Histo_DY3 = (TH1F*) myFile_DY3->Get("hcount");
  DY_events.push_back(Histo_DY3->GetBinContent(2));
  myFile_DY3->Close();
  TFile * myFile_DY4 = new TFile("DY4JetsToLL.root");
  TH1F * Histo_DY4 = (TH1F*) myFile_DY4->Get("hcount");
  DY_events.push_back(Histo_DY4->GetBinContent(2));
  myFile_DY4->Close();
  return DY_events ;
}

void PUWeight(std::vector<double>& puwt) {

  puwt.clear();
  TFile* PUData = TFile::Open("dataMoriondPU.root");
  TH1F* HistoPUData = (TH1F*) PUData->Get("pileup");
  HistoPUData->Scale(1.0/HistoPUData->Integral());
  int nbData = HistoPUData->GetNbinsX();

  TFile* PUMC= TFile::Open("mcMoriondPU.root");
  TH1F* HistoPUMC= (TH1F *) PUMC->Get("pileup");
  HistoPUMC->Scale(1.0/HistoPUMC->Integral());
  int nbMC = HistoPUMC->GetNbinsX();

  for (int k=1; k<=std::max(nbData,nbMC); ++k) {
    double wt(1.0);
    if (k <= nbData && k <= nbMC) {
      double PUMC_   = HistoPUMC->GetBinContent(k);
      double PUData_ = HistoPUData->GetBinContent(k);
      if (PUMC_ ==0)
	std::cout << "PUMC_ is zero!!! for bin " << k << "\n";
      else
	wt = PUData_/PUMC_;
    }
    puwt.push_back(wt);
  }
  PUData->Close(); PUMC->Close();
}

double XSection(std::string OutName) {
  ////////////////////////////////////////////////////////////////
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns
  ////////////////////////////////////////////////////////////////
    
  //    https://docs.google.com/spreadsheets/d/1rWM3AlFKO8IJVaeoQkWZYWwSvicQ1QCXYSzH74QyZqE/edit?alt=json#gid=398123591
    
  if (OutName.find("WJetsToLNu") != std::string::npos) return 50690;   // As we have large cut at Skim, this one is not needed
  else if (OutName.find("W1JetsToLNu") != std::string::npos) return 9644.5 ;
  else if (OutName.find("W2JetsToLNu") != std::string::npos) return 3144.5 ;
  else if (OutName.find("W3JetsToLNu") != std::string::npos) return 954.8 ;
  else if (OutName.find("W4JetsToLNu") != std::string::npos) return 485.6 ;

  else if (OutName.find("WJetsToLNu_HT-70to100") != std::string::npos) return 0;
  else if (OutName.find("WJetsToLNu_HT-100to200") != std::string::npos) return 1343;
  else if (OutName.find("WJetsToLNu_HT-200to400") != std::string::npos) return 359.6;
  else if (OutName.find("WJetsToLNu_HT-400to600") != std::string::npos) return 48.85;
  else if (OutName.find("WJetsToLNu_HT-600to800") != std::string::npos) return 12.05;
  else if (OutName.find("WJetsToLNu_HT-800to1200") != std::string::npos) return 5.501;
  else if (OutName.find("WJetsToLNu_HT-1200to2500") != std::string::npos) return 1.329;
  else if (OutName.find("WJetsToLNu_HT-2500toInf") != std::string::npos) return 0.03216;
    
  else if (OutName.find("DYJetsToLL") != std::string::npos) {
    //    std::cout << "DY" << std::endl;
    return 4895; // As we have large cut at Skim, this one is not needed
  }
  else if (OutName.find("DY1JetsToLL") != std::string::npos) return 1012.5;
  else if (OutName.find("DY2JetsToLL") != std::string::npos) return 332.8;
  else if (OutName.find("DY3JetsToLL") != std::string::npos) return 101.8;
  else if (OutName.find("DY4JetsToLL") != std::string::npos) return 54.8;
  
  else if (OutName.find("DYJetsToLL_M-50_HT-70to100") != std::string::npos) return 0;
  else if (OutName.find("DYJetsToLL_M-50_HT-100to200") != std::string::npos) return 148;
  else if (OutName.find("DYJetsToLL_M-50_HT-200to400") != std::string::npos) return 40.94;
  else if (OutName.find("DYJetsToLL_M-50_HT-400to600") != std::string::npos) return 5.497;
  else if (OutName.find("DYJetsToLL_M-50_HT-600to800") != std::string::npos) return 1.354;
  else if (OutName.find("DYJetsToLL_M-50_HT-800to1200") != std::string::npos) return 0.625;
  else if (OutName.find("DYJetsToLL_M-50_HT-1200to2500") != std::string::npos) return 0.151;
  else if (OutName.find("DYJetsToLL_M-50_HT-2500toInf") != std::string::npos) return 0.003647;
    
  //Di-boson
  else if (OutName.find("WW") != std::string::npos) return 115.0;
  else if (OutName.find("WZ") != std::string::npos) return 47.13;
  else if (OutName.find("ZZ") != std::string::npos) return 16.523;
  else if (OutName.find("zzTo4L") != std::string::npos) return 1.2;
  else if (OutName.find("HH") != std::string::npos) return 0.0334;
    
  //SingleTop
  else if (OutName.find("ST_t-channel_antitop") != std::string::npos) return 80.95;
  else if (OutName.find("ST_t-channel_top") != std::string::npos) return 136.02;
  else if (OutName.find("ST_tW_antitop_5f") != std::string::npos) return 35.6;
  else if (OutName.find("ST_tW_top_5f") != std::string::npos) return 35.6;
    
  else if (OutName.find("TT") != std::string::npos) return (831.76);

  else if (OutName.find("Codex") != std::string::npos ) return      1.0;
    
  else if (OutName.find("QCD_Pt-20toInf_MuEnrichedPt15") != std::string::npos) return     720648000  * 0.00042 ;
  else if (OutName.find("QCD") != std::string::npos) {
    //    std::cout << "QCD" << std::endl;
    return 720648000 * 0.00042 ;
  }
  else {
    std::cout<<"\n\n*********\nNot Listed in XSection menu !!!! Watch cout    "<<OutName<< "\n\n*********\n";
    return 1;
  }
}

int  main(int argc, char** argv) {

  std::vector<std::string> input;
  for (int f = 1; f < argc; f++) {
    input.push_back(*(argv + f));
    std::cout << "\n INPUT NAME IS:   " << input[f - 1] << "\n";
  }

  bool Muon = (input[0]=="Muon") ? true : false;
  std::vector<double> W_events  = W_EvenetMultiplicity();
  std::vector<double> DY_events = DY_EvenetMultiplicity();
  std::vector<double> puwt;
  PUWeight(puwt);

  for (unsigned int k=1; k<input.size(); k++ ) {
    ControlPlot c1(input[k],Muon);
    c1.Loop(puwt, W_events, DY_events);
    c1.SavePlot();
  }
}

ControlPlot::ControlPlot(std::string& fname, bool muon) : fname_(fname), 
							  f_Double_(0),fout_(0),
							  HistoTot_(0),h_mass_(0),
							  Muon_(muon),njets(0) {
  if      (fname_ == "DY1JetsToLL") njets =1;
  else if (fname_ == "DY2JetsToLL") njets =2;
  else if (fname_ == "DY3JetsToLL") njets =3;
  else if (fname_ == "DY4JetsToLL") njets =4;
  else if (fname_ == "W1JetsToLNu") njets =1;
  else if (fname_ == "W2JetsToLNu") njets =2;
  else if (fname_ == "W3JetsToLNu") njets =3;
  else if (fname_ == "W4JetsToLNu") njets =4;
  InitTree();
  if(fname_ == "DY1JetsToLL" || fname_ == "WJetsToLNu") njets = num_gen_jets;
  cout << "njets=  " << njets << endl;
}

void ControlPlot::InitTree() {

  f_Double_ = TFile::Open(fname_.c_str());
  std::cout << "\n  Now is running on ------->   " 
	    << std::string(f_Double_->GetName()) << "\n";

  HistoTot_ = (TH1F*) f_Double_->Get("hcount");
  Run_Tree_ = (TTree*) f_Double_->Get("EventTree");    

  /////////////////////////   General Info
  Run_Tree_->SetBranchAddress("isData", &isData);
  Run_Tree_->SetBranchAddress("run", &run);
  Run_Tree_->SetBranchAddress("lumis", &lumis);
  Run_Tree_->SetBranchAddress("event", &event);
  Run_Tree_->SetBranchAddress("genWeight",&genWeight);
  Run_Tree_->SetBranchAddress("HLTEleMuX", &HLTEleMuX);
  Run_Tree_->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled);
  Run_Tree_->SetBranchAddress("puTrue", &puTrue);        
    
  /////////////////////////   Tau Info
  Run_Tree_->SetBranchAddress("nTau", &nTau);
  Run_Tree_->SetBranchAddress("tauPt"  ,&tauPt);
  Run_Tree_->SetBranchAddress("tauEta"  ,&tauEta);
  Run_Tree_->SetBranchAddress("tauPhi"  ,&tauPhi);
  Run_Tree_->SetBranchAddress("tauMass"  ,&tauMass);
  Run_Tree_->SetBranchAddress("tauCharge"  ,&tauCharge);
  Run_Tree_->SetBranchAddress("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding);
  Run_Tree_->SetBranchAddress("tauByTightMuonRejection3", &tauByTightMuonRejection3);
  Run_Tree_->SetBranchAddress("tauByLooseMuonRejection3", &tauByLooseMuonRejection3);
  Run_Tree_->SetBranchAddress("tauByMVA6TightElectronRejection"  ,&tauByMVA6TightElectronRejection);
  Run_Tree_->SetBranchAddress("tauByLooseCombinedIsolationDeltaBetaCorr3Hits",&tauByLooseCombinedIsolationDeltaBetaCorr3Hits);
  Run_Tree_->SetBranchAddress("tauByMediumCombinedIsolationDeltaBetaCorr3Hits",&tauByMediumCombinedIsolationDeltaBetaCorr3Hits);
  Run_Tree_->SetBranchAddress("tauByMVA6LooseElectronRejection", &tauByMVA6LooseElectronRejection);
  Run_Tree_->SetBranchAddress("tauDxy",&tauDxy);
  Run_Tree_->SetBranchAddress("tauByMediumIsolationMVArun2v1DBoldDMwLT", &tauByMediumIsolationMVArun2v1DBoldDMwLT);
  Run_Tree_->SetBranchAddress("tauByLooseIsolationMVArun2v1DBoldDMwLT", &tauByLooseIsolationMVArun2v1DBoldDMwLT);
  /////////////////////////   Mu Info
  Run_Tree_->SetBranchAddress("nMu", &nMu);
  Run_Tree_->SetBranchAddress("eleEn", &eleEn);
  Run_Tree_->SetBranchAddress("muEn", &muEn);
  Run_Tree_->SetBranchAddress("tauEnergy", &tauEnergy);
  Run_Tree_->SetBranchAddress("muPt"  ,&muPt);
  Run_Tree_->SetBranchAddress("muEta"  ,&muEta);
  Run_Tree_->SetBranchAddress("muPhi"  ,&muPhi);
  Run_Tree_->SetBranchAddress("muIsoTrk", &muIsoTrk);
  Run_Tree_->SetBranchAddress("muCharge",&muCharge);
  //Run_Tree_->SetBranchAddress("muIsMediumID",&muIsMediumID);
  //Run_Tree_->SetBranchAddress("muIsLooseID",&muIsLooseID);
  Run_Tree_->SetBranchAddress("muPFChIso", &muPFChIso);
  Run_Tree_->SetBranchAddress("muPFPhoIso", &muPFPhoIso);
  Run_Tree_->SetBranchAddress("muPFNeuIso", &muPFNeuIso);
  Run_Tree_->SetBranchAddress("muPFPUIso", &muPFPUIso);
  Run_Tree_->SetBranchAddress("muD0",&muD0);
  Run_Tree_->SetBranchAddress("muDz",&muDz);
  Run_Tree_->SetBranchAddress("muIDbit", &muIDbit);
  
  /////////////////////////   Ele Info
  Run_Tree_->SetBranchAddress("nEle", &nEle);
  Run_Tree_->SetBranchAddress("elePt"  ,&elePt);
  Run_Tree_->SetBranchAddress("eleEta"  ,&eleEta);
  Run_Tree_->SetBranchAddress("elePhi"  ,&elePhi);
  Run_Tree_->SetBranchAddress("elePFChIso", &elePFChIso);
  Run_Tree_->SetBranchAddress("eleIDMVA", &eleIDMVA);
  Run_Tree_->SetBranchAddress("eleCharge",&eleCharge);
  Run_Tree_->SetBranchAddress("eleSCEta",&eleSCEta);
  Run_Tree_->SetBranchAddress("elePFChIso", &elePFChIso);
  Run_Tree_->SetBranchAddress("elePFPhoIso", &elePFPhoIso);
  Run_Tree_->SetBranchAddress("elePFNeuIso", &elePFNeuIso);
  Run_Tree_->SetBranchAddress("elePFPUIso", &elePFPUIso);
  Run_Tree_->SetBranchAddress("eleIDMVANonTrg", &eleIDMVANonTrg);
  Run_Tree_->SetBranchAddress("eleD0",&eleD0);
  Run_Tree_->SetBranchAddress("eleDz",&eleDz);
  Run_Tree_->SetBranchAddress("eleMissHits", &eleMissHits);
  Run_Tree_->SetBranchAddress("eleConvVeto", &eleConvVeto);
  /////////////////////////   Jet Info
  Run_Tree_->SetBranchAddress("nJet",&nJet);
  Run_Tree_->SetBranchAddress("jetPt",&jetPt);
  Run_Tree_->SetBranchAddress("jetEta",&jetEta);
  Run_Tree_->SetBranchAddress("jetPhi",&jetPhi);
  Run_Tree_->SetBranchAddress("jetEn",&jetEn);
  //Run_Tree_->SetBranchAddress("jetpfCombinedInclusiveSecondaryVertexV2BJetTags",&jetpfCombinedInclusiveSecondaryVertexV2BJetTags);
  Run_Tree_->SetBranchAddress("jetCSV2BJetTags",&jetCSV2BJetTags);
  Run_Tree_->SetBranchAddress("mcMomPID" ,&mcMomPID);
  Run_Tree_->SetBranchAddress("mcGMomPID" ,&mcGMomPID);
  Run_Tree_->SetBranchAddress("nMC" ,&nMC);
  Run_Tree_->SetBranchAddress("num_gen_jets", &num_gen_jets);
  Run_Tree_->SetBranchAddress("mcPID", &mcPID);
  Run_Tree_->SetBranchAddress("mcPt", &mcPt);
  Run_Tree_->SetBranchAddress("mcEta", &mcEta);
  Run_Tree_->SetBranchAddress("mcPhi", &mcPhi);
  Run_Tree_->SetBranchAddress("mcE", &mcE);
  Run_Tree_->SetBranchAddress("mcMomPt", &mcMomPt);
  Run_Tree_->SetBranchAddress("mcMomEta", &mcMomEta);
  Run_Tree_->SetBranchAddress("mcMomPhi", &mcMomPhi);
  Run_Tree_->SetBranchAddress("mcMomMass", &mcMomMass);
  Run_Tree_->SetBranchAddress("mcHadronPt", &mcHadronPt);
  Run_Tree_->SetBranchAddress("mcHadronEta", &mcHadronEta);
  Run_Tree_->SetBranchAddress("mcHadronPhi", &mcHadronPhi);
  Run_Tree_->SetBranchAddress("mcHadronMass", &mcHadronMass);
  Run_Tree_->SetBranchAddress("mcHadronE", &mcHadronE);
  Run_Tree_->SetBranchAddress("mcHadronGMomPID", &mcHadronGMomPID);
  Run_Tree_->SetBranchAddress("mcHadronMomPID", &mcHadronMomPID);
  Run_Tree_->SetBranchAddress("mcHadronMomPt", &mcHadronMomPt);
  Run_Tree_->SetBranchAddress("mcHadronMomMass", &mcHadronMomMass);
  Run_Tree_->SetBranchAddress("mcHadronMomEta", &mcHadronMomEta);
  Run_Tree_->SetBranchAddress("mcHadronMomPhi", &mcHadronMomPhi);
  Run_Tree_->SetBranchAddress("mcStatusFlag", &mcStatusFlag);
  Run_Tree_->SetBranchAddress("mcMass", &mcMass);
    
  /////////////////////////   MET Info
  Run_Tree_->SetBranchAddress("pfMET",&pfMET);
  Run_Tree_->SetBranchAddress("pfMETPhi",&pfMETPhi);

}

void ControlPlot::Loop(std::vector<double>& puwt, std::vector<double> W_events,
		       std::vector<double>& DY_events) {

  string Charge[2] = {"opposite_sign","same_sign"};
  string totcharge[2] = {"totCharge==0","totCharge!=0"};
  string taumucharge[2] = {"totCharge==0_and_mucharge==0","totCharge==0_and_mucharge_not_equal_to_zero_i.e_signal_region"};
  string tauCharge_[2] = {"opposite_sign_tau","same_sign_tau"};
  string muCharge_[2] = {"opposite_sign_muon","same_sign_muon"};

  std::string output = fname_;
  size_t root = output.find(".root");
  output.erase(root);
  output += "_1controlPlot.root";
  std::cout << "\n\n\n OUTPUT NAME IS:    " << output << std::endl;
  fout_    = TFile::Open(output.c_str(), "RECREATE");
  h_mass_  = new TH1F("mass","",20,70,110);
  myMap1_.clear();

  Int_t nentries_wtn = (Int_t) Run_Tree_->GetEntries();
  std::cout <<"nentries_wtn====" << nentries_wtn << "\n";

  for (Int_t i = 0; i < nentries_wtn; i++) {

    //    for ( Int_t i = 0; i < 10000; i++) { 
    cout << "i== " << i << endl;
    Run_Tree_->GetEntry(i);
    cout << "i== " << i << endl;
    if (i % 1000 == 0) 
      std::cout << "\r  Processed events: " << std::setw(8) << i << " of "
		<< std::setw(8) << nentries_wtn << std::endl;
    cout << "i== " << i << endl;
    bool PassTrigger = (Muon_) ? (HLTEleMuX >> 14 & 1) : (HLTEleMuX >> 34 & 1);
    if (!PassTrigger) continue;                                       
      
    int count_bjet(0);
    for (int ijet= 0 ; ijet < nJet ; ijet++) {
      if (jetPt->at(ijet) > 20 && fabs(jetEta->at(ijet)) < 2.5 && 
	  jetCSV2BJetTags->at(ijet) > 0.679) { count_bjet +=1; break; }
    }
    if (count_bjet > 0) continue; 
      
    float LumiWeight(1), GetGenWeight(1), PUWeight(1);
    if (!isData) {
      if (HistoTot_) 
	LumiWeight = weightCalc(HistoTot_, fname_, njets, W_events, DY_events);
      GetGenWeight = genWeight;
      unsigned int puNU = (unsigned int)(puTrue->at(0)*10);
      if (puNU < puwt.size()) PUWeight = puwt[puNU];
    }
     
    double weight = GetGenWeight*PUWeight*LumiWeight;

    ///////////////////////////////////////////////
    //Important Analysis Loop Will Happen Here!!!//
    ///////////////////////////////////////////////
      
    bool firstPart(true);
    std::vector<int> vec_muele, vec_tau;
    if (Muon_) {
      for  (int imu=0 ; imu < nMu; imu++){
	bool mupt = (firstPart) ? (muPt->at(imu) >18) : (muPt->at(imu) >9);
	UShort_t id = (muIDbit->at(imu) >> 1 & 1);
	float IsoMu=muPFChIso->at(imu)/muPt->at(imu);
	if ( (muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu) )  > 0.0)
	  IsoMu= ( muPFChIso->at(imu) + muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu))/muPt->at(imu);
	  
	if(mupt && id && fabs(muEta->at(imu)) <2.4 && IsoMu <0.3) {
	  firstPart = false;
	  vec_muele.push_back(imu);
	}
      }
    } else {
      for (int iele = 0; iele < nEle; ++iele) {
	bool elept = (firstPart) ? (elePt->at(iele) >18) : (elePt->at(iele) >15);
	float IsoEle=elePFChIso->at(iele)/elePt->at(iele);
	if ( (elePFNeuIso->at(iele) + elePFPhoIso->at(iele) - 0.5* elePFPUIso->at(iele))  > 0.0)
	  IsoEle= (elePFChIso->at(iele) + elePFNeuIso->at(iele) + elePFPhoIso->at(iele) - 0.5* elePFPUIso->at(iele))/elePt->at(iele);
	
	bool eleMVAId= false;
	if (fabs (eleSCEta->at(iele)) < 0.8 && eleIDMVA->at(iele) > 0.967083) eleMVAId= true;
	else if (fabs (eleSCEta->at(iele)) >  0.8 &&fabs (eleSCEta->at(iele)) <  1.5 && eleIDMVA->at(iele) > 0.929117) eleMVAId= true;
	else if ( fabs (eleSCEta->at(iele)) >  1.5 && eleIDMVA->at(iele) > 0.726311 ) eleMVAId= true;
	else eleMVAId= false;

	if(elept && eleMVAId && fabs(eleEta->at(iele)) <2.5 && IsoEle <0.1) {
	  firstPart = false;
	  vec_muele.push_back(iele);
	}
      }
    }
       
    for  (int itau=0 ; itau < nTau; itau++) {
      if (tauPt->at(itau) > 20 && fabs(tauEta->at(itau)) < 2.3 && 
	  tauByMVA6LooseElectronRejection->at(itau) !=0 && 
	  tauByTightMuonRejection3->at(itau) !=0 && 
	  taupfTausDiscriminationByDecayModeFinding->at(itau) !=0 && 
	  tauByLooseIsolationMVArun2v1DBoldDMwLT->at(itau) !=0) 
	vec_tau.push_back(itau);
    }
      
    std::string part;
    double pt[2], eta[2], phi[2], ene[2];
    int    chrg[2];
    if (Muon_) {
      part = "#mu_";
      for (int k=0; k<2; ++k) {
	pt[k]  = muPt->at(vec_muele[k]);
	eta[k] = muEta->at(vec_muele[k]);
	phi[k] = muPhi->at(vec_muele[k]);
	ene[k] = muEn->at(vec_muele[k]);
	chrg[k]= muCharge->at(vec_muele[k]);
      }
    } else {
      part = "#ele_";
      for (int k=0; k<2; ++k) {
	pt[k]  = elePt->at(vec_muele[k]);
	eta[k] = eleEta->at(vec_muele[k]);
	phi[k] = elePhi->at(vec_muele[k]);
	ene[k] = eleEn->at(vec_muele[k]);
	chrg[k]= eleCharge->at(vec_muele[k]);
      }
    }
  
    bool charge   =  (chrg[0]*chrg[1] > 0);
    bool samesign = charge ? true : false;
    bool oppositesign = !(samesign);
    bool mucharge[2] = {oppositesign, samesign};
    bool zerojet = vec_tau.size() ==0;
    bool onejet = vec_tau.size() ==1;
    bool twojet = vec_tau.size() ==2;
    bool jet[3] = {zerojet, onejet, twojet};
    std::string Jet[3] = {"zero_tau_jet_region", "one_tau_jet_region", "two_tau_jet_region"};

    TLorentzVector m1,m2,M;
    m1.SetPtEtaPhiE(pt[0],eta[0],phi[0],ene[0]);
    m2.SetPtEtaPhiE(pt[1],eta[1],phi[1],ene[1]);
    M = m1+m2;
    h_mass_->Fill(M.M(),weight);
    
    for (int i=0; i<2; i++) {
      if (!mucharge[i]) continue;
      std::string title = Charge[i]+"_and_no_cut_on_#tau_leg";
      if (M.M() >=20 && M.M() <=200) {
	plotFill("InvariantMass_of_"+part+"pair_with_" +title+"_with_mass_20-200GeV",M.M(),90,20,200,weight);
	plotFill("pt_distribution_of_Leading_"+part+title+"_with_mass_20-200GeV",pt[0],50,0,100,weight);
	plotFill("pt_distribution_of_SubLeading_"+part+title+"_with_mass_20-200GeV",pt[1],50,0,100,weight);
	plotFill("eta_distribution_of_Leading_"+part+title+"_with_mass_20-200GeV",eta[0],20,-2.4,2.4,weight);
	plotFill("eta_distribution_of_SubLeading_"+part+title+"_with_mass_20-200GeV",eta[1],20,-2.4,2.4,weight);
      } 
      
      if (M.M() >=60 && M.M() <=120) {
	plotFill("InvariantMass_of_"+part+"pair_with_" +title+"_with_mass_60-120GeV",M.M(),30,60,120,weight);
	plotFill("pt_distribution_of_Leading_"+part+title+"_with_mass_60-120GeV",pt[0],50,0,100,weight);
	plotFill("pt_distribution_of_SubLeading_"+part+title+"_with_mass_60-120GeV",pt[1],50,0,100,weight);
	plotFill("eta_distribution_of_Leading_"+part+title+"_with_mass_60-120GeV",eta[0],20,-2.4,2.4,weight);
	plotFill("eta_distribution_of_SubLeading_"+part+title+"_with_mass_60-120GeV",eta[1],20,-2.4,2.4,weight);
      }
      for (int j=0; j<3; j++ ){
	if(!jet[j]) continue;
	std::string title = Charge[i]+"_in_"+Jet[j];
	plotFill("InvariantMass_of_"+part+"pair_with_" +title,M.M(),20,70,110,weight);
      }
    }
    if (vec_tau.size() >1) {
      TLorentzVector t1,t2,T,HH;
      t1.SetPtEtaPhiE(tauPt->at(vec_tau[0]),tauEta->at(vec_tau[0]),tauPhi->at(vec_tau[0]),tauEnergy->at(vec_tau[0]));
      t2.SetPtEtaPhiE(tauPt->at(vec_tau[1]),tauEta->at(vec_tau[1]),tauPhi->at(vec_tau[1]),tauEnergy->at(vec_tau[1]));
      T  = t1+t2;
      HH = T+M;
      charge = tauCharge->at(vec_tau[0])*tauCharge->at(vec_tau[1]) >0;
      samesign = charge ? true : false;
      oppositesign = charge ? false : true;
      bool taucharge[2] = {oppositesign, samesign};
      for (int i=0; i<2; i++) {
	if (!mucharge[i]) continue;
	for (int j=0; j<2; j++) {
	  if (!taucharge[j]) continue;
	  std::string title = muCharge_[i]+"_"+tauCharge_[j]+"_region";
	  plotFill("InvariantMass_of_tau_pair_with_"+title,T.M(),6,50,110,weight);
	  plotFill("pt_distribution_of_Leading_#tau_with_"+title,tauPt->at(vec_tau[0]),75,0,150,weight);
	  plotFill("pt_distribution_of_SubLeading_#tau_with_"+title,tauPt->at(vec_tau[1]),75,0,150,weight);
	  plotFill("eta_distribution_of_Leading_#tau_with_"+title,tauEta->at(vec_tau[0]),20,-2.4,2.4,weight);
	  plotFill("eta_distribution_of_SubLeading_#tau_with_"+title,tauEta->at(vec_tau[1]),20,-2.4,2.4,weight);
	}
      }
	
      bool taumuCharge = ((chrg[0]+chrg[1]+tauCharge->at(vec_tau[0])+tauCharge->at(vec_tau[1])) == 0);
      bool taumucharge_[2], mucharge_[2];
      taumucharge_[0] = taumuCharge ? true : false;
      taumucharge_[1] = taumuCharge ? false : true;
      for (int i=0; i<2; i++) {
	if (!taumucharge_[i]) continue;
	plotFill("InvariantMass_of_4_particle_with_"+totcharge[i],HH.M(),7,250,900,weight);
      }
      if (taumucharge_[0]) {
	bool muCharge_ = ((chrg[0]+chrg[1]) == 0);
	mucharge_[0] = muCharge_ ? true : false;
	mucharge_[1] = muCharge_ ? false : true;
	for (int i=0; i<2; i++) {
	  if (!mucharge_[i]) continue;
	  plotFill("InvariantMass_of_4_particle_with_"+taumucharge[i],HH.M(),7,250,900,weight);
	}
      }
      for (int i=0; i<2; i++) {
	for (int j=0; j<2; j++) {
	  bool taumucharge = ((tauCharge->at(vec_tau[i]) + chrg[j]) == 0);
	  std::string title = (taumucharge) ? "opposite_sign" : "same_sign";
	  TLorentzVector t, m, H;
	  t.SetPtEtaPhiE(tauPt->at(vec_tau[i]),tauEta->at(vec_tau[i]),tauPhi->at(vec_tau[i]),tauEnergy->at(vec_tau[i]));
	  m.SetPtEtaPhiE(pt[j],eta[j],phi[j],ene[j]);
	  H = t+m;
	  plotFill("InvariantMass_of_taumu_with_"+title,H.M(),16,40,200,weight);
	}
      }
    }
  }
  
  //end of analysis code
}

void ControlPlot::SavePlot() {
  // close and write histograms/file
  if (fout_) {
    fout_->cd();
    if (h_mass_) h_mass_->Write();
    std::map<string, TH1F*>::const_iterator iMap1 = myMap1_.begin();
    for (; iMap1 != myMap1_.end(); ++iMap1) {
      TH1F* hist = nplot1(iMap1->first);
      if (hist != 0) hist->Write();
    }
    fout_->Close();
  }
  if (f_Double_) f_Double_->Close();
}

void ControlPlot::plotFill(std::string name, float x, int nx, float nxmin, 
			   float nxmax, double weight) {
  if (myMap1_.find(name) == myMap1_.end()) {
    myMap1_[name] = new TH1F(name.c_str(), name.c_str(), nx, nxmin, nxmax);
    myMap1_[name]->SetDefaultSumw2();
  }
  myMap1_[name]->Fill(x,weight);
}

TH1F* ControlPlot::nplot1(std::string name) {
  if (myMap1_.find(name) != myMap1_.end())
    return myMap1_[name];
  else
    return 0;
}

float ControlPlot::weightCalc(TH1F *Histo, std::string& outputName, int njet, 
			      std::vector<double>& W_events, 
			      std::vector<double>& DY_events) {
    
    
    
  //    std::cout<<"--->  Check Name is "<<newOut<<"\n";
    
  float LOtoNLO_DY = 1.177814096;
  float LOtoNLO_W = 1.189386467;
  //    float LOtoNLO_DY = 1.230888662;
  // float LOtoNLO_DY = 1; // Now we boson have pt dependent SF
  //    float LOtoNLO_W = 1.213783784;
  //float LOtoNLO_W = 1;  // Now we boson have pt dependent SF
  //    float luminosity=2154;
  //    float luminosity=    3990;
  //    float luminosity=    6260;
  //    float luminosity=    9235;
  //    float luminosity=    12900;
  //  5.761+2.573+4.248+4.009+3.102+7.540+8.391+0.215
  //78173197+33279413+26691984+27025933+20178544+44581284+46809967+1218674=277958996
  float luminosity=    35867;
  //float luminosity=    35839;
    
    
  size_t isDoubleMu = outputName.find("DoubleMuon");
  if (isDoubleMu != std::string::npos) return 1;
  else if(outputName.find("JetsToLL") != std::string::npos) {
    if (njet == 0) return luminosity*LOtoNLO_DY / (DY_events[0] / XSection("DYJetsToLL"));
    else if (njet == 1) return luminosity*LOtoNLO_DY / (DY_events[1] / XSection("DY1JetsToLL") + DY_events[0] / XSection("DYJetsToLL"));
    else if (njet == 2) return luminosity*LOtoNLO_DY / (DY_events[2] / XSection("DY2JetsToLL") + DY_events[0] / XSection("DYJetsToLL"));
    else if (njet == 3) return luminosity*LOtoNLO_DY / (DY_events[3] / XSection("DY3JetsToLL") + DY_events[0] / XSection("DYJetsToLL"));
    else if (njet == 4) return luminosity*LOtoNLO_DY / (DY_events[4] / XSection("DY4JetsToLL") + DY_events[0] / XSection("DYJetsToLL"));
    else {
      std::cout<<"**********   wooow  ********* There is a problem here\n";
      return 0;
    }
  }

  else if (outputName.find("JetsToLNu") != std::string::npos) {
    if (njet == 0) return luminosity*LOtoNLO_W / (W_events[0] / XSection("WJetsToLNu"));
    else if (njet == 1) return luminosity*LOtoNLO_W / (W_events[1] / XSection("W1JetsToLNu") + W_events[0] / XSection("WJetsToLNu"));
    else if (njet == 2) return luminosity*LOtoNLO_W / (W_events[2] / XSection("W2JetsToLNu") + W_events[0] / XSection("WJetsToLNu"));
    else if (njet == 3) return luminosity*LOtoNLO_W / (W_events[3] / XSection("W3JetsToLNu") + W_events[0] / XSection("WJetsToLNu"));
    else if (njet == 4) return luminosity*LOtoNLO_W / (W_events[4] / XSection("W4JetsToLNu") + W_events[0] / XSection("WJetsToLNu"));
    else {
      std::cout<<"**********   wooow  ********* There is a problem here\n";
      return 0;
    }
  }
  else
    return luminosity * XSection(outputName)*1.0 / Histo->GetBinContent(2);
    
}
