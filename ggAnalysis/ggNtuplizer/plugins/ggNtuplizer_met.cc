#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t metFilters_;
float genMET_;
float genMETPhi_;
float pfMET_;
float pfMETPhi_;
float pfMETsumEt_;
float pfMETmEtSig_;
float pfMETSig_;
float pfMET_T1JERUp_;
float pfMET_T1JERDo_;
float pfMET_T1JESUp_;
float pfMET_T1JESDo_;
float pfMET_T1MESUp_;
float pfMET_T1MESDo_;
float pfMET_T1EESUp_;
float pfMET_T1EESDo_;
float pfMET_T1PESUp_;
float pfMET_T1PESDo_;
float pfMET_T1TESUp_;
float pfMET_T1TESDo_;
float pfMET_T1UESUp_;
float pfMET_T1UESDo_;
float pfMET_T1TxyPhi_;
float pfMET_T1TxyPt_;
float pfMETPhi_T1JESUp_;
float pfMETPhi_T1JESDo_;
float pfMETPhi_T1UESUp_;
float pfMETPhi_T1UESDo_;
double signif_dxx_;
double signif_dyy_;
double signif_dyx_;
double signif_dxy_;

void ggNtuplizer::branchesMET(TTree* tree) {

  if (doGenParticles_) {
    tree->Branch("genMET",      &genMET_);
    tree->Branch("genMETPhi",   &genMETPhi_);
  }
  tree->Branch("metFilters",       &metFilters_);
  tree->Branch("pfMET",            &pfMET_);
  tree->Branch("pfMETPhi",         &pfMETPhi_);
  tree->Branch("pfMETsumEt",       &pfMETsumEt_);
  tree->Branch("pfMETmEtSig",      &pfMETmEtSig_);
  tree->Branch("pfMETSig",         &pfMETSig_);
  tree->Branch("pfMET_T1JERUp",    &pfMET_T1JERUp_);
  tree->Branch("pfMET_T1JERDo",    &pfMET_T1JERDo_);
  tree->Branch("pfMET_T1JESUp",    &pfMET_T1JESUp_);
  tree->Branch("pfMET_T1JESDo",    &pfMET_T1JESDo_);
  /*
  tree->Branch("pfMET_T1MESUp",    &pfMET_T1MESUp_);
  tree->Branch("pfMET_T1MESDo",    &pfMET_T1MESDo_);
  tree->Branch("pfMET_T1EESUp",    &pfMET_T1EESUp_);
  tree->Branch("pfMET_T1EESDo",    &pfMET_T1EESDo_);
  tree->Branch("pfMET_T1PESUp",    &pfMET_T1PESUp_);
  tree->Branch("pfMET_T1PESDo",    &pfMET_T1PESDo_);
  tree->Branch("pfMET_T1TESUp",    &pfMET_T1TESUp_);
  tree->Branch("pfMET_T1TESDo",    &pfMET_T1TESDo_);
  */
  tree->Branch("pfMET_T1UESUp",    &pfMET_T1UESUp_);
  tree->Branch("pfMET_T1UESDo",    &pfMET_T1UESDo_);
  tree->Branch("pfMETPhi_T1JESUp", &pfMETPhi_T1JESUp_);
  tree->Branch("pfMETPhi_T1JESDo", &pfMETPhi_T1JESDo_);
  tree->Branch("pfMETPhi_T1UESUp", &pfMETPhi_T1UESUp_);
  tree->Branch("pfMETPhi_T1UESDo", &pfMETPhi_T1UESDo_);
  tree->Branch("signif_dxx", &signif_dxx_);
  tree->Branch("signif_dyy", &signif_dyy_);
  tree->Branch("signif_dxy", &signif_dxy_);
  tree->Branch("signif_dyx", &signif_dyx_);
}

void ggNtuplizer::fillMET(const edm::Event& e, const edm::EventSetup& es) {

  metFilters_ = 0;

  if (addFilterInfoMINIAOD_) {
    string filterNamesToCheck[9] = {
      "Flag_HBHENoiseFilter",
      "Flag_HBHENoiseIsoFilter", 
      "Flag_globalTightHalo2016Filter",
      "Flag_goodVertices",
      "Flag_eeBadScFilter",
      "Flag_EcalDeadCellTriggerPrimitiveFilter",
      "Flag_badMuons",
      "Flag_duplicateMuons",
      "Flag_noBadMuons"
    };

    edm::Handle<edm::TriggerResults> patFilterResultsHandle;
    e.getByToken(patTrgResultsLabel_, patFilterResultsHandle);
    edm::TriggerResults const& patFilterResults = *patFilterResultsHandle;
    
    auto&& filterNames = e.triggerNames(patFilterResults);

    // === the following lines allow us to find the filters stored in the event ! ===
    /*
    edm::TriggerNames const& triggerNames = e.triggerNames(patFilterResults);
    for ( edm::TriggerNames::Strings::const_iterator triggerName = triggerNames.triggerNames().begin();
	  triggerName != triggerNames.triggerNames().end(); ++triggerName ) {
      int triggerId = triggerNames.triggerIndex(*triggerName);
      if ( triggerId >= 0 && triggerId < (int)triggerNames.size() ) {
	std::string triggerDecision = ( patFilterResultsHandle->accept(triggerId) ) ? "passed" : "failed";
	  
	std::cout << " triggerName = " << (*triggerName) << " " << triggerDecision << std::endl;
      }
    }
    */

    for (unsigned iF = 0; iF < 9; ++iF) {
      unsigned index = filterNames.triggerIndex(filterNamesToCheck[iF]);
      if ( index == filterNames.size() ) 
	LogDebug("METFilters") << filterNamesToCheck[iF] << " is missing, exiting";
      else {
	if ( !patFilterResults.accept(index) ) {
	  if (iF <=5) metFilters_ += pow(2, iF+1);
	  else metFilters_ += pow(2, iF+3);
	}
      }
    }
  }
  
  if (addFilterInfoMINIAOD_) {
    edm::Handle<bool> ifilterbadChCand;
    e.getByToken(BadChCandFilterToken_, ifilterbadChCand);
    bool filterbadChCandidate_ = *ifilterbadChCand;
    
    edm::Handle<bool> ifilterbadPFMuon;
    e.getByToken(BadPFMuonFilterToken_, ifilterbadPFMuon);
    bool filterbadPFMuon_ = *ifilterbadPFMuon;
    
    if ( !filterbadPFMuon_      ) metFilters_ += pow(2, 7);
    if ( !filterbadChCandidate_ ) metFilters_ += pow(2, 8);
  }
  
  edm::Handle<edm::View<pat::MET> > pfMETHandle;
  e.getByToken(pfMETlabel_, pfMETHandle);

  genMET_      = -99;
  genMETPhi_   = -99;
  pfMET_       = -99;
  pfMETPhi_    = -99;
  pfMETsumEt_  = -99;
  pfMETmEtSig_ = -99;
  pfMETSig_    = -99;
  signif_dxx_  = -99;
  signif_dyy_  = -99;
  signif_dxy_  = -99;
  signif_dyx_  = -99;

  if (pfMETHandle.isValid()) {
    const pat::MET *pfMET = 0;
    pfMET = &(pfMETHandle->front());

    pfMET_       = pfMET->et();
    pfMETPhi_    = pfMET->phi();
    pfMETsumEt_  = pfMET->sumEt();
    pfMETmEtSig_ = (pfMET->mEtSig() < 1.e10) ? pfMET->mEtSig() : 0;
    pfMETSig_    = (pfMET->significance() < 1.e10) ? pfMET->significance() : 0;;
    signif_dxx_  =  pfMET->getSignificanceMatrix()(0,0);
    signif_dxy_  =  pfMET->getSignificanceMatrix()(0,1);
    signif_dyy_  =  pfMET->getSignificanceMatrix()(1,1);
    signif_dyx_  =  pfMET->getSignificanceMatrix()(1,0);

    if (!isAOD_) {
      // Type1MET uncertainties =======================================
      pfMET_T1JERUp_ = pfMET->shiftedPt(pat::MET::JetResUp);
      pfMET_T1JERDo_ = pfMET->shiftedPt(pat::MET::JetResDown);
      pfMET_T1JESUp_ = pfMET->shiftedPt(pat::MET::JetEnUp);
      pfMET_T1JESDo_ = pfMET->shiftedPt(pat::MET::JetEnDown);
      /*
      pfMET_T1MESUp_ = pfMET->shiftedPt(pat::MET::MuonEnUp);
      pfMET_T1MESDo_ = pfMET->shiftedPt(pat::MET::MuonEnDown);
      pfMET_T1EESUp_ = pfMET->shiftedPt(pat::MET::ElectronEnUp);
      pfMET_T1EESDo_ = pfMET->shiftedPt(pat::MET::ElectronEnDown);
      pfMET_T1PESUp_ = pfMET->shiftedPt(pat::MET::PhotonEnUp);
      pfMET_T1PESDo_ = pfMET->shiftedPt(pat::MET::PhotonEnDown);
      pfMET_T1TESUp_ = pfMET->shiftedPt(pat::MET::TauEnUp);
      pfMET_T1TESDo_ = pfMET->shiftedPt(pat::MET::TauEnDown);
      */
      pfMET_T1UESUp_ = pfMET->shiftedPt(pat::MET::UnclusteredEnUp);
      pfMET_T1UESDo_ = pfMET->shiftedPt(pat::MET::UnclusteredEnDown);

      pfMETPhi_T1JESUp_ = pfMET->shiftedPhi(pat::MET::JetEnUp);
      pfMETPhi_T1JESDo_ = pfMET->shiftedPhi(pat::MET::JetEnDown);
      pfMETPhi_T1UESUp_ = pfMET->shiftedPhi(pat::MET::UnclusteredEnUp);
      pfMETPhi_T1UESDo_ = pfMET->shiftedPhi(pat::MET::UnclusteredEnDown);
    }

    if (!e.isRealData()) {
      genMET_    = pfMET->genMET()->et();
      genMETPhi_ = pfMET->genMET()->phi();
    }

  } 

}
