//________________________________________________________________
// 
//
//Date 05/10/2016: This is to read geant.root file
//                      -Nihar
//Date 10/26/2016: Added some histograms to compare to
//                 pythia.root file
//                      -Derek
//Date 10/28/2016: Added/updated methods to read in the
//                 trees contained in '*.pythia.root'
//                 files.
//                      -Derek 
//Date 05/10/2017: Added members/methods to specify
//                 pythia file
//                      -Derek
//Date 05/11/2017: Added additional histograms
//                      -Derek
//Date 05/22/2017: Replaced (most) histograms with NTuple.
//                      -Derek
//Date 05/23/2017: Added method to determine charge of pythia
//                 track.
//                      -Derek
//Date 02/02/2018: Converted NTuple to a TTree.
//                      -Derek
//____________________________________________________________


#include <cmath>
#include <cassert>
#include "StJetTreeMcMaker.h"
#include "StChain.h"
#include "StEvent.h"
#include "StEventTypes.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcUtil/filters/StEmcFilter.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TFile.h"
#include <vector>
#include "TVector.h"
#include "TVector3.h"
#include <math.h>
#include "TNtuple.h"
#include <TString.h>
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "StEventUtilities/StuRefMult.hh"
#include "StEventUtilities/StuFtpcRefMult.hh"
#include "StThreeVector.hh"
#include "StEvent/StEnumerations.h"
#include "StMcEvent.hh"
#include "StMcVertex.hh"
#include "StMcTrack.hh"
#include "StMcEventMaker/StMcEventMaker.h"
#include "StMcEventTypes.hh"
#include "StPrimaryVertex.h"
#include "StBTofHeader.h"
//#include "StMuEvent.h"
#include <TSQLServer.h>
#include <TSQLRow.h>
#include <TSQLResult.h>

#include "St_db_Maker/St_db_Maker.h"
#include "StDetectorDbMaker/StDetectorDbTriggerID.h"

#include <emcStatus.h>
#include "tables/St_emcStatus_Table.h"

#include "tables/St_emcPed_Table.h"
//#include "L2Result.h"
#include "StTriggerData.h"
#include "StTriggerData2007.h"
//#include "StDaqLib/TRG/trgStructures.h"
#include "StDaqLib/TRG/trgStructures2007.h"
//#include "StDaqLib/TRG/L2pedResults2006.h"
//#include "StDaqLib/TRG/L2gammaResult2007.h"
#include "StDaqLib/EMC/StEmcDecoder.h"
//#include "L2gammaResult.h"
#include "StEvent/StTriggerIdCollection.h"

#include "StEvent/StTriggerId.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
#include "StMuDSTMaker/COMMON/StMuTofHit.h"
#include "StMuDSTMaker/COMMON/StMuBTofPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StTriggerData.h"

 
#include "GammaJetTrack.h"
#include "GammaJetEvent.h"
#include "GammaJetTower.h"
#include "GammaJetTowerUtil.h"
#include <exception>
#include <vector>

void muEventInfo(StMuEvent&, const Int_t&);
ClassImp(StJetTreeMcMaker)

//------------------------------------------------------------------------------------
void StJetTreeMcMaker::setDbMaker(St_db_Maker* dbMaker)
{
  mDbMaker = dbMaker;
}
//-----------------------------------------------------------------------------------
StJetTreeMcMaker::StJetTreeMcMaker(const char *name,char *dataType):StMaker(name){ 
  

 pTrMatchArr = 0;
 pTriMatchEtaArr =0;
 pTrMatchPhiArr = 0;
 pTrIndexArray = 0;
 pTrEtaArray =0;
 pTrPhiArray =0;

  

//+++++++++++
 muDst = 0;   // NRS: intialization
 muEvent =0;


  outFile="____.root";
}

//_____________________________________________________________________________
StJetTreeMcMaker::~StJetTreeMcMaker(){
  //
}

//_____________________________________________________________________________
// Init - is a first method the top level StChain 
// calls to initialize all its makers 
Int_t StJetTreeMcMaker::Init(){
  // Create Histograms 

 
  mRandom = new TRandom();
  File_output = new TFile(outFile,"RECREATE");

  event = new GammaJetEvent();
  outTree = new TTree("Gfmtodst","Gfmtodst");
  outTree->SetAutoSave(-500000000);  // autosave activated for each 5 MB
  
  eventBranch = outTree->Branch("EventList",&event,1000000);
  
  EvValues = new float[81]; // for event info


  // [10.28.2016, Derek]
  TTree *tree;
  TFile *f = (TFile*) gROOT -> GetListOfFiles() -> FindObject(sPythia.Data());
  if (!f || !f->IsOpen()) {
    f = new TFile(sPythia.Data());
    f -> GetObject("PythiaTree", tree); 
  }

  // [10.28.2016, Derek]
  if (!tree) return -666;
  pyChain   = tree;
  NumEvts   = 0;
  pyCurrent = -1;
  pyChain -> SetMakeClass(1);
  pyChain -> SetBranchAddress("fUniqueID", &fUniqueID, &b_PythiaBranch_fUniqueID);
  pyChain -> SetBranchAddress("fBits", &fBits, &b_PythiaBranch_fBits);
  pyChain -> SetBranchAddress("mRunId", &mRunId, &b_PythiaBranch_mRunId);
  pyChain -> SetBranchAddress("mEventId", &mEventId, &b_PythiaBranch_mEventId);
  pyChain -> SetBranchAddress("mProcessId", &mProcessId, &b_PythiaBranch_mProcessId);
  pyChain -> SetBranchAddress("mTune", &mTune, &b_PythiaBranch_mTune);
  pyChain -> SetBranchAddress("mVertex.fUniqueID", &mVertex_fUniqueID, &b_PythiaBranch_mVertex_fUniqueID);
  pyChain -> SetBranchAddress("mVertex.fBits", &mVertex_fBits, &b_PythiaBranch_mVertex_fBits);
  pyChain -> SetBranchAddress("mVertex.fX", &mVertex_fX, &b_PythiaBranch_mVertex_fX);
  pyChain -> SetBranchAddress("mVertex.fY", &mVertex_fY, &b_PythiaBranch_mVertex_fY);
  pyChain -> SetBranchAddress("mVertex.fZ", &mVertex_fZ, &b_PythiaBranch_mVertex_fZ);
  pyChain -> SetBranchAddress("mS", &mS, &b_PythiaBranch_mS);
  pyChain -> SetBranchAddress("mT", &mT, &b_PythiaBranch_mT);
  pyChain -> SetBranchAddress("mU", &mU, &b_PythiaBranch_mU);
  pyChain -> SetBranchAddress("mPt", &mPt, &b_PythiaBranch_mPt);
  pyChain -> SetBranchAddress("mCosTheta", &mCosTheta, &b_PythiaBranch_mCosTheta);
  pyChain -> SetBranchAddress("mX1", &mX1, &b_PythiaBranch_mX1);
  pyChain -> SetBranchAddress("mX2", &mX2, &b_PythiaBranch_mX2);
  pyChain -> SetBranchAddress("mMstu72", &mMstu72, &b_PythiaBranch_mMstu72);
  pyChain -> SetBranchAddress("mMstu73", &mMstu73, &b_PythiaBranch_mMstu73);
  pyChain -> SetBranchAddress("mMstp111", &mMstp111, &b_PythiaBranch_mMstp111);
  pyChain -> SetBranchAddress("mPartonALL", &mPartonALL, &b_PythiaBranch_mPartonALL);
  pyChain -> SetBranchAddress("mDF1[35]", mDF1, &b_PythiaBranch_mDF1);
  pyChain -> SetBranchAddress("mDF2[35]", mDF2, &b_PythiaBranch_mDF2);
  pyChain -> SetBranchAddress("mF1[2]", mF1, &b_PythiaBranch_mF1);
  pyChain -> SetBranchAddress("mF2[2]", mF2, &b_PythiaBranch_mF2);
  pyChain -> SetBranchAddress("mDF1NNPDF[101]", mDF1NNPDF, &b_PythiaBranch_mDF1NNPDF);
  pyChain -> SetBranchAddress("mDF2NNPDF[101]", mDF2NNPDF, &b_PythiaBranch_mDF2NNPDF);
  pyChain -> SetBranchAddress("mF1NNPDF", &mF1NNPDF, &b_PythiaBranch_mF1NNPDF);
  pyChain -> SetBranchAddress("mF2NNPDF", &mF2NNPDF, &b_PythiaBranch_mF2NNPDF);
  pyChain -> SetBranchAddress("mParticles", &mParticles_, &b_PythiaBranch_mParticles_);
  pyChain -> SetBranchAddress("mParticles.fUniqueID", mParticles_fUniqueID, &b_mParticles_fUniqueID);
  pyChain -> SetBranchAddress("mParticles.fBits", mParticles_fBits, &b_mParticles_fBits);
  pyChain -> SetBranchAddress("mParticles.fLineColor", mParticles_fLineColor, &b_mParticles_fLineColor);
  pyChain -> SetBranchAddress("mParticles.fLineStyle", mParticles_fLineStyle, &b_mParticles_fLineStyle);
  pyChain -> SetBranchAddress("mParticles.fLineWidth", mParticles_fLineWidth, &b_mParticles_fLineWidth);
  pyChain -> SetBranchAddress("mParticles.fPdgCode", mParticles_fPdgCode, &b_mParticles_fPdgCode);
  pyChain -> SetBranchAddress("mParticles.fStatusCode", mParticles_fStatusCode, &b_mParticles_fStatusCode);
  pyChain -> SetBranchAddress("mParticles.fMother[2]", mParticles_fMother, &b_mParticles_fMother);
  pyChain -> SetBranchAddress("mParticles.fDaughter[2]", mParticles_fDaughter, &b_mParticles_fDaughter);
  pyChain -> SetBranchAddress("mParticles.fWeight", mParticles_fWeight, &b_mParticles_fWeight);
  pyChain -> SetBranchAddress("mParticles.fCalcMass", mParticles_fCalcMass, &b_mParticles_fCalcMass);
  pyChain -> SetBranchAddress("mParticles.fPx", mParticles_fPx, &b_mParticles_fPx);
  pyChain -> SetBranchAddress("mParticles.fPy", mParticles_fPy, &b_mParticles_fPy);
  pyChain -> SetBranchAddress("mParticles.fPz", mParticles_fPz, &b_mParticles_fPz);
  pyChain -> SetBranchAddress("mParticles.fE", mParticles_fE, &b_mParticles_fE);
  pyChain -> SetBranchAddress("mParticles.fVx", mParticles_fVx, &b_mParticles_fVx);
  pyChain -> SetBranchAddress("mParticles.fVy", mParticles_fVy, &b_mParticles_fVy);
  pyChain -> SetBranchAddress("mParticles.fVz", mParticles_fVz, &b_mParticles_fVz);
  pyChain -> SetBranchAddress("mParticles.fVt", mParticles_fVt, &b_mParticles_fVt);
  pyChain -> SetBranchAddress("mParticles.fPolarTheta", mParticles_fPolarTheta, &b_mParticles_fPolarTheta);
  pyChain -> SetBranchAddress("mParticles.fPolarPhi", mParticles_fPolarPhi, &b_mParticles_fPolarPhi);


  // [10.26.2016 / 10.28.2016 / 10.31.2016 / 05.11.2017, Derek]
  const int    nMult = 500;
  const double mult1 = 0.;
  const double mult2 = 500.;
  hGeantCounter     = new TH1D("hGeantCounter", "Total multiplicity (Geant)", nMult, mult1, mult2);
  hPythiaCounter    = new TH1D("hPythiaCounter", "Total multiplicity (Pythia)", nMult, mult1, mult2);
  hFinalCounter     = new TH1D("hFinalCounter", "FS multiplicity [daughter1=daughter2=0] (Pythia)", nMult, mult1, mult2);
  hFinalCountS1     = new TH1D("hFinalCountS1", "No. of FS tracks w/ status=1 (Pythia)", nMult, mult1, mult2);
  hFinalCountS2     = new TH1D("hFinalCountS2", "No. of FS tracks w/ status=2 (Pythia)", nMult, mult1, mult2);
  hFinalCountS3     = new TH1D("hFinalCountS3", "No. of FS tracks w/ status=3 (Pythia)", nMult, mult1, mult2);
  hNumOutsideAccept = new TH1D("hNumOutsideAccept", "No. of FS tracks [status=1] outside acceptance (Pythia)", nMult, mult1, mult2);
  hNumInsideAccept  = new TH1D("hNumInsideAccept", "No. of FS tracks [status=1] inside acceptance (Pythia)", nMult, mult1, mult2);
  hPyNumMid         = new TH1D("hPyNumMid", "No. of FS tracks [status=1] with |eta|<1 (Pythia)", nMult, mult1, mult2);
  hGeNumMid         = new TH1D("hGeNumMid", "No. of tracks with |eta|<1 (Geant)", nMult, mult1, mult2);
  hPyMultPos        = new TH1D("hPyMultPos", "No. of (+) FS tracks [status=1] inside acceptance (Pythia)", nMult, mult1, mult2);
  hPyMultNeg        = new TH1D("hPyMultNeg", "No. of (-) FS tracks [status=1] inside acceptace (Pythia)", nMult, mult1, mult2);
  hPyMultNeu        = new TH1D("hPyMultNeu", "No. of (0) FS tracks [status=1] inside acceptance (Pythia)", nMult, mult1, mult2);
  hGeMultPos        = new TH1D("hGeMultPos", "No. of (+) tracks in acceptance (Geant)", nMult, mult1, mult2);
  hGeMultNeg        = new TH1D("hGeMultNeg", "No. of (-) tracks in acceptance (Geant)", nMult, mult1, mult2);
  hGeMultNeu        = new TH1D("hGeMultNeu", "No. of (0) tracks in acceptance (Geant)", nMult, mult1, mult2);
  hGeantCounter     -> Sumw2();
  hPythiaCounter    -> Sumw2();
  hFinalCounter     -> Sumw2();
  hFinalCountS1     -> Sumw2();
  hFinalCountS2     -> Sumw2();
  hFinalCountS3     -> Sumw2();
  hNumOutsideAccept -> Sumw2();
  hNumInsideAccept  -> Sumw2();
  hPyNumMid         -> Sumw2();
  hGeNumMid         -> Sumw2();
  hPyMultPos        -> Sumw2();
  hPyMultNeg        -> Sumw2();
  hPyMultNeu        -> Sumw2();
  hGeMultPos        -> Sumw2();
  hGeMultNeg        -> Sumw2();
  hGeMultNeu        -> Sumw2();

  // [02.02.2018, Derek]
  File_output -> cd();
  _tPythia = new TTree("PyTree", "a tree for pythia");
  _tPythia -> Branch("EventId", &_PyEvtId, "EventId/I");
  _tPythia -> Branch("RunId", &_PyRunId, "RunId/I");
  _tPythia -> Branch("ProcessId", &_PyProId, "ProcessId/I");
  _tPythia -> Branch("NumParticle", &_PyNumPar, "NumParticle/I");
  _tPythia -> Branch("PartonicPt", &_PyPtPart, "PartonicPt/F");
  _tPythia -> Branch("MandelstamS", &_PyS, "MandelstamS/F");
  _tPythia -> Branch("MandelstamT", &_PyT, "MandelstamT/F");
  _tPythia -> Branch("MandelstamU", &_PyU, "MandelstamU/F");
  _tPythia -> Branch("VtxX", &_PyVx, "VtxX/D");
  _tPythia -> Branch("VtxY", &_PyVy, "VtxY/D");
  _tPythia -> Branch("VtxZ", &_PyVz, "VtxZ/D");
  _tPythia -> Branch("ParId", &_PyParId);
  _tPythia -> Branch("ParStatus", &_PyParStatus);
  _tPythia -> Branch("ParKid1", &_PyParKid1);
  _tPythia -> Branch("ParKid2", &_PyParKid2);
  _tPythia -> Branch("ParMom1", &_PyParMom1);
  _tPythia -> Branch("ParMom2", &_PyParMom2);
  _tPythia -> Branch("ParCharge", &_PyParChrg);
  _tPythia -> Branch("ParWeight", &_PyParWeight);
  _tPythia -> Branch("ParEta", &_PyParEta);
  _tPythia -> Branch("ParTheta", &_PyParTheta);
  _tPythia -> Branch("ParPx", &_PyParPx);
  _tPythia -> Branch("ParPy", &_PyParPy);
  _tPythia -> Branch("ParPz", &_PyParPz);
  _tPythia -> Branch("ParPt", &_PyParPt);
  _tPythia -> Branch("ParE", &_PyParE);
  _tPythia -> Branch("ParMass", &_PyParMass);
  _tPythia -> SetAutoSave(-500000000);
  return StMaker::Init();

}

//_____________________________________________________________________________
/// Make - this method is called in loop for each mStEvent
Int_t StJetTreeMcMaker::Make()
{
  //doMuEvent();  
  //__________
  doMCEvent();

  return kStOk;
}
//_______________________________________
Int_t StJetTreeMcMaker::doMuEvent()
{

  muDst =  (StMuDst*) GetInputDS("MuDst");
  if(!muDst){
    cout << "No MuDst ....." << endl;
    return kStOK;
  }
  
  muEvent = (StMuEvent*) muDst->event();

  cout<<" I am Here in MuDSt...."<<endl;
  Int_t refmult = muEvent->refMult();

  Int_t RunId   = muEvent->runId();
  Int_t eventId   = muEvent->eventId();
  
  //Int_t nPrimaryVertex = 0;
  //  nPrimaryVertex = muDst->numberOfPrimaryVertices();

  float Vz= 0., Vx=0., Vy=0., Vxy=0.;

  StThreeVectorF Vtx=muEvent->primaryVertexPosition();

  Vz = Vtx.z();
  Vx = Vtx.x();
  Vy = Vtx.y();

  Vxy = TMath::Sqrt( Vx*Vx + Vy*Vy );


  cout<<"MuEvent:  Refmult= "<<refmult<<"  RunId= "<<RunId<<" evntId= "<<eventId<<" V x y z: "<<Vx<<" "<<Vy<<"  "<<Vz<<endl;
  
  

  return kStOk;
}

///____________________________________________
Int_t StJetTreeMcMaker::doMCEvent()
{  


  mcEvent =  (StMcEvent*) GetDataSet("StMcEvent");       

  Int_t event_counter =0;

  event->ResetEvent();
  event_counter++;
  GammaJetTrack TrackInfo;  
  //___________________


  // [10.28.2016, Derek]
  Int_t    mcEvntNumber = mcEvent -> eventNumber();
  Int_t    mcRunNumber  = mcEvent->runNumber();
  Long64_t numBytes     = pyChain -> GetEntry(NumEvts);
  if (numBytes < 0) {
    cerr << "OH NO!!!!! Something weird at event " << NumEvts << "!"
         << endl;
    return -666;
  }
  NumEvts++;



  ///______________________________________
  const StPtrVecMcTrack& mcTracks = mcEvent->primaryVertex()->daughters();
  StMcTrackConstIterator mcTrkIter = mcTracks.begin();
  Int_t track_counter = 0;
  Int_t trkIndex_counter=0;

  // [05.22.2017 / 05.23.2017, Derek]
  Int_t nGeMid = 0;
  Int_t nGePos = 0;
  Int_t nGeNeg = 0;
  Int_t nGeNeu = 0;

  // geant track loop
  for ( ; mcTrkIter != mcTracks.end(); ++mcTrkIter) {
    StMcTrack* track = *mcTrkIter;

    Double_t mcTrk_px = track->momentum().x();
    Double_t mcTrk_py = track->momentum().y();
    Double_t mcTrk_pz = track->momentum().z();
    Double_t mcTrk_pt = sqrt(mcTrk_px*mcTrk_px + mcTrk_py*mcTrk_py);
    Double_t mcTrk_eng = track->energy();
    Double_t mcTrk_eta = track->pseudoRapidity();
    Int_t    mcTrk_pdgId = track->pdgId();
    Int_t    mcTrk_geantId = track->geantId();

    Int_t mcTrk_chrg = -9999;
    Int_t mcTrk_dEdX = -9999;
    if( track->particleDefinition()){
      mcTrk_chrg = track->particleDefinition()->charge();
      mcTrk_dEdX = track->particleDefinition()->pdgEncoding();
    }
    else{
      cout << "Particle with no encoding " << endl;
      mcTrk_chrg = 0.;
      mcTrk_dEdX = 0.;

    }
    Int_t mcTrk_NoOfFitHits = track->tpcHits().size();
    

    TrackInfo.Clear();
    TrackInfo.SetnHitsFit(0);
    TrackInfo.SetnHitsPoss(0);
    TrackInfo.SetTrackFlag(0);
    TrackInfo.SetPdgId(mcTrk_pdgId);
    TrackInfo.SetGeantId(mcTrk_geantId);
    TrackInfo.SetpZ(mcTrk_pz);
    TrackInfo.SetpY(mcTrk_py);
    TrackInfo.SetpX(mcTrk_px);
    TrackInfo.SetpT(mcTrk_pt);
    TrackInfo.SetdEdx(mcTrk_dEdX);
    TrackInfo.SetCharge(mcTrk_chrg);
    TrackInfo.SetTOFBeta(0);
    TrackInfo.SetEta(mcTrk_eta);
    TrackInfo.SetPhi(0);
    TrackInfo.SetnSigElectron(0);
    TrackInfo.SetnSigPion(0);
    TrackInfo.SetnSigKaon(0);
    TrackInfo.SetnSigProton(0);
    TrackInfo.SetDCAg(0);
    TrackInfo.SetnHits(0);
    TrackInfo.SetdEdxHits(0);
    TrackInfo.SetFirstZPoint(0);
    TrackInfo.SetLastZPoint(0);
    TrackInfo.SetTOFSigElectron(0);
    TrackInfo.SetTOFSigPion(0);
    TrackInfo.SetTOFSigKaon(0);
    TrackInfo.SetTOFSigProton(0);
    TrackInfo.SetPathLength(0);
    TrackInfo.SettimeOfflight(0);
    TrackInfo.SettrkIndex(trkIndex_counter);

    event->AddTrack(&TrackInfo,track_counter);
    track_counter++;

    // [05.22.2017 / 05.23.2017, Derek]
    if (abs(mcTrk_eta) < 1.)
      ++nGeMid;
    if (abs(mcTrk_eta) < 1.8) {
      if (mcTrk_chrg == 1.)
        ++nGePos;
      if (mcTrk_chrg == -1.)
        ++nGeNeg;
      if (mcTrk_chrg == 0.)
        ++nGeNeu;
    }

    Double_t mcTrk_p     = sqrt(mcTrk_pt*mcTrk_pt + mcTrk_pz*mcTrk_pz);
    Double_t mcTrk_theta = asin(mcTrk_pt / mcTrk_p);


  }  // end Geant track loop

  // [10.26.2016 / 05.22.2017, Derek]
  hGeantCounter -> Fill(track_counter);
  hGeNumMid     -> Fill(nGeMid);
  hGeMultPos    -> Fill(nGePos);
  hGeMultNeg    -> Fill(nGeNeg);
  hGeMultNeu    -> Fill(nGeNeu);


  // reset pythia members [02.02.2018, Derek]
  _PyEvtId  = -1;
  _PyRunId  = -1;
  _PyProId  = -1;
  _PyNumPar = -1;
  _PyPtPart = -1.;
  _PyS      = -1.;
  _PyT      = -1.;
  _PyU      = -1.;
  _PyVx     = -1000.;
  _PyVy     = -1000.;
  _PyVz     = -1000.;
  _PyParId.clear();
  _PyParStatus.clear();
  _PyParKid1.clear();
  _PyParKid2.clear();
  _PyParMom1.clear();
  _PyParMom2.clear();
  _PyParChrg.clear();
  _PyParWeight.clear();
  _PyParEta.clear();
  _PyParTheta.clear();
  _PyParPx.clear();
  _PyParPy.clear();
  _PyParPz.clear();
  _PyParPt.clear();
  _PyParE.clear();
  _PyParMass.clear();


  // Pythia track loop [05.22.2017, Derek]
  int nFinal = 0;
  int nS1    = 0;
  int nS2    = 0;
  int nS3    = 0;
  int nOut   = 0;
  int nIn    = 0;
  int nMid   = 0;
  int nPos   = 0;
  int nNeg   = 0;
  int nNeu   = 0;
  for (int j = 0; j < mParticles_; j++) {

    int    d1     = mParticles_fDaughter[j][0];
    int    d2     = mParticles_fDaughter[j][1];
    int    status = mParticles_fStatusCode[j];
    int    pID    = mParticles_fPdgCode[j];
    float  chrg   = GetCharge(pID);
    double Px     = mParticles_fPx[j];
    double Py     = mParticles_fPy[j];
    double Pz     = mParticles_fPz[j];
    double Pt     = sqrt(Px*Px + Py*Py);
    double P      = sqrt(Pt*Pt + Pz*Pz);
    double theta  = asin(Pt / P);
    double eta    = -1. * log(tan(theta / 2.));
    if (Pz < 0) eta = -1. * eta;

    // additional members [02.02.2018, Derek]
    Int_t    m1     = mParticles_fMother[j][0];
    Int_t    m2     = mParticles_fMother[j][1];
    Float_t  weight = mParticles_fWeight[j]; 
    Double_t mass   = mParticles_fCalcMass[j];
    Double_t energy = mParticles_fE[j];

    // [05.22.2017, Derek]
    if ((d1 == 0) && (d2 == 0)) {
      ++nFinal;
      if (status == 3)
        ++nS3;
      if (status == 2)
        ++nS2;
      if (status == 1) {
        ++nS1;
        if (abs(eta) < 1.)
          ++nMid;
        if (abs(eta) < 1.8) {
          ++nIn;
          if (chrg == 1.)
            ++nPos;
          if (chrg == -1.)
            ++nNeg;
          if (chrg == 0.)
            ++nNeu;
        }
        if (abs(eta) > 1.8)
          ++nOut;
      }
    }

    // fill pythia particle members [02.02.2018, Derek]
    _PyParId.push_back(pID);
    _PyParStatus.push_back(status);
    _PyParKid1.push_back(d1);
    _PyParKid2.push_back(d2);
    _PyParMom1.push_back(m1);
    _PyParMom2.push_back(m2);
    _PyParChrg.push_back(chrg);
    _PyParWeight.push_back(weight);
    _PyParEta.push_back(eta);
    _PyParTheta.push_back(theta);
    _PyParPx.push_back(Px);
    _PyParPy.push_back(Py);
    _PyParPz.push_back(Pz);
    _PyParPt.push_back(Pt);
    _PyParE.push_back(energy);
    _PyParMass.push_back(mass);

  }  // end Pythia track loop

  // fill pythia event members [02.02.2018, Derek]
  _PyEvtId  = mcEvntNumber;
  _PyRunId  = mRunId;
  _PyProId  = mProcessId;
  _PyNumPar = mParticles_;
  _PyPtPart = mPt;
  _PyS      = mS;
  _PyT      = mT;
  _PyU      = mU;
  _PyVx     = mVertex_fX;
  _PyVy     = mVertex_fY;
  _PyVz     = mVertex_fZ;


  // [10.28.2016 / 10.31.2016 / 05.11.2017 / 05.22.2017 / 05.23.2017, Derek]
  hPythiaCounter    -> Fill(mParticles_);
  hFinalCounter     -> Fill(nFinal);
  hFinalCountS1     -> Fill(nS1);
  hFinalCountS2     -> Fill(nS2);
  hFinalCountS3     -> Fill(nS3);
  hNumOutsideAccept -> Fill(nOut);
  hNumInsideAccept  -> Fill(nIn);
  hPyNumMid         -> Fill(nMid);
  hPyMultPos        -> Fill(nPos);
  hPyMultNeg        -> Fill(nNeg);
  hPyMultNeu        -> Fill(nNeu);


  Int_t mcEvntRefmult = mcEvent->eventGeneratorFinalStateTracks();
  Int_t mcEvntNPrimTrk = mcEvent->numberOfPrimaryTracks();
  Double_t mcEvntReactPlane = mcEvent->phiReactionPlane();
  Double_t mcEvntVrtX = mcEvent->primaryVertex()->position().x();
  Double_t mcEvntVrtY = mcEvent->primaryVertex()->position().y();
  Double_t mcEvntVrtZ = mcEvent->primaryVertex()->position().z();
  
  
  cout<<"Event Info "<<mcEvntNumber<<"  "<<mcRunNumber<<"   "<<mcEvntRefmult<<"  "<<mcEvntNPrimTrk<<"   "<<mcEvntReactPlane<<"  "<<mcEvntVrtX<<"  "<<mcEvntVrtY<<"  "<<mcEvntVrtZ<<endl;

  EvValues[0]   = mRunId;  // [02.02.2018, Derek]
  EvValues[1]   = mcEvntNumber;
  EvValues[2]   = 0;  //Trigger Id check getTrgId()
  EvValues[3]   = 0;
  EvValues[4]   = mcEvntNPrimTrk;
  EvValues[5]   = mcEvntRefmult;
  EvValues[6]   = 0;
  EvValues[7]   = mcEvntVrtX;
  EvValues[8]   = mcEvntVrtY;
  EvValues[9]   = mcEvntVrtZ;
  EvValues[10]   = 0;
  EvValues[11]   = 0;      //zdcConincidenceRate
  EvValues[12]   = 0;   //bbcCoincidenceRate
  EvValues[13]   = 0;  //backgroundRate
  EvValues[14]   = 0;  //bbcBlueBackgroundRate
  EvValues[15]   = 0; //bbcYellowBackgroundRate
  EvValues[16]   =0;
  EvValues[17]   =0;
  EvValues[18]   =0 ; //bTOFTrayMultiplicity (need to add)
  EvValues[19]   =0;
  EvValues[20]   = 0;
  EvValues[80]   = 0;
  event->SetEventAttributes(EvValues);

  // [10.28.2016, Derek]
  if (NumEvts != mcEvntNumber) {
    cerr << "WHOAH!! Event numbers don't match up!\n"
         << "  (Pythia no., Geant no.) = (" << NumEvts << ", "
         << mcEvntNumber <<")"
         << endl;
  }

  // fill pythia and geant tree [02.02.2018, Derek]
  _tPythia -> Fill();
  outTree  -> Fill(); 
  return kStOk;

}

//-----------------------------------------------------------------------------------------------------
Int_t StJetTreeMcMaker::Finish()
{
  // [05.22.2017 / 02.02.2018, Derek] 
  File_output       -> cd();
  outTree           -> Write();
  _tPythia          -> Write();
  hGeantCounter     -> Write();
  hPythiaCounter    -> Write();
  hFinalCounter     -> Write();
  hFinalCountS1     -> Write();
  hFinalCountS2     -> Write();
  hFinalCountS3     -> Write();
  hNumOutsideAccept -> Write();
  hNumInsideAccept  -> Write();
  hPyNumMid         -> Write();
  hGeNumMid         -> Write();
  hPyMultPos        -> Write();
  hPyMultNeg        -> Write();
  hPyMultNeu        -> Write();
  hGeMultPos        -> Write();
  hGeMultNeg        -> Write();
  hGeMultNeu        -> Write();
  File_output       -> Close();
  return kStOk;
}

//-----------------------------------------------------------------------------------------------------void
void  StJetTreeMcMaker::Clear(Option_t *opt)
{
  delete [] pTrMatchArr; pTrMatchArr = 0;
  delete []  pTriMatchEtaArr; pTriMatchEtaArr =0;
  delete []  pTrMatchPhiArr; pTrMatchPhiArr = 0;
  delete []  pTrIndexArray; pTrIndexArray = 0;
  delete []  pTrEtaArray; pTrEtaArray =0;
  delete []  pTrPhiArray; pTrPhiArray =0;


  StMaker::Clear();

  // return kStOk;
 }

// [05.10.2017, Derek]
void StJetTreeMcMaker::SetPythiaFile(char *py) {

  TString sPath(py);
  sPythia = sPath;

}  // end 'SetPythiaFile(char*)'



// [05.23.2017, Derek]
Float_t StJetTreeMcMaker::GetCharge(const Int_t pID) {

  Int_t   absID  = abs(pID);
  Float_t charge = -2.;
  switch (absID) {
    case 11:
      if (pID == -11)
        charge = 1.;
      else
        charge = -1.;
      break;
    case 13:
      if (pID == -13)
        charge = 1.;
      else
        charge = -1.;
      break;
    case 211:
      if (pID == 211)
        charge = 1.;
      else
        charge = -1.;
      break;
    case 321:
      if (pID == 321)
        charge = 1.;
      else
        charge = -1.;
      break;
    case 2212:
      if (pID == 2212)
        charge = 1.;
      else
        charge = -1.;
      break;
    case 3112:
      if (pID == -3112)
        charge = 1.;
      else
        charge = -1.;
      break;
    case 3222:
      if (pID == 3222)
        charge = 1.;
      else
        charge = -1.;
      break;
    case 3312:
      if (pID == -3312)
        charge = 1.;
      else
        charge = -1.;
      break;
    case 3334:
      if (pID == -3334)
        charge = 1.;
      else
        charge = -1.;
      break;
    default:
      charge = 0.;
      break;
  }
  return charge;

}  // end 'GetCharge(Int_t)'

// End ------------------------------------------------------------------------
