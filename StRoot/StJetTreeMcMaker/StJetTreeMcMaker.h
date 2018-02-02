//________________________________________________________________
//This is the Only clustering algorithim for Gamma-hadron correlation
//analysis. This class uses BEMC, BSMD, and TPC information to
//find HighTower info, pi0/photon discrimination info and hadrons
//Tracks projection info to project on to BEMC.
//
//
// STAR published paper: Phys. Rev. C 82, 034909 (2010)
// based on this clustering algorithm.
// Contact: A. Hamed and Saskia Mioduszewski
//
// Date 07/2014: TPC tracks from StMuTrack -  Nihar r. sahoo
//             : globalDca has been implemented
//
//Date  12/2015: "StJetTreeMcMaker" inherited from "StThirdMaker"
//              It makes "JetTree" for Full Jet Reconstruction
//              for Gamma+Jet analysis.      - Nihar r. Sahoo
//Date  06/2016: all events, tracks and tower attributes are added and 
//                                           - Nihar
//Date 10/26/2016: Added some histograms to compare to
//                 pythia.root file
//                                           -Derek
//Date 10/28/2016: Added/updated methods to read in the trees
//                 contained in '*.pythia.root' files.
//                                           -Derek 
//Date 05/10/2017: Added members/methods to specify pythia file
//                                           -Derek
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


#ifndef STAR_StJetTreeMcMaker
#define STAR_StJetTreeMcMaker

#ifndef StMaker_H
#include "StMaker.h"
#include <string>
#include <vector>


#include "StBichsel/Bichsel.h"
#include "StBichsel/dEdxParameterization.h"


#endif


class StEvent;
class StMcEvent;
class TString;
class Bichsel;
class TNtuple;
class TFile;
class TH1F;
class TH2F;
class StPrimaryVertex;

class StTriggerData;
class StTriggerData2007;

class StEmcCluster;
class StEmcCollection;
class StEmcDetector;
class StEmcModule;
class StEmcModuleHitCollection;
class StEmcRawHit;
class StEmcFilter;
//class StThreeVectorF;
class TRandom;
class St_db_Maker;
class StEmcDecoder;
class emcPed_st;
class StMcVertex;


class StMuDstMaker;
//class vector<float>;
//_______Class for Tree

class GammaJetEvent;
class GammaJetTrack;
class GammaJetTower;
class GammaJetTowerUtil;
class TClonesArray;
class TLorentzVector;

class StMuDst;
class StMuEvent;
class StMuTrack;
class StEmcPosition;
class StBemcTables; //v3.14


#define TowerHVchangeMax 400

const  Int_t max_pTracks = 10000;
const  Int_t max_gTracks = 10000;
const  Int_t kMaxmParticles = 10000;



class StJetTreeMcMaker : public StMaker 
{
 private: //BetheBloch      	mBB;
  TRandom*     mRandom;  

  St_db_Maker* mDbMaker;

  StEmcFilter*    mEmcFilter;
  StEvent           *mStEvent; 
  StMcEvent           *mcEvent;
  StEmcDecoder* mEmcDecoder;
  StMuDstMaker* mMuDstMaker;
  StMuEvent *muEvent;
  StMuDst* muDst;


  // [10.28.2016, Derek]
  TTree          *pyChain;   // pointer to the analyzed TTree or TChain
  Int_t           pyCurrent; // current Tree number in a TChain
  Int_t           NumEvts;   // no. of events processed

   // Declaration of leaf types [10.28.2016, Derek]
   //StPythiaEvent   *PythiaBranch;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           mRunId;
   Int_t           mEventId;
   Int_t           mProcessId;
   Int_t           mTune;
   UInt_t          mVertex_fUniqueID;
   UInt_t          mVertex_fBits;
   Double_t        mVertex_fX;
   Double_t        mVertex_fY;
   Double_t        mVertex_fZ;
   Float_t         mS;
   Float_t         mT;
   Float_t         mU;
   Float_t         mPt;
   Float_t         mCosTheta;
   Float_t         mX1;
   Float_t         mX2;
   Int_t           mMstu72;
   Int_t           mMstu73;
   Int_t           mMstp111;
   Float_t         mPartonALL;
   Float_t         mDF1[35];
   Float_t         mDF2[35];
   Float_t         mF1[2];
   Float_t         mF2[2];
   Float_t         mDF1NNPDF[101];
   Float_t         mDF2NNPDF[101];
   Float_t         mF1NNPDF;
   Float_t         mF2NNPDF;
   Int_t           mParticles_;
   UInt_t          mParticles_fUniqueID[kMaxmParticles];   //[mParticles_]
   UInt_t          mParticles_fBits[kMaxmParticles];   //[mParticles_]
   Short_t         mParticles_fLineColor[kMaxmParticles];   //[mParticles_]
   Short_t         mParticles_fLineStyle[kMaxmParticles];   //[mParticles_]
   Short_t         mParticles_fLineWidth[kMaxmParticles];   //[mParticles_]
   Int_t           mParticles_fPdgCode[kMaxmParticles];   //[mParticles_]
   Int_t           mParticles_fStatusCode[kMaxmParticles];   //[mParticles_]
   Int_t           mParticles_fMother[kMaxmParticles][2];   //[mParticles_]
   Int_t           mParticles_fDaughter[kMaxmParticles][2];   //[mParticles_]
   Float_t         mParticles_fWeight[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fCalcMass[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPx[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPy[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPz[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fE[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVx[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVy[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVz[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVt[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPolarTheta[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPolarPhi[kMaxmParticles];   //[mParticles_]

   // List of branches [10.28.2016, Derek]
   TBranch        *b_PythiaBranch_fUniqueID;   //!
   TBranch        *b_PythiaBranch_fBits;   //!
   TBranch        *b_PythiaBranch_mRunId;   //!
   TBranch        *b_PythiaBranch_mEventId;   //!
   TBranch        *b_PythiaBranch_mProcessId;   //!
   TBranch        *b_PythiaBranch_mTune;   //!
   TBranch        *b_PythiaBranch_mVertex_fUniqueID;   //!
   TBranch        *b_PythiaBranch_mVertex_fBits;   //!
   TBranch        *b_PythiaBranch_mVertex_fX;   //!
   TBranch        *b_PythiaBranch_mVertex_fY;   //!
   TBranch        *b_PythiaBranch_mVertex_fZ;   //!
   TBranch        *b_PythiaBranch_mS;   //!
   TBranch        *b_PythiaBranch_mT;   //!
   TBranch        *b_PythiaBranch_mU;   //!
   TBranch        *b_PythiaBranch_mPt;   //!
   TBranch        *b_PythiaBranch_mCosTheta;   //!
   TBranch        *b_PythiaBranch_mX1;   //!
   TBranch        *b_PythiaBranch_mX2;   //!
   TBranch        *b_PythiaBranch_mMstu72;   //!
   TBranch        *b_PythiaBranch_mMstu73;   //!
   TBranch        *b_PythiaBranch_mMstp111;   //!
   TBranch        *b_PythiaBranch_mPartonALL;   //!
   TBranch        *b_PythiaBranch_mDF1;   //!
   TBranch        *b_PythiaBranch_mDF2;   //!
   TBranch        *b_PythiaBranch_mF1;   //!
   TBranch        *b_PythiaBranch_mF2;   //!
   TBranch        *b_PythiaBranch_mDF1NNPDF;   //!
   TBranch        *b_PythiaBranch_mDF2NNPDF;   //!
   TBranch        *b_PythiaBranch_mF1NNPDF;   //!
   TBranch        *b_PythiaBranch_mF2NNPDF;   //!
   TBranch        *b_PythiaBranch_mParticles_;   //!
   TBranch        *b_mParticles_fUniqueID;   //!
   TBranch        *b_mParticles_fBits;   //!
   TBranch        *b_mParticles_fLineColor;   //!
   TBranch        *b_mParticles_fLineStyle;   //!
   TBranch        *b_mParticles_fLineWidth;   //!
   TBranch        *b_mParticles_fPdgCode;   //!
   TBranch        *b_mParticles_fStatusCode;   //!
   TBranch        *b_mParticles_fMother;   //!
   TBranch        *b_mParticles_fDaughter;   //!
   TBranch        *b_mParticles_fWeight;   //!
   TBranch        *b_mParticles_fCalcMass;   //!
   TBranch        *b_mParticles_fPx;   //!
   TBranch        *b_mParticles_fPy;   //!
   TBranch        *b_mParticles_fPz;   //!
   TBranch        *b_mParticles_fE;   //!
   TBranch        *b_mParticles_fVx;   //!
   TBranch        *b_mParticles_fVy;   //!
   TBranch        *b_mParticles_fVz;   //!
   TBranch        *b_mParticles_fVt;   //!
   TBranch        *b_mParticles_fPolarTheta;   //!
   TBranch        *b_mParticles_fPolarPhi;   //!

 
 public:
  void          setPrint(Bool_t);
  StJetTreeMcMaker(const char *name="ThirdMaker", char* dataType = "");
  virtual       ~StJetTreeMcMaker();

  Int_t doMCEvent();
  Int_t  doMuEvent();
  virtual Int_t Finish();
  void     Clear(Option_t *option="");
  virtual Int_t Init();
  virtual Int_t Make();
  void SetFileName( char* name){outFile=name;};
  void SetPythiaFile(char *py);
 
  Int_t ResizeHit(); 
  void  setDbMaker(St_db_Maker*);
  Bichsel* m_dEdxParameterization;
  StEmcFilter* getEmcFilter() { return mEmcFilter; };

  Bool_t Check_hot_Tower(Int_t TwrId);
  Bool_t Check_hot_EtaStrip(Int_t EtaStpId);

  // [05.23.2017, Derek]
  Float_t GetCharge(const Int_t pID);


  const StMuTrack* ptrack;  
  StEmcPosition *mPosition;

  StBemcTables    *mTables;


  
  ///////////////////////

  TH1F *mHistBg;
  TH1D *hTotalEbemc;
  TH1D *hEbemcTwr;
  // [10.26.2016, Derek]
  TH1D *hGeantCounter;
  // [10.28.2016, Derek]
  TH1D *hPythiaCounter;
  TH1D *hFinalCounter;
  TH1D *hFinalCountS1;
  TH1D *hFinalCountS2;
  TH1D *hFinalCountS3;
  // [10.31.2016, Derek]
  TH1D *hNumOutsideAccept;
  // [05.11.2017, Derek]
  TH1D *hNumInsideAccept;
  // [05.22.2017, Derek]
  TH1D *hPyNumMid;
  TH1D *hGeNumMid;
  // [05.23.2017, Derek]
  TH1D *hPyMultPos;
  TH1D *hPyMultNeg;
  TH1D *hPyMultNeu;
  TH1D *hGeMultPos;
  TH1D *hGeMultNeg;
  TH1D *hGeMultNeu;

  // for pythia tree[02.02.2018, Derek]
  TTree   *_tPythia;

  Int_t    _PyEvtId;
  Int_t    _PyRunId;
  Int_t    _PyProId;
  Int_t    _PyNumPar;
  Float_t  _PyPtPart;
  Float_t  _PyS;
  Float_t  _PyT;
  Float_t  _PyU;
  Double_t _PyVx;
  Double_t _PyVy;
  Double_t _PyVz;

  vector<Int_t>    _PyParId;
  vector<Int_t>    _PyParStatus;
  vector<Int_t>    _PyParKid1;
  vector<Int_t>    _PyParKid2;
  vector<Int_t>    _PyParMom1;
  vector<Int_t>    _PyParMom2;
  vector<Float_t>  _PyParChrg;
  vector<Float_t>  _PyParWeight;
  vector<Double_t> _PyParEta;
  vector<Double_t> _PyParTheta;
  vector<Double_t> _PyParPx;
  vector<Double_t> _PyParPy;
  vector<Double_t> _PyParPz;
  vector<Double_t> _PyParPt;
  vector<Double_t> _PyParE;
  vector<Double_t> _PyParMass;

  // [05.10.2017, Derek]
  TString sPythia;
  TFile   *fPythia;

  TH1F* hya;
  int  mEventCounter;
  TFile* File_output;
  char* outFile;
  
  
  ///________
  
  TBranch *eventBranch;
  
  Float_t *EvValues;
  
  TTree *outTree;
  GammaJetEvent *event;

  
  
  //_____________for matching tracks
  Int_t *pTrMatchArr;//!
  Float_t *pTriMatchEtaArr;//!
  Float_t *pTrMatchPhiArr;//!
  Int_t *pTrIndexArray;//!
  Float_t *pTrEtaArray;//!
  Float_t *pTrPhiArray;//!
  
  
  ClassDef(StJetTreeMcMaker, 1)   //StAF chain virtual base class for Makers
};






#endif
