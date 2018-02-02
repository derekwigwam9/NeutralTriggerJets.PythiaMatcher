// 'MatchEvents.C'
// 10.27.2016
//
// Compares the events located in a '*.geant.root' file against the events
// located in a '*.pythia.root' file (or even against a '*.mudst.file' if
// you want).
//
// Usage Note: This has been modified.  The macro assumes that
//   the provided Pythia and MuDst file are in the same directory as the
//   provided geant file.  Provide the FULL NAME AND PATH for
//   the geant file, but provide ONLY THE NAME of the corresponding
//   Pythia and MuDst files.  [Derek, 02.02.2018]

const Int_t  Nevt       = 2000;
const Int_t  Nmu        = 1;
const char  *PythiaFile = "pt7_9_10162033_1.pythia.root";
const char  *MuDstFile  = "pt7_9_10162033_1.MuDst.root";
const char  *GeantFile  = "/projecta/projectdirs/starprod/embedding/production2009_200GeV/Jet_pp200_2009.elz14/SL11d_embed/10162033/pt7_9_10162033_1.geant.root";
const char  *OutFile    = "test.root";

void MatchEvents(const Int_t nEvt=Nevt, const Int_t nMu=Nmu, const char *geFile=GeantFile, const char *muFile=MuDstFile, const char *pyFile=PythiaFile, const char *oFile=OutFile) {

  gROOT   -> Macro("LoadLogger.C");
  gROOT   -> Macro("loadMuDst.C");
  gSystem -> Load("StarMagField.so");
  gSystem -> Load("StMagF");
  gSystem -> Load("StDetectorDbMaker");
  gSystem -> Load("StTpcDb");
  gSystem -> Load("St_db_Maker");
  gSystem -> Load("StDbUtilities");
  gSystem -> Load("StMcEvent");
  gSystem -> Load("StMcEventMaker");
  gSystem -> Load("StDaqLib");
  gSystem -> Load("StEmcRawMaker");
  gSystem -> Load("StEmcADCtoEMaker");
  gSystem -> Load("StEpcMaker");
  gSystem -> Load("StTriggerUtilities");
  gSystem -> Load("StDbBroker");
  gSystem -> Load("libgeometry_Tables");
  gSystem -> Load("StEEmcUtil");
  gSystem -> Load("StEEmcDbMaker");
  gSystem -> Load("StPreEclMaker");
  gSystem -> Load("StEpcMaker");
  gSystem -> Load("StEmcTriggerMaker");
  gSystem -> Load("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Embedding/Run9pp/PythiaMatching/.sl64_gcc482/lib/libStJetTreeMcMaker.so");
  gSystem -> Load("StJetTreeMcMaker");


  // determine pythia and mudst path [Derek, 02.02.2018]
  TString geant(geFile);
  TString pythia(pyFile);
  TString mudst(muFile);
  TString pPath("");
  TString uPath("");
  Ssiz_t  nPath = geant.Last('/');
  Ssiz_t  nName = nPath + 1;
  pPath.Append(geant.Data(), nName);
  pPath.Append(pythia.Data());
  uPath.Append(geant.Data(), nName);
  uPath.Append(mudst.Data());
  cout << "I/O Info: geant  = '" << geant.Data() << "'\n"
       << "          pythia = '" << pPath.Data() << "'\n"
       << "          mudst  = '" << uPath.Data() << "'"
       << endl;


  StChain *chain = new StChain("StChain");
	
  // I/O maker
  StIOMaker *ioMaker = new StIOMaker;
  ioMaker -> SetFile(geant);
  ioMaker -> SetIOMode("r");
  ioMaker -> SetBranch("*", 0, "0");            // Deactivate all branches
  ioMaker -> SetBranch("geantBranch", 0, "r");  // Activate geant Branch

  // StMcEvent maker
  StMcEventMaker *mcEventMaker = new StMcEventMaker;
  mcEventMaker -> doPrintEventInfo  = false;
  mcEventMaker -> doPrintMemoryInfo = false;

  // MuDst maker
  StMuDstMaker     *muDstMaker = new StMuDstMaker(0, 0, "", uPath.Data(), "", nMu);
  StJetTreeMcMaker *tMcMaker   = new StJetTreeMcMaker("tMcMaker");
  tMcMaker -> SetFileName(oFile);
  tMcMaker -> SetPythiaFile(pPath.Data());


  StMemStat memory;
  memory.PrintMem(NULL);
	
  chain -> Init();
  cout << "chain initialized" << endl;
	
  TStopwatch total;
  TStopwatch timer;
	
  int i=0;
  while((i < nEvt) && (chain->Make() == kStOk)) {
    if (i % 100000 == 0) {
      cout << "done with event " << i
           << "\tcpu: "<< timer.CpuTime() << "\treal: " << timer.RealTime()
           << "\tratio: " << (timer.CpuTime() / timer.RealTime());
      timer.Start();
      memory.PrintMem(NULL);
    }
    i++;
    chain -> Clear();
  }
	
  chain -> ls(3);
  chain -> Finish();

  // announce end
  cout << "\tcpu: " << total.CpuTime() << "\treal: " << total.RealTime()
       << "\tratio: " << (total.CpuTime() / total.RealTime())
       << endl;

  cout << "\n"
       << "-------------\n"
       << "(-: Done :-) \n"
       << "-------------\n"
       << endl;

}

// End ------------------------------------------------------------------------
