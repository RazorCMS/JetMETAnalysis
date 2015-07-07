//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 31 12:29:31 2014 by ROOT version 5.34/18
// from TChain promptanaTree/tree/
//////////////////////////////////////////////////////////

#ifndef rootNtupleClass_h
#define rootNtupleClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;

class rootNtupleClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Bool_t          isData;
   vector<bool>    *vertexisValid;
   Double_t        time;
   vector<double>  *ak5PFJetChargedEmEnergyFraction;
   vector<double>  *ak5PFJetChargedHadronEnergyFraction;
   vector<double>  *ak5PFJetChargedMuEnergyFraction;
   vector<double>  *ak5PFJetChargedMultiplicity;
   vector<double>  *ak5PFJetEnergy;
   vector<double>  *ak5PFJetEta;
   vector<double>  *ak5PFJetMuonMultiplicity;
   vector<double>  *ak5PFJetNeutralEmEnergyFraction;
   vector<double>  *ak5PFJetNeutralHadronEnergyFraction;
   vector<double>  *ak5PFJetNeutralMultiplicity;
   vector<double>  *ak5PFJetPhi;
   vector<double>  *ak5PFJetpT;
   vector<double>  *ak5PFJetscaleL2L3;
   vector<double>  *CaloTowersEcalTime;
   vector<double>  *CaloTowersEmEt;
   vector<double>  *CaloTowersEmEtVtx0;
   vector<double>  *CaloTowersEta;
   vector<double>  *CaloTowersHadEt;
   vector<double>  *CaloTowersHadEtVtx0;
   vector<double>  *CaloTowersHcalTime;
   vector<double>  *CaloTowersOuterEt;
   vector<double>  *CaloTowersPhi;
   vector<double>  *calometEmEtFraction;
   vector<double>  *calometEmEtInEB;
   vector<double>  *calometEmEtInEE;
   vector<double>  *calometEmEtInHF;
   vector<double>  *calometEtFractionHadronic;
   vector<double>  *calometHadEtInHB;
   vector<double>  *calometHadEtInHE;
   vector<double>  *calometHadEtInHF;
   vector<double>  *calometHadEtInHO;
   vector<double>  *calometMETInmHF;
   vector<double>  *calometMETInpHF;
   vector<double>  *calometMETPhiInmHF;
   vector<double>  *calometMETPhiInpHF;
   vector<double>  *calometMETSigCornell;
   vector<double>  *calometMaxEtInEmTowers;
   vector<double>  *calometMaxEtInHadTowers;
   vector<double>  *calometPhi;
   vector<double>  *calometPt;
   vector<double>  *calometPx;
   vector<double>  *calometPy;
   vector<double>  *calometSETInmHF;
   vector<double>  *calometSETInpHF;
   vector<double>  *calometSumEt;
   vector<double>  *pfmetPhi;
   vector<double>  *pfmetPt;
   vector<double>  *pfmetPx;
   vector<double>  *pfmetPy;
   vector<double>  *pfmetSumEt;
   vector<double>  *tracksChi2;
   vector<double>  *tracksDXY;
   vector<double>  *tracksDXYError;
   vector<double>  *tracksEta;
   vector<double>  *tracksEtaError;
   vector<double>  *tracksNDOF;
   vector<double>  *tracksPhi;
   vector<double>  *tracksPhiError;
   vector<double>  *tracksPt;
   vector<double>  *tracksPtError;
   vector<double>  *vertexChi2;
   vector<double>  *vertexNDF;
   vector<double>  *vertexSumPt;
   vector<double>  *vertexSumPtW5;
   vector<double>  *vertexX;
   vector<double>  *vertexXErr;
   vector<double>  *vertexY;
   vector<double>  *vertexYErr;
   vector<double>  *vertexZ;
   vector<double>  *vertexZErr;
   vector<float>   *ECALnoiseECALenergy11x11;
   vector<float>   *ECALnoiseECALenergy13x13;
   vector<float>   *ECALnoiseECALenergy3x3;
   vector<float>   *ECALnoiseECALenergy5x5;
   vector<float>   *ECALnoiseECALenergy7x7;
   vector<float>   *ECALnoiseECALenergy9x9;
   vector<float>   *ECALnoiseECalEBRechiEnergy;
   vector<float>   *ECALnoiseECalEBRechiEta;
   vector<float>   *ECALnoiseECalEBRechiPhi;
   vector<float>   *ECALnoiseECalEBRechiSwissCross;
   vector<float>   *ECALnoiseECalEBRechiTime;
   vector<float>   *ECALnoiseECalEBRechiiEta;
   vector<float>   *ECALnoiseECalEBRechiiPhi;
   vector<float>   *ECALnoiseECalEBSeedChi2;
   vector<float>   *ECALnoiseECalEBSeedE3x3;
   vector<float>   *ECALnoiseECalEBSeedEBottom;
   vector<float>   *ECALnoiseECalEBSeedELeft;
   vector<float>   *ECALnoiseECalEBSeedERight;
   vector<float>   *ECALnoiseECalEBSeedETop;
   vector<float>   *ECALnoiseECalEBSeedEnergy;
   vector<float>   *ECALnoiseECalEBSeedEta;
   vector<float>   *ECALnoiseECalEBSeedPhi;
   vector<float>   *ECALnoiseECalEBSeedTime;
   vector<float>   *ECALnoiseECalEBSeediEta;
   vector<float>   *ECALnoiseECalEBSeediPhi;
   vector<float>   *ECALnoiseHCALenergy3x3;
   vector<float>   *ECALnoiseHCALenergy5x5;
   vector<float>   *ECALnoiseHCALenergy7x7;
   vector<float>   *ECALnoiseHCALenergy9x9;
   vector<float>   *ECALnoiseHCALenergyUp;
   vector<int>     *ak5PFJetNConstituents;
   vector<int>     *ak5PFJetNJets;
   vector<int>     *CaloTowersIeta;
   vector<int>     *CaloTowersIphi;
   vector<int>     *CaloTowersTowerStatusWord;
   vector<int>     *ECALnoiseECalEBSeedRecoFlag;
   vector<int>     *tracksAlgorithm;
   vector<int>     *tracksNumberOfValidHits;
   vector<int>     *tracksNumberOfValidPixelHits;
   vector<int>     *tracksNumberOfValidStripHits;
   vector<int>     *tracksQuality;
   vector<int>     *vertexNTracks;
   vector<int>     *vertexNTracksW5;
   UInt_t          bunch;
   UInt_t          event;
   UInt_t          ls;
   UInt_t          orbit;
   UInt_t          run;

   // List of branches
   TBranch        *b_isData;   //!
   TBranch        *b_vertexisValid;   //!
   TBranch        *b_time;   //!
   TBranch        *b_ak5PFJetChargedEmEnergyFraction;   //!
   TBranch        *b_ak5PFJetChargedHadronEnergyFraction;   //!
   TBranch        *b_ak5PFJetChargedMuEnergyFraction;   //!
   TBranch        *b_ak5PFJetChargedMultiplicity;   //!
   TBranch        *b_ak5PFJetEnergy;   //!
   TBranch        *b_ak5PFJetEta;   //!
   TBranch        *b_ak5PFJetMuonMultiplicity;   //!
   TBranch        *b_ak5PFJetNeutralEmEnergyFraction;   //!
   TBranch        *b_ak5PFJetNeutralHadronEnergyFraction;   //!
   TBranch        *b_ak5PFJetNeutralMultiplicity;   //!
   TBranch        *b_ak5PFJetPhi;   //!
   TBranch        *b_ak5PFJetpT;   //!
   TBranch        *b_ak5PFJetscaleL2L3;   //!
   TBranch        *b_CaloTowersEcalTime;   //!
   TBranch        *b_CaloTowersEmEt;   //!
   TBranch        *b_CaloTowersEmEtVtx0;   //!
   TBranch        *b_CaloTowersEta;   //!
   TBranch        *b_CaloTowersHadEt;   //!
   TBranch        *b_CaloTowersHadEtVtx0;   //!
   TBranch        *b_CaloTowersHcalTime;   //!
   TBranch        *b_CaloTowersOuterEt;   //!
   TBranch        *b_CaloTowersPhi;   //!
   TBranch        *b_calometEmEtFraction;   //!
   TBranch        *b_calometEmEtInEB;   //!
   TBranch        *b_calometEmEtInEE;   //!
   TBranch        *b_calometEmEtInHF;   //!
   TBranch        *b_calometEtFractionHadronic;   //!
   TBranch        *b_calometHadEtInHB;   //!
   TBranch        *b_calometHadEtInHE;   //!
   TBranch        *b_calometHadEtInHF;   //!
   TBranch        *b_calometHadEtInHO;   //!
   TBranch        *b_calometMETInmHF;   //!
   TBranch        *b_calometMETInpHF;   //!
   TBranch        *b_calometMETPhiInmHF;   //!
   TBranch        *b_calometMETPhiInpHF;   //!
   TBranch        *b_calometMETSigCornell;   //!
   TBranch        *b_calometMaxEtInEmTowers;   //!
   TBranch        *b_calometMaxEtInHadTowers;   //!
   TBranch        *b_calometPhi;   //!
   TBranch        *b_calometPt;   //!
   TBranch        *b_calometPx;   //!
   TBranch        *b_calometPy;   //!
   TBranch        *b_calometSETInmHF;   //!
   TBranch        *b_calometSETInpHF;   //!
   TBranch        *b_calometSumEt;   //!
   TBranch        *b_pfmetPhi;   //!
   TBranch        *b_pfmetPt;   //!
   TBranch        *b_pfmetPx;   //!
   TBranch        *b_pfmetPy;   //!
   TBranch        *b_pfmetSumEt;   //!
   TBranch        *b_tracksChi2;   //!
   TBranch        *b_tracksDXY;   //!
   TBranch        *b_tracksDXYError;   //!
   TBranch        *b_tracksEta;   //!
   TBranch        *b_tracksEtaError;   //!
   TBranch        *b_tracksNDOF;   //!
   TBranch        *b_tracksPhi;   //!
   TBranch        *b_tracksPhiError;   //!
   TBranch        *b_tracksPt;   //!
   TBranch        *b_tracksPtError;   //!
   TBranch        *b_vertexChi2;   //!
   TBranch        *b_vertexNDF;   //!
   TBranch        *b_vertexSumPt;   //!
   TBranch        *b_vertexSumPtW5;   //!
   TBranch        *b_vertexX;   //!
   TBranch        *b_vertexXErr;   //!
   TBranch        *b_vertexY;   //!
   TBranch        *b_vertexYErr;   //!
   TBranch        *b_vertexZ;   //!
   TBranch        *b_vertexZErr;   //!
   TBranch        *b_ECALnoiseECALenergy11x11;   //!
   TBranch        *b_ECALnoiseECALenergy13x13;   //!
   TBranch        *b_ECALnoiseECALenergy3x3;   //!
   TBranch        *b_ECALnoiseECALenergy5x5;   //!
   TBranch        *b_ECALnoiseECALenergy7x7;   //!
   TBranch        *b_ECALnoiseECALenergy9x9;   //!
   TBranch        *b_ECALnoiseECalEBRechiEnergy;   //!
   TBranch        *b_ECALnoiseECalEBRechiEta;   //!
   TBranch        *b_ECALnoiseECalEBRechiPhi;   //!
   TBranch        *b_ECALnoiseECalEBRechiSwissCross;   //!
   TBranch        *b_ECALnoiseECalEBRechiTime;   //!
   TBranch        *b_ECALnoiseECalEBRechiiEta;   //!
   TBranch        *b_ECALnoiseECalEBRechiiPhi;   //!
   TBranch        *b_ECALnoiseECalEBSeedChi2;   //!
   TBranch        *b_ECALnoiseECalEBSeedE3x3;   //!
   TBranch        *b_ECALnoiseECalEBSeedEBottom;   //!
   TBranch        *b_ECALnoiseECalEBSeedELeft;   //!
   TBranch        *b_ECALnoiseECalEBSeedERight;   //!
   TBranch        *b_ECALnoiseECalEBSeedETop;   //!
   TBranch        *b_ECALnoiseECalEBSeedEnergy;   //!
   TBranch        *b_ECALnoiseECalEBSeedEta;   //!
   TBranch        *b_ECALnoiseECalEBSeedPhi;   //!
   TBranch        *b_ECALnoiseECalEBSeedTime;   //!
   TBranch        *b_ECALnoiseECalEBSeediEta;   //!
   TBranch        *b_ECALnoiseECalEBSeediPhi;   //!
   TBranch        *b_ECALnoiseHCALenergy3x3;   //!
   TBranch        *b_ECALnoiseHCALenergy5x5;   //!
   TBranch        *b_ECALnoiseHCALenergy7x7;   //!
   TBranch        *b_ECALnoiseHCALenergy9x9;   //!
   TBranch        *b_ECALnoiseHCALenergyUp;   //!
   TBranch        *b_ak5PFJetNConstituents;   //!
   TBranch        *b_ak5PFJetNJets;   //!
   TBranch        *b_CaloTowersIeta;   //!
   TBranch        *b_CaloTowersIphi;   //!
   TBranch        *b_CaloTowersTowerStatusWord;   //!
   TBranch        *b_ECALnoiseECalEBSeedRecoFlag;   //!
   TBranch        *b_tracksAlgorithm;   //!
   TBranch        *b_tracksNumberOfValidHits;   //!
   TBranch        *b_tracksNumberOfValidPixelHits;   //!
   TBranch        *b_tracksNumberOfValidStripHits;   //!
   TBranch        *b_tracksQuality;   //!
   TBranch        *b_vertexNTracks;   //!
   TBranch        *b_vertexNTracksW5;   //!
   TBranch        *b_bunch;   //!
   TBranch        *b_event;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_run;   //!

   rootNtupleClass(TTree *tree=0);
   virtual ~rootNtupleClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef rootNtupleClass_cxx
rootNtupleClass::rootNtupleClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("promptanaTree/tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("promptanaTree/tree","");
      chain->Add("test_nocorr.root/promptanaTree/tree");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

rootNtupleClass::~rootNtupleClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t rootNtupleClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t rootNtupleClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void rootNtupleClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vertexisValid = 0;
   ak5PFJetChargedEmEnergyFraction = 0;
   ak5PFJetChargedHadronEnergyFraction = 0;
   ak5PFJetChargedMuEnergyFraction = 0;
   ak5PFJetChargedMultiplicity = 0;
   ak5PFJetEnergy = 0;
   ak5PFJetEta = 0;
   ak5PFJetMuonMultiplicity = 0;
   ak5PFJetNeutralEmEnergyFraction = 0;
   ak5PFJetNeutralHadronEnergyFraction = 0;
   ak5PFJetNeutralMultiplicity = 0;
   ak5PFJetPhi = 0;
   ak5PFJetpT = 0;
   ak5PFJetscaleL2L3 = 0;
   CaloTowersEcalTime = 0;
   CaloTowersEmEt = 0;
   CaloTowersEmEtVtx0 = 0;
   CaloTowersEta = 0;
   CaloTowersHadEt = 0;
   CaloTowersHadEtVtx0 = 0;
   CaloTowersHcalTime = 0;
   CaloTowersOuterEt = 0;
   CaloTowersPhi = 0;
   calometEmEtFraction = 0;
   calometEmEtInEB = 0;
   calometEmEtInEE = 0;
   calometEmEtInHF = 0;
   calometEtFractionHadronic = 0;
   calometHadEtInHB = 0;
   calometHadEtInHE = 0;
   calometHadEtInHF = 0;
   calometHadEtInHO = 0;
   calometMETInmHF = 0;
   calometMETInpHF = 0;
   calometMETPhiInmHF = 0;
   calometMETPhiInpHF = 0;
   calometMETSigCornell = 0;
   calometMaxEtInEmTowers = 0;
   calometMaxEtInHadTowers = 0;
   calometPhi = 0;
   calometPt = 0;
   calometPx = 0;
   calometPy = 0;
   calometSETInmHF = 0;
   calometSETInpHF = 0;
   calometSumEt = 0;
   pfmetPhi = 0;
   pfmetPt = 0;
   pfmetPx = 0;
   pfmetPy = 0;
   pfmetSumEt = 0;
   tracksChi2 = 0;
   tracksDXY = 0;
   tracksDXYError = 0;
   tracksEta = 0;
   tracksEtaError = 0;
   tracksNDOF = 0;
   tracksPhi = 0;
   tracksPhiError = 0;
   tracksPt = 0;
   tracksPtError = 0;
   vertexChi2 = 0;
   vertexNDF = 0;
   vertexSumPt = 0;
   vertexSumPtW5 = 0;
   vertexX = 0;
   vertexXErr = 0;
   vertexY = 0;
   vertexYErr = 0;
   vertexZ = 0;
   vertexZErr = 0;
   ECALnoiseECALenergy11x11 = 0;
   ECALnoiseECALenergy13x13 = 0;
   ECALnoiseECALenergy3x3 = 0;
   ECALnoiseECALenergy5x5 = 0;
   ECALnoiseECALenergy7x7 = 0;
   ECALnoiseECALenergy9x9 = 0;
   ECALnoiseECalEBRechiEnergy = 0;
   ECALnoiseECalEBRechiEta = 0;
   ECALnoiseECalEBRechiPhi = 0;
   ECALnoiseECalEBRechiSwissCross = 0;
   ECALnoiseECalEBRechiTime = 0;
   ECALnoiseECalEBRechiiEta = 0;
   ECALnoiseECalEBRechiiPhi = 0;
   ECALnoiseECalEBSeedChi2 = 0;
   ECALnoiseECalEBSeedE3x3 = 0;
   ECALnoiseECalEBSeedEBottom = 0;
   ECALnoiseECalEBSeedELeft = 0;
   ECALnoiseECalEBSeedERight = 0;
   ECALnoiseECalEBSeedETop = 0;
   ECALnoiseECalEBSeedEnergy = 0;
   ECALnoiseECalEBSeedEta = 0;
   ECALnoiseECalEBSeedPhi = 0;
   ECALnoiseECalEBSeedTime = 0;
   ECALnoiseECalEBSeediEta = 0;
   ECALnoiseECalEBSeediPhi = 0;
   ECALnoiseHCALenergy3x3 = 0;
   ECALnoiseHCALenergy5x5 = 0;
   ECALnoiseHCALenergy7x7 = 0;
   ECALnoiseHCALenergy9x9 = 0;
   ECALnoiseHCALenergyUp = 0;
   ak5PFJetNConstituents = 0;
   ak5PFJetNJets = 0;
   CaloTowersIeta = 0;
   CaloTowersIphi = 0;
   CaloTowersTowerStatusWord = 0;
   ECALnoiseECalEBSeedRecoFlag = 0;
   tracksAlgorithm = 0;
   tracksNumberOfValidHits = 0;
   tracksNumberOfValidPixelHits = 0;
   tracksNumberOfValidStripHits = 0;
   tracksQuality = 0;
   vertexNTracks = 0;
   vertexNTracksW5 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("vertexisValid", &vertexisValid, &b_vertexisValid);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("ak5PFJetChargedEmEnergyFraction", &ak5PFJetChargedEmEnergyFraction, &b_ak5PFJetChargedEmEnergyFraction);
   fChain->SetBranchAddress("ak5PFJetChargedHadronEnergyFraction", &ak5PFJetChargedHadronEnergyFraction, &b_ak5PFJetChargedHadronEnergyFraction);
   fChain->SetBranchAddress("ak5PFJetChargedMuEnergyFraction", &ak5PFJetChargedMuEnergyFraction, &b_ak5PFJetChargedMuEnergyFraction);
   fChain->SetBranchAddress("ak5PFJetChargedMultiplicity", &ak5PFJetChargedMultiplicity, &b_ak5PFJetChargedMultiplicity);
   fChain->SetBranchAddress("ak5PFJetEnergy", &ak5PFJetEnergy, &b_ak5PFJetEnergy);
   fChain->SetBranchAddress("ak5PFJetEta", &ak5PFJetEta, &b_ak5PFJetEta);
   fChain->SetBranchAddress("ak5PFJetMuonMultiplicity", &ak5PFJetMuonMultiplicity, &b_ak5PFJetMuonMultiplicity);
   fChain->SetBranchAddress("ak5PFJetNeutralEmEnergyFraction", &ak5PFJetNeutralEmEnergyFraction, &b_ak5PFJetNeutralEmEnergyFraction);
   fChain->SetBranchAddress("ak5PFJetNeutralHadronEnergyFraction", &ak5PFJetNeutralHadronEnergyFraction, &b_ak5PFJetNeutralHadronEnergyFraction);
   fChain->SetBranchAddress("ak5PFJetNeutralMultiplicity", &ak5PFJetNeutralMultiplicity, &b_ak5PFJetNeutralMultiplicity);
   fChain->SetBranchAddress("ak5PFJetPhi", &ak5PFJetPhi, &b_ak5PFJetPhi);
   fChain->SetBranchAddress("ak5PFJetpT", &ak5PFJetpT, &b_ak5PFJetpT);
   fChain->SetBranchAddress("ak5PFJetscaleL2L3", &ak5PFJetscaleL2L3, &b_ak5PFJetscaleL2L3);
   fChain->SetBranchAddress("CaloTowersEcalTime", &CaloTowersEcalTime, &b_CaloTowersEcalTime);
   fChain->SetBranchAddress("CaloTowersEmEt", &CaloTowersEmEt, &b_CaloTowersEmEt);
   fChain->SetBranchAddress("CaloTowersEmEtVtx0", &CaloTowersEmEtVtx0, &b_CaloTowersEmEtVtx0);
   fChain->SetBranchAddress("CaloTowersEta", &CaloTowersEta, &b_CaloTowersEta);
   fChain->SetBranchAddress("CaloTowersHadEt", &CaloTowersHadEt, &b_CaloTowersHadEt);
   fChain->SetBranchAddress("CaloTowersHadEtVtx0", &CaloTowersHadEtVtx0, &b_CaloTowersHadEtVtx0);
   fChain->SetBranchAddress("CaloTowersHcalTime", &CaloTowersHcalTime, &b_CaloTowersHcalTime);
   fChain->SetBranchAddress("CaloTowersOuterEt", &CaloTowersOuterEt, &b_CaloTowersOuterEt);
   fChain->SetBranchAddress("CaloTowersPhi", &CaloTowersPhi, &b_CaloTowersPhi);
   fChain->SetBranchAddress("calometEmEtFraction", &calometEmEtFraction, &b_calometEmEtFraction);
   fChain->SetBranchAddress("calometEmEtInEB", &calometEmEtInEB, &b_calometEmEtInEB);
   fChain->SetBranchAddress("calometEmEtInEE", &calometEmEtInEE, &b_calometEmEtInEE);
   fChain->SetBranchAddress("calometEmEtInHF", &calometEmEtInHF, &b_calometEmEtInHF);
   fChain->SetBranchAddress("calometEtFractionHadronic", &calometEtFractionHadronic, &b_calometEtFractionHadronic);
   fChain->SetBranchAddress("calometHadEtInHB", &calometHadEtInHB, &b_calometHadEtInHB);
   fChain->SetBranchAddress("calometHadEtInHE", &calometHadEtInHE, &b_calometHadEtInHE);
   fChain->SetBranchAddress("calometHadEtInHF", &calometHadEtInHF, &b_calometHadEtInHF);
   fChain->SetBranchAddress("calometHadEtInHO", &calometHadEtInHO, &b_calometHadEtInHO);
   fChain->SetBranchAddress("calometMETInmHF", &calometMETInmHF, &b_calometMETInmHF);
   fChain->SetBranchAddress("calometMETInpHF", &calometMETInpHF, &b_calometMETInpHF);
   fChain->SetBranchAddress("calometMETPhiInmHF", &calometMETPhiInmHF, &b_calometMETPhiInmHF);
   fChain->SetBranchAddress("calometMETPhiInpHF", &calometMETPhiInpHF, &b_calometMETPhiInpHF);
   fChain->SetBranchAddress("calometMETSigCornell", &calometMETSigCornell, &b_calometMETSigCornell);
   fChain->SetBranchAddress("calometMaxEtInEmTowers", &calometMaxEtInEmTowers, &b_calometMaxEtInEmTowers);
   fChain->SetBranchAddress("calometMaxEtInHadTowers", &calometMaxEtInHadTowers, &b_calometMaxEtInHadTowers);
   fChain->SetBranchAddress("calometPhi", &calometPhi, &b_calometPhi);
   fChain->SetBranchAddress("calometPt", &calometPt, &b_calometPt);
   fChain->SetBranchAddress("calometPx", &calometPx, &b_calometPx);
   fChain->SetBranchAddress("calometPy", &calometPy, &b_calometPy);
   fChain->SetBranchAddress("calometSETInmHF", &calometSETInmHF, &b_calometSETInmHF);
   fChain->SetBranchAddress("calometSETInpHF", &calometSETInpHF, &b_calometSETInpHF);
   fChain->SetBranchAddress("calometSumEt", &calometSumEt, &b_calometSumEt);
   fChain->SetBranchAddress("pfmetPhi", &pfmetPhi, &b_pfmetPhi);
   fChain->SetBranchAddress("pfmetPt", &pfmetPt, &b_pfmetPt);
   fChain->SetBranchAddress("pfmetPx", &pfmetPx, &b_pfmetPx);
   fChain->SetBranchAddress("pfmetPy", &pfmetPy, &b_pfmetPy);
   fChain->SetBranchAddress("pfmetSumEt", &pfmetSumEt, &b_pfmetSumEt);
   fChain->SetBranchAddress("tracksChi2", &tracksChi2, &b_tracksChi2);
   fChain->SetBranchAddress("tracksDXY", &tracksDXY, &b_tracksDXY);
   fChain->SetBranchAddress("tracksDXYError", &tracksDXYError, &b_tracksDXYError);
   fChain->SetBranchAddress("tracksEta", &tracksEta, &b_tracksEta);
   fChain->SetBranchAddress("tracksEtaError", &tracksEtaError, &b_tracksEtaError);
   fChain->SetBranchAddress("tracksNDOF", &tracksNDOF, &b_tracksNDOF);
   fChain->SetBranchAddress("tracksPhi", &tracksPhi, &b_tracksPhi);
   fChain->SetBranchAddress("tracksPhiError", &tracksPhiError, &b_tracksPhiError);
   fChain->SetBranchAddress("tracksPt", &tracksPt, &b_tracksPt);
   fChain->SetBranchAddress("tracksPtError", &tracksPtError, &b_tracksPtError);
   fChain->SetBranchAddress("vertexChi2", &vertexChi2, &b_vertexChi2);
   fChain->SetBranchAddress("vertexNDF", &vertexNDF, &b_vertexNDF);
   fChain->SetBranchAddress("vertexSumPt", &vertexSumPt, &b_vertexSumPt);
   fChain->SetBranchAddress("vertexSumPtW5", &vertexSumPtW5, &b_vertexSumPtW5);
   fChain->SetBranchAddress("vertexX", &vertexX, &b_vertexX);
   fChain->SetBranchAddress("vertexXErr", &vertexXErr, &b_vertexXErr);
   fChain->SetBranchAddress("vertexY", &vertexY, &b_vertexY);
   fChain->SetBranchAddress("vertexYErr", &vertexYErr, &b_vertexYErr);
   fChain->SetBranchAddress("vertexZ", &vertexZ, &b_vertexZ);
   fChain->SetBranchAddress("vertexZErr", &vertexZErr, &b_vertexZErr);
   fChain->SetBranchAddress("ECALnoiseECALenergy11x11", &ECALnoiseECALenergy11x11, &b_ECALnoiseECALenergy11x11);
   fChain->SetBranchAddress("ECALnoiseECALenergy13x13", &ECALnoiseECALenergy13x13, &b_ECALnoiseECALenergy13x13);
   fChain->SetBranchAddress("ECALnoiseECALenergy3x3", &ECALnoiseECALenergy3x3, &b_ECALnoiseECALenergy3x3);
   fChain->SetBranchAddress("ECALnoiseECALenergy5x5", &ECALnoiseECALenergy5x5, &b_ECALnoiseECALenergy5x5);
   fChain->SetBranchAddress("ECALnoiseECALenergy7x7", &ECALnoiseECALenergy7x7, &b_ECALnoiseECALenergy7x7);
   fChain->SetBranchAddress("ECALnoiseECALenergy9x9", &ECALnoiseECALenergy9x9, &b_ECALnoiseECALenergy9x9);
   fChain->SetBranchAddress("ECALnoiseECalEBRechiEnergy", &ECALnoiseECalEBRechiEnergy, &b_ECALnoiseECalEBRechiEnergy);
   fChain->SetBranchAddress("ECALnoiseECalEBRechiEta", &ECALnoiseECalEBRechiEta, &b_ECALnoiseECalEBRechiEta);
   fChain->SetBranchAddress("ECALnoiseECalEBRechiPhi", &ECALnoiseECalEBRechiPhi, &b_ECALnoiseECalEBRechiPhi);
   fChain->SetBranchAddress("ECALnoiseECalEBRechiSwissCross", &ECALnoiseECalEBRechiSwissCross, &b_ECALnoiseECalEBRechiSwissCross);
   fChain->SetBranchAddress("ECALnoiseECalEBRechiTime", &ECALnoiseECalEBRechiTime, &b_ECALnoiseECalEBRechiTime);
   fChain->SetBranchAddress("ECALnoiseECalEBRechiiEta", &ECALnoiseECalEBRechiiEta, &b_ECALnoiseECalEBRechiiEta);
   fChain->SetBranchAddress("ECALnoiseECalEBRechiiPhi", &ECALnoiseECalEBRechiiPhi, &b_ECALnoiseECalEBRechiiPhi);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedChi2", &ECALnoiseECalEBSeedChi2, &b_ECALnoiseECalEBSeedChi2);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedE3x3", &ECALnoiseECalEBSeedE3x3, &b_ECALnoiseECalEBSeedE3x3);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedEBottom", &ECALnoiseECalEBSeedEBottom, &b_ECALnoiseECalEBSeedEBottom);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedELeft", &ECALnoiseECalEBSeedELeft, &b_ECALnoiseECalEBSeedELeft);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedERight", &ECALnoiseECalEBSeedERight, &b_ECALnoiseECalEBSeedERight);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedETop", &ECALnoiseECalEBSeedETop, &b_ECALnoiseECalEBSeedETop);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedEnergy", &ECALnoiseECalEBSeedEnergy, &b_ECALnoiseECalEBSeedEnergy);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedEta", &ECALnoiseECalEBSeedEta, &b_ECALnoiseECalEBSeedEta);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedPhi", &ECALnoiseECalEBSeedPhi, &b_ECALnoiseECalEBSeedPhi);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedTime", &ECALnoiseECalEBSeedTime, &b_ECALnoiseECalEBSeedTime);
   fChain->SetBranchAddress("ECALnoiseECalEBSeediEta", &ECALnoiseECalEBSeediEta, &b_ECALnoiseECalEBSeediEta);
   fChain->SetBranchAddress("ECALnoiseECalEBSeediPhi", &ECALnoiseECalEBSeediPhi, &b_ECALnoiseECalEBSeediPhi);
   fChain->SetBranchAddress("ECALnoiseHCALenergy3x3", &ECALnoiseHCALenergy3x3, &b_ECALnoiseHCALenergy3x3);
   fChain->SetBranchAddress("ECALnoiseHCALenergy5x5", &ECALnoiseHCALenergy5x5, &b_ECALnoiseHCALenergy5x5);
   fChain->SetBranchAddress("ECALnoiseHCALenergy7x7", &ECALnoiseHCALenergy7x7, &b_ECALnoiseHCALenergy7x7);
   fChain->SetBranchAddress("ECALnoiseHCALenergy9x9", &ECALnoiseHCALenergy9x9, &b_ECALnoiseHCALenergy9x9);
   fChain->SetBranchAddress("ECALnoiseHCALenergyUp", &ECALnoiseHCALenergyUp, &b_ECALnoiseHCALenergyUp);
   fChain->SetBranchAddress("ak5PFJetNConstituents", &ak5PFJetNConstituents, &b_ak5PFJetNConstituents);
   fChain->SetBranchAddress("ak5PFJetNJets", &ak5PFJetNJets, &b_ak5PFJetNJets);
   fChain->SetBranchAddress("CaloTowersIeta", &CaloTowersIeta, &b_CaloTowersIeta);
   fChain->SetBranchAddress("CaloTowersIphi", &CaloTowersIphi, &b_CaloTowersIphi);
   fChain->SetBranchAddress("CaloTowersTowerStatusWord", &CaloTowersTowerStatusWord, &b_CaloTowersTowerStatusWord);
   fChain->SetBranchAddress("ECALnoiseECalEBSeedRecoFlag", &ECALnoiseECalEBSeedRecoFlag, &b_ECALnoiseECalEBSeedRecoFlag);
   fChain->SetBranchAddress("tracksAlgorithm", &tracksAlgorithm, &b_tracksAlgorithm);
   fChain->SetBranchAddress("tracksNumberOfValidHits", &tracksNumberOfValidHits, &b_tracksNumberOfValidHits);
   fChain->SetBranchAddress("tracksNumberOfValidPixelHits", &tracksNumberOfValidPixelHits, &b_tracksNumberOfValidPixelHits);
   fChain->SetBranchAddress("tracksNumberOfValidStripHits", &tracksNumberOfValidStripHits, &b_tracksNumberOfValidStripHits);
   fChain->SetBranchAddress("tracksQuality", &tracksQuality, &b_tracksQuality);
   fChain->SetBranchAddress("vertexNTracks", &vertexNTracks, &b_vertexNTracks);
   fChain->SetBranchAddress("vertexNTracksW5", &vertexNTracksW5, &b_vertexNTracksW5);
   fChain->SetBranchAddress("bunch", &bunch, &b_bunch);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("run", &run, &b_run);
   Notify();
}

Bool_t rootNtupleClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void rootNtupleClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t rootNtupleClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef rootNtupleClass_cxx
