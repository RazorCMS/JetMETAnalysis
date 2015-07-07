#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TLine.h>

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{
  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}

double DeltaPhi(double jetPhi1,double jetPhi2){
  double deltaphi=fabs(jetPhi1-jetPhi2);
  if(deltaphi>M_PI){
    deltaphi=2*M_PI-deltaphi;
  }
  return deltaphi;
}

//first proposal - give only index
bool JetIdloose(double ak5JetJIDresEMF,double ak5JetJIDfHPD,int ak5JetJIDn90Hits, double ak5JetEta){
  bool jetidresEMF=true;
  bool jetidloose=false;
  double fhpdmax = 0.98;
  double n90hitsmin =1;
  double emf_min = 0.01;

  if(fabs(ak5JetEta)<2.6 && ak5JetJIDresEMF<emf_min) jetidresEMF=false;
  if(jetidresEMF && ak5JetJIDfHPD<fhpdmax && ak5JetJIDn90Hits>n90hitsmin) {
    jetidloose=true;
  }
  return jetidloose;
}

//originally floats in the PFJet class, but filled as double in the ntuples
//bool PFJetIDloose(float CHF, float NHF, float CEF, float NEF,double eta){
bool pass_PFJetIDloose(double CHF, double NHF, double CEF, double NEF, double NCJP, int NJC, double eta){
  bool pass_result = false;
  bool jetidCHF=true;

  if(fabs(eta)<2.4 && CHF<=0.0 && NJC <= 1 && NCJP <= 0.){jetidCHF=false;}
  if( jetidCHF &&(NHF<1.0)&&(CEF<1.00)&&(NEF<1.00) )
    {
      pass_result = true;
    }
  
  if(fabs(eta)>2.4)
    if ( NHF<1.0 && NEF<1.0 && NCJP > 1.)
      pass_result = true;
  
  return pass_result;
}

void analysisClass::Loop()
{
  std::cout << "analysisClass::Loop() begins" <<std::endl; 

  if (fChain == 0) return;

  //////////book histos here
  int   Nbins_METSumET = 400;
  float Max_METSumET = 400;

  int   Nbins_SumET = 500;
  float Max_SumET = 500;

  int   Nbins_Phi = 50;
  float Max_Phi = 3.15;


  // Vertex histograms
  TH1F *h_Vx   = new TH1F ("h_Vx","h_Vx",600, -15, 15);
  TH1F *h_Vy   = new TH1F ("h_Vy","h_Vy",600, -15, 15);
  
  h_Vx->Sumw2();
  h_Vy->Sumw2();
  
  // calo MET and calo MET default
  TH1F *h_calometPt   = new TH1F ("h_calometPt","h_calometPt",Nbins_METSumET, 0, Max_METSumET);
  TH1F *h_caloSumet   = new TH1F ("h_caloSumet","h_caloSumet",Nbins_SumET, 0, Max_SumET);
  TH1F *h_calometPtDflt   = new TH1F ("h_calometPtDflt","h_calometPtDflt",Nbins_METSumET,0,Max_METSumET);

  TH1F *h_calometPxy   = new TH1F ("h_calometPxy","h_calometPxy",Nbins_METSumET/4,-Max_METSumET/8,Max_METSumET/8);
  TH1F *h_calometPx   = new TH1F ("h_calometPx","h_calometPx",Nbins_METSumET,-Max_METSumET/2,Max_METSumET/2);
  TH1F *h_calometPy   = new TH1F ("h_calometPy","h_calometPy",Nbins_METSumET,-Max_METSumET/2,Max_METSumET/2);

  h_calometPt->Sumw2();
  h_calometPtDflt->Sumw2();
  h_calometPxy->Sumw2();
  h_calometPx->Sumw2();
  h_calometPy->Sumw2();
  h_caloSumet->Sumw2();

  TH2F *h2_DMET_VS_METDflt = new TH2F ("h2_DMET_VS_METDflt","h2_DMET_VS_METDflt", Nbins_METSumET,0,Max_METSumET, 200, -100., 100.);

  //MPT
  TH1F *h_nGoodTracks = new TH1F ("h_nGoodTracks","h_nGoodTracks",500,0,500);
  TH1F *h_nValidHits = new TH1F ("h_nValidHits","h_nValidHits",100,0,100);

  TH1F *h_trackPt     = new TH1F ("h_trackPt","h_trackPt",2*Nbins_METSumET,0,Max_METSumET);
  TH1F *h_MPT   = new TH1F ("h_MPT","h_MPT",Nbins_METSumET,0,Max_METSumET);
  TH1F *h_SumPT   = new TH1F ("h_SumPT","h_SumPT",Nbins_METSumET,0,Max_METSumET);
  TH1F *h_MPTPhi  = new TH1F ("h_MPTPhi","h_MPTPhi",Nbins_Phi,0,2*Max_Phi);
  TH1F *h_MPx   = new TH1F ("h_MPx","h_MPx",Nbins_METSumET,-Max_METSumET/2,Max_METSumET/2);
  TH1F *h_MPy   = new TH1F ("h_MPy","h_MPy",Nbins_METSumET,-Max_METSumET/2,Max_METSumET/2);

  h_nGoodTracks->Sumw2();
  h_nValidHits->Sumw2();
  h_MPT->Sumw2();
  h_trackPt->Sumw2();
  h_SumPT->Sumw2();
  h_MPTPhi->Sumw2();
  h_MPx->Sumw2(); 
  h_MPy->Sumw2(); 

  // ECAL spikes plots
  TH1F *h_S4OS1_EBEE   = new TH1F ("h_S4OS1_EBEE","h_S4OS1_EBEE",200,-0.5,1.5);
  h_S4OS1_EBEE->Sumw2();
  TH2F *h2_ECalSeedET_Vs_S4_calometoverS1= new TH2F ("h2_ECalSeedET_Vs_S4_calometoverS1","h2_ECalSeedET_Vs_S4_calometoverS1", 1000,-0.5,1.5 , 1000, 0., 1000.);
  TH2F *h2_ECalSeedET_Vs_S4_calomet = new TH2F ("h2_ECalSeedET_Vs_S4_calomet","h2_ECalSeedET_Vs_S4_calomet", 1000,0,1000 , 1000, 0., 1000.);
  TH2F *h2_ECalSpikePhi_calometPhi = new TH2F ("h2_ECalSpikePhi_calometPhi","h2_ECalSpikePhi_calometPhi", Nbins_Phi,-Max_Phi,Max_Phi, Nbins_Phi,-Max_Phi,Max_Phi);
  TH2F *h2_ECalSpikePhi_pfmetPhi = new TH2F ("h2_ECalSpikePhi_pfmetPhi","h2_ECalSpikePhi_pfmetPhi", Nbins_Phi,-Max_Phi,Max_Phi, Nbins_Phi,-Max_Phi,Max_Phi);
  TH2F *h2_ECalSeedET_Vs_S4_pfmet = new TH2F ("h2_ECalSeedET_Vs_S4_pfmet","h2_ECalSeedET_Vs_S4_pfmet", 1000,0,1000 , 1000, 0., 1000.);
  TH2F *h2_ECalSeedET_Vs_S4_calomet_098 = new TH2F ("h2_ECalSeedET_Vs_S4_calomet_098","h2_ECalSeedET_Vs_S4_calomet_098", 1000,0,1000 , 1000, 0., 1000.);
  TH2F *h2_ECalSeedET_Vs_S4_calomet_0999 = new TH2F ("h2_ECalSeedET_Vs_S4_calomet_0999","h2_ECalSeedET_Vs_S4_calomet_0999", 1000,0,1000 , 1000, 0., 1000.);

  TH2F *h2_ECalSeedET_Vs_calometClean = new TH2F ("h2_ECalSeedET_Vs_calometClean","h2_ECalSeedET_Vs_calometClean", 1000,0,1000 , 1000, 0., 1000.);
  TH2F *h2_calomet_Vs_calometClean = new TH2F ("h2_calomet_Vs_calometClean","h2_calomet_Vs_calometClean; Cleaned MET [GeV]; Default MET [GeV]", 1000,0,1000 , 1000, 0., 1000.);

  //calomet in dijets (loose)
  TH1F *h_dijetLoose_jet1Eta   = new TH1F ("h_dijetLoose_jet1Eta","h_dijetLoose_jet1Eta",100, -3.5,3.5);
  TH1F *h_dijetLoose_jet2Eta   = new TH1F ("h_dijetLoose_jet2Eta","h_dijetLoose_jet2Eta",100, -3.5,3.5);

  TH1F *h_dijetLoose_calometPt   = new TH1F ("h_dijetLoose_calometPt","h_dijetLoose_calometPt",Nbins_METSumET/4, 0, Max_METSumET/4);
  TH1F *h_dijetLoose_calometPxy   = new TH1F ("h_dijetLoose_calometPxy","h_dijetLoose_calometPxy",Nbins_METSumET/4, -Max_METSumET/8, Max_METSumET/8);
  TH1F *h_dijetLoose_caloSumet   = new TH1F ("h_dijetLoose_caloSumet","h_dijetLoose_caloSumet",Nbins_SumET, 0, Max_SumET);

  h_dijetLoose_jet1Eta->Sumw2();
  h_dijetLoose_jet2Eta->Sumw2();

  h_dijetLoose_calometPt->Sumw2();
  h_dijetLoose_calometPxy->Sumw2();
  h_dijetLoose_caloSumet->Sumw2();

  //Vertex
  TH1F *h_AllVertexZ    = new TH1F ("h_AllVertexZ","h_AllVertexZ",100,-100,100);
  TH1F *h_AllVertexChi2 = new TH1F ("h_AllVertexChi2","h_AllVertexChi",100,0,100);
  TH1F *h_AllVertexNDOF = new TH1F ("h_AllVertexNDOF","h_AllVertexNDOF",50,0,50);
  TH1F *h_AllVertexChi2_0_NDOF = new TH1F ("h_AllVertexChi2_0_NDOF","h_AllVertexChi2_0_NDOF",200,0,40);
  TH1F *h_AllVertexNtrk = new TH1F ("h_AllVertexNtrk","h_AllVertexNtrk",50,0,50);
  TH1F *h_AllNVertex    = new TH1F ("h_AllNVertex","h_AllNVertex",50,0,50);

  h_AllVertexZ->Sumw2();
  h_AllVertexChi2->Sumw2(); 
  h_AllVertexNDOF->Sumw2(); 
  h_AllVertexChi2_0_NDOF->Sumw2();
  h_AllVertexNtrk->Sumw2(); 
  h_AllNVertex->Sumw2();

  /////////initialize variables

  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl; 

  Long64_t nb = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) 
    //for (Long64_t jentry=0; jentry<2000;jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      //      if(jentry>10000) break;
      nb = fChain->GetEntry(jentry); 
      
      if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl; 
      ////////////////////// User's code starts here ///////////////////////
            
      //#####################
      //## Trigger selection
      //#####################

      int pass_BPTX              = 0;
      int pass_BSC_BeamHaloVeto  = 1;
      int pass_BSC_MB                 = 0;
      int pass_HBHENoiseFilterResult  = 0;
      
      //## pass_BPTX - Two beams crossing at CMS (only Data)
      if(isData==1)
	{
	  // if(L1TechBits->at(0)==1)
	    pass_BPTX = 1;
	}
      else if(isData==0)
	pass_BPTX = 1;
      
      //## pass_BSC_BeamHaloVeto - Veto on BSC Beam Halo Triggers firing
      if(isData==1)
	{
	  // if( L1TechBits->at(36) == 1 || L1TechBits->at(37) == 1 || L1TechBits->at(38) == 1 || L1TechBits->at(39) == 1 )
	    pass_BSC_BeamHaloVeto = 1;
	}
      else if(isData == 0)
	pass_BSC_BeamHaloVeto = 1;

      //## pass_BSC_MB - BSC MinBias triggers firing (both Data and MC)
      if(isData==1)
	//  	if ( HLTResults->at(6) == 1 )
	  pass_BSC_MB = 1;
      
      if(isData==0)
	  pass_BSC_MB = 1;
      
      //HBHE noise
      if(isData==1)
	//        if(HBHENoiseFilterResult==1)
          pass_HBHENoiseFilterResult = 1;
      if(isData==0)
        pass_HBHENoiseFilterResult = 1;
      //############################
      //## Calculate Reco Quantities 
      //############################

      //pass_GoodVertex 
      //https://twiki.cern.ch/twiki/bin/viewauth/CMS/TRKPromptFeedBack#Event_and_track_selection_recipe
      int pass_GoodVertex = 1;

      if(vertexZ->size() == 0) pass_GoodVertex = 0;
      for (int ii=0; ii<vertexZ->size(); ii++)
	if( vertexChi2->at(ii) != 0. && vertexNDF->at(ii) != 0 && vertexNDF->at(ii) > 4 && fabs(vertexZ->at(ii)) <= 15. )
	  {
	    pass_GoodVertex = 1;
	    break;
	  }

      //=================================================================

      // Set the evaluation of the cuts to false and clear the variable values and filled status
      resetCuts();

      // Set the value of the variableNames listed in the cutFile to their current value
      fillVariableWithValue("pass_BPTX", pass_BPTX);
      fillVariableWithValue("pass_BSC_BeamHaloVeto", pass_BSC_BeamHaloVeto);
      fillVariableWithValue("pass_BSC_MB", pass_BSC_MB);
      fillVariableWithValue("pass_GoodVertex", pass_GoodVertex);
      fillVariableWithValue("pass_HBHENoiseFilterResult", pass_HBHENoiseFilterResult);

      // Evaluate cuts (but do not apply them)
      evaluateCuts();

      //###########################
      //## Start filling histograms
      //###########################

      if( passedCut("all") )
	{
	  //Vertex
	  h_AllNVertex->Fill(vertexZ->size());
	  for (int ii=0; ii<vertexZ->size(); ii++)
	    {
	      if(vertexNTracksW5->at(ii)==0)
		continue;

	      h_AllVertexZ->Fill(vertexZ->at(ii));
	      h_AllVertexChi2->Fill(vertexChi2->at(ii));
	      h_AllVertexNDOF->Fill(vertexNDF->at(ii));
	      h_AllVertexNtrk->Fill(vertexNTracks->at(ii));
	      if(vertexNDF->at(ii)!=0)
		h_AllVertexChi2_0_NDOF->Fill( vertexChi2->at(ii) / vertexNDF->at(ii) );
	      
	    }
	}

      if( passedCut("all") )
	{

	  for (int ii=0; ii<vertexZ->size(); ii++)
	    {
	      h_Vx->Fill(vertexX->at(ii));
	      h_Vy->Fill(vertexY->at(ii));
	    }
	  
	  for (int ii=0; ii<ECALnoiseECalEBRechiEnergy->size(); ii++)
	    {
	      if(abs(int(ECALnoiseECalEBRechiiEta->at(ii)))==85)
	  	{
	  	  h_S4OS1_EBEE->Fill(ECALnoiseECalEBRechiSwissCross->at(ii));
	  	  h2_ECalSeedET_Vs_S4_calometoverS1->Fill( ECALnoiseECalEBRechiSwissCross->at(ii), ECALnoiseECalEBRechiEnergy->at(ii) );
	  	  h2_ECalSeedET_Vs_S4_calomet->Fill(calometPt->at(0), ECALnoiseECalEBRechiEnergy->at(ii) );
	  	  h2_ECalSeedET_Vs_S4_pfmet->Fill(pfmetPt->at(0), ECALnoiseECalEBRechiEnergy->at(ii) );
		  
	  	  h2_ECalSpikePhi_calometPhi->Fill(calometPhi->at(0), ECALnoiseECalEBRechiPhi->at(ii) );
	  	  h2_ECalSpikePhi_pfmetPhi->Fill(pfmetPhi->at(0), ECALnoiseECalEBRechiPhi->at(ii) );

	  	  if( ECALnoiseECalEBRechiSwissCross->at(ii) < 0.98)
	  	    h2_ECalSeedET_Vs_S4_calomet_098->Fill(calometPt->at(0), ECALnoiseECalEBRechiEnergy->at(ii));
	  	  if( ECALnoiseECalEBRechiSwissCross->at(ii) < 0.99)
	  	    h2_ECalSeedET_Vs_S4_calomet_0999->Fill(calometPt->at(0), ECALnoiseECalEBRechiEnergy->at(ii));		  
	  	}
	    }

	 	  
	  //## Reconstructed MPT from tracks
	  int    nGoodTracks=0;
	  double MPx=0.;
	  double MPy=0.;
	  double MPT=0.;
	  double SumPT=0.;
	  double MPTPhi=0.;
	  
	  for (int ii=0; ii<tracksPt->size(); ii++)
	    {
	      int trackFlags = tracksQuality->at(ii);
	      int highPurityFlag = 3;
	      if( ( trackFlags & 1 << highPurityFlag) > 0)
		if(tracksPt->at(ii) > 1. )
		  {
		    MPx += -1.*(tracksPt->at(ii)*cos(tracksPhi->at(ii)));
		    MPy += -1.*(tracksPt->at(ii)*sin(tracksPhi->at(ii)));
		    SumPT+=tracksPt->at(ii);
		    h_trackPt->Fill( tracksPt->at(ii) );
		    h_nValidHits->Fill( tracksNumberOfValidHits->at(ii) );
		    nGoodTracks++;
		  }	  
	    }
	  
	  TVector2 *mpt2 = new TVector2(MPx,MPy);
	  MPT    = mpt2->Mod();
	  MPTPhi = mpt2->Phi();
	  
	  //MPT
	  h_nGoodTracks->Fill(nGoodTracks);	       
	  if(nGoodTracks>0)
	    {
	      h_MPT->Fill( MPT );
	      h_SumPT->Fill( SumPT);
	      h_MPTPhi->Fill( MPTPhi );
	      
	      h_MPx->Fill( MPx );
	      h_MPy->Fill( MPy );
	      
	    }

	  //#########################
	  //## inclusive MET
	  //#########################
	  
	  //CaloMET
	    {
	      h_calometPt->Fill( calometPt->at(0) );
	      
	      h_calometPxy->Fill( calometPx->at(0) );
	      h_calometPxy->Fill( calometPy->at(0) );
	      
	      h_caloSumet->Fill( calometSumEt->at(0) ); 
	      
	      h_calometPx->Fill( calometPx->at(0) ); 
	      h_calometPy->Fill( calometPy->at(0) ); 
	    }
	  
    }//-------------- passed cuts "0"
  ////////////////////// User's code ends here ///////////////////////
} // End loop over events

  //////////write histos 

  //## 1D histograms
h_Vx->Write();
h_Vy->Write();

//MPT
h_nGoodTracks->Write();
h_nValidHits->Write();
h_MPT->Write();
h_MPx->Write();
h_MPy->Write();
h_SumPT->Write();
h_MPTPhi->Write();
h_trackPt->Write();

//calomet
h_dijetLoose_jet1Eta->Write();
h_dijetLoose_jet2Eta->Write();

h_calometPt->Write(); 
h_calometPtDflt->Write();
h2_DMET_VS_METDflt->Write();

h_calometPxy->Write(); 
h_caloSumet->Write(); 

h_calometPx->Write(); 
h_calometPy->Write(); 
  
//Dijets (loose)
h_dijetLoose_calometPt->Write(); 
h_dijetLoose_calometPxy->Write(); 
h_dijetLoose_caloSumet->Write(); 
 
//Vertex
h_AllVertexZ->Write();
h_AllVertexChi2->Write(); 
h_AllVertexNDOF->Write(); 
h_AllVertexChi2_0_NDOF->Write();
h_AllVertexNtrk->Write(); 
h_AllNVertex->Write();

  h_S4OS1_EBEE->Write();
  h2_ECalSeedET_Vs_S4_calometoverS1->Write();
  h2_ECalSeedET_Vs_S4_calomet->Write();
  h2_ECalSeedET_Vs_S4_pfmet->Write();
  h2_ECalSpikePhi_calometPhi->Write();
  h2_ECalSpikePhi_pfmetPhi->Write();
  h2_ECalSeedET_Vs_S4_calomet_098->Write();
  h2_ECalSeedET_Vs_S4_calomet_0999->Write();

  h2_ECalSeedET_Vs_calometClean->Write();
  h2_calomet_Vs_calometClean->Write();

//## 2D histograms

std::cout << "analysisClass::Loop() ends" <<std::endl; 
}
