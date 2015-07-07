#include "JetMETAnalysis/PromptAnalysis/interface/PromptAna_ECALspikes.h"
#include "DataFormats/EcalDigi/interface/EcalMGPASample.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"
#include "FWCore/Framework/interface/Event.h"

using namespace std;
using namespace edm;

PromptAna_ECALspikes::PromptAna_ECALspikes(const edm::ParameterSet& iConfig) 
  : inputTag(iConfig.getParameter<edm::InputTag>("InputTag"))
    , prefix  (iConfig.getParameter<std::string>  ("Prefix"  ))
    , suffix  (iConfig.getParameter<std::string>  ("Suffix"  ))
    // , ebDigisColl(iConfig.getParameter<edm::InputTag>  ("EBDigisColl"  ))
{
  produces <std::vector<float> >  ( prefix + "ECalEBSeedEta"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBSeedPhi"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBSeediEta"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBSeediPhi"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBSeedEnergy"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBSeedE3x3"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBSeedERight"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBSeedELeft"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBSeedETop"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBSeedEBottom"  + suffix );

  produces <std::vector<float> >  ( prefix + "ECalEBRechiEnergy"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBRechiTime"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBRechiSwissCross"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBRechiiEta"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBRechiiPhi"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBRechiEta"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBRechiPhi"  + suffix );

  // produces <std::vector<float> >  ( prefix + "ECalEBDigiEnergy"  + suffix );
  // produces <std::vector<int> >    ( prefix + "ECalEBDigiiEta"  + suffix );
  // produces <std::vector<int> >    ( prefix + "ECalEBDigiiPhi"  + suffix );

  produces <std::vector<float> >  ( prefix + "ECalEBSeedChi2"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECalEBSeedTime"  + suffix );
  produces <std::vector<int> >    ( prefix + "ECalEBSeedRecoFlag"  + suffix );
  //  produces <std::vector<float> >  ( prefix + "ECalEBUnCalSeedChi2"  + suffix );

  produces <std::vector<float> >  ( prefix + "HCALenergyUp"   + suffix );
  produces <std::vector<float> >  ( prefix + "HCALenergy3x3"  + suffix );
  produces <std::vector<float> >  ( prefix + "HCALenergy5x5"  + suffix );
  produces <std::vector<float> >  ( prefix + "HCALenergy7x7"  + suffix );
  produces <std::vector<float> >  ( prefix + "HCALenergy9x9"  + suffix );

  produces <std::vector<float> >  ( prefix + "ECALenergy3x3"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECALenergy5x5"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECALenergy7x7"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECALenergy9x9"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECALenergy11x11"  + suffix );
  produces <std::vector<float> >  ( prefix + "ECALenergy13x13"  + suffix );
}

void PromptAna_ECALspikes::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  std::auto_ptr<std::vector<float> >   ecalebseedeta      ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebseedphi      ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebseedieta      ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebseediphi      ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebseedenergy   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebseeden3x3    ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebseedenright  ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebseedenleft   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebseedentop    ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebseedenbottom ( new std::vector<float>()  ) ;

  std::auto_ptr<std::vector<float> >   ecalebrechitenergy ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebrechittime   ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebrechitswissX ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebrechitieta ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebrechitiphi ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebrechiteta ( new std::vector<float>()  ) ;
  std::auto_ptr<std::vector<float> >   ecalebrechitphi ( new std::vector<float>()  ) ;

  // std::auto_ptr<std::vector<float> >   ecalebdigienergy ( new std::vector<float>()  ) ;
  // std::auto_ptr<std::vector<int> >     ecalebdigiieta   ( new std::vector<int>()  ) ;
  // std::auto_ptr<std::vector<int> >     ecalebdigiiphi   ( new std::vector<int>()  ) ;

  std::auto_ptr<std::vector<float> >   ecalebseedchi2      ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<float> >   ecalebseedtime      ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<int> >     ecalebseedflag      ( new std::vector<int>()   ) ;
  //  std::auto_ptr<std::vector<float> >   ecalebuncalseedchi2 ( new std::vector<float>() ) ;

  std::auto_ptr<std::vector<float> >   hcalenergyup        ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<float> >   hcalenergy3x3       ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<float> >   hcalenergy5x5       ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<float> >   hcalenergy7x7       ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<float> >   hcalenergy9x9       ( new std::vector<float>() ) ;

  std::auto_ptr<std::vector<float> >   ecalenergy3x3       ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<float> >   ecalenergy5x5       ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<float> >   ecalenergy7x7       ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<float> >   ecalenergy9x9       ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<float> >   ecalenergy11x11     ( new std::vector<float>() ) ;
  std::auto_ptr<std::vector<float> >   ecalenergy13x13     ( new std::vector<float>() ) ;

  //Get the ECAL SuperClusters
  edm::Handle<reco::SuperClusterCollection> pHybridSuperClusters;
  iEvent.getByLabel("hybridSuperClusters", pHybridSuperClusters);
  //iEvent.getByLabel("correctedHybridSuperClusters", pHybridSuperClusters);
  const reco::SuperClusterCollection* hybridSuperClusters = pHybridSuperClusters.product();
  
  //Get ECAL Rechit collections
  edm::Handle<EcalRecHitCollection> barrelEcalRecHitsH;
   // edm::Handle<EcalUncalibratedRecHitCollection> barrelEcalUnCalRecHitsH;

  iEvent.getByLabel("ecalRecHit","EcalRecHitsEB", barrelEcalRecHitsH);  
  //  iEvent.getByLabel("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEB", barrelEcalUnCalRecHitsH);
  
  allEBRecHits = barrelEcalRecHitsH.product();
  //  allEBUnCalRecHits = barrelEcalUnCalRecHitsH.product();
  
  //Create a CaloRecHitMetaCollection
  edm::Handle< HBHERecHitCollection > HcalRechitsH ;
  iEvent.getByLabel("hbhereco","",HcalRechitsH);

  allHcalRecHits = HcalRechitsH.product();

  // get the  calo topology  from the event setup:
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  const CaloTopology *topology = pTopology.product();
  
  // get ECAL barrel topology
  edm::ESHandle<CaloTopology> theCaloTopology;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopology); 
  const CaloSubdetectorTopology* theEBTopology   = theCaloTopology->getSubdetectorTopology(DetId::Ecal,EcalBarrel);

  // get HcalTopology
  edm::ESHandle<HcalTopology> htopo;
  iSetup.get<IdealGeometryRecord>().get(htopo);
  const HcalTopology* theHBHETopology = htopo.product();

  //Get CaloGeometry
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloSubdetectorGeometry *geometry_eb = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  const CaloSubdetectorGeometry *geometry_hcal = geoHandle->getSubdetectorGeometry(DetId::Hcal, 4);

  // Get APD digis for ECAL spike studies
  // edm::Handle<EBDigiCollection>  EBdigis;
  // iEvent.getByLabel(ebDigisColl, EBdigis);

  // Barrel
  // Get calorimetry
  edm::ESHandle<CaloGeometry> calo;
  iSetup.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry * theCaloGeometry;
  theCaloGeometry = (const CaloGeometry*)calo.product();

  const CaloSubdetectorGeometry * geom;
    
  geom = theCaloGeometry->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
  for(EBRecHitCollection::const_iterator recHit = allEBRecHits->begin();
      recHit!= allEBRecHits->end();
      recHit++)
    {
      EBDetId detId(recHit->id());
      const CaloCellGeometry* cell = geom->getGeometry(detId);

      if (recHit->energy()>4.)
	if(cell != 0) 
	  {
	    CaloNavigator<DetId> cursorE = CaloNavigator<DetId>(recHit->detid(), theEBTopology );

	    float s4 = 0;
	    float swissX = 0.;

	    cursorE.home();
	    float e1 = recHitEnergyECAL( *cursorE, allEBRecHits );

	    cursorE.offsetBy( 1, 0 );
	    s4 += recHitEnergyECAL( *cursorE, allEBRecHits ) ;
	    
	    cursorE.home();
	    cursorE.offsetBy( -1, 0 );
	    s4 += recHitEnergyECAL( *cursorE, allEBRecHits ) ;

	    cursorE.home();
	    cursorE.offsetBy( 0, 1 );	   
	    s4 += recHitEnergyECAL( *cursorE, allEBRecHits ) ;
	    
	    cursorE.home();
	    cursorE.offsetBy( 0, -1 );
	    s4 += recHitEnergyECAL( *cursorE, allEBRecHits ) ;

	    if ( e1 == 0 ) 
	      swissX = 0;
	    else 
	      swissX = 1 - s4 / e1;
	    
	    ecalebrechitenergy ->push_back ( recHit->energy() );
	    ecalebrechittime   ->push_back ( recHit->time() );
	    ecalebrechitswissX ->push_back ( swissX );
	    ecalebrechitieta   ->push_back ( detId.ieta() );
	    ecalebrechitiphi   ->push_back ( detId.iphi() );
	    	    
	    GlobalPoint posi = cell->getPosition();
	    
	    ///////////////////////////////////////////////////////
	    // std::map<int, double, std::less<int> > gainConv_;
	    // double barrelADCtoGeV_;
	    // double endcapADCtoGeV_;
	    // gainConv_[1] = 1.;
	    // gainConv_[2] = 2.;
	    // gainConv_[3] = 12.;
	    // gainConv_[0] = 12.;   // saturated channels
	    // barrelADCtoGeV_ = 0.035;
	    // endcapADCtoGeV_ = 0.06;
	        
	    // std::vector<double> ebAnalogSignal ;
	    // std::vector<double> ebADCCounts ;
	    // std::vector<double> ebADCGains ;
	    // ebAnalogSignal.reserve(EBDataFrame::MAXSAMPLES);
	    // ebADCCounts.reserve(EBDataFrame::MAXSAMPLES);
	    // ebADCGains.reserve(EBDataFrame::MAXSAMPLES);
	    
	    // float denergy = -9999.;
	    // int dieta     = -9999;
	    // int diphi     = -9999;
	    
	    // for ( EBDigiCollection::const_iterator digiItr= EBdigis->begin(); digiItr != EBdigis->end(); ++digiItr )   // Loop on EB crystals
	      // {		    
	      // 	EBDataFrame ebdf( *digiItr );		
	      // 	EBDetId ebid = ebdf.id();
		
	      // 	if(ebid.iphi()!=detId.iphi() || ebid.ieta()!=detId.ieta())
	      // 	  continue;
		
	      // 	double Emax = 0. ;
	      // 	int Pmax = 0 ;
	      // 	double pedestalPreSample = 0.;
	      // 	double pedestalPreSampleAnalog = 0.;
		
	      // 	for (int sample = 0 ; sample < ebdf.size(); ++sample) 
	      // 	  {
	      // 	    ebAnalogSignal[sample] = 0.;
	      // 	    ebADCCounts[sample] = 0.;
	      // 	    ebADCGains[sample] = -1.;
	      // 	  }
		
	      // 	for (int sample = 0 ; sample < ebdf.size(); ++sample) 
	      // 	  {	     
	      // 	    EcalMGPASample mySample = ebdf[sample];
		    
	      // 	    ebADCCounts[sample] = (mySample.adc()) ;
	      // 	    ebADCGains[sample]  = (mySample.gainId()) ;
	      // 	    ebAnalogSignal[sample] = (ebADCCounts[sample]*gainConv_[(int)ebADCGains[sample]]*barrelADCtoGeV_);
	      // 	    if (Emax < ebAnalogSignal[sample] ) 
	      // 	      {
	      // 		Emax = ebAnalogSignal[sample] ;
	      // 		Pmax = sample ;
	      // 	      }
	      // 	    if ( sample < 3 ) 
	      // 	      {
	      // 		pedestalPreSample += ebADCCounts[sample] ;
	      // 		pedestalPreSampleAnalog += ebADCCounts[sample]*gainConv_[(int)ebADCGains[sample]]*barrelADCtoGeV_ ;
	      // 	      }
	      // 	  }
		
	      // 	pedestalPreSample /= 3. ; 
	      // 	pedestalPreSampleAnalog /= 3. ; 
	      // 	double Erec = Emax - pedestalPreSampleAnalog*gainConv_[(int)ebADCGains[Pmax]];
		
	      // 	denergy = Erec;
	      // 	dieta   =  ebid.ieta();
	      // 	diphi   =  ebid.iphi();
	      // }
	    ///////////////////////////////////////////////////////	    
	    // ecalebdigienergy  ->push_back ( denergy );
	    // ecalebdigiieta    ->push_back ( dieta );
	    // ecalebdigiiphi    ->push_back ( diphi );
	    
	    ecalebrechiteta   ->push_back ( posi.eta() );
	    ecalebrechitphi   ->push_back ( posi.phi() );
	  }
    }

  // loop over the super clusters in EB
  for(reco::SuperClusterCollection::const_iterator aClus = hybridSuperClusters->begin();
      aClus != hybridSuperClusters->end(); aClus++)
    {
      //Get vector of rechit DetIDs
      std::vector< std::pair<DetId, float> > v_id = aClus->hitsAndFractions();

      //Find the DetID of the most energetic crystal
      DetId seedId; 
      float emax = -1.; 
      for ( size_t i = 0; i < v_id.size(); ++i ) 
	{	  
	  float en = recHitEnergyECAL( v_id[i].first, allEBRecHits);
	  if( emax < en)  { emax = en; seedId = v_id[i].first; }
	}
      
      //Get (Eta,Phi) of the seed crystal
      const CaloCellGeometry *this_cell = geometry_eb->getGeometry(seedId);
      GlobalPoint posi = this_cell->getPosition();
      float seedeta    =  posi.eta();
      float seedphi    = posi.phi();

      EBDetId ebId( seedId );
      float seedieta    =  ebId.ieta();
      float seediphi    = ebId.iphi();
      
      //Get Chi2, time etc.
      float seedchi2      = 0.;
      float seedtime      = 0.;
      int seedflag        = 0.;
      //      float uncalseedchi2 = 0.;

      EcalRecHitCollection::const_iterator it = allEBRecHits->find( seedId );
      //      EcalUncalibratedRecHitCollection::const_iterator itUnCal = allEBUnCalRecHits->find( seedId );
      
      if ( it != allEBRecHits->end() ) 
	{
	  seedchi2 = (*it).chi2();
	  seedtime = (*it).time();
	  seedflag = (int)(*it).recoFlag();
	}
      
//       if( itUnCal != allEBUnCalRecHits->end() )
// 	{
// 	  uncalseedchi2    = (*itUnCal).chi2();
// 	}
      
      //////////////////////////////////////////////
      /////ECAL energies around ECAL seeds//////////
      //////////////////////////////////////////////

      float ECALenergy3x3   = 0.;
      float ECALenergy5x5   = 0.;
      float ECALenergy7x7   = 0.;
      float ECALenergy9x9   = 0.;
      float ECALenergy11x11 = 0.;
      float ECALenergy13x13 = 0.;
      
      //Sum up energy around spike crystal in ECAL
      CaloNavigator<DetId> cursorE = CaloNavigator<DetId>(seedId, theEBTopology );
      
      //In a 3x3 matrix around the spike
      for ( int ii = -1; ii <= 1; ++ii ) for ( int jj = -1; jj <= 1; ++jj ) 
	  {
	    cursorE.home();
	    cursorE.offsetBy( ii, jj );
	    if( recHitEnergyECAL( *cursorE, allEBRecHits ) > 0.3 )
	      ECALenergy3x3 += recHitEnergyECAL( *cursorE, allEBRecHits );
	  }

      //In a 5x5 matrix around the spike
      for ( int ii = -2; ii <= 2; ++ii ) for ( int jj = -2; jj <= 2; ++jj ) 
	  {
	    cursorE.home();
	    cursorE.offsetBy( ii, jj );
	    if( recHitEnergyECAL( *cursorE, allEBRecHits ) > 0.3 )
	      ECALenergy5x5 += recHitEnergyECAL( *cursorE, allEBRecHits );
	  }

      //In a 7x7 matrix around the spike
      for ( int ii = -3; ii <= 3; ++ii ) for ( int jj = -3; jj <= 3; ++jj ) 
	  {
	    cursorE.home();
	    cursorE.offsetBy( ii, jj );
	    if( recHitEnergyECAL( *cursorE, allEBRecHits ) > 0.3 )
	      ECALenergy7x7 += recHitEnergyECAL( *cursorE, allEBRecHits );
	  }

      //In a 9x9 matrix around the spike
      for ( int ii = -4; ii <= 4; ++ii ) for ( int jj = -4; jj <= 4; ++jj ) 
	  {
	    cursorE.home();
	    cursorE.offsetBy( ii, jj );
	    if( recHitEnergyECAL( *cursorE, allEBRecHits ) > 0.3 )
	      ECALenergy9x9 += recHitEnergyECAL( *cursorE, allEBRecHits );
	  }

      //In a 11x11 matrix around the spike
      for ( int ii = -5; ii <= 5; ++ii ) for ( int jj = -5; jj <= 5; ++jj ) 
	  {
	    cursorE.home();
	    cursorE.offsetBy( ii, jj );
	    if( recHitEnergyECAL( *cursorE, allEBRecHits ) > 0.3 )
	      ECALenergy11x11 += recHitEnergyECAL( *cursorE, allEBRecHits );
	  }

      //In a 13x13 matrix around the spike
      for ( int ii = -6; ii <= 6; ++ii ) for ( int jj = -6; jj <= 6; ++jj ) 
	  {
	    cursorE.home();
	    cursorE.offsetBy( ii, jj );
	    if( recHitEnergyECAL( *cursorE, allEBRecHits ) > 0.3 )
	      ECALenergy13x13 += recHitEnergyECAL( *cursorE, allEBRecHits );
	  }
      
      
      /////////////////////////////////////////////////
      ////////HCAL energies around ECAL seeds//////////
      /////////////////////////////////////////////////

      DetId hcalDetId ;
      hcalDetId = geometry_hcal->getClosestCell(posi);
      
      float HCALenergyUp  = 0.;
      float HCALenergy3x3 = 0.;
      float HCALenergy5x5 = 0.;
      float HCALenergy7x7 = 0.;
      float HCALenergy9x9 = 0.;

      //energy in HCAL tower above the spike
      HBHERecHitCollection::const_iterator it1 = allHcalRecHits->find( hcalDetId );
      if ( it1 != allHcalRecHits->end() )
	if( (*it1).energy()>1.0 )
	  HCALenergyUp = (*it1).energy();

      //Sum up energy around spike crystal in ECAL
      CaloNavigator<DetId> cursor = CaloNavigator<DetId>(hcalDetId, theHBHETopology );
      
      //In a 3x3 matrix around the spike
      for ( int ii = -1; ii <= 1; ++ii ) for ( int jj = -1; jj <= 1; ++jj ) 
	  {
	    cursor.home();
	    cursor.offsetBy( ii, jj );
	    if( recHitEnergyHCAL( *cursor, allHcalRecHits ) > 1.0 )
	      HCALenergy3x3 += recHitEnergyHCAL( *cursor, allHcalRecHits );
	  }

      //In a 5x5 matrix around the spike
      for ( int ii = -2; ii <= 2; ++ii ) for ( int jj = -2; jj <= 2; ++jj ) 
	  {
	    cursor.home();
	    cursor.offsetBy( ii, jj );
	    if( recHitEnergyHCAL( *cursor, allHcalRecHits ) > 1.0 )
	      HCALenergy5x5 += recHitEnergyHCAL( *cursor, allHcalRecHits );
	  }

      //In a 7x7 matrix around the spike
      for ( int ii = -3; ii <= 3; ++ii ) for ( int jj = -3; jj <= 3; ++jj ) 
	  {
	    cursor.home();
	    cursor.offsetBy( ii, jj );
	    if( recHitEnergyHCAL( *cursor, allHcalRecHits ) > 1.0 )
	      HCALenergy7x7 += recHitEnergyHCAL( *cursor, allHcalRecHits );
	  }

      //In a 9x9 matrix around the spike
      for ( int ii = -4; ii <= 4; ++ii ) for ( int jj = -4; jj <= 4; ++jj ) 
	  {
	    cursor.home();
	    cursor.offsetBy( ii, jj );
	    if( recHitEnergyHCAL( *cursor, allHcalRecHits ) > 1.0 )
	      HCALenergy9x9 += recHitEnergyHCAL( *cursor, allHcalRecHits );
	  }

      //store the values for each seeds
      ecalebseedeta       ->push_back ( seedeta );
      ecalebseedphi       ->push_back ( seedphi );
      ecalebseedieta       ->push_back ( seedieta );
      ecalebseediphi       ->push_back ( seediphi );
      ecalebseedenergy    ->push_back ( emax );

      ecalebseedchi2      ->push_back ( seedchi2 );
      ecalebseedtime      ->push_back ( seedtime );
      ecalebseedflag      ->push_back ( seedflag );
      //      ecalebuncalseedchi2 ->push_back ( uncalseedchi2 );

      //store the values for each super cluster, to identify spikes
      ecalebseeden3x3    ->push_back ( EcalClusterTools::e3x3 ( *aClus , allEBRecHits, &(*topology) ) );
      ecalebseedenright  ->push_back ( EcalClusterTools::eRight ( *aClus , allEBRecHits, &(*topology) ) );
      ecalebseedenleft   ->push_back ( EcalClusterTools::eLeft ( *aClus , allEBRecHits, &(*topology) ) );
      ecalebseedentop    ->push_back ( EcalClusterTools::eTop ( *aClus , allEBRecHits, &(*topology) ) );
      ecalebseedenbottom ->push_back ( EcalClusterTools::eBottom ( *aClus , allEBRecHits, &(*topology) ) );
      
      //HCAL energies
      hcalenergyup        ->push_back ( HCALenergyUp );
      hcalenergy3x3       ->push_back ( HCALenergy3x3);
      hcalenergy5x5       ->push_back ( HCALenergy5x5);
      hcalenergy7x7       ->push_back ( HCALenergy7x7);
      hcalenergy9x9       ->push_back ( HCALenergy9x9);
      
      ecalenergy3x3       ->push_back ( ECALenergy3x3);
      ecalenergy5x5       ->push_back ( ECALenergy5x5);
      ecalenergy7x7       ->push_back ( ECALenergy7x7);
      ecalenergy9x9       ->push_back ( ECALenergy9x9);
      ecalenergy11x11     ->push_back ( ECALenergy11x11);
      ecalenergy13x13     ->push_back ( ECALenergy13x13);
    }
  
  iEvent.put( ecalebseedeta          ,  prefix + "ECalEBSeedEta"  + suffix );
  iEvent.put( ecalebseedphi          ,  prefix + "ECalEBSeedPhi"  + suffix );
  iEvent.put( ecalebseedieta          ,  prefix + "ECalEBSeediEta"  + suffix );
  iEvent.put( ecalebseediphi          ,  prefix + "ECalEBSeediPhi"  + suffix );
  iEvent.put( ecalebseedenergy       ,  prefix + "ECalEBSeedEnergy"  + suffix );

  iEvent.put( ecalebseeden3x3        ,  prefix + "ECalEBSeedE3x3"  + suffix );
  iEvent.put( ecalebseedenright      ,  prefix + "ECalEBSeedERight"  + suffix );
  iEvent.put( ecalebseedenleft       ,  prefix + "ECalEBSeedELeft"  + suffix );
  iEvent.put( ecalebseedentop        ,  prefix + "ECalEBSeedETop"  + suffix );
  iEvent.put( ecalebseedenbottom     ,  prefix + "ECalEBSeedEBottom"  + suffix );

  iEvent.put( ecalebrechitenergy         ,  prefix + "ECalEBRechiEnergy"  + suffix );
  iEvent.put( ecalebrechittime          ,  prefix + "ECalEBRechiTime"  + suffix );
  iEvent.put( ecalebrechitswissX         ,  prefix + "ECalEBRechiSwissCross"  + suffix );
  iEvent.put( ecalebrechitieta         ,  prefix + "ECalEBRechiiEta"  + suffix );
  iEvent.put( ecalebrechitiphi         ,  prefix + "ECalEBRechiiPhi"  + suffix );
  iEvent.put( ecalebrechiteta         ,  prefix + "ECalEBRechiEta"  + suffix );
  iEvent.put( ecalebrechitphi         ,  prefix + "ECalEBRechiPhi"  + suffix );

  // iEvent.put( ecalebdigienergy         ,  prefix + "ECalEBDigiEnergy"  + suffix );
  // iEvent.put( ecalebdigiieta           ,  prefix + "ECalEBDigiiEta"  + suffix );
  // iEvent.put( ecalebdigiiphi           ,  prefix + "ECalEBDigiiPhi"  + suffix );

  iEvent.put( ecalebseedchi2         ,  prefix + "ECalEBSeedChi2"  + suffix );
  iEvent.put( ecalebseedtime         ,  prefix + "ECalEBSeedTime"  + suffix );
  iEvent.put( ecalebseedflag         ,  prefix + "ECalEBSeedRecoFlag"  + suffix );
  //  iEvent.put( ecalebuncalseedchi2    ,  prefix + "ECalEBUnCalSeedChi2"  + suffix );

  iEvent.put( hcalenergyup           ,  prefix + "HCALenergyUp"  + suffix );
  iEvent.put( hcalenergy3x3          ,  prefix + "HCALenergy3x3"  + suffix );
  iEvent.put( hcalenergy5x5          ,  prefix + "HCALenergy5x5"  + suffix );
  iEvent.put( hcalenergy7x7          ,  prefix + "HCALenergy7x7"  + suffix );
  iEvent.put( hcalenergy9x9          ,  prefix + "HCALenergy9x9"  + suffix );

  iEvent.put( ecalenergy3x3          ,  prefix + "ECALenergy3x3"  + suffix );
  iEvent.put( ecalenergy5x5          ,  prefix + "ECALenergy5x5"  + suffix );
  iEvent.put( ecalenergy7x7          ,  prefix + "ECALenergy7x7"  + suffix );
  iEvent.put( ecalenergy9x9          ,  prefix + "ECALenergy9x9"  + suffix );
  iEvent.put( ecalenergy11x11        ,  prefix + "ECALenergy11x11"  + suffix );
  iEvent.put( ecalenergy13x13        ,  prefix + "ECALenergy13x13"  + suffix );
}

float PromptAna_ECALspikes::recHitEnergyECAL(DetId id, const EcalRecHitCollection *recHits)
{
  if ( id == DetId(0) ) {
    return 0;
  } else {
    EcalRecHitCollection::const_iterator it = recHits->find( id );
    if ( it != recHits->end() ) {
      return (*it).energy();
     } else {
       return 0;
     }
   }
   return 0;
 }

float PromptAna_ECALspikes::recHitEnergyHCAL(DetId id, const HBHERecHitCollection *recHits)
{
  if ( id == DetId(0) ) {
    return 0;
  } else {
    HBHERecHitCollection::const_iterator it = recHits->find( id );
    if ( it != recHits->end() ) {
      return (*it).energy();
     } else {
       return 0;
     }
   }
   return 0;
 }
