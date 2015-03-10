#include "SeedClusteringFitter.h"

SeedClusteringFitter::SeedClusteringFitter():TrackFitter(0){

}

SeedClusteringFitter::SeedClusteringFitter(int nb):TrackFitter(nb){

  //Layer infos
  m_nLayer = 6;

  m_pMeanLayerRadius = new float[m_nLayer];
  m_pMeanLayerRadius[0] = 22.59;
  m_pMeanLayerRadius[1] = 35.57;
  m_pMeanLayerRadius[2] = 50.75;
  m_pMeanLayerRadius[3] = 68.18;
  m_pMeanLayerRadius[4] = 88.61;
  m_pMeanLayerRadius[5] = 107.67;  


  //Slopes processing
  float step = (MAX_SLOPE - MIN_SLOPE) / (N_SLOPES-1);
  float currentSlope = MIN_SLOPE;
  while(currentSlope <= MAX_SLOPE){
    m_vSlope.push_back(currentSlope);
    currentSlope += step;
  }

  unsigned int bin_mask_height = int(ceil(BINMASK_PHI_MAX * BINMASK_PHI_RES));

  //Ordinate processing
  step = (BINMASK_PHI_MAX) / (bin_mask_height-1);
  float currentIntercept = 0.0;
  while(currentIntercept <= BINMASK_PHI_MAX){
    m_vIntercept.push_back(currentIntercept);
    currentIntercept += step;
  }

  //Binary mask allocation
  m_pppBinMask = new bool ** [m_vIntercept.size()];
  for (unsigned int i = 0; i < m_vIntercept.size(); i++){
    m_pppBinMask[i] = new bool * [m_nLayer];
    for (unsigned int j = 0; j < m_nLayer; j++){
      m_pppBinMask[i][j] = new bool[2];
    }
  }
//cout<<"BINMASK ALLOCATION OK, "<<"binMaskHeight = "<<m_vIntercept.size()<<endl;

  //Seed mask allocation
  m_ppSlopeSeedMask = new int * [m_vSlope.size()];
  for (unsigned int i = 0; i < m_vSlope.size(); i++){
    m_ppSlopeSeedMask[i] = new int [m_nLayer];
  }

  for (unsigned int slopeIndex = 0; slopeIndex < m_vSlope.size(); slopeIndex++){

    float slope_value = m_vSlope[slopeIndex];
    
    //Determination of the corresponding pixel for each row and for each mask layer
    for (unsigned int layerIndex = 0; layerIndex < m_nLayer; layerIndex++){
      int maskIndex = int(trunc(slope_value * m_pMeanLayerRadius[layerIndex] * BINMASK_PHI_RES));
      m_ppSlopeSeedMask[slopeIndex][layerIndex] = maskIndex;
    }
  }

//cout<<"SEEDS PRE-PROCESSING OK, "<<m_vSlope.size()<<" Seeds pre-processed"<<endl;
}

SeedClusteringFitter::~SeedClusteringFitter(){

  //BinMasj
  for (unsigned int i = 0; i < m_vIntercept.size(); i++){
    for (unsigned int j =0; j < m_nLayer; j++){
      delete m_pppBinMask[i][j];
    }
    delete m_pppBinMask[i];
  }
  delete m_pppBinMask;

  m_pppBinMask = NULL;

  for (unsigned int i = 0; i < m_vSlope.size(); i++)
    delete m_ppSlopeSeedMask[i];

  delete m_ppSlopeSeedMask;

  m_ppSlopeSeedMask = NULL;

  delete m_pMeanLayerRadius;
  
  m_pMeanLayerRadius = NULL;

}

void SeedClusteringFitter::initialize(){

}

void SeedClusteringFitter::mergePatterns(){
  //cout<<"Merging of patterns not implemented"<<endl;
}

void SeedClusteringFitter::mergeTracks(){
  unsigned int index = 0;
  vector<Track*>::iterator it = tracks.begin();
  while(it!=tracks.end()){
    Track* newTrack = *it;
    bool found = false;
    for(unsigned int i=0;i<index;i++){
      Track* ref = tracks[i];
      float dpt,dphi,dz,deta;
      dpt = fabs(newTrack->getCurve()-ref->getCurve());
      dphi = fabs(newTrack->getPhi0()-ref->getPhi0());
      dz = fabs(newTrack->getZ0()-ref->getZ0());
      deta = fabs(newTrack->getEta0()-ref->getEta0());
      found = (deta<0.02) &&
	(dphi<0.005) &&
	(dpt<0.1) &&
	(dz<0.3);
      if(found)
	break;
    }
    if(found)
      tracks.erase(it);
    else{
      index++;
      it++;
    }
  }
}

void SeedClusteringFitter::LinearLeastSquareRegress (std::vector<std::pair <float, float> > * pvXY, std::pair <float, float> * pResultAB){

  float Xmoy = 0.0;
  float Ymoy = 0.0;
  for(unsigned int hitIndex=0;hitIndex < pvXY->size();hitIndex++){
    Xmoy+= pvXY->at(hitIndex).first;
    Ymoy+= pvXY->at(hitIndex).second;
  }

  Xmoy/= pvXY->size();
  Ymoy/= pvXY->size();
      
  float numeratorSum = 0.0;
  float denominatorSum = 0.0;

  for(unsigned int hitIndex=0;hitIndex < pvXY->size();hitIndex++){
    float X = pvXY->at(hitIndex).first;
    float Y = pvXY->at(hitIndex).second;
   
    numeratorSum+= (X-Xmoy)*(Y-Ymoy);
    denominatorSum+= (X-Xmoy)*(X-Xmoy);
  }

  pResultAB->first = numeratorSum/denominatorSum;
  pResultAB->second = Ymoy - (pResultAB->first * Xmoy);

}

void SeedClusteringFitter::fit(vector<Hit*> hits){

  //Barrel sector check
  for (unsigned int i = 0; i<hits.size();i++){
    if (hits[i]->getLayer() > 10 || hits[i]->getLayer() < 5){
      cout<<"Endcap and Hybrid sectors not implemented yet"<<endl;
      return;      
    }
  }

//cout<<"START OF THE SeedClusteringFit ("<<hits.size()<<" hits)"<<endl;

  float sector_phi_start_value = sec_phi - 0.4;

  std::vector< std::pair<float, float> > vRotatedXYCoordinates;

  //Rotation of all the XY coordinates in the range [0 ; PI/2]
  for (unsigned int i = 0; i<hits.size(); i++){
    vRotatedXYCoordinates.push_back(std::pair<float, float> (hits[i]->getX() * cos(-sector_phi_start_value) - hits[i]->getY() * sin(-sector_phi_start_value),
                                                             hits[i]->getX() * sin(-sector_phi_start_value) + hits[i]->getY() * cos(-sector_phi_start_value)));
  }

  //Seed generation threshold (min number of stubs in a row)
  unsigned int seeding_threshold;

  if(hits.size() > 250) {
    seeding_threshold = 6;
  }
  else if (hits.size() > 140){
    seeding_threshold = 5;
  }
  else if (hits.size() > 60){
    seeding_threshold = 4;
  }
  else {
    seeding_threshold = 3;
  }

  //Binary Mask Reset
  for (unsigned int i = 0; i < m_vIntercept.size(); i++){
    for (unsigned int j = 0; j < m_nLayer; j++){
      for (unsigned int k =0; k < 2; k++){
        m_pppBinMask[i][j][k] = false;
      }
    }
  }


  //Binary Mask Fill
  for (unsigned int hitIndex=0; hitIndex<hits.size(); hitIndex++){
      
    float X = vRotatedXYCoordinates[hitIndex].first;
    float Y = vRotatedXYCoordinates[hitIndex].second;

    float R = sqrt(X*X + Y*Y);
    float Phi = asin(Y/R);

    if (int(round(Phi * BINMASK_PHI_RES - 0.5)) < 0 || (unsigned int)(round(Phi * BINMASK_PHI_RES - 0.5))>=m_vIntercept.size())
      cout<<"MASK DOES'NT COVER ALL THE PHI RANGE OF THE SECTOR !!!"<<endl;

    unsigned int layerIndex = hits[hitIndex]->getLayer() - 5;
    unsigned int ordinateIndex = int(round(Phi * BINMASK_PHI_RES - 0.5));

    if (hits[hitIndex]->getLayer() <= 7 || abs(hits[hitIndex]->getBend()) <= BM_BEND_DC_THRESHOLD){
      //DC zone
      m_pppBinMask[ordinateIndex][layerIndex][0] = true;
      m_pppBinMask[ordinateIndex][layerIndex][1] = true;
    }
    else if (hits[hitIndex]->getBend() > BM_BEND_DC_THRESHOLD){
      //Really positive bend
      m_pppBinMask[ordinateIndex][layerIndex][0] = true;
    }
    else {
      //Really negative bend
      m_pppBinMask[ordinateIndex][layerIndex][1] = true;
    }
  }

//cout<<"MASK FILL OK"<<endl;

  std::vector <float> vCandidateSlope;
  std::vector <float> vCandidateIntercept;

  //Seeds generation (check rows on bin mask, if numStub in row > seeding_threshold, record slope and intercept of the row)
  for (unsigned int interceptIndex=0; interceptIndex<m_vIntercept.size(); interceptIndex++){
    for (unsigned int slopeIndex=0; slopeIndex<m_vSlope.size(); slopeIndex++){

      int min_max_ordinate_index_value = m_ppSlopeSeedMask[slopeIndex][m_nLayer-1] + interceptIndex;
      
      if (min_max_ordinate_index_value < int(m_vIntercept.size()) && min_max_ordinate_index_value >= 0){

        unsigned int bendIndex = 0;
        unsigned int cptStubsOnSeed = 0;

        if(m_vSlope[slopeIndex]>= 0){
          bendIndex = 1;
        }
        else{
          bendIndex = 0;
        }

        for (unsigned int layerIndex=0; layerIndex < m_nLayer; layerIndex++){

          int ordinate_index_value = m_ppSlopeSeedMask[slopeIndex][layerIndex] + interceptIndex;

	  if(m_pppBinMask[ordinate_index_value][layerIndex][bendIndex] == true){
            cptStubsOnSeed++;
          }
        }

        if (cptStubsOnSeed >= seeding_threshold) {
          vCandidateSlope.push_back(m_vSlope[slopeIndex]);
          vCandidateIntercept.push_back(m_vIntercept[interceptIndex]);
        }

      }
    }
  }

//cout<<"CANDIDATE SEEDS GENERATION OK, "<<vCandidateSlope.size()<<" Seeds generated"<<endl;

  //Creation of the Track candidates
  set<int> activated_layers;
  std::vector<Hit*> vCandidateHits;

  std::vector < std::pair<float,float> > vCandidatesRPHI;
  std::vector < std::pair<float,float> > vInternalLayersRZ;

  for (unsigned int candidateIndex=0; candidateIndex<vCandidateSlope.size(); candidateIndex++){
    vCandidateHits.clear();
    activated_layers.clear();

    for (unsigned int hitIndex=0; hitIndex<hits.size(); hitIndex++){

      float X = vRotatedXYCoordinates[hitIndex].first;
      float Y = vRotatedXYCoordinates[hitIndex].second;

      float R = sqrt(X*X + Y*Y);
      float Phi = asin(Y/R);

      float A = vCandidateSlope[candidateIndex];
      float B = vCandidateIntercept[candidateIndex];

      if (abs(A * R + B - Phi) <= ACCUMULATION_THRESHOLD){

        if ( (hits[hitIndex]->getLayer() <= 7) || 
           (abs(hits[hitIndex]->getBend()) <= TC_BEND_DC_THRESHOLD ) ||
           ( A <= 0 && hits[hitIndex]->getBend() > TC_BEND_DC_THRESHOLD) ||
           ( A >= 0 && hits[hitIndex]->getBend() < -TC_BEND_DC_THRESHOLD) ){
         
          activated_layers.insert(hits[hitIndex]->getLayer());
          vCandidateHits.push_back(hits[hitIndex]);
        }
      }
    }

    //If the Track candidate have hits on 5 or more different layers
    if (activated_layers.size() >= 5){
      vCandidatesRPHI.clear();
      vInternalLayersRZ.clear();

      float pt = 0.0;
      float phi = 0.0;
      float eta = 0.0;
      float z0 = 0.0;     

      for(unsigned int hitIndex=0;hitIndex < vCandidateHits.size();hitIndex++){

        float X = vCandidateHits[hitIndex]->getX() * cos(- sector_phi_start_value) - vCandidateHits[hitIndex]->getY() * sin(- sector_phi_start_value);
        float Y = vCandidateHits[hitIndex]->getX() * sin(- sector_phi_start_value) + vCandidateHits[hitIndex]->getY() * cos(- sector_phi_start_value);

        float R = sqrt(X*X + Y*Y);

        float Phi = atan(Y/X);

        vCandidatesRPHI.push_back( std::pair<float, float> ( R, Phi ) );
        
        //If the hit is on a Z precise Layer
        if (vCandidateHits[hitIndex]->getLayer() <= 7){

          float Z = vCandidateHits[hitIndex]->getZ();

          vInternalLayersRZ.push_back( std::pair<float, float> ( R, Z ) );
        }
      }
      
      std::pair<float, float> resultABwithRPHI(0.0,0.0);
      std::pair<float, float> resultABwithRZ(0.0,0.0);

      LinearLeastSquareRegress (&vInternalLayersRZ, &resultABwithRZ);
      LinearLeastSquareRegress (&vCandidatesRPHI, &resultABwithRPHI);

      //RPHI slope is used to define track_pt
      pt = abs( (3.8 * 0.3) / (2 * resultABwithRPHI.first * 100) );

      //Impossible for the system to precisely define higher pt
      if (pt > 300.0)
        pt = 300.0;
      
      //RPHI intercept is used to define track_phi
      phi = resultABwithRPHI.second + sector_phi_start_value;

      //phi is contained in the range ]-PI;PI]
      if (phi > M_PI){
        phi = phi - (2 * M_PI);
      }
      else if (phi < -M_PI){
        phi = phi + (2 * M_PI);
      }

      //RZ slope is used to define track_eta
      float theta = atan(resultABwithRZ.first) + (M_PI/2); //angle with the beam axis
      eta = log( tan( theta / 2 ) );


      //RZ intercept is used to define track_z0
      z0 =  resultABwithRZ.second;

      if (abs(z0)<= 30.0){
        Track* fit_track = new Track();
        fit_track->setCurve(pt);
        fit_track->setPhi0(phi);
        fit_track->setEta0(eta);
        fit_track->setZ0(z0);
      
        for(unsigned int hitIndex=0;hitIndex < vCandidateHits.size();hitIndex++){
          fit_track->addStubIndex(vCandidateHits[hitIndex]->getID());
        }
        tracks.push_back(fit_track);
      }
    }
  }

}

void SeedClusteringFitter::fit(){

  vector<Hit*> activatedHits;

  //////// Get the list of unique stubs from all the patterns ///////////
  set<int> ids;
  int total=0;
  for(unsigned int i=0;i<patterns.size();i++){
    vector<Hit*> allHits = patterns[i]->getHits();
    total+=allHits.size();
    for(unsigned int j=0;j<allHits.size();j++){
      pair<set<int>::iterator,bool> result = ids.insert(allHits[j]->getID());
      if(result.second==true)
	activatedHits.push_back(allHits[j]);
    }
  }
  fit(activatedHits);
}

TrackFitter* SeedClusteringFitter::clone(){
  SeedClusteringFitter* fit = new SeedClusteringFitter(nb_layers);
  fit->setPhiRotation(sec_phi);
  fit->setSectorID(sector_id);
  return fit;
}

