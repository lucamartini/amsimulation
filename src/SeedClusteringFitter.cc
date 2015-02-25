#include "SeedClusteringFitter.h"

SeedClusteringFitter::SeedClusteringFitter():TrackFitter(0){

}

SeedClusteringFitter::SeedClusteringFitter(int nb):TrackFitter(nb){

  m_nLayer = 6;

  m_pMeanLayerRadius = new float[m_nLayer];
 
  m_pMeanLayerRadius[0] = 22.5; 
  m_pMeanLayerRadius[1] = 35.0;
  m_pMeanLayerRadius[2] = 51.5;
  m_pMeanLayerRadius[3] = 68.5;
  m_pMeanLayerRadius[4] = 88.5;
  m_pMeanLayerRadius[5] = 107.5;  

  m_sector_phi_start_value = M_PI - 1.15;

  m_accumulation_threshold = 0.008;

  m_bin_mask_phi_res = 130.0; //in px/rad

  float bin_mask_phi_range = M_PI/2;
  

  float maxSlope = 0.003;
  float minSlope = -maxSlope;

  unsigned int slopes_per_intercept = 65; //an even number is better

  float step = (maxSlope - minSlope) / (slopes_per_intercept-1);
  float currentSlope = minSlope;
  while(currentSlope <= maxSlope){
    m_vSlope.push_back(currentSlope);
    currentSlope += step;
  }

  unsigned int bin_mask_height = int(ceil(bin_mask_phi_range * m_bin_mask_phi_res));

  step = (bin_mask_phi_range) / (bin_mask_height-1);
  float currentIntercept = 0;
  while(currentIntercept <= bin_mask_phi_range){
    m_vIntercept.push_back(currentIntercept);
    currentIntercept += step;
  }

  //Binary mask allocation
  m_ppBinMask = new bool * [m_vIntercept.size()];
  for (unsigned int i = 0; i < m_vIntercept.size(); i++)
    m_ppBinMask[i] = new bool [m_nLayer];

cout<<"BINMASK ALLOCATION OK, "<<"binMaskHeight = "<<m_vIntercept.size()<<endl;

  //Seed mask allocation
  m_ppSlopeSeedMask = new int * [m_vSlope.size()];
  for (unsigned int i = 0; i < m_vSlope.size(); i++)
    m_ppSlopeSeedMask[i] = new int [m_nLayer];


  for (unsigned int slopeIndex = 0; slopeIndex < m_vSlope.size(); slopeIndex++){

    float slope_value = m_vSlope[slopeIndex];
    
    //Determination of the corresponding pixel for each row and for each mask layer
    for (unsigned int layerIndex = 0; layerIndex < m_nLayer; layerIndex++){
      int maskIndex = int(trunc(slope_value * m_pMeanLayerRadius[layerIndex] * m_bin_mask_phi_res));
      m_ppSlopeSeedMask[slopeIndex][layerIndex] = maskIndex;
    }
  }

cout<<"SEEDS PRE-PROCESSING OK, "<<m_vSlope.size()<<" Seeds pre-processed"<<endl;

}


SeedClusteringFitter::~SeedClusteringFitter(){

  for (unsigned int i = 0; i < m_vIntercept.size(); i++)
    delete m_ppBinMask[i];

  delete m_ppBinMask;

  m_ppBinMask = NULL;

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
      
  float numeratorSum =0.0;
  float denominatorSum =0.0;

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
/*
  if(hits.size()>1024){
    cout<<"ERROR : too many stubs for fitting!"<<endl;
    return;
  }
*/

  //Seed generation threshold (min number of stubs in a row)
  unsigned int seeding_threshold;

  if(hits.size() > 300) {
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
  for (unsigned int i = 0; i < m_vIntercept.size(); i++)
    for (unsigned int j = 0; j < m_nLayer; j++)
      m_ppBinMask[i][j] = false;

cout<<"Analysing : "<<hits.size()<<" hits"<<endl;

  //Binary Mask Fill
  for (unsigned int i=0; i<hits.size(); i++){
    float X = hits[i]->getX();
    float Y = hits[i]->getY();

    float R = sqrt(X*X + Y*Y);
    float PHI = asin(abs(Y)/R);
    
    int LAYER = hits[i]->getLayer();

    if (X < 0 && Y > 0)
      PHI = M_PI - PHI;
    
    if (X > 0 && Y < 0)
      PHI = 2*M_PI - PHI;

    if (X < 0 && Y < 0)
      PHI = M_PI + PHI;

    PHI = PHI - m_sector_phi_start_value;

    m_ppBinMask[int(round(PHI * m_bin_mask_phi_res))][LAYER-5] = true;
  }

cout<<"MASK FILL OK"<<endl;

  std::vector <float> vCandidateSlope;
  std::vector <float> vCandidateIntercept;

  //Seeds generation (check rows on bin mask, if numStub in row > seeding_threshold, record slope and intercept of the row)

  for (unsigned int interceptIndex=0; interceptIndex<m_vIntercept.size(); interceptIndex++){
    for (unsigned int slopeIndex=0; slopeIndex<m_vSlope.size(); slopeIndex++){

      int min_max_ordinate_index_value = m_ppSlopeSeedMask[slopeIndex][m_nLayer-1] + interceptIndex;
      
      if (min_max_ordinate_index_value < int(m_vIntercept.size()) && min_max_ordinate_index_value >= 0){
        unsigned int cptStubsOnSeed = 0;

        for (unsigned int layerIndex=0; layerIndex < m_nLayer; layerIndex++){

          int ordinate_index_value = m_ppSlopeSeedMask[slopeIndex][layerIndex] + interceptIndex;

	  if(m_ppBinMask[ordinate_index_value][layerIndex] == true){
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

cout<<"CANDIDATE SEEDS GENERATION OK, "<<vCandidateSlope.size()<<" Seeds generated"<<endl;

  //Creation of the Track candidates
  set<int> activated_layers;
  std::vector<Hit*> vCandidateHits;

  std::vector < std::pair<float,float> > vCandidatesXY;
  std::vector < std::pair<float,float> > vInternalLayersRZ;

  for (unsigned int candidateIndex=0; candidateIndex<vCandidateSlope.size(); candidateIndex++){
    vCandidateHits.clear();

    for (unsigned int hitIndex=0; hitIndex<hits.size(); hitIndex++){
      float X = hits[hitIndex]->getX();
      float Y = hits[hitIndex]->getY();
      float R = sqrt(X*X + Y*Y);
      float PHI = asin(abs(Y)/R);

      if (X < 0 && Y > 0)
        PHI = M_PI - PHI;
    
      if (X > 0 && Y < 0)
        PHI = 2*M_PI - PHI;

      if (X < 0 && Y < 0)
        PHI = M_PI + PHI;

      PHI = PHI - m_sector_phi_start_value;

      float A = vCandidateSlope[candidateIndex];
      float B = vCandidateIntercept[candidateIndex];

      if (abs(A * R + B - PHI) <= m_accumulation_threshold){ 
        activated_layers.insert(hits[hitIndex]->getLayer());
        vCandidateHits.push_back(hits[hitIndex]);
      }
    }

    //If the Track candidate have hits on 5 or more different layers
    if (activated_layers.size() >= 5){
      vCandidatesXY.clear();
      vInternalLayersRZ.clear();

      float pt = 0.0;
      float phi = 0.0;
      float eta = 0.0;
      float z0 = 0.0;     

      for(unsigned int hitIndex=0;hitIndex < vCandidateHits.size();hitIndex++){

        float X = vCandidateHits[hitIndex]->getX();
        float Y = vCandidateHits[hitIndex]->getY();
        float Z = vCandidateHits[hitIndex]->getZ();

        vCandidatesXY.push_back( std::pair<float, float> ( X/(X*X + Y*Y), Y/(X*X + Y*Y) ) );
        
        //If the hit is on a Z precise Layer
        if (vCandidateHits[hitIndex]->getLayer() <= 7){      
          vInternalLayersRZ.push_back( std::pair<float, float> ( sqrt(X*X + Y*Y), Z ) );
        }
      }
      
      std::pair<float, float> resultABwithXY(0.0,0.0);
      std::pair<float, float> resultABwithRZ(0.0,0.0);

      LinearLeastSquareRegress (&vInternalLayersRZ, &resultABwithRZ);
      LinearLeastSquareRegress (&vCandidatesXY, &resultABwithXY);


      pt = 3.8 * 0.3 / 200 * sqrt(resultABwithXY.first*resultABwithXY.first + 1) / abs(resultABwithXY.second);
      phi = atan(resultABwithXY.first);

      eta = -log(tan(atan(resultABwithRZ.first)/2));
      z0 =  resultABwithRZ.second;
      
      if (pt > 200.0)
        pt = 200.0;

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

