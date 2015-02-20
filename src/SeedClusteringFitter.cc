#include "SeedClusteringFitter.h"

SeedClusteringFitter::SeedClusteringFitter():TrackFitter(0){

}

SeedClusteringFitter::SeedClusteringFitter(int nb):TrackFitter(nb){

}

SeedClusteringFitter::~SeedClusteringFitter(){
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

void SeedClusteringFitter::fit(vector<Hit*> hits){
  if(hits.size()>1024){
    cout<<"ERROR : too many stubs for fitting!"<<endl;
    return;
  }

  cout<<"SeedClusteringFitter::fit()"<<endl;
  for(unsigned int i=0;i<hits.size();i++){
    cout<<*(hits[i])<<endl;
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
