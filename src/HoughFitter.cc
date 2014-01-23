#include "HoughFitter.h"

HoughFitter::HoughFitter():HoughFitter(0){
}

HoughFitter::HoughFitter(int nb):TrackFitter(nb){
  cuts = new HoughCut();
  ch = new ComputerHough(cuts);
  h_x= new float[1024];
  h_y= new float[1024];
  h_z= new float[1024];
  h_layer= new unsigned int[1024];
  ch->DefaultCuts();
}

HoughFitter::~HoughFitter(){
  delete ch;
  delete cuts;
  delete [] h_x;
  delete [] h_y;
  delete [] h_z;
  delete [] h_layer;
}

void HoughFitter::initialize(){

}

void HoughFitter::mergePatterns(){
  //cout<<"Merging of patterns not implemented"<<endl;
}

void HoughFitter::mergeTracks(){
  //cout<<"Merging of Tracks not implemented"<<endl;
}

void HoughFitter::fit(vector<Hit*> hits){
  if(hits.size()>1024){
    cout<<"ERROR : too many stubs for fitting!"<<endl;
    return;
  }

  memset(h_x,0,1024*sizeof(float));
  memset(h_y,0,1024*sizeof(float));
  memset(h_z,0,1024*sizeof(float));
  memset(h_layer,0,1024*sizeof(unsigned int));
 
  //Collect the stubs
  for(unsigned int j=0;j<hits.size();j++){
    h_x[j]=hits[j]->getX();
    h_y[j]=hits[j]->getY();
    h_z[j]=hits[j]->getZ();
    h_layer[j]=hits[j]->getLayer();
  }

  //Actual computing
  ch->ComputeOneShot(sector_id,hits.size(),h_x,h_y,h_z,h_layer);
  std::vector<mctrack_t> &v=ch->getCandidates();

  //format the results and feed the tracks vector
  for (std::vector<mctrack_t>::iterator it=v.begin();it!=v.end();it++){
    Track* fit_track = new Track((*it).pt, 0, (*it).phi, (*it).eta, (*it).z0);
    tracks.push_back(fit_track);
  }
  
}

void HoughFitter::fit(){

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

TrackFitter* HoughFitter::clone(){
  HoughFitter* fit = new HoughFitter(nb_layers);
  fit->setPhiRotation(sec_phi);
  fit->setSectorID(sector_id);
  return fit;
}
