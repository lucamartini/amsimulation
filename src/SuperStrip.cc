#include "SuperStrip.h"

SuperStrip::SuperStrip(int s){
  hit=false;
  hit_by_negative_bend=false;
  hit_by_positive_bend=false;
  size=s;
}

SuperStrip::~SuperStrip(){
  clear();
}

short SuperStrip::getSize(){
  return size;
}

bool SuperStrip::isHit(){
  return hit;
}

bool SuperStrip::isHit(bool bend){
  if(bend)
    return hit_by_positive_bend;
  else
    return hit_by_negative_bend;
}

vector<Hit*>& SuperStrip::getHits(){
  return hits;
}

void SuperStrip::clear(){
  hit=false;
  hit_by_negative_bend=false;
  hit_by_positive_bend=false;
  for(unsigned int i=0;i<hits.size();i++){
    delete(hits[i]);
  }
  hits.clear();
}

void SuperStrip::touch(const Hit* h){
  hit=true;
  if(h->getBend())
    hit_by_positive_bend=true;
  else
    hit_by_negative_bend=true;

  Hit* copy = new Hit(*h);
  hits.push_back(copy);
}
