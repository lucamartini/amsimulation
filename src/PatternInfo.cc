#include "PatternInfo.h"

PatternInfo::PatternInfo(){
  nb_tracks = 0;
  minEta = 3;
  maxEta = -3;
  minZ0 = 20;
  maxZ0 = -20;
  minPT = 200;
  maxPT = -10;
}

PatternInfo::PatternInfo(const PatternInfo& pi){
  nb_tracks = pi.nb_tracks;
  minEta = pi.minEta;
  maxEta = pi.maxEta;
  minZ0 = pi.minZ0;
  maxZ0 = pi.maxZ0;
  minPT = pi.minPT;
  maxPT = pi.maxPT;
}

void PatternInfo::addTrack(float eta, float z0, float pt){
  if(eta>maxEta)
    maxEta=eta;
  if(eta<minEta)
    minEta=eta;

  if(z0>maxZ0)
    maxZ0=z0;
  if(z0<minZ0)
    minZ0=z0;

  if(pt>maxPT)
    maxPT=pt;
  if(pt<minPT)
    minPT=pt;

  nb_tracks++;
}

void PatternInfo::merge(const PatternInfo& pi){
  if(pi.minEta<minEta)
    minEta=pi.minEta;
  if(pi.maxEta>maxEta)
    maxEta=pi.maxEta;

  if(pi.minZ0<minZ0)
    minZ0=pi.minZ0;
  if(pi.maxZ0>maxZ0)
    maxZ0=pi.maxZ0;

  if(pi.minPT<minPT)
    minPT=pi.minPT;
  if(pi.maxPT>maxPT)
    maxPT=pi.maxPT;

  nb_tracks+=pi.getNbTracks();
  /*
  if(nb_tracks>10){
    cout<<"Nombre de traces : "<<nb_tracks<<endl;
    cout<<"sumPT : "<<sumPT<<" sumPTSquare : "<<sumPTSquare<<endl;
    cout<<"diviseur : "<<(nb*(nb-1))<<endl;
    cout<<"val : "<<(nb*sumPTSquare-(sumPT*sumPT))<<endl;
    cout<<"total : "<<(nb*sumPTSquare-(sumPT*sumPT))/(nb*(nb-1))<<endl;
    cout<<"resultat : "<<sqrt((nb*sumPTSquare-(sumPT*sumPT))/(nb*(nb-1)))<<endl;
    cout<<"Eta moyen : "<<averageEta<<endl;
    cout<<"Z0 moyen : "<<averageZ0<<endl;
    cout<<"Phi moyen : "<<averagePhi<<endl;
    cout<<"PT moyen : "<<averagePT<<endl;
    cout<<"Ecart Eta : "<<sEta<<endl;
    cout<<"Ecart Z0 : "<<sZ0<<endl;
    cout<<"Ecart PHI : "<<sPhi<<endl;
    cout<<"Ecart PT : "<<sPT<<endl; 
    cout<<endl;
  }
  */
}

float PatternInfo::getMinEta() const{
  return minEta;
}

float PatternInfo::getMaxEta() const{
  return maxEta;
}

float PatternInfo::getMinZ0() const{
  return minZ0;
}

float PatternInfo::getMaxZ0() const{
  return maxZ0;
}

float PatternInfo::getMinPT() const{
  return minPT;
}

float PatternInfo::getMaxPT() const{
  return maxPT;
}

int PatternInfo::getNbTracks() const{
  return nb_tracks;
}
