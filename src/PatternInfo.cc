#include "PatternInfo.h"

PatternInfo::PatternInfo(){
  nb_tracks = 0;
  averageEta = 0;
  averageZ0 = 0;
  averagePhi = 0;
  averagePT = 0;
  sEta = 0;
  sumEta = 0;
  sumEtaSquare = 0;
  sZ0 = 0;
  sumZ0 = 0;
  sumZ0Square = 0;
  sPhi = 0;
  sumPhi = 0;
  sumPhiSquare = 0;
  sPT = 0;
  sumPT = 0;
  sumPTSquare = 0;
}

PatternInfo::PatternInfo(const PatternInfo& pi){
  nb_tracks = pi.nb_tracks;
  averageEta = pi.averageEta;
  averageZ0 = pi.averageZ0;
  averagePhi = pi.averagePhi;
  averagePT = pi.averagePT;
  sEta = pi.sEta;
  sumEta = pi.sumEta;
  sumEtaSquare = pi.sumEtaSquare;
  sZ0 = pi.sZ0;
  sumZ0 = pi.sumZ0;
  sumZ0Square = pi.sumZ0Square;
  sPhi = pi.sPhi;
  sumPhi = pi.sumPhi;
  sumPhiSquare = pi.sumPhiSquare;
  sPT = pi.sPT;
  sumPT = pi.sumPT;
  sumPTSquare = pi.sumPTSquare;
}

void PatternInfo::addTrack(float eta, float z0, float phi, float pt){
  averageEta = (averageEta*nb_tracks+eta)/(nb_tracks+1);
  averageZ0 = (averageZ0*nb_tracks+z0)/(nb_tracks+1);
  averagePhi = (averagePhi*nb_tracks+phi)/(nb_tracks+1);
  averagePT = (averagePT*nb_tracks+pt)/(nb_tracks+1);
  sumEta+=eta;
  sumEtaSquare+=(eta*eta);
  sEta = sqrt(((nb_tracks+1)*sumEtaSquare-(sumEta*sumEta))/(nb_tracks*(nb_tracks+1)));
  sumZ0+=z0;
  sumZ0Square+=(z0*z0);
  sZ0 = sqrt(((nb_tracks+1)*sumZ0Square-(sumZ0*sumZ0))/(nb_tracks*(nb_tracks+1)));
  sumPhi+=phi;
  sumPhiSquare+=(phi*phi);
  sPhi = sqrt(((nb_tracks+1)*sumPhiSquare-(sumPhi*sumPhi))/(nb_tracks*(nb_tracks+1)));
  sumPT+=pt;
  sumPTSquare+=(pt*pt);
  sPT = sqrt(((nb_tracks+1)*sumPTSquare-(sumPT*sumPT))/(nb_tracks*(nb_tracks+1)));
  nb_tracks++;
}

void PatternInfo::merge(const PatternInfo& pi){
  averageEta = (averageEta*nb_tracks+pi.getAverageEta()*pi.getNbTracks())/(nb_tracks+pi.getNbTracks());
  averageZ0 = (averageZ0*nb_tracks+pi.getAverageZ0()*pi.getNbTracks())/(nb_tracks+pi.getNbTracks());
  averagePhi = (averagePhi*nb_tracks+pi.getAveragePhi()*pi.getNbTracks())/(nb_tracks+pi.getNbTracks());
  averagePT = (averagePT*nb_tracks+pi.getAveragePT()*pi.getNbTracks())/(nb_tracks+pi.getNbTracks());

  int nb = nb_tracks+pi.getNbTracks();

  sumEta+=pi.sumEta;
  sumEtaSquare+=pi.sumEtaSquare;
  sEta = sqrt((nb*sumEtaSquare-(sumEta*sumEta))/(nb*(nb-1)));
  if(sEta!=sEta)
    sEta=0;

  sumZ0+=pi.sumZ0;
  sumZ0Square+=pi.sumZ0Square;
  sZ0 = sqrt((nb*sumZ0Square-(sumZ0*sumZ0))/(nb*(nb-1)));
  if(sZ0!=sZ0)
    sZ0=0;

  sumPhi+=pi.sumPhi;
  sumPhiSquare+=pi.sumPhiSquare;
  sPhi = sqrt((nb*sumPhiSquare-(sumPhi*sumPhi))/(nb*(nb-1)));
  if(sPhi!=sPhi)
    sPhi=0;

  sumPT+=pi.sumPT;
  sumPTSquare+=pi.sumPTSquare;
  sPT = sqrt((nb*sumPTSquare-(sumPT*sumPT))/(nb*(nb-1)));
  if(sPT!=sPT)
    sPT=0;

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

float PatternInfo::getAverageEta() const{
  return averageEta;
}

float PatternInfo::getAverageZ0() const{
  return averageZ0;
}

float PatternInfo::getAveragePT() const{
  return averagePT;
}

float PatternInfo::getAveragePhi() const{
  return averagePhi;
}

float PatternInfo::getSDEta() const{
  return sEta;
}

float PatternInfo::getSDZ0() const{
  return sZ0;
}

float PatternInfo::getSDPT() const{
  return sPT;
}

float PatternInfo::getSDPhi() const{
  return sPhi;
}

int PatternInfo::getNbTracks() const{
  return nb_tracks;
}
