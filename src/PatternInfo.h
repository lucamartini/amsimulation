#ifndef _PATTERNINFO_H_
#define _PATTERNINFO_H_

#include <iostream>
#include <sstream>

#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>

using namespace std;

/**
   \brief Statistical informations about the pattern (number of generating tracks, aveage eta, PT, ...)
**/

class PatternInfo{
 private:
  int nb_tracks;
  float averageEta;
  float averageZ0;
  float averagePhi;
  float averagePT;
  float sEta;
  float sumEta;
  float sumEtaSquare;
  float sZ0;
  float sumZ0;
  float sumZ0Square;
  float sPhi;
  float sumPhi;
  float sumPhiSquare;
  float sPT;
  float sumPT;
  float sumPTSquare;
  
 public:
  PatternInfo();
  PatternInfo(const PatternInfo& pi);
  void addTrack(float eta, float z0, float phi, float pt);
  void merge(const PatternInfo& pi);
  float getAverageEta() const;
  float getAverageZ0() const;
  float getAveragePT() const;
  float getAveragePhi() const;
  float getSDEta() const;
  float getSDZ0() const;
  float getSDPT() const;
  float getSDPhi() const;
  int getNbTracks() const;

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const{
    ar << nb_tracks;
    ar << averageEta;
    ar << averageZ0;
    ar << averagePhi;
    ar << averagePT;

    if(sEta!=sEta ||
       sZ0!=sZ0 ||
       sPT!=sPT
       ){
      cout<<"probleme de NaN!"<<endl;
      float val = 0;
      ar << val;
      ar << val;
      ar << val;
      ar << val;
    }
    else{ 
      ar << sEta;
      ar << sZ0;
      ar << sPhi;
      ar << sPT;
    }

    ar << sumEta;
    ar << sumEtaSquare;
    ar << sumZ0;
    ar << sumZ0Square;
    ar << sumPT;
    ar << sumPTSquare;
  }
  
  template<class Archive> void load(Archive & ar, const unsigned int version){
    ar >> nb_tracks; 
    ar >> averageEta;
    ar >> averageZ0;
    ar >> averagePhi;
    ar >> averagePT;
    ar >> sEta;
    ar >> sZ0;
    ar >> sPhi;
    ar >> sPT;
    ar >> sumEta;
    ar >> sumEtaSquare;
    ar >> sumZ0;
    ar >> sumZ0Square;
    ar >> sumPT;
    ar >> sumPTSquare;
  }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

};
#endif
