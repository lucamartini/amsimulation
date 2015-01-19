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
  float minEta;
  float maxEta;
  float minZ0;
  float maxZ0;
  float minPT;
  float maxPT;
  
 public:
  PatternInfo();
  PatternInfo(const PatternInfo& pi);
  void addTrack(float eta, float z0, float pt);
  void merge(const PatternInfo& pi);
  float getMinEta() const;
  float getMaxEta() const;
  float getMinZ0() const;
  float getMaxZ0() const;
  float getMinPT() const;
  float getMaxPT() const;
  int getNbTracks() const;

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const{
    ar << nb_tracks;
    ar << minEta;
    ar << maxEta;
    ar << minZ0;
    ar << maxZ0;
    ar << minPT;
    ar << maxPT;
  }
  
  template<class Archive> void load(Archive & ar, const unsigned int version){
    ar >> nb_tracks;
    ar >> minEta;
    ar >> maxEta;
    ar >> minZ0;
    ar >> maxZ0;
    ar >> minPT;
    ar >> maxPT;
  }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

};
#endif
