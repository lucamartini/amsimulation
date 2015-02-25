#ifndef _SEEDCLUSTERINGFITTER_H_
#define _SEEDCLUSTERINGFITTER_H_

#include "TrackFitter.h"
#include <math.h>

#include <iomanip>
#include <set>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

/**
   \brief Implementation of a Seed Clustering  fitter
**/
class SeedClusteringFitter:public TrackFitter{

 private:

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const
    {
      cout<<"sauvegarde du SeedClusteringFitter"<<endl;
      ar << boost::serialization::base_object<TrackFitter>(*this);
      ar << sec_phi;
    }
  
  template<class Archive> void load(Archive & ar, const unsigned int version)
    {
      cout<<"chargement de SeedClusteringFitter"<<endl;
      ar >> boost::serialization::base_object<TrackFitter>(*this);
      ar >> sec_phi;
    }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

  std::vector <float> m_vSlope;
  std::vector <float> m_vIntercept;

  bool ** m_ppBinMask;
  int ** m_ppSlopeSeedMask;

  unsigned int m_nLayer;

  float * m_pMeanLayerRadius;

  float m_sector_phi_start_value;
  float m_accumulation_threshold;
  float m_bin_mask_phi_res;

 public:

  /**
     \brief Default constructor   
  **/
  SeedClusteringFitter();
  /**
     \brief Constructor   
     \param nb Layers number
  **/  
  SeedClusteringFitter(int nb);
  ~SeedClusteringFitter();

  void initialize();
  void mergePatterns();
  void mergeTracks();
  void LinearLeastSquareRegress (std::vector<std::pair <float, float> > *, std::pair <float, float> *);
  void fit();
  void fit(vector<Hit*> hits);
  TrackFitter* clone();
};
#endif
