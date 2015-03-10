#ifndef _SEEDCLUSTERINGFITTER_H_
#define _SEEDCLUSTERINGFITTER_H_

#include "TrackFitter.h"
#include <math.h>

#include <iomanip>
#include <set>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

#define N_SLOPES 81			//Number of tested slopes for each intercept of the mask
#define MAX_SLOPE 0.0019		//Max slope (rad/cm) (corresponding to the minimum wanted particle pt) = (3.8 * 0.3)/(2 * minPt * 100)
#define MIN_SLOPE -0.0019		//Min slope (rad/cm) (corresponding to the minimum wanted particle pt) = -(3.8 * 0.3)/(2 * minPt * 100)
#define BINMASK_PHI_RES 128.0   	//BinaryMask ordinate resolution in px/rad
#define BINMASK_PHI_MAX M_PI/2		//BinaryMask maximum phi value (at least equal to the phi range of a sector)
#define ACCUMULATION_THRESHOLD 0.004	//Maximum tolerated error between the seed and the hit coordinates for the acumulation	
#define BM_BEND_DC_THRESHOLD 1.5	//Beyond this limit, we consider the bend of the external layers hits gives the correct sign for the slope (used to reduce the number of seeds generated)
#define TC_BEND_DC_THRESHOLD 1.5	//Beyond this limit, we consider the bend of the external layers hits gives the correct sign for the slope (used to reduce the number stubs in a track candidate)

/**
   \brief Implementation of a Seed Clustering  fitter
**/
class SeedClusteringFitter:public TrackFitter{

 private:

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const
    {
      ar << boost::serialization::base_object<TrackFitter>(*this);
      ar << sec_phi;
    }
  
  template<class Archive> void load(Archive & ar, const unsigned int version)
    {
      ar >> boost::serialization::base_object<TrackFitter>(*this);
      ar >> sec_phi;
    }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

  std::vector <float> m_vSlope;
  std::vector <float> m_vIntercept;

  bool *** m_pppBinMask;
  int ** m_ppSlopeSeedMask;
  unsigned int m_nLayer;
  float * m_pMeanLayerRadius;

  void LinearLeastSquareRegress (std::vector<std::pair <float, float> > *, std::pair <float, float> *);

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
  void fit();
  void fit(vector<Hit*> hits);
  TrackFitter* clone();
};
#endif
