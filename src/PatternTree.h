#ifndef _PATTERNTREE_H_
#define _PATTERNTREE_H_

#include <map>
#include <vector>
#include "PatternTrunk.h"
#include "CMSPatternLayer.h"

using namespace std;
/**
   \brief This class is used to store all the patterns (both low and high resolutions). To quickly find a pattern, PatternTrunk objects are stored in a map using the low resolution pattern value as the map key (string made of the decimal representation of the patternLayers values).
**/
class PatternTree{
 public:
  /**
     \brief Constructor
  **/
  PatternTree();
  ~PatternTree();
 /**
     \brief Add a copy of the pattern to the structure or increment the grade of the already stored pattern
     \param ldp The low definition pattern (will be copied)
     \param fdp The corresponding full definition pattern (will be copied, can be NULL if only one definition level is used)
  **/
  void addPattern(Pattern* ldp, Pattern* fdp);
 /**
     \brief Add a pattern to the structure or increment the grade of the already stored pattern. Update the average Pt of the pattern.
     \param ldp The low definition pattern
     \param fdp The corresponding full definition pattern (can be NULL if only one definition level is used)
     \param new_pt The Pt of the track generating the pattern
  **/
  void addPattern(Pattern* ldp, Pattern* fdp, float new_pt);

 /**
     \brief Get a copy of the full definition patterns
     \return A vector of copies
  **/
  vector<GradedPattern*> getFDPatterns();
 /**
     \brief Get a copy of the low definition patterns
     \return A vector of copies
  **/
  vector<GradedPattern*> getLDPatterns();
  /**
     \brief Get the distribution of PT among the LDPatterns
     \return A vector containing the occurences of each PT (bin is 1).
  **/
  vector<int> getPTHisto();
  /**
     \brief Get the number of LD patterns
     \return The number of LD patterns in the structure
  **/
  int getLDPatternNumber();
  /**
     \brief Get the number of FD patterns
     \return The number of FD patterns in the structure
  **/
  int getFDPatternNumber();
  /**
     \brief Link all patterns to the detector structure
     \param d The detector
     \param sec The ladders in the sector (one vector per layer)
     \param modules The modules in the sector (one vector per ladder)
  **/
  void link(Detector& d, const vector< vector<int> >& sec, const vector<map<int, vector<int> > >& modules);
#ifdef USE_CUDA
  /**
     \brief Link all patterns to the detector structure
     \param p the pattern bank structure on the device
     \param d The detector structure on the device
     \param sec The ladders in the sector (one vector per layer)
     \param modules The modules in the sector (one vector per ladder)
     \param layers The layers IDs
  **/
  void linkCuda(patternBank* p, deviceDetector* d, const vector< vector<int> >& sec, const vector<map<int, vector<int> > >& modules, vector<int> layers);
#endif
  /**
     \brief Returns a vector of copies of the active patterns
     \param active_threshold The minimum number of hit super strips to activate the pattern
     \return A vector containing copies of active patterns
  **/
  void getActivePatterns(int active_threshold, vector<GradedPattern*>& active_patterns);
  /**
     \brief Returns a vector of copies of the active patterns
     \param max_nb_missing_hit The maximum number of non active layers to activate the pattern
     \param active_threshold The minimum number of hit super strips to activate the pattern
     \return A vector containing copies of active patterns
  **/
  void getActivePatternsUsingMissingHit(int max_nb_missing_hit, int active_threshold, vector<GradedPattern*>& active_patterns);
  /**
     \brief Replace all LD patterns with adapatative patterns. All FD patterns are removed.
     \param r The number of DC bits used between FD and LD
  **/
  void computeAdaptativePatterns(short r);
  /**
     \brief Add all LD patterns coming from an other PatternTree
     \param p The PatternTree containing the patterns to add
  **/
  void addPatternsFromTree(PatternTree* p);
  
  /**
     \brief Check if the given pattern is contained in the bank (using DC bits)
     \brief Should only be used with a DC bit activated bank
     \param lp The low definition pattern
     \param hp The high definition version of the pattern
     \return True if the pattern is already in the bank, false otherwise
   **/
  bool checkPattern(Pattern* lp, Pattern* hp);

  /**
     \brief Replace the internal map of patterns with a vector of patterns (reduces memory consumption).
     \brief If a method searching for a pattern is called, we will automatically switch back to a map.
   **/
  void switchToVector();

  /**
     \brief Delete the least used patterns to match the given pattern number
     \param nbPatterns The number of patterns to keep
  **/
  void truncate(int nbPatterns);

 private:
  map<string, PatternTrunk*> patterns;
  vector<PatternTrunk*> v_patterns;

  /**
     \brief Add a pattern and update de DC bits if necessary
     \param ldp The pattern to add
  **/
  void addPatternForMerging(GradedPattern* ldp);
  /**
     \brief Update the internal map and clear the internal vector
  **/
  void switchToMap();

  friend class boost::serialization::access;
 
  template<class Archive> void save(Archive & ar, const unsigned int version) const{
    ar << CMSPatternLayer::MOD_START_BIT;
    ar << CMSPatternLayer::PHI_START_BIT;
    ar << CMSPatternLayer::STRIP_START_BIT;
    ar << CMSPatternLayer::SEG_START_BIT;
    ar << CMSPatternLayer::MOD_MASK;
    ar << CMSPatternLayer::PHI_MASK;
    ar << CMSPatternLayer::STRIP_MASK;
    ar << CMSPatternLayer::SEG_MASK;

    ar << patterns;
  }
  
  template<class Archive> void load(Archive & ar, const unsigned int version){
    if(version>0){//The format of the pattern is contained in the file
      ar >> CMSPatternLayer::MOD_START_BIT;
      ar >> CMSPatternLayer::PHI_START_BIT;
      ar >> CMSPatternLayer::STRIP_START_BIT;
      ar >> CMSPatternLayer::SEG_START_BIT;
      
      ar >> CMSPatternLayer::MOD_MASK;
      ar >> CMSPatternLayer::PHI_MASK;
      ar >> CMSPatternLayer::STRIP_MASK;
      ar >> CMSPatternLayer::SEG_MASK;
    }
    else{//we use the old values for retro compatibility
      CMSPatternLayer::MOD_START_BIT = 11;
      CMSPatternLayer::PHI_START_BIT = 7;
      CMSPatternLayer::STRIP_START_BIT = 1;
      CMSPatternLayer::SEG_START_BIT = 0;
      CMSPatternLayer::MOD_MASK = 0x1F;
      CMSPatternLayer::PHI_MASK = 0xF;
      CMSPatternLayer::STRIP_MASK = 0x3F;
      CMSPatternLayer::SEG_MASK = 0x1;
    }
    ar >> patterns;
  }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()
};
bool comparePatterns(PatternTrunk* p1, PatternTrunk* p2);
BOOST_CLASS_VERSION(PatternTree, 1)
#endif
