#ifndef _CMSPATTERNLAYER_H_
#define _CMSPATTERNLAYER_H_

#include <iostream>
#include <sstream>
#include <bitset>
#include "PatternLayer.h"

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

using namespace std;


/**
   \brief First version of a CMS pattern structure
**/

class CMSPatternLayer : public PatternLayer{
 private:
  static const short MOD_START_BIT = 11;
  static const short PHI_START_BIT = 7;
  static const short STRIP_START_BIT = 1;
  static const short SEG_START_BIT = 0;

  static const short MOD_MASK = 0x1F;
  static const short PHI_MASK = 0xF;
  static const short STRIP_MASK = 0x3F;
  static const short SEG_MASK = 0x1;

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const//const boost::serialization::version_type& version) const 
    {
      ar << boost::serialization::base_object<PatternLayer>(*this);
    }
  
  template<class Archive> void load(Archive & ar, const unsigned int version)
    {
      ar >> boost::serialization::base_object<PatternLayer>(*this);
    }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()

 public:
  CMSPatternLayer();
  CMSPatternLayer* clone();
  vector<SuperStrip*> getSuperStrip(int l, const vector<int>& ladd, const map<int, vector<int> >& modules, Detector& d);
  
  /**
     \brief Set the values in the patternLayer
     \param m The module Z position (0 to 13 for modules 23 to 47)
     \param phi The phi position of the module in the sector (0 to 7)
     \param strip The super strip number
     \param seg The segment in the module (0 or 1)
  **/
  void setValues(short m, short phi, short strip, short seg);
  /**
     \brief Returns a string representation of the PatternLayer
     \return A string describing the PatternLayer
  **/
  string toString();
  /**
     \brief Returns the module's Z position
     \return The module's Z position
  **/
  short getModule();
  /**
     \brief Returns the ladder phi position
     \return The ladder's phi position (0 or 1)
  **/
  short getPhi();
  /**
     \brief Returns the Super strip position
     \return The position of the super strip in the segment
  **/
  short getStrip();
  /**
     \brief Returns the position of the segment in the module
     \return The segment's position in the module (0 or 1)
  **/
  short getSegment();

  /**
     \brief Get the ID of the ladder from the ladder ID in the muon file. Used to change the IDs between the root file and the simulation program (if needed).
     \brief Ladders numbering must start at 0
     \param layerID The layer ID (TIB : 5,6,7 - TOB : 8,9,10 - TEC : 11,12,13,14,15,16,17 and 18,19,20,21,22,23,24)
     \param ladderID The ladder ID in the muon file
     \return The ID to use in the program
  **/
  static int getLadderCode(int layerID, int ladderID);

 /**
     \brief Get the ID of the module from the module ID in the muon file. Used to change the IDs between the root file and the simulation program (if needed).
     \brief This method allows to change the module resolution (you can divide the ID by 2) or the numbering.
     \brief Modules numbering must start at 0
     \param layerID The layer ID (TIB : 5,6,7 - TOB : 8,9,10 - TEC : 11,12,13,14,15,16,17 and 18,19,20,21,22,23,24)
     \param moduleID The module ID in the muon file
     \return The ID to use in the program
  **/
  static int getModuleCode(int layerID, int moduleID);

  /**
     \brief Get the code of the segment in the patternLayer from the segment ID in the muon file
     \brief Segment ID must be 0 or 1
     \param layerID The layer ID (TIB : 5,6,7 - TOB : 8,9,10 - TEC : 11,12,13,14,15,16,17 and 18,19,20,21,22,23,24)
     \param ladderID The ladder ID in the muon file
     \param segmentID The segment ID in the muon file
     \return The number to use in the simulation program
  **/
  static int getSegmentCode(int layerID, int ladderID, int segmentID);

  static int getNbLadders(int layerID);

  static int getNbModules(int layerID, int ladderID);  

};

#endif
