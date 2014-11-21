#ifndef _PATTERNLAYER_H_
#define _PATTERNLAYER_H_

#include <iostream>
#include <bitset>
#include <map>
#include <cmath>
#include <cstring>

#include "Detector.h"

#ifdef USE_CUDA
#include "gpu.h"
#include "omp.h"
#endif

#include <boost/serialization/split_member.hpp>
#include <boost/serialization/map.hpp>

using namespace std;

/**
   \brief Layer part of a pattern
   This class contains the bits for one layer.
   This is an abstract class which should be implemented to be used (see class CMSPatternLayer for an example).
**/

class PatternLayer{

 public:
  /**
     \brief Number of bits per patternLayer
  **/
  static const int LAYER_BITS=15;
  /**
     \brief Number of maximum DC bits
  **/
  static const int DC_BITS=3;

  /**
     \brief Cache between DC bits values and positions
  **/
  static map<string, vector<short> > positions_cache;
  
  /**
     \brief Constructor
  **/
  PatternLayer();
  /**
     \brief Destructor
  **/
  virtual ~PatternLayer(){};
  /**
     \brief Returns a copy of the PatternLayer
     \return A pointer on the copy
  **/
  virtual PatternLayer* clone()=0;

  /**
     \brief Returns the ladder's position in phi
     \return The ladder's phi position
  **/
  virtual short getPhi()=0;

  /**
     \brief Set the values of the Don't Care bits
     \param index The index of the bit
     \param val The value of the bit (between 0 and 3)
  **/
  void setDC(int index, char val);
  /**
     \brief Get a DC bit
     \param index The index of the bit
     \return The value of the DC bit in position index
  **/
  char getDC(int index);
  /**
     \brief Returns the int value of the bitset as a string of 5 characters
     \return A string (ie : "28653" or "00142")
  **/
  string getCode();
  /**
     \brief Get the bitset value as an integer
     \return The value of the bitset (ie 28653 or 142)
  **/
  int getIntValue() const;
  /**
     \brief Change the value contained in the bitset
     \param v The new value as an integer
  **/
  void setIntValue(int v);
  /**
     \brief Retrieve the SuperStrip objects corresponding to the PatternLayer from the Detector structure
     \param l The layer of the PatternLayer (starting from 0)
     \param ladd The ladders of this layer for the current sector
     \param modules The modules in the current sector
     \param d The detector structure
     \return A list of SuperStrip*. If no DC bits are used we have only one value.
  **/
  virtual vector<SuperStrip*> getSuperStrip(int l, const vector<int>& ladd, const map<int, vector<int> >& modules, Detector& d)=0;
#ifdef USE_CUDA
  /**
     \brief Retrieve the SuperStrip objects corresponding to the PatternLayer from the Detector structure
     \param l The layer of the PatternLayer (starting from 0)
     \param ladd The ladders of this layer for the current sector
     \param modules The modules in the current sector
     \param layerID The ID of the current layer
     \param v : pointer on an allocated array of 8 integers used for the results (index of the superstrip in the detector).
  **/
  virtual void getSuperStripCuda(int l, const vector<int>& ladd, const map<int, vector<int> >& modules, int layerID, unsigned int* v)=0;
#endif
 /**
     \brief Allows to display a PatternLayer as a string
  **/
  virtual string toString()=0;

 /**
     \brief Returns the Super strip position
     \return The position of the super strip
  **/
  virtual short getStrip()=0;
  /**
     \brief Returns the Super strip encoded value (may be different from the position if using gray code)
     \return The encoded value of the super strip
  **/
  virtual short getStripCode()=0;
  /**
     \brief Get the list of positions from the DC bits
  **/
  vector<short> getPositionsFromDC();

  /**
     \brief Get the number of DC bits used
     \return The number of active DC bits for this PatternLayer
  **/
  int getDCBitsNumber();

  /**
     \brief Set the value of the stub PT for this layer
     \param pt The stub PT value
  **/
  void setPT(float pt);
  /**
     \brief Update the value of the stub PT for this layer according to the value of an other pattern
     \param pt The additional stub PT
  **/
  void updatePT(signed char pt);
  /**
     \brief Get the stub PT value for this layer
     \return -1 if uncertain, 0 if negative, 1 if positive
  **/
  signed char getPT();
  /**
     \brief Check if the PatternLayer is using the stub PT's information
     \return True if used, False otherwise
  **/
  bool getPTUsage();
  /**
     \brief Check if the PatternLayer is a fake one (used on layers not crossed by the track)
     \return True if the PatternLayer is a placeholder
  **/
  virtual bool isFake()=0;

 private:
  /**
     Get the list of positions from the DC bits
  **/
  void getPositionsFromDC(vector<char> dc, vector<short>& positions);

 protected:
  /**
     The value of the pattern layer (15 bits)
   **/
  bitset<LAYER_BITS> bits;
  /**  
      4 possible values for a DC bit:
        - 0 : 0
        - 1 : 1
	- 2 : X (Don't Care)
	- 3 : UNUSED
  **/
  char dc_bits[DC_BITS];
  
  /**
     Stub PT value :
       - 0 if all tracks leading to this pattern have a clearly negative PT
       - 1 if all tracks leading to this pattern have a clearly positive PT
       - -1 if all tracks leading to this pattern do not have the same PT sign or uncertain values
  **/
  signed char ptValue;

  /**
     Flag on the usage of the stub PT information : if true, the last DC bits is taken by this information
  **/
  bool useStubPT;

  friend class boost::serialization::access;
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const//const boost::serialization::version_type& version) const 
    {
      unsigned long i(this->getIntValue());
      ar << i; 
      for(int j=0;j<DC_BITS;j++){
	ar << dc_bits[j];
      }
      ar << useStubPT;
    }
  
  template<class Archive> void load(Archive & ar, const unsigned int version)
    {
      int i;
      ar >> i;
      setIntValue(i); 
      for(int j=0;j<DC_BITS;j++){
	ar >> dc_bits[j];
      }
      if(version>0){
	ar >> useStubPT;
	if(useStubPT){
	  char val = dc_bits[getDCBitsNumber()];
	  if(val!=2)
	    ptValue=val;
	}
      }
    }
  
  BOOST_SERIALIZATION_SPLIT_MEMBER()
  
};
BOOST_CLASS_VERSION(PatternLayer, 1)
#endif
