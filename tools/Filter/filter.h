#ifndef FILTER_H
#define FILTER_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>
#include <stdio.h>  
#include <stdlib.h> 

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <fstream>
#include <string>
#include <sstream> 

///////////////////////////////////
//
//
// Base class for bank root files filtering 
//
// The role of this class is to select events belonging to a given 
// sector in a muons bank file. The result is a skimmed root file with exactly
// the same format, but containing tracks belonging to the sectors 
// (with more than 4 good stubs in the sector, to be precise)
//
// filename    : the name and directory of the input ROOT file to filter
// secfilename : the name and directory of the csv file containing the sectors definition
// outfile     : the name of the output ROOT file containing the filtered events
//           
// secid       : the sector number
// hit_lim     : the minimum number of layer/disks with stubs in the sector
//
// Info about the sector definition:
//
//  http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.HLLHCData
//
///////////////////////////////////

using namespace std;

class filter
{
 public:

  filter(std::string filename, std::string secfilename, 
	 std::string outfile, int secid, int hit_lim);

  void   do_filter(int secid, int hit_lim);    
  void   initTuple(std::string test,std::string out);
  bool   convert(std::string sectorfilename); 
  void   reset();
    
 private:

  TFile  *m_outfile;
  TChain *m_L1TT;

  TTree  *m_efftree;

  std::vector< std::vector<int> >   m_modules; 

  int m_sec_mult;

  int m_stub;
  int mf_stub;

  std::vector<int>   m_stub_modid;
  std::vector<float> m_stub_strip;
  std::vector<float> m_stub_ptGEN;  
  std::vector<float> m_stub_etaGEN;

  std::vector<int>   *pm_stub_modid;
  std::vector<float> *pm_stub_strip;
  std::vector<float> *pm_stub_ptGEN;  
  std::vector<float> *pm_stub_etaGEN;

  std::vector<int>   *mf_stub_modid;
  std::vector<float> *mf_stub_strip;
  std::vector<float> *mf_stub_ptGEN;  
  std::vector<float> *mf_stub_etaGEN;
};

#endif

