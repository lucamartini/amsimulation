#include "UnitTest.h"

#define BOOST_TEST_MODULE ALL_TESTS tests
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/program_options.hpp>
#include <boost/progress.hpp>

#include <iostream>

#include "PatternTree.h"
#include "PatternGenerator.h"
#include "PatternFinder.h"
#include "SectorTree.h"
#include "Detector.h"
#include "PrincipalTrackFitter.h"
#include "PrincipalFitGenerator.h"

#ifndef __APPLE__
BOOST_CLASS_EXPORT_IMPLEMENT(CMSPatternLayer) 
BOOST_CLASS_EXPORT_IMPLEMENT(PrincipalTrackFitter) 
#endif

using namespace std;

BOOST_AUTO_TEST_CASE( CMSPatternLayer_constructor_test )
{
  const int NB_LAYERS = 6;
  const int SEGMENT = 1;
  const int MODULE = 2;
  const int LADDER = 3;
  const int SSTRIP = 4;

  Pattern p(NB_LAYERS);

  for(int i=0;i<NB_LAYERS;i++){
    CMSPatternLayer test;
    test.setValues(MODULE,LADDER,SSTRIP,SEGMENT);
    p.setLayerStrip(i, &test);
  }

  BOOST_CHECK_EQUAL( p.getNbLayers() , NB_LAYERS );

  for(int i=0;i<NB_LAYERS;i++){
    CMSPatternLayer* pl = (CMSPatternLayer*)p.getLayerStrip(i);
    BOOST_CHECK_EQUAL(pl->getModule(),MODULE);
    BOOST_CHECK_EQUAL(pl->getPhi(),LADDER);
    BOOST_CHECK_EQUAL(pl->getStrip(),SSTRIP);
    BOOST_CHECK_EQUAL(pl->getSegment(),SEGMENT);
  }
}


BOOST_AUTO_TEST_CASE( bank_compatibility_test ){

  const int BANK_SIZE = 139922;
  int* pattern100[6];
  int patternLayer0[4] = {0,0,0,0};
  int patternLayer1[4] = {0,0,1,6};
  int patternLayer2[4] = {1,0,2,6};
  int patternLayer3[4] = {2,1,1,1};
  int patternLayer4[4] = {2,0,3,4};
  int patternLayer5[4] = {3,0,3,0};
  
  pattern100[0] = patternLayer0;
  pattern100[1] = patternLayer1;
  pattern100[2] = patternLayer2;
  pattern100[3] = patternLayer3;
  pattern100[4] = patternLayer4;
  pattern100[5] = patternLayer5;

  SectorTree st;
  {
    std::ifstream ifs("./test_data/test_bank.pbk");
    boost::iostreams::filtering_stream<boost::iostreams::input> f;
    f.push(boost::iostreams::gzip_decompressor());
    //we try to read a compressed file
    try { 
      f.push(ifs);
      boost::archive::text_iarchive ia(f);
      ia >> st;
    }
    catch (boost::iostreams::gzip_error& e) {
      if(e.error()==4){//file is not compressed->read it without decompression
	std::ifstream new_ifs("./test_data/test_bank.pbk");
	boost::archive::text_iarchive ia(new_ifs);
	ia >> st;
      }
    }
  }
  
  vector<Sector*> sectors = st.getAllSectors();

  BOOST_CHECK_EQUAL((int)sectors.size(),1);

  for(unsigned int i=0;i<sectors.size();i++){
    Sector* mySector = sectors[i];
    vector<GradedPattern*> patterns = mySector->getPatternTree()->getLDPatterns();

    BOOST_CHECK_EQUAL((int)patterns.size(),BANK_SIZE);

    Pattern* p = patterns[100];
    for(int k=0;k<p->getNbLayers();k++){
      CMSPatternLayer* mp = (CMSPatternLayer*)p->getLayerStrip(k);
      BOOST_CHECK_EQUAL(mp->getModule(),pattern100[k][0]);
      BOOST_CHECK_EQUAL(mp->getSegment(),pattern100[k][1]);
      BOOST_CHECK_EQUAL(mp->getPhi(),pattern100[k][2]);
      BOOST_CHECK_EQUAL(mp->getStrip(),pattern100[k][3]);
    }
  }
  
}

BOOST_AUTO_TEST_CASE( pattern_finding_test ){
  SectorTree st;
  {
    std::ifstream ifs("./test_data/test_bank.pbk");
    boost::iostreams::filtering_stream<boost::iostreams::input> f;
    f.push(boost::iostreams::gzip_decompressor());
    //we try to read a compressed file
    try { 
      f.push(ifs);
      boost::archive::text_iarchive ia(f);
      ia >> st;
    }
    catch (boost::iostreams::gzip_error& e) {
      if(e.error()==4){//file is not compressed->read it without decompression
	std::ifstream new_ifs("./test_data/test_bank.pbk");
	boost::archive::text_iarchive ia(new_ifs);
	ia >> st;
      }
    }
  }

  ///////////////////////////////////////////////////////////////
  // If we don't have a fitter -> create a Hough default one
  vector<Sector*> sectors = st.getAllSectors();
  for(unsigned int i=0;i<sectors.size();i++){
    if(sectors[i]->getFitter()==NULL){
      //TrackFitter* fitter = new KarimakiTrackFitter(sectors[i]->getNbLayers());
      TrackFitter* fitter = new HoughFitter(sectors[i]->getNbLayers());
      sectors[i]->setFitter(fitter);
      sectors[i]->updateFitterPhiRotation();
    }
  }
  ///////////////////////////////////////////////////////////////

  PatternFinder pf(5, &st,  "test_data/PU4T_620SLHC7_light.root",  "test_data/output.root");
  {
    pf.useMissingHitThreshold(1);
    int stop = 10;
    pf.find(0, stop);
  }
}
