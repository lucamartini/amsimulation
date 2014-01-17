#ifndef _GPUPOOLER_H_
#define _GPUPOOLER_H_


#include <sstream>
#include <iostream>
#include <fstream>
#include "FileEventProxy.h"
#include "gpu.h"
#include "PatternFinder.h"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


#define STUBLEN 9
#define PATTERN_STUBLEN 5
#define MAX_NB_STUBS 20000

class GPUPooler{
 private:
  deviceDetector d_detector;
  patternBank d_pb;
  deviceStubs d_stub1;
  deviceStubs d_stub2;
  deviceStubs* d_stubs[2];

  deviceParameters d_param;

  char* buf;
  int32_t* ibuf;
  float* vbuf;
  FileEventProxy* Raw_proxy;

  char** cuda_hits;
  bool** active_stubs;
  int* cuda_nb_hits;
  vector<Hit*> **hits;

  int eventID[2];
  int sectorID;

  FileEventProxy* pattern_proxy;
  char* pattern_buf;
  int32_t* pattern_ibuf;
  float* pattern_vbuf;
  
  SectorTree sectors;

  PatternFinder* pf;

  cudaStream_t** streams;

 public:
  GPUPooler(string sectorFilename, string inputDirectory, string outputDirectory, int patternThreshold);
  ~GPUPooler();
  void loadEvent(string fileName, int stream);
  void saveEvent(string fileName, int stream);
  void sendEventToDevice(int stream);
  void getEventFromDevice(int stream);
  void computeEvent(int stream);
  void loopForEvents(int waitingTime, int timeout);
};
#endif
