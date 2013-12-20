#ifndef _GPU_H_
#define _GPU_H_

#include "gpu_struct.h"

const int NB_LAYER = 16;
const int NB_LADDER = 76;
const int NB_MODULE = 80;
const int NB_SEGMENT = 2;
const int SIZE_SEGMENT = 64;
const int SIZE_MODULE = NB_SEGMENT*SIZE_SEGMENT;
const int SIZE_LADDER = NB_MODULE*SIZE_MODULE;
const int SIZE_LAYER = NB_LADDER*SIZE_LADDER;

const int MAX_NB_STUBS_PER_SSTRIPS = 8;

const int PATTERN_LAYERS = 8;
const int PATTERN_SSTRIPS = 8;
const int PATTERN_SIZE = PATTERN_LAYERS*PATTERN_SSTRIPS;
const unsigned int PATTERN_UNUSED = 4294967295;

const int CUDA_STUB_SIZE = 5;
const int CUDA_MAX_NB_STUBS = 5000;

const int CUDA_NB_CORES = 2304;

const int cuda_layer_index[23]={-1,-1,-1,-1,-1,0,1,2,3,4,5,6,7,8,9,10,-1,-1,11,12,13,14,15};


extern "C"{
//  detector* d_detector;

  void allocateDetector(deviceDetector* d);
  void resetDetector(deviceDetector* d);
  void allocateBank(patternBank* p, int nbPatterns);
  void freeDetector(deviceDetector* d);
  void allocateStubs(deviceStubs* s);
  void freeStubs(deviceStubs* s);
  void freeBank(patternBank* p);
  void cudaSetLink(patternBank* p, int index, unsigned int* vals);
  void cudaSetNbPatterns(patternBank* p, int nb);
  void cudaCopyStubs(char* stubs, deviceStubs* d_stubs, int nb);
  void cudaGetActiveStubs(bool* active_stubs, deviceStubs* d_stubs, int nb);
  void getHitsArray(deviceDetector* det, int* list, int nb);
  void cudaShowBank(patternBank*p);
  void cudaShowStubs(deviceStubs* s, int nb);
  void printCard();
  void resetCard();
  void initialiseTimer();
  void startTimer();
  float stopTimer();
  void deleteTimer();
  void cudaSetHitsWrapper(deviceStubs* d_stubs, int nbStubs, deviceDetector* d_det);
  int cudaGetActivePatternsWrapper(deviceDetector* detector, patternBank* patterns, deviceStubs* stubs, int threshold, int nbThreads, int nbBlocks);
}
#endif
