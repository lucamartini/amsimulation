#ifndef _GPU_H_
#define _GPU_H_

#include "gpu_struct.h"
#include <cuda_runtime.h>




extern "C"{
//  detector* d_detector;

  void allocateDetector(deviceDetector* d);
  void resetDetector(deviceDetector* d, cudaStream_t* stream=NULL);
  void allocateBank(patternBank* p, int nbPatterns);
  void freeDetector(deviceDetector* d);
  void allocateStubs(deviceStubs* s);
  void freeStubs(deviceStubs* s);
  void freeBank(patternBank* p);
  void allocateParameters(deviceParameters* dp);
  void freeParameters(deviceParameters* dp);
  void cudaSetLink(patternBank* p, int index, unsigned int* vals);
  void cudaSetNbPatterns(patternBank* p, int nb);
  void cudaCopyStubs(char* stubs, deviceStubs* d_stubs, int nb, cudaStream_t* stream=NULL);
  void cudaGetActiveStubs(bool* active_stubs, deviceStubs* d_stubs, int* nb, cudaStream_t* stream=NULL);
  void getHitsArray(deviceDetector* det, int* list, int nb);
  void cudaShowBank(patternBank*p);
  void cudaShowStubs(deviceStubs* s, int nb);
  void printCard();
  void resetCard();
  void initialiseTimer();
  void startTimer();
  float stopTimer();
  void deleteTimer();
  void cudaSetHitsWrapper(deviceStubs* d_stubs, int nbStubs, deviceDetector* d_det, cudaStream_t* stream=NULL);
  int cudaGetActivePatternsWrapper(deviceDetector* detector, patternBank* patterns, deviceStubs* stubs, deviceParameters* params, int threshold, int nbThreads, int nbBlocks, cudaStream_t* stream=NULL);
}
#endif
