#ifndef _GPUSTRUCT_H
#define _GPUSTRUCT_H

const int NB_LAYER = 16;
const int NB_LADDER = 76;
const int NB_MODULE = 80;
const int NB_SEGMENT = 2;
const int SIZE_SEGMENT = 64;
const int SIZE_MODULE = NB_SEGMENT*SIZE_SEGMENT;
const int SIZE_LADDER = NB_MODULE*SIZE_MODULE;
const int SIZE_LAYER = NB_LADDER*SIZE_LADDER;

const int MAX_NB_STUBS_PER_SSTRIPS = 8;

const int PATTERN_LAYERS = 9;
const int PATTERN_SSTRIPS = 8;
const int PATTERN_SIZE = PATTERN_LAYERS*PATTERN_SSTRIPS;
const unsigned int PATTERN_UNUSED = 4294967295;

const int CUDA_STUB_SIZE = 5;
const int CUDA_MAX_NB_STUBS = 5000;

const int CUDA_NB_CORES = 2304;

const int cuda_layer_index[23]={-1,-1,-1,-1,-1,0,1,2,3,4,5,6,7,8,9,10,-1,-1,11,12,13,14,15};

typedef struct {
  bool* sstrips;
  int* stubs;
} deviceDetector;

typedef struct {
  unsigned int* banks;
  int* nb_patterns;
} patternBank;

typedef struct {
  char* stubs;
  bool* active_stubs;
  int* nb_stubs;
} deviceStubs;

typedef struct {
  int* result;
  int* threshold;
  int* iter;
  int* nbPatterns;
} deviceParameters;

#endif
