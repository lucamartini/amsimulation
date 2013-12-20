#ifndef _GPUSTRUCT_H
#define _GPUSTRUCT_H

/*
typedef struct {
//  int nb_stubs;
//  int stubs_index[10];
  bool hit;
} superstrip;

typedef struct {
  superstrip** sstrips;
} segment;

typedef struct {
  segment* seg0;
  segment* seg1;
} module;

typedef struct {
  int nb_modules;
  module** modules;
} ladder;

typedef struct {
  int nb_ladders;
  ladder** ladders;
} layer;
*/

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

#endif
