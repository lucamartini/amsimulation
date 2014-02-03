#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

// includes CUDA
#include <cuda_runtime.h>

// includes, project
#include <helper_cuda.h>
#include <helper_functions.h> // helper functions for SDK examples
#include <device_functions.h> // helper functions for SDK examples

#include "gpu.h"

using namespace std;



void printCard(){
    const int kb = 1024;
    const int mb = kb * kb;
    cout << "NBody.GPU" << endl << "=========" << endl << endl;

    cout << "CUDA version:   v" << CUDART_VERSION << endl;    

    int devCount;
    cudaGetDeviceCount(&devCount);
    wcout << "CUDA Devices: " << endl << endl;

    for(int i = 0; i < devCount; ++i)
    {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);
        cout << i << ": " << props.name << ": " << props.major << "." << props.minor << endl;
        cout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << endl;
        cout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << endl;
        cout << "  Constant memory: " << props.totalConstMem / kb << "kb" << endl;
        cout << "  Block registers: " << props.regsPerBlock << endl << endl;

        cout << "  Warp size:         " << props.warpSize << endl;
        cout << "  Threads per block: " << props.maxThreadsPerBlock << endl;
        cout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1]  << ", " << props.maxThreadsDim[2] << " ]" << endl;
        cout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1]  << ", " << props.maxGridSize[2] << " ]" << endl;
        cout << endl;
    }
}

void allocateDetector(deviceDetector* d){
     if(cudaSuccess != cudaMalloc((void**)&d->sstrips, sizeof(bool)*SIZE_LAYER*NB_LAYER))
        cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
 
     // max of MAX_NB_STUBS_PER_SSTRIPS stubs per superstrip
     if(cudaSuccess != cudaMalloc((void**)&d->stubs, sizeof(int)*SIZE_LAYER*NB_LAYER*MAX_NB_STUBS_PER_SSTRIPS))
        cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;

     resetDetector(d);
}

void resetDetector(deviceDetector* d, cudaStream_t* stream){
  if(stream==NULL){
     if(cudaSuccess != cudaMemset(d->sstrips,0,sizeof(bool)*SIZE_LAYER*NB_LAYER))
        cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;	
     if(cudaSuccess != cudaMemset(d->stubs, -1, sizeof(int)*SIZE_LAYER*NB_LAYER*MAX_NB_STUBS_PER_SSTRIPS))
        cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;	  
  }
  else{
     if(cudaSuccess != cudaMemsetAsync(d->sstrips,0,sizeof(bool)*SIZE_LAYER*NB_LAYER, *stream))
        cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;	
     if(cudaSuccess != cudaMemsetAsync(d->stubs, -1, sizeof(int)*SIZE_LAYER*NB_LAYER*MAX_NB_STUBS_PER_SSTRIPS, *stream))
        cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;	  
  }
}

void freeDetector(deviceDetector* d){
     deleteTimer();
     cudaFree(d->sstrips);
     cudaFree(d->stubs);
}

void allocateBank(patternBank* p, int nbPatterns){
     if(cudaSuccess != cudaMalloc((void**)&p->banks, sizeof(unsigned int)*PATTERN_SIZE*nbPatterns))
       cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
     if(cudaSuccess != cudaMalloc((void**)&p->nb_patterns, sizeof(int)))
       cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
     if(cudaSuccess !=cudaMemset(p->banks,255, sizeof(unsigned int)*PATTERN_SIZE*nbPatterns))
       cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
     if(cudaSuccess !=cudaMemset(p->nb_patterns,nbPatterns, sizeof(int)))
       cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
}

void freeBank(patternBank* p){
     cudaFree(p->banks);
     cudaFree(p->nb_patterns);
}

void allocateStubs(deviceStubs* s){
 if(cudaSuccess != cudaMalloc((void**)&s->stubs, sizeof(char)*CUDA_MAX_NB_STUBS*CUDA_STUB_SIZE))
       cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
 if(cudaSuccess != cudaMalloc((void**)&s->active_stubs, sizeof(bool)*CUDA_MAX_NB_STUBS))
       cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
 if(cudaSuccess != cudaMalloc((void**)&s->nb_stubs, sizeof(int)))
       cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
}

void freeStubs(deviceStubs* s){
     cudaFree(s->stubs);
     cudaFree(s->active_stubs);
     cudaFree(s->nb_stubs);
}

void allocateParameters(deviceParameters* dp){
  if(cudaSuccess != cudaMalloc((void**)&dp->result, sizeof(int)))
    cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;

  if(cudaSuccess != cudaMalloc((void**)&dp->threshold, sizeof(int)))
    cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
 
  if(cudaSuccess != cudaMalloc((void**)&dp->iter, sizeof(int)))
    cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
 
  if(cudaSuccess != cudaMalloc((void**)&dp->nbPatterns, sizeof(int)))
    cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
}

void freeParameters(deviceParameters* dp){
  cudaFree(dp->result);
  cudaFree(dp->threshold);
  cudaFree(dp->iter);
  cudaFree(dp->nbPatterns);
}

void cudaSetLink(patternBank* p, int index, unsigned int* vals){
  if(cudaSuccess != cudaMemcpy(p->banks+index,vals,PATTERN_LAYERS*PATTERN_SSTRIPS*sizeof(unsigned int), cudaMemcpyHostToDevice))
      cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
}

void cudaSetNbPatterns(patternBank* p, int nb){
  if(cudaSuccess != cudaMemcpy(p->nb_patterns,&nb,sizeof(int), cudaMemcpyHostToDevice))
      cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl; 
}

void cudaCopyStubs(char* stubs, deviceStubs* d_stubs, int nb, cudaStream_t* stream){
 if(stream==NULL){
     if(cudaSuccess != cudaMemcpy(d_stubs->stubs,stubs,nb*CUDA_STUB_SIZE*sizeof(char), cudaMemcpyHostToDevice))
          cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
     if(cudaSuccess !=cudaMemset(d_stubs->active_stubs,false, nb*sizeof(bool)))
          cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
 }
 else{
     if(cudaSuccess != cudaMemcpyAsync(d_stubs->stubs,stubs,nb*CUDA_STUB_SIZE*sizeof(char), cudaMemcpyHostToDevice, *stream))
          cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
     if(cudaSuccess !=cudaMemsetAsync(d_stubs->active_stubs,false, nb*sizeof(bool), *stream))
          cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
 }
}

void cudaGetActiveStubs(bool* active_stubs, deviceStubs* d_stubs, int* nb, cudaStream_t* stream){
  if(stream==NULL){
    if(cudaSuccess != cudaMemcpy(active_stubs,d_stubs->active_stubs,(*nb)*sizeof(bool), cudaMemcpyDeviceToHost))
        cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
  }
  else{
    if(cudaSuccess != cudaMemcpyAsync(active_stubs,d_stubs->active_stubs,(*nb)*sizeof(bool), cudaMemcpyDeviceToHost, *stream))
        cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
  }
}

void cudaShowBank(patternBank*p){
     unsigned int val;
     for(int i=0;i<10;i++){
       for(int j=0;j<8;j++){
         for(int k=0;k<8;k++){
           if(cudaSuccess != cudaMemcpy(&val,p->banks+i*PATTERN_SIZE+j*PATTERN_LAYERS+k,sizeof(unsigned int), cudaMemcpyDeviceToHost))
             cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
           cout<<val<<" ";
	 }
	 cout<<endl;
       }
       cout<<endl;
     }
     cout<<endl;
     cout<<endl;
}

void cudaShowStubs(deviceStubs* s, int nb){
     char val;
     for(int i=0;i<nb;i++){
       if(cudaSuccess != cudaMemcpy(&val,s->stubs+i*CUDA_STUB_SIZE+0,sizeof(char), cudaMemcpyDeviceToHost))
         cout<<"cudaShowStubs : error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
       cout<<"layer "<<(int)val<<" ";
       if(cudaSuccess != cudaMemcpy(&val,s->stubs+i*CUDA_STUB_SIZE+1,sizeof(char), cudaMemcpyDeviceToHost))
         cout<<"cudaShowStubs : error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
       cout<<"ladder "<<(int)val<<" ";
       if(cudaSuccess != cudaMemcpy(&val,s->stubs+i*CUDA_STUB_SIZE+2,sizeof(char), cudaMemcpyDeviceToHost))
         cout<<"cudaShowStubs : error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
       cout<<"module "<<(int)val<<" ";
       if(cudaSuccess != cudaMemcpy(&val,s->stubs+i*CUDA_STUB_SIZE+3,sizeof(char), cudaMemcpyDeviceToHost))
         cout<<"cudaShowStubs : error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
       cout<<"segment "<<(int)val<<" ";
       if(cudaSuccess != cudaMemcpy(&val,s->stubs+i*CUDA_STUB_SIZE+4,sizeof(char), cudaMemcpyDeviceToHost))
         cout<<"cudaShowStubs : error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
       cout<<"superstrip "<<(int)val<<" ";

       cout<<endl;
     }
}

void resetCard(){
    if(cudaSuccess != cudaDeviceReset())
      cout<<"Error : "<<cudaGetErrorString(cudaGetLastError())<<endl;
}

StopWatchInterface *theTimer=0;

void initialiseTimer()
{
  theTimer=0;
  sdkCreateTimer(&theTimer);
}

void startTimer()
{
   sdkResetTimer(&theTimer);
   sdkStartTimer(&theTimer);
}

float stopTimer()
{
  sdkStopTimer(&theTimer);
  float t=sdkGetTimerValue(&theTimer);
  printf("Processing time: %f (ms)\n",t );
  return t;
}

void deleteTimer()
{
   sdkDeleteTimer(&theTimer);
}

__global__ void cudaSetHits(char* d_stubs, int nbStubs, bool* d_det, int* d_det_stubs){
  int index = threadIdx.x * CUDA_STUB_SIZE + blockIdx.x * blockDim.x * CUDA_STUB_SIZE;

  if(index/CUDA_STUB_SIZE<nbStubs){
    int detIndex = d_stubs[index]*SIZE_LAYER+d_stubs[index+1]*SIZE_LADDER+d_stubs[index+2]*SIZE_MODULE+d_stubs[index+3]*SIZE_SEGMENT+d_stubs[index+4];
    d_det[detIndex]=1;
    int det_stubs_index=detIndex*MAX_NB_STUBS_PER_SSTRIPS;

    for(int count=0;count<MAX_NB_STUBS_PER_SSTRIPS;count++){
      int old = atomicCAS(&d_det_stubs[det_stubs_index+count],-1,index);
      if(old==-1)
        break;
    }
  }
}

__global__ void cudaGetActivePatterns(bool* detector, int* detector_stubs, bool* active_stubs, unsigned int* patterns, int* threshold, int* nbIter, int* nbMaxPatterns, int* nbActivePatterns){
   int index = blockIdx.x * blockDim.x * PATTERN_SIZE * (*nbIter) + threadIdx.x * PATTERN_SIZE * (*nbIter);
   char score = 0;

   for(int l=0;l<(*nbIter);l++){
     if(index<PATTERN_SIZE*(*nbMaxPatterns)){
       score = 0;
       for(int i=0;i<PATTERN_LAYERS;i++){
         for(int j=0;j<PATTERN_SSTRIPS;j++){
           unsigned int ref = patterns[index+i*PATTERN_SSTRIPS+j];
           if(ref==PATTERN_UNUSED)
             break;
           if(detector[ref]){
             score++;
             break;
           }
         }
       }
       if(score>=(*threshold)){
         atomicAdd(nbActivePatterns,1);
         for(int i=0;i<PATTERN_LAYERS;i++){
           for(int j=0;j<PATTERN_SSTRIPS;j++){
             unsigned int ref = patterns[index+i*PATTERN_SSTRIPS+j];
             if(ref==PATTERN_UNUSED)
               break;
             int stub_index = ref*MAX_NB_STUBS_PER_SSTRIPS; 
	     for(int k=0;k<MAX_NB_STUBS_PER_SSTRIPS;k++){
	       if(detector_stubs[stub_index+k]!=-1){
		 active_stubs[detector_stubs[stub_index+k]/CUDA_STUB_SIZE]=1;
               }
 	     }
           }
         }
       }
     }
     index += PATTERN_SIZE;
   }
   // __syncthreads();
}

void cudaSetHitsWrapper(deviceStubs* d_stubs, int nbStubs, deviceDetector* d_det, cudaStream_t* stream){
  if(nbStubs>0){
    cudaSetHits<<<nbStubs,1,0,((stream==NULL)?0:*stream)>>>(d_stubs->stubs,nbStubs,d_det->sstrips,d_det->stubs);
  }
}

void getHitsArray(deviceDetector* det, int* list, int nb){
   if(cudaSuccess !=cudaMemcpy(list,det->stubs,  sizeof(int)*SIZE_LAYER*NB_LAYER*MAX_NB_STUBS_PER_SSTRIPS,cudaMemcpyDeviceToHost))
       cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
}

int cudaGetActivePatternsWrapper(deviceDetector* detector, patternBank* patterns, deviceStubs* stubs, deviceParameters* params, int threshold, int nbThreads, int nbBlocks, cudaStream_t* stream){

  if(stream==NULL){
    if(cudaSuccess !=cudaMemset(params->result,0, sizeof(int)))
         cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
  }
  else{
    if(cudaSuccess !=cudaMemsetAsync(params->result,0, sizeof(int),*stream))
         cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
  }

    if(stream==NULL)
        cudaGetActivePatterns<<<nbBlocks,nbThreads>>>(detector->sstrips, detector->stubs, stubs->active_stubs, patterns->banks, params->threshold, params->iter, params->nbPatterns, params->result);
    else
        cudaGetActivePatterns<<<nbBlocks,nbThreads,0,*stream>>>(detector->sstrips, detector->stubs, stubs->active_stubs, patterns->banks, params->threshold, params->iter, params->nbPatterns, params->result);
  int res=0;

  if(stream==NULL){
    if(cudaSuccess != cudaMemcpy(&res,params->result,sizeof(int), cudaMemcpyDeviceToHost))
      cout<<"error! "<<cudaGetErrorString(cudaGetLastError())<<endl;
  }

  return res;
}
