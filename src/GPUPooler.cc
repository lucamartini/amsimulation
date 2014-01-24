#include "GPUPooler.h"

GPUPooler::GPUPooler(string sectorFilename, string inputDirectory, string outputDirectory, int patternThreshold){
  resetCard(); 
  cout<<"Loading file "<<sectorFilename<<"..."<<endl;
  std::ifstream ifs(sectorFilename.c_str());
  boost::archive::text_iarchive ia(ifs);
  ia >> sectors;

  Raw_proxy = new FileEventProxy(inputDirectory);
  pattern_proxy = new FileEventProxy(outputDirectory);

  allocateDetector(&d_detector);
  allocateBank(&d_pb,sectors.getAllSectors()[0]->getLDPatternNumber());
  allocateStubs(&d_stub1);
  allocateStubs(&d_stub2);
  allocateParameters(&d_param);
  d_stubs[0]=&d_stub1;
  d_stubs[1]=&d_stub2;
  sectors.getAllSectors()[0]->linkCuda(&d_pb,&d_detector);
  buf = (char*)malloc(MAX_NB_STUBS*STUBLEN*sizeof(int32_t));
  ibuf=( int32_t*) buf;
  vbuf=( float*) buf;
  pattern_buf = (char*)malloc(MAX_NB_STUBS*PATTERN_STUBLEN*sizeof(int32_t));
  pattern_ibuf=( int32_t*) pattern_buf;
  pattern_vbuf=( float*) pattern_buf;

  pf = new PatternFinder(sectors.getSuperStripSize(), patternThreshold, &sectors, "","", &d_pb, &d_detector,&d_param);

  cudaMallocHost((void**)&cuda_nb_hits,2*sizeof(int));

  //cuda_hits = new char*[2];
  cudaMallocHost((void**)&cuda_hits,2*sizeof(char*));
  cudaMallocHost((void**)&cuda_hits[0],5000*sizeof(char));
  cudaMallocHost((void**)&cuda_hits[1],5000*sizeof(char));

  cudaMallocHost((void**)&active_stubs,2*sizeof(bool*));
  cudaMallocHost((void**)&active_stubs[0],5000*sizeof(bool));
  cudaMallocHost((void**)&active_stubs[1],5000*sizeof(bool));

  hits = new vector<Hit*>*[2];
  hits[0] = new vector<Hit*>();
  hits[1] = new vector<Hit*>();

  streams = new cudaStream_t*[2];
  streams[0] = new cudaStream_t();
  streams[1] = new cudaStream_t();
  cudaStreamCreate(streams[0]);
  cudaStreamCreate(streams[1]);

  sectorID = sectors.getAllSectors()[0]->getOfficialID();
} 

GPUPooler::~GPUPooler(){
  delete Raw_proxy;
  delete pattern_proxy;
  delete pf;
  
  cudaFree(cuda_nb_hits);

  cudaFree(cuda_hits[0]);
  cudaFree(cuda_hits[1]);
  cudaFree(cuda_hits);
  //  delete [] cuda_hits;

  cudaFree(active_stubs[0]);
  cudaFree(active_stubs[1]);
  cudaFree(active_stubs);

  cudaStreamDestroy(*streams[0]);
  cudaStreamDestroy(*streams[1]);

  delete streams[0];
  delete streams[1];
  delete [] streams;

  free(pattern_buf);
  free(buf); 
  freeParameters(&d_param);
  freeDetector(&d_detector);
  freeBank(&d_pb);
  freeStubs(&d_stub1);
  freeStubs(&d_stub2);
  resetCard();
}

void GPUPooler::loadEvent(string fileName, int stream){
  //cout<<"loadEvent "<<stream<<endl;
  uint32_t size_buf;
  Raw_proxy->Read(fileName,buf,size_buf);
	  
  unsigned pos = fileName.find("_");
  string evt_str = fileName.substr (pos+1);
  std::istringstream ss( evt_str );
  ss >> eventID[stream];

  int nbStubs = size_buf/STUBLEN/sizeof(int32_t);
  
  cuda_nb_hits[stream] = 0;
  hits[stream]->clear();
  
  for(int i=0;i<nbStubs;i++){
    int index = i*STUBLEN;
    
    int32_t tp = ibuf[index];
    int32_t layer = ibuf[index+1];
    int32_t module = CMSPatternLayer::getModuleCode(layer, ibuf[index+2]);
    int32_t ladder = CMSPatternLayer::getLadderCode(layer, ibuf[index+3]);
    int32_t seg = CMSPatternLayer::getSegmentCode(layer, ladder, ibuf[index+4]);
    int32_t strip = ibuf[index+5];
    float x = vbuf[index+6];
    float y = vbuf[index+7];
    float z = vbuf[index+8];
    
    Hit* h = new Hit(layer,ladder, module, seg, strip, i, tp, 0, 0, 0, 0, x, y, z, 0, 0, 0);
    if(sectors.getSector(*h)!=NULL){
      hits[stream]->push_back(h);
      
      int cuda_idx = cuda_nb_hits[stream]*CUDA_STUB_SIZE;
      cuda_hits[stream][cuda_idx]=cuda_layer_index[layer];
      cuda_hits[stream][cuda_idx+1]=ladder;
      cuda_hits[stream][cuda_idx+2]=module;
      cuda_hits[stream][cuda_idx+3]=seg;
      cuda_hits[stream][cuda_idx+4]=(char)(strip/sectors.getSuperStripSize());
      cuda_nb_hits[stream]++;
      
    }
    else
      delete(h);
  }
}
      
void GPUPooler::saveEvent(string fileName, int stream){
  //cout<<"saveEvent "<<stream<<endl;
  int stubIndex = 0;
  unsigned int ns = 0;
  uint32_t pattern_size_buf=0;
  vector<Hit*> hit_list = *hits[stream];
  for(unsigned int i=0;i<hit_list.size();i++){
    if(active_stubs[stream][i]){
      pattern_ibuf[ns*PATTERN_STUBLEN+0]=hit_list[i]->getLayer();
      pattern_ibuf[ns*PATTERN_STUBLEN+1]=hit_list[i]->getParticuleID();
      pattern_vbuf[ns*PATTERN_STUBLEN+2]=hit_list[i]->getX();
      pattern_vbuf[ns*PATTERN_STUBLEN+3]=hit_list[i]->getY();
      pattern_vbuf[ns*PATTERN_STUBLEN+4]=hit_list[i]->getZ();
      pattern_size_buf+=PATTERN_STUBLEN*sizeof(uint32_t);
      ns++;
      //cout<<*hits[i]<<endl;
      stubIndex++;
    }
  }
  
  //cout<<"event "<<eventID[stream]<<" : nombre de stubs distincts apres patterns : "<<stubIndex<<endl;
  
  std::stringstream s;
  s<<"Event_"<<eventID[stream]<<"_"<<sectorID;
  pattern_proxy->Write(s.str(),pattern_buf,pattern_size_buf);
  Raw_proxy->Erase(fileName);
  
  for(unsigned int i=0;i<hit_list.size();i++){
    delete hit_list[i];
  }	
  hit_list.clear();
}

void GPUPooler::sendEventToDevice(int stream){
  //cout<<"sendEvent "<<stream<<endl;
  cudaCopyStubs(cuda_hits[stream],d_stubs[stream],cuda_nb_hits[stream],streams[stream]); 
}

void GPUPooler::getEventFromDevice(int stream){
  //cout<<"getEvent "<<stream<<endl;
  cudaGetActiveStubs(active_stubs[stream],d_stubs[stream],&cuda_nb_hits[stream],streams[stream]); 
}

void GPUPooler::computeEvent(int stream){
  //cout<<"computeEvent sur "<<streams[stream]<<endl;
  pf->findCuda(cuda_nb_hits[stream],d_stubs[stream],streams[stream]);
}

void GPUPooler::loopForEvents(int waitingTime, int timeout){
  stringstream spat;
  spat<<"Event_*";

  int max_loop = timeout/waitingTime;
  int nb_loop = 0;

  cout<<"Looking for events..."<<endl;

  while(nb_loop<max_loop){
    nb_loop++;
    std::vector<std::string> files;
    Raw_proxy->List(files, spat.str());
    
    if(files.size()==0){
      usleep(waitingTime*1000);
      continue;
    }
    
    nb_loop = 0;
    initialiseTimer();		     
    startTimer();

    std::vector<std::string>::iterator it=files.begin();
    int streamIndex = 0;

    loadEvent(*it,streamIndex);
    sendEventToDevice(streamIndex);
    computeEvent(streamIndex);

    while(it+1!=files.end()){
      loadEvent(*(it+1),!streamIndex);
      sendEventToDevice(!streamIndex);
      cudaStreamSynchronize(*streams[streamIndex]);
      computeEvent(!streamIndex);
      getEventFromDevice(streamIndex);
      cudaStreamSynchronize(*streams[streamIndex]);
      saveEvent(*it,streamIndex);

      it++;
      streamIndex=1-streamIndex;
    }
    cudaDeviceSynchronize();
    getEventFromDevice(streamIndex);
    cudaStreamSynchronize(*streams[streamIndex]);
    saveEvent(*it,streamIndex);

    float time_laps = stopTimer();
    cout<<files.size()<<" events computed in "<<time_laps<<" ms"<<endl;
    cout<<endl;
  }
  cout<<"Timeout reached."<<endl;
}
     
