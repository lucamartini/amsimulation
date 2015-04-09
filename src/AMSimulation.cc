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

#ifdef USE_CUDA
#include "GPUPooler.h"
#include "cuda_profiler_api.h"
#endif

#include <TH1I.h>
#include <TFile.h>

#ifndef __APPLE__
BOOST_CLASS_EXPORT_IMPLEMENT(CMSPatternLayer) 
BOOST_CLASS_EXPORT_IMPLEMENT(PrincipalTrackFitter) 
#endif

using namespace std;

/**
   \mainpage
   \section build_sec Building the project
   \subsection CMSSW
   The project uses the libraries from the CMSSW framework. Once CMSSW is configured (cmsenv command) and from the main directory (above ./src) you only need to enter :
   \code
   make
   \endcode

   If everything goes fine, you should get a binary file called "AMSimulation".
 
   \subsection CUDA
   Some features of the program (pattern recognition) can use a GPU card to accelarate the computing. If you want to use this feature you will need:
   - Root (http://root.cern.ch) installed and configured ($ROOTSYS must be pointing on the installation directory and $ROOTSYS/bin must be in the PATH)
   - Boost (http://www.boost.org/) libraries and header files
   - Cuda SDK and code samples (https://developer.nvidia.com/cuda-downloads)
   - A GPU card compatible with CUDA

   You will need to modify the Makefile by setting the following variables:
   \code
   CUDA_ENABLED=true
   CUDA_ROOTDIR=/path/to/local/installation
   CUDA_EXAMPLEDIR=/path/to/cuda/code/samples/common/
   \endcode
   You can then launch the compilation :
   \code
   make
   \endcode
   To check the executable does support CUDA acceleration, you can try the command
   \code
   ./AMSimulation --help | grep GPU
   \endcode
   which should return :
   \code
   --useGPU              Use the GPU card to accelerate the pattern recognition 
   (needs cuda libraries and a configured GPU card)
   \endcode

   \section use_sec Using the program
   \subsection help Informations on how to use the program
   To get the list of options and parameters :
   \code
   ./AMSimulation --help
   \endcode
   \subsection generate Generating a pattern bank
   To generate a patterns bank, use the command :
   \code
   ./AMSimulation --generateBank
   \endcode
   All options can be stored in a file called amsimulation.cfg, here is an example :
   \code
   # Number of strips in a superstrip {16,32,64,128,256,512,1024}
   ss_size=32
   # Number of DC bits to use [0-3]
   dc_bits=3
   # Minimal PT of tracks used to generate a pattern
   pt_min=2
   # Maximal PT of tracks used to generate a pattern
   pt_max=10
   # Minimal ETA of tracks used to generate a pattern
   eta_min=-0.125
   # Maximal ETA of tracks used to generate a pattern
   eta_max=1.375
   #Maximal number of fake superstrips per pattern (used for hybrid/endcap sectors)
   maxFakeSStrips=0
   # Directory containing root files with single muon/antimuon events (local or RFIO)
   input_directory=rfio:/my/rfio/directory/
   # Output file name
   bank_name=testOutput.pbk
   # Coverage [0-1]
   coverage=0.9
   # Root file containing sectors definitions
   sector_file=sec_test.root
   # Index of the sector to be used
   sector_id=17
   # Layers used
   active_layers=6 7 9 10
   \endcode

   Each option contained in the configuration file can be overwritten via the command line, for example :
   \code
   ./AMSimulation --generateBank --ss_size=64
   \endcode
   will set the value of the ss_size option to 64, whatever is contained in the configuration file.

   If you have created 2 banks for the same sector with 2 different PT range (2-10 GeV and 10-50 GeV for example), you can merge the 2 files into a single one using the command :
   \code
   ./AMSimulation --MergeBanks --inputFile testPT2-10.pbk --secondFile testPT10-50.pbk --outputFile testPT2-10-50.pbk
   \endcode
   \subsection find Finding patterns in events
   To search for patterns in events, enter :
   \code
   ./AMSimulation --findPatterns --inputFile <path to Root File containing events (local or RFIO)> --bankFile <path to your pattern bank file> --outputFile <Root output file> --ss_threshold <minimum number of stubs to activate the pattern> --startEvent <Index of first event to analyse> --stopEvent <Index of last event to analyse>
   \endcode

  The program is using a virtual detector during the search process. The geometry of this detector (layers, ladders, Z modules and strips per segment) is contained in the CMSPatternLayer class but can be overwritten with a detector.cfg file located in the root directory of the program. Here is an example of the syntax :
  \code
  #layerID,nb ladders,nb Z modules,nb strips per segment
  5,16,15,1024
  6,24,15,1024
  7,36,15,1024
  8,48,14,1024
  9,60,14,1024
  10,76,14,1024
  \endcode

If you compiled the program using the cuda libraries, you can add the --useGPU flag to use the GPU device to perform the pattern recognition.

   \subsection merge Merging banks
   If you have created 2 banks for the same sector but with different PT range (for example 2 to 10 GeV and 10 to 100 GeV) you can merge the 2 files into a single one by using the command :   
   \code
   ./AMSimulation --MergeBanks --inputFile <PBK file of the first bank> --secondFile <PBK file of the second bank> --outputFile  <output PBK file name>
   \endcode

   Note that if you are using DC bits, the resulting bank may be smaller than the addition of the 2 original banks.

   \subsection alter Altering banks
   If your bank is using fake superstrips, you can split the PBK file into several files according to the number of fake superstrips in the patterns. You can then use a different threshold for each bank. To split a PBK file :
   \code
   ./AMSimulation --alterBank --maxFS=<max number of fake superstrips in a pattern> --minFS=<min number of fake superstrips in a pattern> --bankFile=<input bank file> --outputFile=<output bank file>
   \endcode

   \subsection view Viewing the content of a pattern bank
   You can display the patterns contained in a patterns bank file using the command:
   \code
   ./AMSimulation --printBank --bankFile <You patterns bank file>
   \endcode
   The interpretation of the DC bits is a bit tricky. Inside the pattern, the superstrip position is stored using gray encoding. If you have 'strip 3 (11X)' you have to convert the decimal value 3 to gray code : '10' then append the DC bits to that value : '1011X' which corresponds to 2 values : '10110' and '10111'. Those two gray values correspond to the decimal values 27 and 26 which are the high resolution indices of the superstrips in the module.


   If you are more interested in the encoding of the patterns you can display the patterns with a decimal representation of each pattern layer :
   \code
   ./AMSimulation --printBankBinary --bankFile <You patterns bank file>
   \endcode

   It should display one pattern per line. The interpretation process is the following : if you have '11014 (00X)' you have to convert the decimal value 11014 to binary : '00101 0110 000011 0'. The first 5 bits are the Z position of the module (00101 is 5 in decimal), then you have the index of the ladder (0110 is 6), the last bit is the segment (0 is 0). 000011 gives you the superstrip position and is encoded using gray code. You have to append the DC bits : 00001100X which corresponds to 2 values : 000011000 and 000011001. The decimal values of these gray encoded binary values are 16 and 17 which are the indices of the high resolution superstrips inside the module.

   You can also display the patterns using a format compatible with AM05 chips : for every pattern you will get one 18 bits value per line corresponding to the value you have to upload in the AM05 chip. DC bits encoding is taken care of. 
   \code
   ./AMSimulation --printBankAM05 --bankFile <Your patterns bank file>
   \endcode

   You can use the --nbActiveLayers option to select the patterns having a specific number of active layers (useful if you want to group patterns by threshold value). In trigger towers using 9 layers, the patterns will be split in 2 groups of 8 layers patterns.

   \author Guillaume Baulieu g.baulieu@ipnl.in2p3.fr
 **/

bool sorting (GradedPattern* p1, GradedPattern* p2) { return (*p2<*p1); }

void getLayers(vector<int> &l){
  cout<<"Enter the layer numbers (separated by spaces) :"<<endl;
  string result;
  getline(cin, result);
  std::istringstream is( result );
  int n;
  while( is >> n ) {
    l.push_back(n);
  }
}

float getPhiMin(){
  cout<<"Enter the minimum PHI0 :"<<endl;
  float result;
  cin>>result;
  if(result<0)
    result=0;
  return result;
}

float getPhiMax(){
  cout<<"Enter the maximum PHI0 :"<<endl;
  float result;
  cin>>result;
  if(result<0)
    result=0;
  return result;
}

float getEtaMin(){
  cout<<"Enter the minimum ETA :"<<endl;
  float result;
  cin>>result;
  if(result<0)
    result=0;
  return result;
}

float getEtaMax(){
  cout<<"Enter the maximum ETA :"<<endl;
  float result;
  cin>>result;
  if(result<0)
    result=0;
  return result;
}

string getSectorDefFileName(){
  cout<<"Enter the sector definition file name :"<<endl;
  string result;
  cin>>result;
  return result;
}

void getSectors(const vector<int> &layers, SectorTree &st){
  bool go=true;
  int sectorNb=1;
  while(go){
    cout<<"Enter the sector num "<<sectorNb<<" : "<<endl;
    Sector s(layers);
    for(unsigned int i=0;i<layers.size();i++){
      cout<<"\tFirst Ladder for layer "<<layers[i]<<" :"<<endl;
      int first;
      int nb;
      cin>>first;
      cout<<"\tNumber of ladders for layer "<<layers[i]<<" :"<<endl;
      cin>>nb;
      s.addLadders(layers[i],first,nb);
      if(!go)
	break;
    }
    if(go){
      st.addSector(s);
      sectorNb++;
    }
  }
}

vector< vector<int> > getRestrictions(const vector<int> &layers){
  vector<vector<int> > res;
  for(unsigned int i=0;i<layers.size();i++){
    vector<int> l;
    cout<<"\tLadder numbers for layer "<<layers[i]<<" (separated by spaces) :"<<endl;
    string result;
    getline(cin, result);
    getline(cin, result);
    std::istringstream is( result );
    int n;
    while( is >> n ) {
      if(n<0){
	break;
      }
      l.push_back(n);
    }
    if(l.size()>0)
      res.push_back(l);
    else
      return res;
  }
  return res;
}

void createAnalysis(SectorTree &st){
  vector<Sector*> list = st.getAllSectors();
  int nbLayers = 0;
  if(list.size()>0)
    nbLayers = list[0]->getNbLayers();
  vector<TH1I*> modulesPlot;
  
  for(int i=0;i<nbLayers;i++){
    ostringstream oss;
    oss<<"Layer "<<i;
    modulesPlot.push_back(new TH1I(oss.str().c_str(),"Module Z position", 14, 0, 14));
  }
  // We put all patterns in the same vector
  vector<GradedPattern*> allPatterns;
  int nbTracks=0;
  for(unsigned int i=0;i<list.size();i++){
    vector<GradedPattern*> patterns = list[i]->getPatternTree()->getLDPatterns();
    for(unsigned int j=0;j<patterns.size();j++){
      nbTracks+=patterns[j]->getGrade();
      allPatterns.push_back(patterns[j]);
    }
  }
  //sorting the patterns
  sort(allPatterns.begin(), allPatterns.end(), sorting);
  
  //float patterns[allPatterns.size()];
  //float tracks[allPatterns.size()];
  //float tracks_nb[allPatterns.size()];
  //float average_pt[allPatterns.size()];
  //int coveredTracks = 0;
  //TH1I* histo = new TH1I("Pattern Eff","Nb Tracks per pattern", allPatterns[0]->getGrade()+10, 0, allPatterns[0]->getGrade()+10);

  //Creates the layer plots
  /*
  for(unsigned int k=0;k<allPatterns.size();k++){
    patterns[k]=k;
    coveredTracks+=allPatterns[k]->getGrade();
    tracks[k]=coveredTracks*100/(float)nbTracks;
    histo->Fill(allPatterns[k]->getGrade());
    average_pt[k]=allPatterns[k]->getAveragePt();
    tracks_nb[k]=allPatterns[k]->getGrade();
    for(int j=0;j<nbLayers;j++){
      CMSPatternLayer* pl = (CMSPatternLayer*)allPatterns[k]->getLayerStrip(j);
      modulesPlot[j]->Fill(pl->getModule());
    }
  }
  histo->SetFillColor(41);
  histo->Write();
  delete histo;
  */
  for(unsigned int k=0;k<list.size();k++){
    vector<int> PT = list[k]->getPatternTree()->getPTHisto();
    TH1I* pt_histo = new TH1I("PT sector "+k,"PT of pattern generating tracks", 110, 0, 110);
    for(int i=0;i<101;i++){
      //cout<<PT[i]<<"-";
      for(int j=0;j<PT[i];j++){
	pt_histo->Fill(i);
      }
    }
    //cout<<endl;
    pt_histo->SetFillColor(41);
    pt_histo->Write();
    delete pt_histo;

    int pattID;
    int out_nbTracks;
    float out_minZ0, out_maxZ0;
    float out_minEta, out_maxEta;
    float out_minPT, out_maxPT;
    TTree *OUT    = new TTree("PatternInfo", "Description of patterns");
    OUT->Branch("pattern_id",    &pattID);
    OUT->Branch("nbTracks",      &out_nbTracks);
    OUT->Branch("minZ0",         &out_minZ0);
    OUT->Branch("maxZ0",         &out_maxZ0);
    OUT->Branch("minEta",        &out_minEta);
    OUT->Branch("maxEta",        &out_maxEta);
    OUT->Branch("minPT",         &out_minPT);
    OUT->Branch("maxPT",         &out_maxPT); 
  
    vector<GradedPattern*> patterns = list[k]->getPatternTree()->getLDPatterns();
    for(unsigned int i=0;i<patterns.size();i++){
      const PatternInfo& pi = patterns[i]->getStatisticalInformations();
      pattID=i;

      out_nbTracks = pi.getNbTracks();

      out_minZ0 = pi.getMinZ0();
      out_maxZ0 = pi.getMaxZ0();

      out_minEta = pi.getMinEta();
      out_maxEta = pi.getMaxEta();

      out_minPT = pi.getMinPT();
      out_maxPT = pi.getMaxPT();

      OUT->Fill();
    }
    OUT->Write();
  }

  /*
  TGraph* nbPatt = new TGraph(allPatterns.size(),patterns,tracks);
  nbPatt->GetXaxis()->SetTitle("Patterns bank size");
  nbPatt->GetYaxis()->SetTitle("Tracks covered (%)");
  nbPatt->SetTitle("Patterns coverage");
  nbPatt->Write();
  delete nbPatt;
  */
  /*
  TGraph* ptFreq = new TGraph(allPatterns.size(),tracks_nb,average_pt);
  ptFreq->GetXaxis()->SetTitle("nb tracks per pattern");
  ptFreq->GetYaxis()->SetTitle("average Pt per pattern");
  ptFreq->SetTitle("Average Pt / pattern size");
  ptFreq->Write();
  delete ptFreq;
  */
  for(int i=0;i<nbLayers;i++){
    modulesPlot[i]->Write();
    delete modulesPlot[i];
  }

  for(unsigned int k=0;k<allPatterns.size();k++){
    delete allPatterns[k];
  }

}

 /**
     \brief Display the layers, ladders and modules hit by the tracks contained in the given file.
     \param fileName The name of the root file
     \param tracker_layers Gives the number of the layers in the tracker
     \param restriction If the vector contains data, only tracks going throught these ladders will be tacken into account
     \param phi_min Minimum value of PHI0 for the selected tracks
     \param phi_max Maximum value of PHI0 for the selected tracks
     \param eta_min Minimum value of etaGEN for the selected tracks
     \param eta_max Maximum value of etaGEN for the selected tracks
  **/
void createFromSimu(string fileName, vector<int> tracker_layers, vector< vector<int> > restriction, float phi_min, float phi_max, float eta_min, float eta_max){

  map<int, set<int> > usedLadders;
  map<int, map<int, set<int> > > usedModules;

  TChain* TT = new TChain("L1TrackTrigger");
  TT->Add(fileName.c_str());

  //--> Signification (et dimension) des variables

  // Stub info (fait a partir de paires de clusters matches)

  //static const int      m_stub_MAX    = 10000;     // Nombre maximal de stubs
  
  int m_stub;
  vector<int>           m_stub_layer;  // Layer du stub (5 a 10 pour les 6 layers qui nous interessent)
  vector<int>           m_stub_module; // Position en Z du module contenant le stub
  vector<int>           m_stub_ladder; // Position en PHI du module contenant le stub
  vector<int>           m_stub_seg;    // Segment du module contenant le stub
  vector<int>           m_stub_strip;  // Strip du cluster interne du stub
  vector<float>         m_stub_pxGEN;  // pxGEN de la particule originelle
  vector<float>         m_stub_pyGEN;  // pyGEN de la particule originelle
  vector<float>         m_stub_etaGEN;  // etaGEN de la particule originelle
 
  vector<int>           *p_m_stub_layer = &m_stub_layer;
  vector<int>           *p_m_stub_module = &m_stub_module;
  vector<int>           *p_m_stub_ladder = &m_stub_ladder;
  vector<int>           *p_m_stub_seg = &m_stub_seg;
  vector<int>           *p_m_stub_strip = &m_stub_strip;
  vector<float>         *p_m_stub_pxGEN = &m_stub_pxGEN;
  vector<float>         *p_m_stub_pyGEN = &m_stub_pyGEN;
  vector<float>         *p_m_stub_etaGEN = &m_stub_etaGEN;
  
  TT->SetBranchAddress("STUB_n",         &m_stub);
  TT->SetBranchAddress("STUB_layer",     &p_m_stub_layer);
  TT->SetBranchAddress("STUB_module",    &p_m_stub_module);
  TT->SetBranchAddress("STUB_ladder",    &p_m_stub_ladder);
  TT->SetBranchAddress("STUB_seg",       &p_m_stub_seg);
  TT->SetBranchAddress("STUB_strip",     &p_m_stub_strip);
  TT->SetBranchAddress("STUB_pxGEN",     &p_m_stub_pxGEN);
  TT->SetBranchAddress("STUB_pyGEN",     &p_m_stub_pyGEN);
  TT->SetBranchAddress("STUB_etaGEN",    &p_m_stub_etaGEN);

  int n_entries_TT = TT->GetEntries();

  int nbInLayer=0;

  int minLayer = *(min_element(tracker_layers.begin(),tracker_layers.end()));

  int layers[tracker_layers.size()];
  int ladder_per_layer[tracker_layers.size()];
  int module_per_layer[tracker_layers.size()];
  
  int nbUsedTracks = 0;

  float found_phi_min = 1000;
  float found_phi_max = -1000;

  float phi0=0;

  for(int i=0;i<n_entries_TT;i++){
    TT->GetEntry(i);
    for(unsigned int j=0;j<tracker_layers.size();j++){
      layers[j]=-1;
    }
    for(unsigned int j=0;j<tracker_layers.size();j++){
      ladder_per_layer[j]=-1;
    }
    for(unsigned int j=0;j<tracker_layers.size();j++){
      module_per_layer[j]=-1;
    }
    //check the layers of the stubs
    for(int j=0;j<m_stub;j++){

      float pt_GEN = sqrt(m_stub_pxGEN[j]*m_stub_pxGEN[j]+m_stub_pyGEN[j]*m_stub_pyGEN[j]);

      if(pt_GEN<2){//we only need particules with PT>2
	continue;
      }

      phi0 = atan2(m_stub_pyGEN[j], m_stub_pxGEN[j]);
     
      if(phi0<phi_min){//The stub is coming from a particule outside the considered sector
	continue;
      }

      if(phi0>phi_max){//The stub is coming from a particule outside the considered sector
	continue;
      }

      if(m_stub_etaGEN[j]<eta_min){//The stub is coming from a particule outside the considered sector
	continue;
      }

      if(m_stub_etaGEN[j]>eta_max){//The stub is coming from a particule outside the considered sector
	continue;
      }

      int layer = m_stub_layer[j];
      if(((unsigned int)(layer-minLayer)<tracker_layers.size()) // layer is not above the last considered layer
	 && layers[layer-minLayer]!=-1){//we have 2 stubs on the same layer-> problem
	layers[layer-minLayer]=-1;
	continue;
      }
      if(find(tracker_layers.begin(),tracker_layers.end(), layer)!=tracker_layers.end()){ // is this layer in the layer list?
	layers[layer-minLayer]=j;
	ladder_per_layer[layer-minLayer]=m_stub_ladder[j];
	module_per_layer[layer-minLayer] = m_stub_module[j];
      }
    }
    /**************************************
    Selection on the stubs/layer
    We need at least one stub per layer
    **************************************/
    bool missing_stub = false;
    for(unsigned int j=0;j<tracker_layers.size();j++){
      if(layers[j]==-1){
	missing_stub=true;
	
      }
    }
    if(missing_stub)
      continue;//no stub on each layer -> drop the event    
    nbInLayer++;
    //restriction to some ladders
    bool useTrack = true;
    for(unsigned int j=0;j<restriction.size();j++){
      vector<int> rLayer = restriction[j];
      if(find(rLayer.begin(), rLayer.end(),ladder_per_layer[j])==rLayer.end()){//the track is no going throught the right ladder
	useTrack = false;
	break;
      }
    }
    if(useTrack){
      nbUsedTracks++;
      if(phi0<found_phi_min)
	found_phi_min=phi0;
      if(phi0>found_phi_max)
	found_phi_max=phi0;
      for(unsigned int j=0;j<tracker_layers.size();j++){
	int layer_id=tracker_layers[j];

	usedLadders[layer_id].insert(CMSPatternLayer::getLadderCode(layer_id, ladder_per_layer[j]));
	usedModules[layer_id][ladder_per_layer[j]].insert(CMSPatternLayer::getModuleCode(layer_id, module_per_layer[j]));

      }
    }
    if(nbUsedTracks>100000)
      break;
  }

  cout<<"Nb Events : "<<n_entries_TT<<endl;
  cout<<"Nb Events with stubs on all layers : "<<nbInLayer<<endl;
  cout<<"PHI min : "<<found_phi_min<<endl;
  cout<<"PHI max : "<<found_phi_max<<endl;

  for(map<int, set<int> >::const_iterator it_layer=usedLadders.begin();it_layer!=usedLadders.end();it_layer++){
    cout<<"Layer "<<it_layer->first<<" : "<<endl;
    for(set<int>::const_iterator it_lad=it_layer->second.begin();it_lad!=it_layer->second.end();it_lad++){
      cout<<"    "<<*it_lad<<" : ";
      set<int> modules = usedModules[it_layer->first][*it_lad];
      for(set<int>::const_iterator it_mod=modules.begin();it_mod!=modules.end();it_mod++){
	cout<<*it_mod<<" ";
      }
      cout<<endl;
    }
    cout<<endl;
  }

  delete TT;
}

// Get the first ladder and number of ladders from a ladders list
void getOrderData(vector<int> list, int* first, int* nb){
  sort(list.begin(),list.end());
  *nb = (int)list.size();
  if((*nb)>0){
    int index = list.size()-1;
    *first = list[index];//greatest value
    
    //we decrease the greatest value and stop as soon as there is a gap
    while(index>0 && list[index-1]==(*first)-1){
      index--;
      (*first)--;
    }
  }
}

// Get the first module and number of modules from a modules list
void getOrderDataForModules(vector<int> list, int maxNbModules, int* first, int* nb){
  sort(list.begin(),list.end());

  //Some sectors have holes... If we have one big gap and others are just 1 missing module->we fill these
  int maxGap=0;
  vector<int> newList;
  for(unsigned int i=0;i<list.size()-1;i++){
    newList.push_back(list[i]);
    int dif=list[i+1]-list[i];
    if(dif>maxGap)
      maxGap=dif;
    if(dif==2){
      newList.push_back(list[i]+1);
    }
  }
  newList.push_back(list[list.size()-1]);
  int dif=maxNbModules-list[list.size()-1]+list[0];
  if(dif>maxGap)
    maxGap=dif;
  if(dif==2){
    newList.push_back(list[list.size()-1]+1);
  }
  
  if(maxGap>2)
    list=newList;

  *nb = (int)list.size();
  if((*nb)>0){
    int index = list.size()-1;
    *first = list[index];//greatest value
    
    //we decrease the greatest value and stop as soon as there is a gap
    while(index>0 && list[index-1]==(*first)-1){
      index--;
      (*first)--;
    }
  }
}

void createSectorFromRootFile(SectorTree* st, string fileName, vector<int> layers, int sector_id){

  Sector s(layers);
  TChain* tree = NULL;
  vector<int> modules;
  vector<float> coords;

  if(fileName.substr(fileName.length()-4,fileName.length()).compare(".csv")==0){
    ifstream is(fileName.c_str());
    string line;
    int line_index=0;
    if (is.is_open()){
      while ( is.good() ){
	getline (is,line);
	if(line_index==sector_id+1){ // line describing the sector we want
	  std::stringstream ss(line);
	  std::string item;
	  int column_index=0;
	  while (std::getline(ss, item, ',')) {
	    if(column_index>1){   
	      int number=0;
	      std::istringstream ss( item );
	      ss >> number;
	      if(number!=0){
		//cout<<number<<endl;
		modules.push_back(number);
	      }
	    }
	    column_index++;
	  }
	}
	line_index++;
      }
      is.close();
    }
  }
  else{
    tree = new TChain("Sectors");
    tree->Add(fileName.c_str());
    vector< vector<int> > m_mod_tot;  // Secteurs dans le endcap
    vector< vector<float> > m_coords;
    vector< vector<int> > * p_m_mod_tot =  &m_mod_tot;
    vector< vector<float> > * p_m_coords = &m_coords;
    tree->SetBranchAddress("sectors",   &p_m_mod_tot);
    tree->SetBranchAddress("sectors_coord",    &p_m_coords);
    
    tree->GetEntry(0);
  
    modules = m_mod_tot[sector_id];
    cout<<"\tmodule number : "<<modules.size()<<endl;

    coords = m_coords[sector_id];
  }

  map<int, vector<int> > ladders_from_layer;
  map<int, map<int, vector<int> > > modules_from_ladder;

  for(unsigned int i=0;i<modules.size();i++){

    ostringstream oss;
    oss<<std::setfill('0');
    oss<<setw(6)<<modules[i];
    
    int layer,ladder,module;
    istringstream ss_layer(oss.str().substr(0,2));
    ss_layer>>layer;
    istringstream ss_ladder(oss.str().substr(2,2));
    ss_ladder>>ladder;
    ladder = CMSPatternLayer::getLadderCode(layer,ladder);

    ////// TMP FIX FOR TKLAYOUT NUMBERING (10 JUL 2013)
    
    if(layer<11){
      int tmp_nb_ladders = CMSPatternLayer::getNbLadders(layer);
      ladder = (ladder+tmp_nb_ladders*1/4) % tmp_nb_ladders;
    } 

    //////////////////////////////////////////////////

    istringstream ss_module(oss.str().substr(4,2));
    ss_module>>module;
    module = CMSPatternLayer::getModuleCode(layer,module);

    vector<int> tmp_ladders = ladders_from_layer[layer];
    if(find(tmp_ladders.begin(),tmp_ladders.end(),ladder)==tmp_ladders.end()){
      ladders_from_layer[layer].push_back(ladder);
    }

    vector<int> tmp_modules = modules_from_ladder[layer][ladder];
    if(find(tmp_modules.begin(),tmp_modules.end(),module)==tmp_modules.end()){
      modules_from_ladder[layer][ladder].push_back(module);
    }
  }

  for(unsigned int i=0;i<layers.size();i++){
    vector<int> tmp_ladders = ladders_from_layer[layers[i]];
    int first=0, nb=0;
    getOrderData(tmp_ladders, &first, &nb);
    s.addLadders(layers[i], first, nb);
    for(unsigned int k=0;k<tmp_ladders.size();k++){
      vector<int> tmp_modules = modules_from_ladder[layers[i]][tmp_ladders[k]];
      getOrderDataForModules(tmp_modules, CMSPatternLayer::getNbModules(layers[i], tmp_ladders[k]), &first, &nb);
      s.addModules(layers[i], tmp_ladders[k], first, nb);
    }
  }

  if(tree!=NULL){
    delete tree;
    cout<<"\tPHI min : "<<coords[0]<<endl;  
    cout<<"\tPHI max : "<<coords[1]<<endl;  
    cout<<"\tETA min : "<<coords[2]<<endl;  
    cout<<"\tETA max : "<<coords[3]<<endl;  
    cout<<endl;
  }

  s.setOfficialID(sector_id);

  st->addSector(s);
}

int main(int av, char** ac){
  namespace po = boost::program_options;
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("generateBank", "Generates a pattern bank from root simulation file (needs --ss_size --dc_bits --pt_min --pt_max --eta_min --eta_max --coverage --input_directory --bank_name --sector_file --sector_id --active_layers)")
    ("testSectors", "Get the tracks sectors")
    ("MergeSectors", "Merge 2 root files having same events but different sectors (needs --inputFile --secondFile and --outputFile)")
    ("MergeBanks", "Merge 2 bank files having only 1 sector (needs --inputFile --secondFile and --outputFile)")
    ("buildFitParams", "Computes the Fit parameters for the given bank using tracks from the given directory (needs --bankFile, --input_directory and --outputFile)")
    ("findPatterns", "Search for patterns in an event file (needs --ss_threshold --inputFile, --bankFile, --outputFile, --startEvent and --stopEvent)")
#ifdef USE_CUDA
    ("useGPU", "Use the GPU card to accelerate the pattern recognition (needs cuda libraries and a configured GPU card)")
#endif
    ("printBank", "Display all patterns from a bank (needs --bankFile)")
    ("printBankBinary", "Display all patterns from a bank using a decimal representation of the binary values (needs --bankFile)")
    ("printBankAM05", "Display all patterns from a bank with the format used for the AM05 chip (needs --bankFile and optionaly --nbActiveLayers)")
    ("testCode", "Dev tests")
    ("analyseBank", "Creates histograms from a pattern bank file (needs --bankFile and --outputFile)")
    ("stubsToSuperstrips", "Display each stub of an event file as a superstrip (used for tests) (needs --inputFile, --bankFile, --startEvent, --stopEvent)")
    ("alterBank", "Creates a new bank from an existing one, the existing bank is not modified (used with --bankFile and --outputFile and --truncate or --minFS and/or --maxFS)")
    ("inputFile", po::value<string>(), "The file to analyse")
    ("secondFile", po::value<string>(), "Second file to merge")
    ("bankFile", po::value<string>(), "The patterns bank file to use")
    ("outputFile", po::value<string>(), "The root output file")
    ("ss_threshold", po::value<int>(), "The minimum number of hit superstrips to activate a pattern")
    ("ss_missingHits", po::value<int>(), "The maximum number of non activated layers to activate a pattern. --ss_threshold is used as a mandatory minimum value.")
    ("startEvent", po::value<int>(), "The first event index")
    ("stopEvent", po::value<int>(), "The last event index")
    ("decode", po::value<int>(), "Decode the given super strip")
    ("ss_size", po::value<string>(), "Number of strips in a super strip {16,32,64,128,256,512,1024}. Either one value for all detector or one value per layer")
    ("dc_bits", po::value<int>(), "Number of used DC bits [0-3]")
    ("pt_min", po::value<int>(), "Only tracks having a greater PT will be used to generate a pattern")
    ("pt_max", po::value<int>(), "Only tracks having a smaller PT will be used to generate a pattern")
    ("eta_min", po::value<float>(), "Only tracks having a greater ETA will be used to generate a pattern")
    ("eta_max", po::value<float>(), "Only tracks having a smaller ETA will be used to generate a pattern")
    ("maxFakeSStrips", po::value<int>(), "The maximum number of fake superstrips used in a pattern (0=desactivated)")
    ("coverage", po::value<float>(), "Value to reach to stop the process [0-1]")
    ("input_directory", po::value<string>(), "The directory containing the single particule root files (local or RFIO)")
    ("sector_file", po::value<string>(), "The file (.root or .csv) containing the sectors definition")
    ("sector_id", po::value<int>(), "The index of the sector to use in the sector file. In a CSV file the first sector has index 0.")
    ("active_layers", po::value<string>(), "The layers to use in the sector (8 at most). If a layer ID is prefixed with '+' it will never contain a fake superstrip, if it's prefixed with '-' it will always contain a fake superstrip.")
    ("bank_name", po::value<string>(), "The bank file name")    
    ("minFS", po::value<int>(), "Used with --alterBank : only patterns with at least minFS fake stubs will be kept in the new bank")
    ("maxFS", po::value<int>(), "Used with --alterBank : only patterns with at most maxFS fake stubs will be kept in the new bank")
    ("truncate", po::value<int>(), "Used with --alterBank : gives the number of patterns to keep in the new bank, starting with the most used ones.")
    ("nbActiveLayers", po::value<int>(), "Used with --printBankAM05 : only patterns with this exact number of active layers will be printed")
    ;
     
  po::variables_map vm;
  std::ifstream in_file( "amsimulation.cfg" ); 
  po::store(po::parse_command_line(av, ac, desc), vm);
  po::store(po::parse_config_file(in_file, desc), vm);
  po::notify(vm);    
  
  if (vm.count("help")) {
    cout << desc << "\n";

    return 1;
  }

  if (vm.count("analyseBank")) {
    SectorTree st;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
      //Decompression
      boost::iostreams::filtering_stream<boost::iostreams::input> f;
      f.push(boost::iostreams::gzip_decompressor());
      try { 
	f.push(ifs);
	boost::archive::text_iarchive ia(f);
	ia >> st;
      }
      catch (boost::iostreams::gzip_error& e) {
	if(e.error()==4){//file is not compressed->read it without decompression
	  std::ifstream new_ifs(vm["bankFile"].as<string>().c_str());
	  boost::archive::text_iarchive ia(new_ifs);
	  ia >> st;
	}
      }
    }
    TFile f( vm["outputFile"].as<string>().c_str(), "recreate");
    createAnalysis(st);
  } else if (vm.count("generateBank")) {
    
    vector<int> layers;
    SectorTree st;
    string stripSizes="";
    int dcBits=0;
    string partDirName="";
    string bankFileName="";
    string rootFileName="";
    string activeLayers = "";
    vector<int> active_layers;
    vector<int> forced_layers;
    vector<int> desactivated_layers;
    float threshold=0;
    float min=0;
    float max=0;
    float minEta=0;
    float maxEta=0;
    int maxNbFake=0;
    int sector_tklayout_id=0;
    map<int,pair<float,float> > eta = CMSPatternLayer::getLayerDefInEta();
    
    try{
      stripSizes=vm["ss_size"].as<string>();
      cout<<"Superstrip sizes : "<<stripSizes<<endl;
      dcBits=vm["dc_bits"].as<int>();
      cout<<"DC bits number : "<<dcBits<<endl;
      min=vm["pt_min"].as<int>();
      cout<<"PT min : "<<min<<endl;
      max=vm["pt_max"].as<int>();
      cout<<"PT max : "<<max<<endl;
      minEta=vm["eta_min"].as<float>();
      cout<<"ETA min : "<<minEta<<endl;
      maxEta=vm["eta_max"].as<float>();
      cout<<"ETA max : "<<maxEta<<endl;
      maxNbFake=vm["maxFakeSStrips"].as<int>();
      cout<<"Max number of fake superstrips : "<<maxNbFake<<endl;
      cout<<"Coverage : "<<vm["coverage"].as<float>()<<endl;
      partDirName=vm["input_directory"].as<string>();
      cout<<"Using particules from "<<partDirName<<endl;
      bankFileName=vm["bank_name"].as<string>();
      cout<<"Output file name is "<<bankFileName<<endl;
      activeLayers=vm["active_layers"].as<string>();
      cout<<"Using layers "<<activeLayers<<endl;
      sector_tklayout_id=vm["sector_id"].as<int>();
      cout<<"Using sector "<<sector_tklayout_id<<" from "<<vm["sector_file"].as<string>()<<endl;

      //Get the active/forced/inactive layers
      //layersWithoutSigns : list of all layers

      //forceLayers : layers prefixed with a '+' : no fake stub will be allowed
      //if there is no stub on this layer -> there will be no pattern

      //inactiveLayers : layers prefixed with a '-' : only fake stubs on this layer
      //even if there is a stub it will be replaced with a fake one

      string layersWithoutSigns = activeLayers;
      string forceLayers="";
      string inactiveLayers="";
      size_t found = layersWithoutSigns.find("+");
      while (found!=string::npos){
	layersWithoutSigns.erase(found,1);//remove the '+'
	size_t endIndex = layersWithoutSigns.find(" ", found);
	if(endIndex!=string::npos)
	  endIndex=endIndex-found+1;
	forceLayers.append(layersWithoutSigns.substr(found,endIndex));//add the layer number and the following space
	found = layersWithoutSigns.find("+"); // search for the next '+'
      }
      found = layersWithoutSigns.find("-");
      while (found!=string::npos){
	layersWithoutSigns.erase(found,1);
	size_t endIndex = layersWithoutSigns.find(" ", found);
	if(endIndex!=string::npos)
	  endIndex=endIndex-found+1;
	inactiveLayers.append(layersWithoutSigns.substr(found,endIndex));
	found = layersWithoutSigns.find("-");
      }

      std::istringstream is( layersWithoutSigns );
      int n;
      while( is >> n ) {
	active_layers.push_back(n);
      }
      std::istringstream fl( forceLayers );
      while( fl >> n ) {
	forced_layers.push_back(n);
      }
      std::istringstream il( inactiveLayers );
      while( il >> n ) {
	desactivated_layers.push_back(n);
      }

      //remove the force_layers from the eta list -> no fake stub will be added
      for(unsigned int i=0;i<forced_layers.size();i++){
	map<int,pair<float,float> >::iterator it = eta.find(forced_layers[i]);
	eta.erase(it);
      }

      //change the eta definition of the desactivated layers
      //they will not be reachable and a fake stub will be created
      for(unsigned int i=0;i<desactivated_layers.size();i++){
	eta[desactivated_layers[i]].first=10;
	eta[desactivated_layers[i]].second=10;
      }

      vector<int> layer_superstrip_size;
      std::istringstream strips( stripSizes );
      while( strips >> n ) {
	layer_superstrip_size.push_back(n);
      }
      if(layer_superstrip_size.size()==1){
	st.setSuperStripSize(layer_superstrip_size[0]);
      }
      else if(layer_superstrip_size.size()==active_layers.size()){//one size per layer
	for(unsigned int ln=0;ln<active_layers.size();ln++){
	  st.setSuperStripSize(layer_superstrip_size[ln],active_layers[ln]);
	}
      }
      else{
	st.setSuperStripSize(-1);
      }
      

      size_t end_index = bankFileName.find(".pbk");
      if(end_index==string::npos)
	end_index=bankFileName.length()-4;
      rootFileName = bankFileName.substr(0,end_index)+"_report.root";
      threshold=vm["coverage"].as<float>();
      createSectorFromRootFile(&st,vm["sector_file"].as<string>(), active_layers, sector_tklayout_id);
    }
    catch(boost::bad_any_cast e){
      cout<<"At least one option is missing! Please check : "<<endl;
      cout<<desc<<endl;
      return -1;
    }
    
    vector<Sector*> list = st.getAllSectors();
    cout<<"Sector :"<<endl;
    for(unsigned int i=0;i<list.size();i++){
      cout<<*list[i];
      cout<<endl;
    }

    if(active_layers.size()>9){
      cout<<"ERROR : your sector contains "<<active_layers.size()<<" layers : maximum number of layers is 9!"<<endl;
      return -1;
    }
    
    PatternGenerator pg;
    pg.setLayers(active_layers);
    pg.setInactiveLayers(desactivated_layers);
    pg.setParticuleDirName(partDirName);
    pg.setMinPT(min);
    pg.setMaxPT(max);
    pg.setMinEta(minEta);
    pg.setMaxEta(maxEta);
    pg.setMaxFakeSuperstrips(maxNbFake);
    TFile f(rootFileName.c_str(), "recreate");
    pg.setVariableResolution(dcBits);
    pg.generate(&st, 100000, threshold, eta);


    if(pg.getVariableResolutionState()>0){
      cout<<"LD Patterns : "<<st.getLDPatternNumber()<<endl;
    }
  
    cout<<"Saving SectorTree...";
    {
      const SectorTree& ref = st;
      std::ofstream ofs(bankFileName.c_str());
      // Compression part
      boost::iostreams::filtering_stream<boost::iostreams::output> f;
      f.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
      f.push(ofs);
      //
      boost::archive::text_oarchive oa(f);
      oa << ref;
      cout<<"done."<<endl;
    }
    
  }
  else if(vm.count("testSectors")) {
    vector<int> layers;
    getLayers(layers);
    float phi0_min = getPhiMin();
    float phi0_max = getPhiMax();
    float eta_min = getEtaMin();
    float eta_max = getEtaMax();
    vector< vector<int> > restriction = getRestrictions(layers);
    string fn = getSectorDefFileName();
    createFromSimu(fn, layers, restriction, phi0_min, phi0_max, eta_min, eta_max);
  }
  else if(vm.count("decode")) {
    CMSPatternLayer p;
    int val = vm["decode"].as<int>();
    p.setIntValue(val);
    cout<<p.toString()<<endl;
  }
  else if(vm.count("stubsToSuperstrips")) {
    SectorTree st;
    cout<<"Loading pattern bank..."<<endl;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
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
	  std::ifstream new_ifs(vm["bankFile"].as<string>().c_str());
	  boost::archive::text_iarchive ia(new_ifs);
	  ia >> st;
	}
      }
    }
    PatternFinder pf(0, &st,  vm["inputFile"].as<string>().c_str(),  "none");
    {
      int start = vm["startEvent"].as<int>();
      int stop = vm["stopEvent"].as<int>();
      pf.displayEventsSuperstrips(start, stop);
    }
  }
  else if(vm.count("findPatterns")) {
    SectorTree st;
    cout<<"Loading pattern bank..."<<endl;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
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
	  std::ifstream new_ifs(vm["bankFile"].as<string>().c_str());
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

#ifdef USE_CUDA
    if(vm.count("useGPU")){
      deviceDetector d_detector;
      patternBank d_pb;
      deviceStubs d_stubs;
      deviceParameters d_param;
      resetCard();
      allocateDetector(&d_detector);
      allocateBank(&d_pb,st.getAllSectors()[0]->getLDPatternNumber());
      allocateStubs(&d_stubs);
      allocateParameters(&d_param);

      st.getAllSectors()[0]->linkCuda(&d_pb,&d_detector);

      PatternFinder pf(vm["ss_threshold"].as<int>(), &st,  vm["inputFile"].as<string>().c_str(),  vm["outputFile"].as<string>().c_str(), 
		       &d_pb, &d_detector, &d_param); 
      {
	boost::progress_timer t;
	int start = vm["startEvent"].as<int>();
	int stop = vm["stopEvent"].as<int>();
	pf.findCuda(start, stop, &d_stubs);
	cout<<"Time used to analyse "<<stop-start+1<<" events : "<<endl;
      }

      freeParameters(&d_param);
      freeDetector(&d_detector);
      freeBank(&d_pb);
      freeStubs(&d_stubs);
      resetCard();
    }
    else{
#endif
      int nbMissingHit=0;
      int threshold=0;
      if(vm.count("ss_missingHits")){
	  nbMissingHit=vm["ss_missingHits"].as<int>();
	  threshold=vm["ss_threshold"].as<int>();
      }
      else{
	nbMissingHit=-1;
	threshold=vm["ss_threshold"].as<int>();
      }
      PatternFinder pf(threshold, &st,  vm["inputFile"].as<string>().c_str(),  vm["outputFile"].as<string>().c_str());
      {
	boost::progress_timer t;
	int start = vm["startEvent"].as<int>();
	int stop = vm["stopEvent"].as<int>();

	if(vm.count("ss_missingHits")){
	    pf.useMissingHitThreshold(nbMissingHit);
	}
	pf.find(start, stop);
	cout<<"Time used to analyse "<<stop-start+1<<" events : "<<endl;
      }
#ifdef USE_CUDA
    }
#endif
  }
  else if(vm.count("buildFitParams")) {
    SectorTree st;
    cout<<"Loading pattern bank..."<<endl;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
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
	  std::ifstream new_ifs(vm["bankFile"].as<string>().c_str());
	  boost::archive::text_iarchive ia(new_ifs);
	  ia >> st;
	}
      }
    }

    map<int,pair<float,float> > eta_limits = CMSPatternLayer::getLayerDefInEta();

    PrincipalFitGenerator pfg(vm["input_directory"].as<string>().c_str(), &st);
    pfg.generate(eta_limits, 2, 100, 0, 0.87);
    
    cout<<"Saving SectorTree...";
    {
      const SectorTree& ref = st;
      std::ofstream ofs(vm["outputFile"].as<string>().c_str());
      // Compression part
      boost::iostreams::filtering_stream<boost::iostreams::output> f;
      f.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
      f.push(ofs);
      //
      boost::archive::text_oarchive oa(f);
      oa << ref;
      cout<<"done."<<endl;
    }
    
  }
  else if(vm.count("printBank")) {
    SectorTree st;
    cout<<"Loading pattern bank..."<<endl;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
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
	  std::ifstream new_ifs(vm["bankFile"].as<string>().c_str());
	  boost::archive::text_iarchive ia(new_ifs);
	  ia >> st;
	}
      }
    }
    vector<Sector*> sectors = st.getAllSectors();
    for(unsigned int i=0;i<sectors.size();i++){
      Sector* mySector = sectors[i];
      vector<GradedPattern*> patterns = mySector->getPatternTree()->getLDPatterns();
      for(unsigned int j=0;j<patterns.size();j++){
	Pattern* p = patterns[j];
	for(int k=0;k<p->getNbLayers();k++){
	  PatternLayer* mp = p->getLayerStrip(k);
	  cout<<((CMSPatternLayer*)mp)->toString()<<" - ";
	}
	cout<<endl;
      }
    }
  }
  else if(vm.count("printBankBinary")) {
    SectorTree st;
    cout<<"Loading pattern bank..."<<endl;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
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
	  std::ifstream new_ifs(vm["bankFile"].as<string>().c_str());
	  boost::archive::text_iarchive ia(new_ifs);
	  ia >> st;
	}
      }
    }
    vector<Sector*> sectors = st.getAllSectors();
    for(unsigned int i=0;i<sectors.size();i++){
      Sector* mySector = sectors[i];
      vector<GradedPattern*> patterns = mySector->getPatternTree()->getLDPatterns();
      for(unsigned int j=0;j<patterns.size();j++){
	Pattern* p = patterns[j];
	for(int k=0;k<p->getNbLayers();k++){
	  PatternLayer* mp = p->getLayerStrip(k);
	  cout<<((CMSPatternLayer*)mp)->toStringBinary()<<" - ";
	}
	cout<<endl;
      }
    }
  } 
  else if(vm.count("printBankAM05")) {
    SectorTree st;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
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
	  std::ifstream new_ifs(vm["bankFile"].as<string>().c_str());
	  boost::archive::text_iarchive ia(new_ifs);
	  ia >> st;
	}
      }
    }
    vector<Sector*> sectors = st.getAllSectors();
    for(unsigned int i=0;i<sectors.size();i++){
      Sector* mySector = sectors[i];
      int sector_id = mySector->getOfficialID();
      vector<GradedPattern*> patterns = mySector->getPatternTree()->getLDPatterns();
      vector<int> layers = mySector->getLayersID();
      int expected_active_layers = -1;
      bool hybrid_sector = layers.size()>8;//we have more layers than input buses
      bool first_pass=true;

      int biggestID = -1;//last layer of the sector
      int biggestBarrelID = -1;//last barrel layer of the sector
      int biggestBarrelIndex = -1;//index of last barrel layer of the sector
      for(unsigned int j=0;j<layers.size();j++){
	biggestID=layers[j];
	if(layers[j]<11){
	  biggestBarrelID=layers[j];
	  biggestBarrelIndex=j;
	}
      }

      if(vm.count("nbActiveLayers"))
	expected_active_layers = vm["nbActiveLayers"].as<int>(); // we only want patterns with this number of active layers (ie layers without fake superstrips)

      cout<<"** Bank of patterns formatted for AM05 chip. Patterns are displayed with 1 line per layer (from innermost to outermost layer) and separated with an empty line."<<endl;
      cout<<"** Each line contains the decimal representation of the 18 bits value to store in the AM chip followed by the number of DC bits configured on that layer."<<endl;
      cout<<"** The format of the 16 bits superstrips to send to the AM chip is (all positions are relative to the trigger tower): "<<endl;
      int total_nb_bits = PatternLayer::getSizeFromMask(CMSPatternLayer::MOD_MASK)+
	PatternLayer::getSizeFromMask(CMSPatternLayer::PHI_MASK)+
	PatternLayer::getSizeFromMask(CMSPatternLayer::SEG_MASK)+
	PatternLayer::getSizeFromMask(CMSPatternLayer::STRIP_MASK);

      bool overflow = false;
      if(total_nb_bits+patterns[0]->getLayerStrip(0)->getDCBitsNumber()>18){//format+DC bit too large
	if(PatternLayer::getSizeFromMask(CMSPatternLayer::MOD_MASK)==5){//can we take of a module bit?
	  overflow=true;
	  cout<<"** \t1 bit set at 0"<<endl;
	  cout<<"** \t4 bits for the position of the module on the ladder"<<endl;
	}
	else{//not enough room...
	  cout<<"** Your format can not fit in the 18 bits of the AM05 chip!!! ("<<total_nb_bits+patterns[0]->getLayerStrip(0)->getDCBitsNumber()<<" bits needed)"<<endl;
	}
      }
      else{//we have enough room
	if(!overflow){//module is not already treated
	  cout<<"** \t"<<PatternLayer::getSizeFromMask(CMSPatternLayer::MOD_MASK)<<" bits for the position of the module on the ladder"<<endl;
	}
	cout<<"** \t"<<PatternLayer::getSizeFromMask(CMSPatternLayer::PHI_MASK)<<" bits for the position of the ladder on the layer"<<endl;
	cout<<"** \t"<<PatternLayer::getSizeFromMask(CMSPatternLayer::SEG_MASK)<<" bit for the position of the segment on the module";
	if(CMSPatternLayer::INNER_LAYER_SEG_DIVIDE==2 && CMSPatternLayer::OUTER_LAYER_SEG_DIVIDE==2) 
	  cout<<" (forced to 0 on barrel layers)"<<endl;
	else if(CMSPatternLayer::INNER_LAYER_SEG_DIVIDE==2){
	  cout<<" (forced to 0 on the 3 innermost barrel layers)"<<endl;
	}
	else if(CMSPatternLayer::OUTER_LAYER_SEG_DIVIDE==2){
	  cout<<" (forced to 0 on the 3 outermost barrel layers)"<<endl;
	}
	cout<<"** \t"<<PatternLayer::getSizeFromMask(CMSPatternLayer::STRIP_MASK)<<" bits for the position of the superstrip on the segment (stub's strip index divided by the superstrip size (see below) ENCODED IN GRAY CODE)"<<endl;
      }
      cout<<"**"<<endl;
      while(true){
	cout<<"** The 8 input buses are used for the following layers (CMS IDs - superstrips size between parenthesis) : ";
	for(unsigned int j=0;j<layers.size();j++){
	  if(hybrid_sector && first_pass && layers[j]==biggestID)
	    continue;
	  if(hybrid_sector && !first_pass && layers[j]==biggestBarrelID)
	    continue;
	  cout<<layers[j]<<" ("<<st.getSuperStripSize(layers[j])<<") - ";
	}
	for(unsigned int j=layers.size();j<8;j++){
	  cout<<"Unused - ";
	}
	cout<<endl;
	cout<<"**"<<endl;

	if(sector_id!=-1)
	  cout<<"********* PATTERNS FOR SECTOR "<<sector_id<<" **********"<<endl;
	else
	  cout<<"********** PATTERNS **********"<<endl;
	cout<<endl;

	for(unsigned int j=0;j<patterns.size();j++){
	  Pattern* p = patterns[j];

	  int nb_active_layers = p->getNbLayers()-p->getNbFakeSuperstrips();
	  if(expected_active_layers!=-1 && expected_active_layers!=nb_active_layers) // the pattern does not have the expected number of active layers
	    continue;

	  if(hybrid_sector && first_pass && !p->getLayerStrip(8)->isFake())
	    continue;

	  if(hybrid_sector && !first_pass && p->getLayerStrip(8)->isFake())
	    continue;

	  for(int k=0;k<p->getNbLayers();k++){

	    if(hybrid_sector && first_pass && k==8)
	      continue;
	    if(hybrid_sector && !first_pass && k==biggestBarrelIndex)
	      continue;

	    PatternLayer* mp = p->getLayerStrip(k);
	    cout<<((CMSPatternLayer*)mp)->toAM05Format()<<endl;
	  }
	  //unused layers set to 15 (2 DC with "cannot be activated" value)
	  for(int k=p->getNbLayers();k<8;k++){
	    cout<<"15 2"<<endl;
	  }
	  cout<<endl;
	}
	cout<<"********** END OF PATTERNS **********"<<endl;
	if(!hybrid_sector)
	  break;
	else if(!first_pass)
	  break;
	first_pass=false;
      }
    }
  }
  else if(vm.count("alterBank")) {
    SectorTree st;
    SectorTree st2;
    SectorTree* save=NULL;
    int minFS=-1;
    int maxFS=-1;
    int newNbPatterns = -1;
    if(vm.count("minFS"))
      minFS = vm["minFS"].as<int>();
    if(vm.count("maxFS"))
      maxFS = vm["maxFS"].as<int>();
    if(vm.count("truncate"))
      newNbPatterns = vm["truncate"].as<int>();

    if(minFS<0 && maxFS<0 && newNbPatterns<0){
      cout<<"Missing parameter : you need to set minFS, maxFS or truncate!"<<endl;
      return -1;
    }

    if(maxFS<0)
      maxFS=1000;

    cout<<"Loading pattern bank..."<<endl;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
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
	  std::ifstream new_ifs(vm["bankFile"].as<string>().c_str());
	  boost::archive::text_iarchive ia(new_ifs);
	  ia >> st;
	}
      }
    }

    if(newNbPatterns>0){
      cout<<"loaded "<<st.getAllSectors()[0]->getLDPatternNumber()<<" patterns for sector "<<st.getAllSectors()[0]->getOfficialID()<<endl;
      st.getAllSectors()[0]->getPatternTree()->truncate(newNbPatterns);
      save=&st;
    }
    else{
      vector<int> size_on_layers = st.getSuperStripSizeLayers();
      for(unsigned int i=0;i<size_on_layers.size();i++){
	st2.setSuperStripSize(st.getSuperStripSize(size_on_layers[i]), size_on_layers[i]);
      }
      
      cout<<"Altering bank..."<<endl;
      vector<Sector*> sectors = st.getAllSectors();
      for(unsigned int i=0;i<sectors.size();i++){
	Sector* mySector = sectors[i];
	st2.addSector(*mySector);
	Sector* newSector = st2.getAllSectors()[i];
	vector<GradedPattern*> patterns = mySector->getPatternTree()->getLDPatterns();
	for(unsigned int j=0;j<patterns.size();j++){
	  GradedPattern* p = patterns[j];
	  int nbFS = p->getNbFakeSuperstrips();
	  if(nbFS<=maxFS && nbFS>=minFS){
	    //add the pattern
	    newSector->getPatternTree()->addPattern(p,NULL,p->getAveragePt());
	  }
	}
	cout<<"Sector "<<mySector->getOfficialID()<<" :\n\tinput bank : "<<mySector->getPatternTree()->getLDPatternNumber()<<" patterns\n\toutput bank : "<<newSector->getPatternTree()->getLDPatternNumber()<<" patterns."<<endl;
      }
      save=&st2;
    }

    cout<<"Saving new bank in "<<vm["outputFile"].as<string>().c_str()<<"..."<<endl;
    {
      const SectorTree& ref = *save;
      std::ofstream ofs(vm["outputFile"].as<string>().c_str());
      // Compression part
      boost::iostreams::filtering_stream<boost::iostreams::output> f;
      f.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
      f.push(ofs);
      //
      boost::archive::text_oarchive oa(f);
      oa << ref;
    }
  }
  else if(vm.count("MergeSectors")) {
    PatternFinder::mergeFiles(vm["outputFile"].as<string>().c_str(), vm["inputFile"].as<string>().c_str(), vm["secondFile"].as<string>().c_str());
  }
  else if(vm.count("MergeBanks")) {
    SectorTree st1;
    cout<<"Loading pattern bank from "<<vm["inputFile"].as<string>().c_str()<<"..."<<endl;
    {
      std::ifstream ifs(vm["inputFile"].as<string>().c_str());
      //Decompression
      boost::iostreams::filtering_stream<boost::iostreams::input> f;
      f.push(boost::iostreams::gzip_decompressor());
      //we try to read a compressed file
      try { 
	f.push(ifs);
	boost::archive::text_iarchive ia(f);
	ia >> st1;
      }
      catch (boost::iostreams::gzip_error& e) {
	if(e.error()==4){//file is not compressed->read it without decompression
	  std::ifstream new_ifs(vm["inputFile"].as<string>().c_str());
	  boost::archive::text_iarchive ia(new_ifs);
	  ia >> st1;
	}
      }
    }
    vector<Sector*> list1 = st1.getAllSectors();
    unsigned int nbSectors1 = list1.size();
    if(nbSectors1>1){
      cout<<"You can only merge banks containing 1 sector ("<<nbSectors1<<" found)"<<endl;
      return -1;
    }

    int nbPatterns1 = list1[0]->getPatternTree()->getLDPatternNumber();
    cout<<nbPatterns1<<" patterns found."<<endl;
    SectorTree st2;
    cout<<"Loading pattern bank from "<<vm["secondFile"].as<string>().c_str()<<"..."<<endl;
    {
      std::ifstream ifs(vm["secondFile"].as<string>().c_str());
      boost::iostreams::filtering_stream<boost::iostreams::input> f;
      f.push(boost::iostreams::gzip_decompressor());
      //we try to read a compressed file
      try { 
	f.push(ifs);
	boost::archive::text_iarchive ia(f);
	ia >> st2;
      }
      catch (boost::iostreams::gzip_error& e) {
	if(e.error()==4){//file is not compressed->read it without decompression
	  std::ifstream new_ifs(vm["secondFile"].as<string>().c_str());
	  boost::archive::text_iarchive ia(new_ifs);
	  ia >> st2;
	}
      }
    } 
    vector<Sector*> list2 = st2.getAllSectors();
    unsigned int nbSectors2 = list2.size();
    if(nbSectors2>1){
      cout<<"You can only merge banks containing 1 sector ("<<nbSectors2<<" found)"<<endl;
      return -1;
    }
    vector<int> layers_1 = st1.getSuperStripSizeLayers();
    vector<int> layers_2 = st2.getSuperStripSizeLayers();
    for(unsigned int i=0;i<layers_1.size();i++){
      if(st1.getSuperStripSize(layers_1[i])!=st2.getSuperStripSize(layers_1[i])){
	cout<<"You can only merge banks using the same superstrip size (on layer "<<layers_1[i]<<" : "<<st1.getSuperStripSize(layers_1[i])<<" and "<<st2.getSuperStripSize(layers_1[i])<<" found)"<<endl;
	return -1;
      }
    }
    for(unsigned int i=0;i<layers_2.size();i++){
      if(st1.getSuperStripSize(layers_2[i])!=st2.getSuperStripSize(layers_2[i])){
	cout<<"You can only merge banks using the same superstrip size (on layer "<<layers_2[i]<<" : "<<st1.getSuperStripSize(layers_2[i])<<" and "<<st2.getSuperStripSize(layers_2[i])<<" found)"<<endl;
	return -1;
      }
    }
    if(list1[0]->getNbLayers()!=list2[0]->getNbLayers()){
      cout<<"You can only merge banks using the same number of layers ("<<list1[0]->getNbLayers()<<" and "<<list2[0]->getNbLayers()<<" found)"<<endl;
      return -1;
    }
    int nbPatterns2 = list2[0]->getPatternTree()->getLDPatternNumber();
    cout<<nbPatterns2<<" patterns found."<<endl;

    cout<<"Merging banks..."<<endl;
    if(nbPatterns1>nbPatterns2){
      list1[0]->getPatternTree()->addPatternsFromTree(list2[0]->getPatternTree());
      cout<<"-> "<<list1[0]->getPatternTree()->getLDPatternNumber()<<" patterns."<<endl;
      cout<<"Saving new bank in "<<vm["outputFile"].as<string>().c_str()<<"..."<<endl;
      {
	const SectorTree& ref = st1;
	std::ofstream ofs(vm["outputFile"].as<string>().c_str());
	// Compression part
	boost::iostreams::filtering_stream<boost::iostreams::output> f;
	f.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
	f.push(ofs);
	//
	boost::archive::text_oarchive oa(f);
	oa << ref;
      }
    }
    else{
      list2[0]->getPatternTree()->addPatternsFromTree(list1[0]->getPatternTree());
      cout<<"-> "<<list2[0]->getPatternTree()->getLDPatternNumber()<<" patterns."<<endl;
      cout<<"Saving new bank in "<<vm["outputFile"].as<string>().c_str()<<"..."<<endl;
      {
	const SectorTree& ref = st2;
	std::ofstream ofs(vm["outputFile"].as<string>().c_str());
	// Compression part
	boost::iostreams::filtering_stream<boost::iostreams::output> f;
	f.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
	f.push(ofs);
	//
	boost::archive::text_oarchive oa(f);
	oa << ref;
      }
    }
  }
  else if(vm.count("testCode")) {
#ifdef USE_CUDA
    //cuProfilerStart();
    GPUPooler *gp = new GPUPooler("612_SLHC6_MUBANK_lowmidhig_sec26_ss32_cov40.pbk", "/dev/shm/RawData","/dev/shm/Pattern",5);
    gp->loopForEvents(100,60000);
    delete gp;
    //cuProfilerStop();
#else
    string result;
    {
      SectorTree sTest;
      cout<<"Loading pattern bank..."<<endl;
      {
	string file = "/home/infor/baulieu/private/cms/slc6/CMSSW_6_2_0_SLHC15/src/amsimulation_dev/test3_Fountain_NoSeg_2DC_lowmidhig.pbk";
	std::ifstream ifs(file.c_str());
	boost::iostreams::filtering_stream<boost::iostreams::input> f;
	f.push(boost::iostreams::gzip_decompressor());
	//we try to read a compressed file
	try { 
	  f.push(ifs);
	  boost::archive::text_iarchive ia(f);
	  ia >> sTest;
	}
	catch (boost::iostreams::gzip_error& e) {
	  if(e.error()==4){//file is not compressed->read it without decompression
	    std::ifstream new_ifs(file.c_str());
	    boost::archive::text_iarchive ia(new_ifs);
	    ia >> sTest;
	  }
	}
      }
      cout<<"Sector :"<<endl;
      cout<<*(sTest.getAllSectors()[0])<<endl;
      cout<<"loaded "<<sTest.getAllSectors()[0]->getLDPatternNumber()<<" patterns for sector "<<sTest.getAllSectors()[0]->getOfficialID()<<endl;
      sTest.getAllSectors()[0]->getPatternTree()->truncate(1000000);

      cout<<"Saving new bank in testOut.pbk..."<<endl;
      {
	const SectorTree& ref = sTest;
	std::ofstream ofs("./testOut.pbk");
	// Compression part
	boost::iostreams::filtering_stream<boost::iostreams::output> f;
	f.push(boost::iostreams::gzip_compressor(boost::iostreams::gzip_params(boost::iostreams::gzip::best_compression)));
	f.push(ofs);
	//
	boost::archive::text_oarchive oa(f);
	oa << ref;
      }
    }
#endif
  }
}
