#include <fstream>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
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


#include <TH1I.h>
#include <TFile.h>


using namespace std;

/**
   \mainpage
   \section build_sec Building the project
   In order to build the project, you will need:
   - Root (http://root.cern.ch) installed and configured ($ROOTSYS must be pointing on the installation directory and $ROOTSYS/bin must be in the PATH)
   - Boost (http://www.boost.org/) libraries and header files
   
   Then from the main directory (above ./src) you only need to type :
   \code
   make
   \endcode

   If everything goes fine, you should get a binary file called "AMSimulation".
 
   \section use_sec Using the program
   \subsection generate Generating a pattern bank

   A text user interface is available but not always up to date (does not support the Z for the sectors yet). The easiest way to generate banks is to write a small code creating the sector and running the generation. For an example look at the AMSimulation.cc file in the testCode section. This section is run with the command :
   \code
   ./AMSimulation --testCode
   \endcode

   If you want to use the text user interface (not up to date) to generate a pattern bank file from simulation files (muons + anti muons), enter:
   \code
   ./AMSimulation --generateBank
   \endcode
   You will be asked the following informations :
   \code
   Enter the layer numbers (separated by spaces) :
   \endcode
   You need to enter the layer's ID numbers, for example : "8 9 10" for the 3 outermost layers
   \code
   Enter the sector num 1 :
   \endcode
   You have to give the sectors for which you want to generate the patterns. You need at least one sector. Enter the number of the ladders on the first layer, separated with spaces, then press the enter key. Enter the ladders on the second layer and so on... Once the first sector is entered, you will be asked for the second one. Proceed as previously until you have entered all your sectors. When you are done, enter '-1' as the ladder number of the first layer of the new sector, this will stop the sectors edition.
   \code
   Enter the super strip size :
   \endcode
   This is the number of strips contained in a superstrip. The value can be 8, 16, 32, 64, 128, 256, 512, 1024.
   \code
   Enter the number of DC bits :
   \endcode
   The number of DC bits used for adaptative patterns. This ranges from 0 (no adaptative patterns) to 3.
   \code
   Enter the minimum PT :
   \endcode
   Tracks having a PT bellow this threshold will not be used
   \code
   Enter the maximum PT :
   \endcode
   Tracks having a PT above this threshold will not be used
   \code
   Enter the muon files directory name :
   \endcode
   The directory containing the root files with PG muons and anti muons (local or RFIO).
   \code
   Enter the patterns bank file name :
   \endcode
   The name of the file that will contain the patterns
   \code
   Enter the root output file name :
   \endcode
   The name of the Root file that will contain the generation histograms
   \code
   Enter the bank growth threshold :
   \endcode
   The threshold value to stop the generation (the value is in the range [0-1], for example 0.95). This is the percentage of recognized tracks with the current patterns bank.

   \subsection find Finding patterns in events
   To search for patterns in events, enter :
   \code
   ./AMSimulation --findPatterns --inputFile <path to Root File containing events (local or RFIO)> --bankFile <path to your pattern bank file> --outputFile <Root output file> --ss_threshold <minimum number of stubs to activate the pattern> --startEvent <Index of first event to analyse> --stopEvent <Index of last event to analyse>
   \endcode

  The program is using a virtual detector during the search process. The geometry of this detector (layers, ladders, Z modules and strips per segment) is contained in the detector.cfg file located in the root directory of the program. Here is an example of the syntax :
  \code
  #layerID,nb ladders,nb Z modules,nb strips per segment
  5,16,15,1024
  6,24,15,1024
  7,36,15,1024
  8,48,14,1024
  9,60,14,1024
  10,76,14,1024
  \endcode

   \subsection merge Merging result files
   To merge root files containing results from patterns recognition, you can use the hadd binary distributed with Root. This will work to merge results concerning the same sectors but different events.
   
   If you need to merge files containing the same events but different sectors you can use (NOT UP TO DATE!!):
   \code
   ./AMSimulation --MergeSectors --inputFile  <Root file for events A to X in sector 1> --secondFile <Root file for events A to X in sector 2> --outputFile <Resulting Root file for events A to X in sector 1 & 2>
   \endcode

   \subsection view Viewing the content of a pattern bank
   You can display the patterns contained in a patterns bank file using the command :
   \code
   ./AMSimulation --printBank --bankFile <You patterns bank file>
   \endcode

   It should display one pattern per line.


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

int getSuperStripSize(){
  cout<<"Enter the super strip size :"<<endl;
  int result;
  cin>>result;
  return result;
}

int getDCBitsNumber(){
  cout<<"Enter the number of DC bits :"<<endl;
  int result;
  cin>>result;
  if(result<0)
    result=0;
  if(result>3)
    result=3;
  return result;
}

float getThreshold(){
  cout<<"Enter the bank growth threshold :"<<endl;
  float result;
  cin>>result;
  return result;
}

float getMinPT(){
  cout<<"Enter the minimum PT :"<<endl;
  float result;
  cin>>result;
  if(result<0)
    result=0;
  return result;
}

float getMaxPT(){
  cout<<"Enter the maximum PT :"<<endl;
  float result;
  cin>>result;
  if(result<0)
    result=0;
  return result;
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

string getParticuleDirName(){
  cout<<"Enter the muon files directory name :"<<endl;
  string result;
  cin>>result;
  return result;
}

string getPatternBankFileName(){
  cout<<"Enter the patterns bank file name :"<<endl;
  string result;
  cin>>result;
  return result;
}

string getSectorDefFileName(){
  cout<<"Enter the sector definition file name :"<<endl;
  string result;
  cin>>result;
  return result;
}

string getRootOutputFileName(){
  cout<<"Enter the root output file name :"<<endl;
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

int main(int av, char** ac){
  namespace po = boost::program_options;
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("analyseBank", "Creates histograms from a pattern bank file")
    ("inputFile", po::value<string>(), "The file to analyse")
    ("secondFile", po::value<string>(), "Second file to merge")
    ("bankFile", po::value<string>(), "The patterns bank file to use")
    ("outputFile", po::value<string>(), "The root output file")
    ("ss_threshold", po::value<int>(), "The minimum number of hit superstrips to activate a pattern")
    ("startEvent", po::value<int>(), "The first event index")
    ("stopEvent", po::value<int>(), "The last event index")
    ("decode", po::value<int>(), "Decode the given super strip")
    ("generateBank", "Generates a pattern bank from root simulation file")
    ("testSectors", "Get the tracks sectors")
    ("MergeSectors", "Merge 2 root files having same events but different sectors (needs --inputFile --secondFile and --outputFile)")
    ("MergeBanks", "Merge 2 bank files having only 1 sector (needs --inputFile --secondFile and --outputFile)")
    ("buildFitParams", "Computes the Fit parameters for the given bank using tracks from the given directory (needs --bankFile, --inputFile and --outputFile)")
    ("findPatterns", "Search for patterns in an event file (needs --ss_threshold --inputFile, --bankFile, --outputFile, --startEvent and --stopEvent)")
    ("printBank", "Display all patterns from a bank (needs --bankFile)")
    ("testCode", "Dev tests")
    ;
     
  po::variables_map vm;
  po::store(po::parse_command_line(av, ac, desc), vm);
  po::notify(vm);    
  
  if (vm.count("help")) {
    cout << desc << "\n";

    return 1;
  }

  if (vm.count("analyseBank")) {
    SectorTree st;
    {
      std::ifstream ifs(vm["inputFile"].as<string>().c_str());
      boost::archive::text_iarchive ia(ifs);
      ia >> st;
    }
    TFile f( vm["outputFile"].as<string>().c_str(), "recreate");
    createAnalysis(st);
  } else if (vm.count("generateBank")) {
    vector<int> layers;
    SectorTree st;
    int stripSize;
    int dcBits;
    string partDirName;
    string bankFileName;
    string rootFileName;
    float threshold;
    float min;
    float max;
    float minEta;
    float maxEta;
    map<int,pair<float,float> > eta;
    getLayers(layers);
    getSectors(layers,st);
    stripSize=getSuperStripSize();
    dcBits=getDCBitsNumber();
    min=getMinPT();
    max=getMaxPT();
    minEta=getEtaMin();
    maxEta=getEtaMax();
    partDirName=getParticuleDirName();
    bankFileName=getPatternBankFileName();
    rootFileName=getRootOutputFileName();
    threshold=getThreshold();
    PatternGenerator pg(stripSize);//Super strip size
    pg.setLayers(layers);
    pg.setParticuleDirName(partDirName);
    pg.setMinPT(min);
    pg.setMaxPT(max);
    pg.setMinEta(minEta);
    pg.setMaxEta(maxEta);
    TFile f(rootFileName.c_str(), "recreate");
    pg.setVariableResolution(dcBits);
    pg.generate(&st, 40000, threshold, eta);


    if(pg.getVariableResolutionState()>0){
      cout<<"HD Patterns : "<<st.getFDPatternNumber()<<endl;
      cout<<"LD Patterns : "<<st.getLDPatternNumber()<<endl;
    }
  
    cout<<"Saving SectorTree...";
    {
      const SectorTree& ref = st;
      std::ofstream ofs(bankFileName.c_str());
      boost::archive::text_oarchive oa(ofs);
      oa << ref;
      cout<<"done."<<endl;
    }

    vector<Sector*> list = st.getAllSectors();
    int nb_patterns = 0;
    cout<<"************************************************"<<endl;
    cout<<"We have "<<list.size()<<" sectors :"<<endl;
    for(unsigned int i=0;i<list.size();i++){
      cout<<*list[i];
      int nb = list[i]->getLDPatternNumber();
      cout<<nb<<" generated patterns"<<endl;
      cout<<endl;
      nb_patterns+=nb;
    }

    cout<<"Total number of patterns : "<<nb_patterns<<endl;
    cout<<"************************************************"<<endl;
    
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
  else if(vm.count("findPatterns")) {
    SectorTree st;
    cout<<"Loading pattern bank..."<<endl;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
      boost::archive::text_iarchive ia(ifs);
      ia >> st;
    }

    ///////////////////////////////////////////////////////////////
    // If we don't have a fitter -> create a Karimaki default one
    vector<Sector*> sectors = st.getAllSectors();
    for(unsigned int i=0;i<sectors.size();i++){
      if(sectors[i]->getFitter()==NULL){
	TrackFitter* fitter = new KarimakiTrackFitter(sectors[i]->getNbLayers());
	sectors[i]->setFitter(fitter);
	sectors[i]->updateFitterPhiRotation();
      }
    }
    ///////////////////////////////////////////////////////////////

    PatternFinder pf(st.getSuperStripSize(), vm["ss_threshold"].as<int>(), &st,  vm["inputFile"].as<string>().c_str(),  vm["outputFile"].as<string>().c_str());
    {
      boost::progress_timer t;
      int start = vm["startEvent"].as<int>();
      int stop = vm["stopEvent"].as<int>();
      vector<Sector*> pattern_list2 = pf.find(start, stop);
      cout<<"Time used to analyse "<<stop-start+1<<" events : "<<endl;
    }
  }
  else if(vm.count("buildFitParams")) {
    SectorTree st;
    cout<<"Loading pattern bank..."<<endl;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
      boost::archive::text_iarchive ia(ifs);
      ia >> st;
    }

    map<int,pair<float,float> > eta_limits;// eta values for which each layer does exist
    eta_limits[6]=pair<float,float>(-1.69,1.69);
    eta_limits[7]=pair<float,float>(-1.41,1.41);
    eta_limits[8]=pair<float,float>(-1.19,1.19);
    eta_limits[9]=pair<float,float>(-1.02,1.02);
    eta_limits[10]=pair<float,float>(-0.87,0.87);
    eta_limits[11]=pair<float,float>(1.12,2.19);
    eta_limits[12]=pair<float,float>(1.19,2.19);
    eta_limits[13]=pair<float,float>(1.28,2.19);
    eta_limits[14]=pair<float,float>(1.35,2.19);
    eta_limits[15]=pair<float,float>(1.43,2.19);

    PrincipalFitGenerator pfg(vm["inputFile"].as<string>().c_str(), &st);
    pfg.generate(eta_limits, 2, 100, 0, 0.87);
    
    cout<<"Saving SectorTree...";
    {
      const SectorTree& ref = st;
      std::ofstream ofs(vm["outputFile"].as<string>().c_str());
      boost::archive::text_oarchive oa(ofs);
      oa << ref;
      cout<<"done."<<endl;
    }
    
  }
  else if(vm.count("printBank")) {
    SectorTree st;
    cout<<"Loading pattern bank..."<<endl;
    {
      std::ifstream ifs(vm["bankFile"].as<string>().c_str());
      boost::archive::text_iarchive ia(ifs);
      ia >> st;
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
    /*
    vector<Sector*> sectors = st.getAllSectors();
    for(unsigned int i=0;i<sectors.size();i++){
      set<string> combinaisons;
      Sector* mySector = sectors[i];
      vector<GradedPattern*> patterns = mySector->getPatternTree()->getLDPatterns();
      cout<<patterns.size()<<" patterns in the bank"<<endl;
      for(unsigned int j=0;j<patterns.size();j++){
	ostringstream oss;
	Pattern* p = patterns[j];
	for(int k=0;k<p->getNbLayers()-3;k++){
	  PatternLayer* mp = p->getLayerStrip(k);
	  oss<<((CMSPatternLayer*)mp)->getPhi()<<"/"<<((CMSPatternLayer*)mp)->getModule()<<" - ";
	}
	combinaisons.insert(oss.str());
      }
      cout<<"Nb combinaisons : "<<combinaisons.size()<<endl;
      std::set<string>::iterator it;
      for(it=combinaisons.begin();it!=combinaisons.end();it++){
	cout<< *it <<endl;
      }
    }
    */
  }
  else if(vm.count("MergeSectors")) {
    PatternFinder::mergeFiles(vm["outputFile"].as<string>().c_str(), vm["inputFile"].as<string>().c_str(), vm["secondFile"].as<string>().c_str());
  }
  else if(vm.count("MergeBanks")) {
    SectorTree st1;
    cout<<"Loading pattern bank from "<<vm["inputFile"].as<string>().c_str()<<"..."<<endl;
    {
      std::ifstream ifs(vm["inputFile"].as<string>().c_str());
      boost::archive::text_iarchive ia(ifs);
      ia >> st1;
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
      boost::archive::text_iarchive ia(ifs);
      ia >> st2;
    } 
    vector<Sector*> list2 = st2.getAllSectors();
    unsigned int nbSectors2 = list2.size();
    if(nbSectors2>1){
      cout<<"You can only merge banks containing 1 sector ("<<nbSectors2<<" found)"<<endl;
      return -1;
    }
    if(st1.getSuperStripSize()!=st2.getSuperStripSize()){
      cout<<"You can only merge banks using the same superstrip size ("<<st1.getSuperStripSize()<<" and "<<st2.getSuperStripSize()<<" found)"<<endl;
      return -1;
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
	boost::archive::text_oarchive oa(ofs);
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
	boost::archive::text_oarchive oa(ofs);
	oa << ref;
      }
    }
  }
  else if(vm.count("testCode")) {
    //cout<<"Nothing to be done"<<endl;
    /*
    cout<<"Taille CMSPatternLayer : "<<sizeof(CMSPatternLayer)*8<<" bits"<<endl;
    cout<<"Taille bitset<15> : "<<sizeof(bitset<15>)*8<<" bits"<<endl;
    cout<<"Taille char : "<<sizeof(char)*8<<" bits"<<endl;
    cout<<"Taille  vector< PatternLayer* > : "<<sizeof( vector< PatternLayer* >)*8<<" bits"<<endl;
    cout<<"Taille char[3] : "<<sizeof(char[3])*8<<" bits"<<endl;
    cout<<"Taille GradedPattern : "<<sizeof(GradedPattern)*8<<" bits"<<endl;
    cout<<"Taille PatternTrunk : "<<sizeof(PatternTrunk)*8<<" bits"<<endl;
    cout<<"Taille PatternTree : "<<sizeof(PatternTree)*8<<" bits"<<endl;
    cout<<"Taille SuperStrip : "<<sizeof(SuperStrip)*8<<" bits"<<endl;
    cout<<"Taille Pattern : "<<sizeof(Pattern)*8<<" bits"<<endl;
    cout<<"Taille Hit : "<<sizeof(Hit)*8<<" bits"<<endl;
    cout<<"Taille short : "<<sizeof(short)*8<<" bits"<<endl;
    cout<<"Taille int : "<<sizeof(int)*8<<" bits"<<endl;
    cout<<"Taille float : "<<sizeof(float)*8<<" bits"<<endl;
    cout<<"Taille string : "<<sizeof(string)*8<<" bits"<<endl;
    cout<<"Taille string \"1234567891234567\" : "<<sizeof("1234567891234567")*8<<" bits"<<endl;

    CMSPatternLayer tt;
    tt.setValues(14,7,42,1);
    cout<<tt.toString()<<endl;
    */

    string result;
    {
      /*
      SectorTree sTest;
      cout<<"Loading pattern bank..."<<endl;
      {
	std::ifstream ifs("/home/infor/baulieu/private/cms/AMSimulation/bank_GIO_5L_32ss_2DC_PT2-100_10.pbk");
	boost::archive::text_iarchive ia(ifs);
	ia >> sTest;
      }
      cout<<"Sector :"<<endl;
      cout<<*(sTest.getAllSectors()[0])<<endl;
      cout<<"loaded "<<sTest.getAllSectors()[0]->getLDPatternNumber()<<" patterns"<<endl;
      
      //SectorTree newST;
      //newST.setSuperStripSize(16);
      //Sector newSector(*(st.getAllSectors()[0]));
      //Sector newSector;
      //newST.addSector(newSector);
      //PatternTree* newPT = newST.getAllSectors()[0]->getPatternTree();
      //vector<GradedPattern*> list = st.getAllSectors()[0]->getPatternTree()->getLDPatterns();
      
      //struct sysinfo info;
      //sysinfo(&info);
      //cout<<"RAM free : "<<info.freeram<<endl;
      */

      /*
      cout<<"saving pattern bank..."<<endl;
      {
	const SectorTree& ref = newST;
	cout<<"trying to save "<<newST.getAllSectors()[0]->getLDPatternNumber()<<" patterns"<<endl;
	std::ofstream ofs("/home/infor/baulieu/private/cms/AMSimulation/bank16_3L_78_DCXXX_70_.pbk");
	boost::archive::text_oarchive oa(ofs);
	oa << ref;
      } 
      */

      vector<int> layers;
      SectorTree st;
      int stripSize;
      int dcBits;
      string partDirName;
      string bankFileName;
      string rootFileName;
      float threshold;
      float min;
      float max;
      float minEta;
      float maxEta;
      
      /*      
      //ENDCAP

      layers.push_back(11);
      layers.push_back(12);
      layers.push_back(13);
      layers.push_back(14);
      layers.push_back(15);

      Sector s(layers);
      s.addLadder(11,1,10);
      s.addModules(11,1,1,3);
      s.addModules(11,2,2,4);
      s.addModules(11,3,1,4);
      s.addModules(11,4,1,5);
      s.addModules(11,5,1,5);
      s.addModules(11,6,1,5);
      s.addModules(11,7,1,6);
      s.addModules(11,8,2,6);
      s.addModules(11,9,2,8);
      s.addModules(11,10,3,9);

      s.addLadder(12,1,11);
      s.addModules(12,1,1,3);
      s.addModules(12,2,1,4);
      s.addModules(12,3,1,4);
      s.addModules(12,4,1,5);
      s.addModules(12,5,1,5);
      s.addModules(12,6,1,5);
      s.addModules(12,7,2,6);
      s.addModules(12,8,2,7);
      s.addModules(12,9,2,8);
      s.addModules(12,10,2,10);
      s.addModules(12,11,3,11);

      s.addLadder(13,2,11);
      s.addModules(13,2,1,4);
      s.addModules(13,3,1,4);
      s.addModules(13,4,1,5);
      s.addModules(13,5,1,5);
      s.addModules(13,6,1,6);
      s.addModules(13,7,2,6);
      s.addModules(13,8,2,7);
      s.addModules(13,9,2,8);
      s.addModules(13,10,2,10);
      s.addModules(13,11,2,11);
      s.addModules(13,12,3,12);

      s.addLadder(14,3,11);
      s.addModules(14,3,1,4);
      s.addModules(14,4,1,5);
      s.addModules(14,5,1,5);
      s.addModules(14,6,1,6);
      s.addModules(14,7,1,6);
      s.addModules(14,8,2,7);
      s.addModules(14,9,2,8);
      s.addModules(14,10,2,10);
      s.addModules(14,11,2,11);
      s.addModules(14,12,2,12);
      s.addModules(14,13,2,14);

      s.addLadder(15,11, 11);
      s.addModules(15,3,1,4);
      s.addModules(15,4,1,5);
      s.addModules(15,5,1,5);
      s.addModules(15,6,2,6);
      s.addModules(15,7,2,7);
      s.addModules(15,8,2,8);
      s.addModules(15,9,2,10);
      s.addModules(15,10,2,11);
      s.addModules(15,11,2,12);
      s.addModules(15,12,2,14);
      s.addModules(15,13,2,15);

      */
      
      //BARREL
      layers.push_back(5);
      layers.push_back(6);
      layers.push_back(7);
      layers.push_back(8);
      layers.push_back(9);
      layers.push_back(10);

      Sector s(layers);
      
      s.addLadders(5,2,1);
      s.addModules(5,2,13,8);
      
      s.addLadders(6,2,3);
      s.addModules(6,2,12,8);
      s.addModules(6,3,12,8);
      s.addModules(6,4,12,8);

      s.addLadders(7,3,3);
      s.addModules(7,3,11,10);
      s.addModules(7,4,11,10);
      s.addModules(7,5,11,10);
      
      s.addLadders(8,4,5);
      s.addModules(8,4,10,11);
      s.addModules(8,5,10,11);
      s.addModules(8,6,10,11);
      s.addModules(8,7,10,11);
      s.addModules(8,8,10,11);
      
      s.addLadders(9,4,8);
      s.addModules(9,4,10,13);
      s.addModules(9,5,10,13);
      s.addModules(9,6,10,13);
      s.addModules(9,7,10,13);
      s.addModules(9,8,10,13);
      s.addModules(9,9,10,13);
      s.addModules(9,10,10,13);
      s.addModules(9,11,10,13);

      s.addLadders(10,4,11);
      s.addModules(10,4,10,14);
      s.addModules(10,5,10,14);
      s.addModules(10,6,10,14);
      s.addModules(10,7,10,14);
      s.addModules(10,8,10,14);
      s.addModules(10,9,10,14);
      s.addModules(10,10,10,14);
      s.addModules(10,11,10,14);
      s.addModules(10,12,10,14);
      s.addModules(10,13,10,14);
      s.addModules(10,14,10,14);

      cout<<s<<endl;
      
      /*
      //BARREL+ENDCAP
      layers.push_back(6);
      layers.push_back(7);
      layers.push_back(8);
      layers.push_back(9);
      layers.push_back(11);
      layers.push_back(12);

      Sector s(layers);

      s.addLadder(6,4);
      s.addModules(6,4,4,13);

      s.addLadder(7,4);
      s.addModules(7,4,5,16);
      s.addLadder(7,5);
      s.addModules(7,5,5,16);
      s.addLadder(7,6);
      s.addModules(7,6,5,16);
      
      s.addLadder(8,5);
      s.addModules(8,5,6,12);
      s.addLadder(8,6);
      s.addModules(8,6,6,12);
      s.addLadder(8,7);
      s.addModules(8,7,5,12);
      s.addLadder(8,8);
      s.addModules(8,8,6,12);
      s.addLadder(8,9);
      s.addModules(8,9,6,12);
      
      s.addLadder(9,5);
      s.addModules(9,5,8,12);
      s.addLadder(9,6);
      s.addModules(9,6,8,12);
      s.addLadder(9,7);
      s.addModules(9,7,7,12);
      s.addLadder(9,8);
      s.addModules(9,8,7,12);
      s.addLadder(9,9);
      s.addModules(9,9,7,12);
      s.addLadder(9,10);
      s.addModules(9,10,8,12);
      s.addLadder(9,11);
      s.addModules(9,11,8,12);

      s.addLadder(11,9);
      s.addModules(11,9,2,8);
      s.addLadder(11,10);
      s.addModules(11,10,3,9);
      s.addLadder(11,11);
      s.addModules(11,11,2,12);
      s.addModule(11,11,60);
      s.addLadder(11,12);
      s.addModules(11,12,3,11);
      s.addLadder(11,13);
      s.addModules(11,13,2,14);
      s.addLadder(11,14);
      s.addModules(11,14,3,13);

      s.addLadder(12,10);
      s.addModules(12,10,2,10);
      s.addLadder(12,11);
      s.addModules(12,11,2,12);
      s.addLadder(12,12);
      s.addModules(12,12,2,12);
      s.addLadder(12,13);
      s.addModules(12,13,2,14);
      s.addModule(12,13,69);
      s.addLadder(12,14);
      s.addModules(12,14,2,16);
      */
      //TrackFitter* fitter = new PrincipalTrackFitter(s.getNbLayers(),1000);
      //s.setFitter(fitter);

      st.addSector(s);

      map<int,pair<float,float> > eta_limits;// eta values for which each layer does exist
      eta_limits[5]=pair<float,float>(-1.69,1.69);
      eta_limits[6]=pair<float,float>(-1.69,1.69);
      eta_limits[7]=pair<float,float>(-1.41,1.41);
      eta_limits[8]=pair<float,float>(-1.19,1.19);
      eta_limits[9]=pair<float,float>(-1.02,1.02);
      eta_limits[10]=pair<float,float>(-0.87,0.87);
      eta_limits[11]=pair<float,float>(1.12,2.1);
      eta_limits[12]=pair<float,float>(1.19,2.25);
      eta_limits[13]=pair<float,float>(1.28,2.4);
      eta_limits[14]=pair<float,float>(1.35,2.5);
      eta_limits[15]=pair<float,float>(1.43,2.5);
            
      stripSize=32;
      dcBits=3;
      min=50;
      max=100;
      minEta=0;
      maxEta=0.87;
      //partDirName="/gridgroup/cms/viret/SLHC/MuBank/";
      //partDirName="/scratch/baulieu/bank/";
      //partDirName="/Users/baulieu/Documents/dev/CMS/AMSimulation/data/1301/bank/";
      //partDirName="rfio:/dpm/in2p3.fr/home/cms/data/store/user/gbaulieu/SLHC/bankgen_1212/";
      //partDirName="/home/infor/baulieu/private/cms/amsimulation/testFiles/";
      partDirName="rfio:/dpm/in2p3.fr/home/cms/data/store/user/sviret/SLHC/GEN/MUBANK_61/";
      bankFileName="bank_6L_BARREL_091_138_32ss_3DC_PT50-100_90.pbk";
      rootFileName="analysis_6L_BARREL_091_138_32ss_3DC_PT50-100_90.root";
      threshold=0.9;

      PatternGenerator pg(stripSize);//Super strip size
      pg.setLayers(layers);
      pg.setMinPT(min);
      pg.setMaxPT(max);
      pg.setMinEta(minEta);
      pg.setMaxEta(maxEta);
      pg.setParticuleDirName(partDirName);
      TFile f(rootFileName.c_str(), "recreate");
      pg.setVariableResolution(dcBits);
      pg.generate(&st, 40000, threshold, eta_limits);
      
      
      if(pg.getVariableResolutionState()>0){
	cout<<"HD Patterns : "<<st.getFDPatternNumber()<<endl;
	cout<<"LD Patterns : "<<st.getLDPatternNumber()<<endl;
      }
      
      cout<<"Saving SectorTree...";
      {
	const SectorTree& ref = st;
	std::ofstream ofs(bankFileName.c_str());
	boost::archive::text_oarchive oa(ofs);
	oa << ref;
	cout<<"done."<<endl;
      }
      
    }
  }
}
