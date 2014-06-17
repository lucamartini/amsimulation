#ifndef _PRINCIPALFITGENERATOR_H_
#define _PRINCIPALFITGENERATOR_H_

#include <string>
#include <TChain.h>
#include <TSystemDirectory.h>
#include "SectorTree.h"
#include "PrincipalTrackFitter.h"


/**
   \brief Class used to generate the parameters of a Principal Components Analysis Fitter
**/
class PrincipalFitGenerator{

 private:
  string inputDirectory;
  SectorTree* st;

  // Containers to load the TTree branches
  int m_stub;
  vector<int>                   m_stub_modid;  // Layer/ladder/module du stub
  vector<float>                 m_stub_ptGEN;  // Pt de la particule initiale (en GeV/c)
  vector<float>                 m_stub_phiGEN;  // Phi0 de la particule initiale (en GeV/c)
  vector<float>                 m_stub_etaGEN;  // Eta de la particule initiale
  vector<float>                 m_stub_x;      // x coordinate of the hit
  vector<float>                 m_stub_y;      // y coordinate of the hit
  vector<float>                 m_stub_z;      // z coordinate of the hit
  vector<int>                   m_stub_pdg;
  vector<float>                 m_stub_x0;
  vector<float>                 m_stub_y0;
  vector<float>                 m_stub_z0;

  

  vector<int>                   *p_m_stub_modid;
  vector<float>                 *p_m_stub_ptGEN;
  vector<float>                 *p_m_stub_phiGEN;
  vector<float>                 *p_m_stub_etaGEN;
  vector<float>                 *p_m_stub_x;
  vector<float>                 *p_m_stub_y;
  vector<float>                 *p_m_stub_z;
  vector<int>                   *p_m_stub_pdg;
  vector<float>                 *p_m_stub_x0;
  vector<float>                 *p_m_stub_y0;
  vector<float>                 *p_m_stub_z0;

  TChain* createTChain();
  void generatePrincipal(map<int,pair<float,float> > eta_limits, float min_pt, float max_pt, float min_eta, float max_eta);
  void generateMultiDim(map<int,pair<float,float> > eta_limits, float min_pt, float max_pt, float min_eta, float max_eta);
  bool selectTrack(vector<int> &tracker_layers, int* layers, vector<int> &ladder_per_layer, vector<int> &module_per_layer,
		   map<int,pair<float,float> > &eta_limits, float min_pt, float max_pt, float min_eta, float max_eta);

 public:
  /**
     \brief Constructor of the class
     \param f The name of the input directory : directory containing root files with single muons events, used to compute the parameters.
     \param s The SectorTree containing the sector structure and the patterns.
  **/
  PrincipalFitGenerator(string f, SectorTree *s);
  /**
     \brief Method used to compute the parameters used for the PCA
     \param eta_limits Range in eta of each layer
     \param min_pt Minimum PT of the tracks used for the computing
     \param max_pt Maximum PT of the tracks used for the computing
     \param min_eta Minimum ETA of the tracks used for the computing
     \param max_eta Maximum ETA of the tracks used for the computing
   **/
  void generate(map<int,pair<float,float> > eta_limits, float min_pt, float max_pt, float min_eta, float max_eta);
};

#endif
