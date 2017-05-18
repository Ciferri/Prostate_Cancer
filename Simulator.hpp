/**
 * @(#) Simulator.hpp
 */

#ifndef DEF_Simulator
#define DEF_Simulator

#include "Model.hpp"
#include "Treatment.hpp"
#include "SimMeth.hpp"
#include <iostream>
#include <fstream>

typedef std::vector<double> OutputData;
typedef std::vector<OutputData *> OutputDataList;

class Simulator{	
public:
  Simulator();
  Simulator(Model *model, double currentTime, double DT);
  ~Simulator();
  void setModel(Model *model);
  //void simulate( double currentTime, double DT );
  void start(double duration);
  void stop();
  
private:
  OutputDataList *m_outList;
  SimMeth *m_simMeth;
  std::ofstream m_outFile;
  double m_DT;	
  double m_currentTime;
  Model *m_model;
  Treatment *m_treatment;
};

#endif
