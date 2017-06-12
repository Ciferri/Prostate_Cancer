/**
 * @file Simulator.hpp
 * @brief
 * @author Alfredo Hernandez
 * @author Carlos Sosa Marrero
 * @date 05.19.17 
 */

#ifndef DEF_SIMULATOR
#define DEF_SIMULATOR

#include "Model.hpp"
#include "SimMeth.hpp"

#include <iostream>
#include <fstream>

typedef std::vector<double> OutputData;
typedef std::vector<OutputData *> OutputDataList;

class Simulator{	
public:
  Simulator();
  Simulator(Model *model, const double DT);
  ~Simulator();
  void initSim();
  void setModel(Model *model);
  void simulate(const double currentTime, const double simTime);
  void stop();
  
private:
  OutputDataList *m_outList;
  SimMeth *m_simMeth;
  std::ofstream m_outFile;
  double m_DT;	
  double m_currentTime;
  Model *m_model;
};

#endif
