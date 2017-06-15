/**
 * @file Simulator.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17 
 */

#ifndef DEF_SIMULATOR
#define DEF_SIMULATOR

#include "Model.hpp"
#include "SimMeth.hpp"
#include <string>
#include <iostream>
#include <fstream>

typedef std::vector<double> OutputData;
typedef std::vector<OutputData *> OutputDataList;

class Simulator{	
public:
  Simulator();
  Simulator(Model *model, const double DT, const std::string nFOut);
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
