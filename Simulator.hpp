/**
 * @file Simulator.hpp
 * @brief
 * @author Alfredo Hernandez
 * @author Carlos Sosa Marrero
 * @date 05.19.17 
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
  Simulator(Model *const model, const double currentTime,
	    const double DT);
  ~Simulator();
  void setModel(Model *const model);
  //void simulate( double currentTime, double DT );
  void start(double const simTime);
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
