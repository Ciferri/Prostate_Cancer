/**
 * @file rootSimulator.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#ifndef DEF_ROOTSIMULATOR
#define DEF_ROOTSIMULATOR

#include "coupler.hpp"
#include "Model.hpp"
#include "Simulator.hpp"

class RootSimulator{
public:
  RootSimulator();
  RootSimulator(Model *coupler, const double DT1, const double DT2);
  ~RootSimulator();
  void initSim();
  void simulate(const double currentTime, const double simTime);
  void stop();
private:
  double m_DT1, m_DT2;
  double m_currentTime;
  Model *m_coupler;
  Simulator *m_sim1, *m_sim2;
};
#endif
