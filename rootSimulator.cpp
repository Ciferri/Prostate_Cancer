/**
 * @file rootSimulator.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#include "rootSimulator.hpp"

using namespace std;

RootSimulator::RootSimulator(){
}


RootSimulator::RootSimulator(Model *coupler, const double DT1,
			     const double DT2){
  m_coupler = coupler;
  m_DT1 = DT1;
  m_DT2 = DT2;
  m_currentTime = 0.0;
  m_sim1 = new Simulator(((Coupler *)m_coupler)->getModel1(), DT1,
			 "out.dat");
  m_sim2 = new Simulator(((Coupler *)m_coupler)->getModel2(), DT2,
			 "out2.dat");
}


void RootSimulator::initSim(){
  m_sim1->initSim();
  m_sim2->initSim();
}


void RootSimulator::simulate(const double currentTime,
			     const double simTime){
  int numIter;
  
  numIter = simTime / m_DT1;
  for(int j(0); j < numIter; j++){
    m_sim2->simulate(m_currentTime, 500);
    m_coupler->updateModel();
    m_sim1->simulate(m_currentTime, m_DT1);
    m_coupler->updateModel();
    m_currentTime += m_DT1;
  }
}
  
void RootSimulator::stop(){
  m_sim1->stop();
  m_sim2->stop();
}
