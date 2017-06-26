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


RootSimulator::RootSimulator(Coupler *coupler, const double DT1,
			     const double DT2){
  m_coupler = coupler;
  m_DT1 = DT1;
  m_DT2 = DT2;
  m_currentTime = 0.0;
  m_sim1 = new Simulator(m_coupler->getModel1(), DT1);
  m_sim2 = new Simulator(m_coupler->getModel2(), DT2);
}


void RootSimulator::initSim(){
  m_sim1->initSim();
  m_sim2->initSim();
}


void RootSimulator::simulate(const double currentTime,
			     const double simTime){
  int numIter;
  ofstream fTum("tum.dat");
  ofstream fTumDens("tumDens.dat");
  ofstream fKilledCells("killedCells.dat");
  ofstream fState("out.dat");
  ofstream fTimer("timer.dat");
  ofstream fPO2("po2.dat");
  
  numIter = simTime / m_DT1;
  for(int j(0); j < numIter; j++){
    //m_sim2->simulate(m_currentTime, 3600);
    //m_coupler->updateModel();
    m_sim1->simulate(m_currentTime, m_DT1);
    cout << m_currentTime << endl;
    fTum << m_currentTime << " " << ((ProstTissue *)m_coupler->
				     getModel1())->getNumTum()
	 << endl;
    fTumDens << m_currentTime << " " << m_coupler->getModel1()
      ->getOut()->at(0) << endl;
    fKilledCells << m_currentTime << " " << m_coupler->getModel1()
      ->getOut()->at(0) / m_coupler->getModel1()->getParam()->at(0)
		 << endl;
    for(int i(0); i < m_coupler->getModel1()->getNumComp(); i++){
      fState << m_coupler->getModel1()->getComp()->at(i)->getOut()
	->at(0) << "\t";
      /*fTimer << m_coupler->getModel1()->getComp()->at(i)->getParam()
	->at(0) << "\t";
	fPO2 << m_coupler->getModel2()->getComp()->at(i)->getOut()
	->at(0) << "\t";*/
    }
    fState << endl;
    /*fTimer << endl;
      fPO2 << endl;*/
    
    m_coupler->updateModel();
    m_currentTime += m_DT1;
  }
  fTumDens.close();
  fKilledCells.close();
  fState.close();
  fTimer.close();
  fPO2.close();
}
  
void RootSimulator::stop(){
  m_sim1->stop();
  m_sim2->stop();
}
