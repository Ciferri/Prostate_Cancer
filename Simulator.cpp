/*Simulator.cpp*/

#include <stdlib.h>

#include "Simulator.hpp"

using namespace std;

Simulator::Simulator(){
  m_model=0;
  m_currentTime = 0.0;
  m_DT = 1; //h
  m_outList = new OutputDataList();
  m_simMeth = new SimMeth(euler);
  m_outFile.open("out.dat");
}


Simulator::Simulator (Model* model, double currentTime, double DT){
  m_model = model;
  m_currentTime = currentTime;
  m_DT = DT;
  m_outList = new OutputDataList();
  m_simMeth = new SimMeth (euler);
  m_outFile.open("out.dat");
}


Simulator::~Simulator(){
  delete m_model;
  delete m_outList;
}


void Simulator::setModel(Model *model){
  m_model = model;
}


void Simulator::stop(){	
}


void Simulator::start(double simTime){
  int index(0);
  int toto(0);
  //Creation of each element of outList
  for (int i=0;i<m_model->getNumOut();i++){
    //vector of numOut pointers on globalTime/DT samples
    m_outList->push_back(new OutputData());

  }
  //Initialization of the state of every cell
  m_model->initModel(m_DT);
  
  if(m_model->getNumComp()==0){
    for(int i=0;i<m_model->getNumOut();i++){
      (m_outList->at(i))->push_back((m_model->getOut())->at(i));
      //Copy of the outputs in outputList
      m_outFile<<(*(m_outList->at(i))).at(toto)<<'\t';
    }           
  }
  else {
    for(int k=0;k<m_model->getNumComp();k++){
      for(int i=0;i<(m_model->getComp()->at(k))->getNumOut();
	  i++){
	index = (m_model->getComp()->at(k))->getNumOut()*k+i;
	(m_outList->at(index))->
	  push_back((m_model->getComp()->at(k)->getOut())->
		    at(i));
	//Copy of the outputs in outputList
	m_outFile<<(*(m_outList->at(index))).at(toto)<<'\t';
      }
    }
  }
  m_outFile<<endl ;
  toto++;


  //It does nothing for the moment
  m_model->startModel();
    
  // Simulation loop
  int numIter;
  numIter = simTime/m_DT;
  
  for(int j=0;j<numIter;j++) {
    m_currentTime = j*m_DT;
    //Update of the state of every cell composing the tissue
    m_model->updateModel(m_currentTime, m_DT);
    m_model->calcModelOut();
             
    if ((j%6)==0) {
      //cout << "Simulator time = " << m_currentTime << endl;
      if(m_model->getNumComp()==0){
	for(int i=0;i<m_model->getNumOut();i++){
	  (m_outList->at(i))->push_back((m_model->getOut())->at(i));
	  //Copy of the outputs in outputList
	  m_outFile<<(*(m_outList->at(i))).at(toto)<<'\t';
	}           
      }
      else {
	for(int k=0;k<m_model->getNumComp();k++){
	  for(int i=0;i<(m_model->getComp()->at(k))->getNumOut();
	      i++){
	    index = (m_model->getComp()->at(k))->getNumOut()*k+i;
	    (m_outList->at(index))->
	      push_back((m_model->getComp()->at(k)->getOut())->
			at(i));
	    //Copy of the outputs in outputList
	    m_outFile<<(*(m_outList->at(index))).at(toto)<<'\t';
	  }
	}
      }
      m_outFile<<endl ;
      toto++;
    }
  }
  m_model->terminateModel();
  m_outFile.close();
}
