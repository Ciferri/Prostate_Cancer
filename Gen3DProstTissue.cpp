/**
 * @file Gen3DProstTissue.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.19.17 
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "Gen3DProstTissue.hpp"

using namespace std;

Gen3DProstTissue::Gen3DProstTissue() :
  Model(DESS, 0, 0, 0, 2, TISSUEROW*TISSUECOL*TISSUELAYER){
  int k(0);
  double inputTimer;

  m_treatment = 0;
    
  //Creation of the cells composing the tissue model
  for(int i=0;i<TISSUEROW;i++){
    for(int j=0;j<TISSUECOL;j++){
      for(int l=0;l<TISSUELAYER;l++){
	m_comp->at(k) = new ProstCell(this);
	m_tissue[i][j][l] = m_comp->at(k);
	m_comp->at(k)->calcModelOut();
	m_numOut += (m_comp->at(k))->getNumOut();
	k++;
      }
    }
  }
}


Gen3DProstTissue::Gen3DProstTissue(const string nFInPO2,
				   const string nFInTum,
				   const string nFInVes,
				   Treatment *const treatment) :
  Model(DESS, 0, 0, 0, 2, TISSUEROW*TISSUECOL*TISSUELAYER){
  int k(0);
  int selInitPhase;
  double doubTime;
  double inputPO2, inputTimer, inputTum, inputVes;
  ifstream fInPO2(nFInPO2.c_str());
  ifstream fInTum(nFInTum.c_str());
  ifstream fInVes(nFInVes.c_str());
  m_treatment = treatment;
  
  //Creation of the cells composing the tissue model
  for(int i=0;i<TISSUEROW;i++){
    for(int j=0;j<TISSUECOL;j++){
      for(int l=0;l<TISSUELAYER;l++){
	m_comp->at(k) = new ProstCell(this);
    	m_tissue[i][j][l] = m_comp->at(k);
	m_comp->at(k)->calcModelOut();
	m_numOut += (m_comp->at(k))->getNumOut();

	k++;
      }
    }
  }

  for(int i=0;i<TISSUEROW;i++){
    for(int j=0;j<TISSUECOL;j++){
      for(int l=0;l<TISSUELAYER;l++){
	for(int ii=-1;ii<=1;ii++){
	  for(int jj=-1;jj<=1;jj++){
	    for(int ll=-1;ll<=1;ll++){
	      if(ii!=0 || jj!=0 || ll!=0){
		if(i+ii>=0 && j+jj>=0 && l+ll>=0 && i+ii<TISSUEROW
		   && j+jj<TISSUECOL && l+ll<TISSUELAYER){
		  ((ProstCell *)m_tissue[i][j][l])
		    ->addToEdge(((ProstCell *)m_tissue[i+ii]
				 [j+jj][l+ll]));
		}
	      }
	    }
	  }
	}   
      }
    }
  }
  
  //Initialization of the PO2 and cells state
  
  if(fInPO2.is_open() == 0){
    cout << "An error occurred while opening initial PO2 data file"
	 << endl;
  }
  else if(fInTum.is_open() == 0){
    cout << "An error occurred while opening initial tumor" <<
      "data file" << endl;
  }
  else if(fInVes.is_open() == 0){
    cout << "An error occurred while opening initial vessel" <<
      "data file" << endl;
  }
  else{
    srand(time(NULL));
    for(int k(0); k < m_numComp; k++){
      if(fInPO2 >> inputPO2){
	((ProstCell *)m_comp->at(k))->setInPO2(inputPO2);
      }
      else{
	cout << "Insufficient data in PO2 file"<<endl;
	break;
      }
      if(fInTum >> inputTum){
	((ProstCell *)m_comp->at(k))->setInTum(inputTum);
      }
      else{
	cout << "Insufficient data in tumor file" << endl;
	break;
      }
      if(fInVes >> inputVes){
	if(inputTum && inputVes){
	  cout << "Conflict between initial data. Cell "<< k <<
	    " is both tumor and vessel" << endl;
	  break;
	}
	else{
	  ((ProstCell *)m_comp->at(k))->setInVes(inputVes);
	}
      }
      else{
	cout << "Insufficient data in vessel file" << endl;
	break;
      }
      
      doubTime = ((ProstCell *)m_comp->at(k))->getDoubTime();
      selInitPhase = rand() % 200;
      if(selInitPhase < 120){
	inputTimer = rand() % (int)(0.55 * doubTime);
      }
      else if(selInitPhase < 170){
	inputTimer = (int)(0.55 * doubTime) +
	  rand() % (int)(0.2 * doubTime);
      }
      else if(selInitPhase < 185){
	inputTimer = (int)(0.75 * doubTime) +
	  rand() % (int)(0.15 * doubTime);
      }
      else{
	inputTimer = (int)(0.9 * doubTime) +
	  rand() % (int)(0.1 * doubTime);
      }
      ((ProstCell *)m_comp->at(k))->setInTimer(inputTimer);
      (m_comp->at(k))->initModel();
      ((ProstCell *)m_comp->at(k))->setInTum(0.0);
      ((ProstCell *)m_comp->at(k))->setInVes(0.0);
    }
    fInPO2.close();
    fInTum.close();
    fInVes.close();
  }
}


Gen3DProstTissue::~Gen3DProstTissue(){
}


int Gen3DProstTissue::calcModelOut(){
  for(int k(0); k < m_numComp; k++){
    (m_comp->at(k))->calcModelOut();
  }
  return 0;
}


int Gen3DProstTissue::initModel(const double DT){
  PAR_INIT_NUM_TUM = getNumTum();
  PAR_NUM_SESSION = 0.0;
  m_flag = 0;
  cout << "Total number of cells = " << m_numComp << endl;
  cout << "Initial number of cells at G1 = " << getNumG1() << endl;
  cout << "Initial number of cells at S = " << getNumS() << endl;
  cout << "Initial number of cells at G2 = "<< getNumG2() << endl;
  cout << "Initial number of cells at M = "<< getNumM() << endl;
  cout << "Initial number of living cells = "
       << getNumAlive() << endl;
  cout << "Initial number of tumor cells = " << getNumTum() << endl;
  cout << "Initial number of vessels = " << getNumVes() << endl;
  cout << "Initial number of dead cells = " << getNumDead() << endl;
  cout << "---------------------------------------------" << endl;
  return 0;
}


//It does nothing for the moment
int Gen3DProstTissue::startModel(){
  for (int k(0) ;k < m_numComp; k++){
    (m_comp->at(k))->startModel();
  }
  return 0;
}


int Gen3DProstTissue::terminateModel(){
  for(int k(0); k < m_numComp; k++){
    (m_comp->at(k))->terminateModel();
  }
  cout << "---------------------------------------------" << endl;
  cout << "Final number of cells at G1 = " << getNumG1() << endl;
  cout << "Final number of cells at S = " << getNumS() << endl;
  cout << "Final number of cells at G2 = " << getNumG2() << endl;
  cout << "Final number of cells at M = " << getNumM() << endl;
  cout << "Final number of living cells = "
       << getNumAlive() << endl;
  cout << "Final number of tumor cells = " << getNumTum() << endl;
  cout << "Final number of vessels = " << getNumVes() << endl;
  cout << "Final number of dead cells = "<< getNumDead() << endl;
  
  return 0;
}


int Gen3DProstTissue::updateModel(const double currentTime,
				  const double DT){
  int i;
  for(int k(0); k < m_numComp; k++){
    (m_comp->at(k))->updateModel(currentTime, DT);
  }

  if(fmod(currentTime, m_treatment->getInterval()) == 0){
      i = currentTime / m_treatment->getInterval();
      if((m_treatment->getSchedule()).at(i)){
	PAR_NUM_SESSION += 1.0;
      }
  }
  
  if(getNumTum() / PAR_INIT_NUM_TUM < 0.5 && m_flag == 0){
      cout << "Total dose needed to kill 50% of tumor cells = " <<
	PAR_NUM_SESSION * m_treatment->getFraction() << endl;
      m_flag++;
    }
    else if(getNumTum() / PAR_INIT_NUM_TUM < 0.2 && m_flag == 1){
      cout << "Total dose needed to kill 80% of tumor cells = " <<
	PAR_NUM_SESSION * m_treatment->getFraction() << endl;
      m_flag++;
    }
    else if(getNumTum() / PAR_INIT_NUM_TUM < 0.1 && m_flag == 2){
      cout << "Total dose needed to kill 90% of tumor cells = " <<
	PAR_NUM_SESSION * m_treatment->getFraction() << endl;
      m_flag++;
    }
    else if(getNumTum() / PAR_INIT_NUM_TUM < 0.05 && m_flag == 3){
      cout << "Total dose needed to kill 95% of tumor cells = " <<
	PAR_NUM_SESSION * m_treatment->getFraction() << endl;
      m_flag++;
    }
    else if(getNumTum() / PAR_INIT_NUM_TUM < 0.01 && m_flag == 4){
      cout << "Total dose needed to kill 99% of tumor cells = " <<
	PAR_NUM_SESSION * m_treatment->getFraction() << endl;
      m_flag++;
    }
    else if(getNumTum() / PAR_INIT_NUM_TUM < 0.001 && m_flag == 5){
      cout << "Total dose needed to kill 99.9% of tumor cells = " <<
	PAR_NUM_SESSION * m_treatment->getFraction() << endl;
      m_flag++;
    }
  return 0;
}


int Gen3DProstTissue::getNumAlive() const{
  int count(0);
  for(int k(0); k < m_numComp; k++){
    if(((ProstCell *)m_comp->at(k))->getAlive()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumDead() const{
  int count(0);
  for(int k(0); k < m_numComp; k++){
    if(((ProstCell *)m_comp->at(k))->getDead()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumG1() const{
  int count(0);
  for(int k(0); k < m_numComp; k++){
    if(((ProstCell *)m_comp->at(k))->getG1()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumG2() const{
  int count(0);
  for(int k(0); k < m_numComp; k++){
    if(((ProstCell *)m_comp->at(k))->getG2()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumM() const{
  int count(0);
  for(int k(0); k < m_numComp; k++){
    if(((ProstCell *)m_comp->at(k))->getM()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumS() const{
  int count(0);
  for(int k(0); k < m_numComp; k++){
    if(((ProstCell *)m_comp->at(k))->getS()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumTum() const{
  int count(0);
  for(int k(0); k < m_numComp; k++){
    if(((ProstCell *)m_comp->at(k))->getTum()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumVes() const{
  int count(0);
  for(int k(0); k < m_numComp; k++){
    if(((ProstCell *)m_comp->at(k))->getVes()){
      count++;
    }
  }
  return count;
}


Treatment *Gen3DProstTissue::getTreatment() const{
  return m_treatment;
}








