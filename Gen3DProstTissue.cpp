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

Gen3DProstTissue::Gen3DProstTissue(const int nrow, const int ncol,
				   const int nlayer) :
  Model(DESS, 0, 0, 0, 2, nrow * ncol * nlayer){
  m_nrow   = nrow;
  m_ncol   = ncol;
  m_nlayer = nlayer;
  
  m_treatment = 0;
    
  //Creation of the cells composing the tissue model
  for(int k(0); k < m_numComp; k++){
    m_comp->at(k) = new ProstCell(this);
    m_numOut += (m_comp->at(k))->getNumOut();
  }
}


Gen3DProstTissue::Gen3DProstTissue(const int nrow, const int ncol,
				   const int nlayer,
				   const string nFInPO2,
				   const string nFInTum,
				   const string nFInVes,
				   Treatment *const treatment) :
  Model(DESS, 0, 0, 0, 2, nrow * ncol * nlayer){
  int selInitPhase;
  int m;
  int x, y, z;
  int *coordxyz;
  double doubTime;
  double inputPO2, inputTimer, inputTum, inputVes;
  ifstream fInPO2(nFInPO2.c_str());
  ifstream fInTum(nFInTum.c_str());
  ifstream fInVes(nFInVes.c_str());

  m_nrow   = nrow;
  m_ncol   = ncol;
  m_nlayer = nlayer;
  
  m_treatment = treatment;
  
  //Creation of the cells composing the tissue model
  for(int k(0); k < m_numComp; k++){
    m_comp->at(k) = new ProstCell(this);
    m_numOut += (m_comp->at(k))->getNumOut();
  }
  //Definition of each cell's edge
  for(int k(0); k < m_numComp; k++){
    coordxyz = kToXyz(k);
    x = coordxyz[0];
    y = coordxyz[1];
    z = coordxyz[2];
    for(int i(-1); i <= 1; i++){
      for(int j(-1); j <= 1; j++){
	for(int l(-1); l <= 1; l++){
	  if(i!=0 || j!=0 || l!=0){
	    m = xyzTok(x + i, y + j, z +l);
	    if(m != -1){
	      ((ProstCell *)m_comp->at(k))
		->addToEdge(((ProstCell *)m_comp->at(m)));
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
  cout << "Initial tumor density: " <<
    (double)getNumTum() / (double)m_numComp * 100.0 << endl;
  cout << "Initial number of vessels = " << getNumVes() << endl;
  cout << "Initial vascular density: " <<
    (double)getNumVes() / (double)m_numComp * 100.0 << endl;
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


int *Gen3DProstTissue::kToXyz(const int k) const{
  int *tab= new int[3];
  if(k > -1 && k < m_nrow * m_ncol * m_nlayer){
    tab[0] = (k % (m_nrow * m_ncol)) / m_ncol;
    tab[1] = (k % (m_nrow * m_ncol)) % m_ncol;
    tab[2] = k / (m_nrow * m_ncol);
  }
  else{
    tab[0] = -1;
    tab[1] = -1;
    tab[2] = -1;
  }
  return tab;
}


int Gen3DProstTissue::xyzTok(const int x, const int y,
			     const int z) const{
  if(x > -1 && x < m_nrow && y > -1 && y < m_ncol && z > -1 &&
     z < m_nlayer){
    return x * m_ncol + y + z * m_nrow * m_ncol;
  }
  else{
    return -1;
  }
}






