
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
  Model(DESS, 0, 0, 0, 1, nrow * ncol * nlayer){
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
  Model(DESS, 0, 0, 0, 1, nrow * ncol * nlayer){
  int selInitPhase;
  double doubTime;
  double inputPO2, inputTimer, inputTum, inputVes;
  vector<Model **> map2D;
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
  for(int i(0); i < m_nrow * m_nlayer; i++){
    map2D.push_back(&(m_comp->at(i * m_ncol)));
  }
  for(int l(0); l < m_nlayer; l++){
    m_map.push_back(&(map2D[l * m_nrow]));
  }
  
  //Definition of each cell's edge
  for(int l(0); l < m_nlayer; l++){
    for(int i(0); i < m_nrow; i++){
      for(int j(0); j < m_ncol; j++){
	for(int ll(-1); ll <= 1; ll++){
	  for(int ii(-1); ii <= 1; ii++){
	    for(int jj(-1); jj <= 1; jj++){
	      if(ii != 0 || jj != 0 || ll != 0){
		if(l + ll >= 0 && l + ll < m_nlayer && i + ii >= 0
		   && i + ii < m_nrow && j + jj >= 0 &&
		   j + jj < m_ncol){
		  ((ProstCell *)m_map[l][i][j])
		    ->addToEdge(((ProstCell *)m_map[l + ll][i + ii]
				 [j + jj]));
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
  for(int k(0); k < m_numComp; k++){
    (m_comp->at(k))->updateModel(currentTime, DT);
  }
  
  if(m_treatment){
    int print(0);
    double tumSurv;
    tumSurv = getNumTum() / PAR_INIT_NUM_TUM;
    if(tumSurv < 0.5){
      print = 1;
    }
    if(tumSurv < 0.2){
      print = 2;
    }
    if(tumSurv < 0.1){
      print = 3;
    }
    if(tumSurv < 0.05){
      print = 4;
    }
    if(tumSurv < 0.01){
      print = 5;
    }
    if(tumSurv < 0.001){
      print = 6;
    }
    if(print > m_flag){
      printNeededDose();
      m_flag++;
    }
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


void Gen3DProstTissue::printNeededDose() const{
  string perc;
  
  switch(m_flag){
  case 0:
    perc = "50";
    break;
  case 1:
    perc = "80";
    break;
  case 2:
    perc = "90";
    break;
  case 3:
    perc = "95";
    break;
  case 4:
    perc = "99";
    break;
  case 5:
    perc = "99.9";
    break;
  }
 
  cout << "Total dose needed to kill " << perc <<
    "% of tumor cells = " <<
    ((ProstCell*)m_comp->at(0))->getAccDose() << endl;
}


