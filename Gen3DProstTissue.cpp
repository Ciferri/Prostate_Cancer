/**
 * @file Gen3DProstTissue.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.19.17 
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

#include "Gen3DProstTissue.hpp"

using namespace std;

Gen3DProstTissue::Gen3DProstTissue(const int nrow, const int ncol,
				   const int nlayer) :
  Model(DESS, 0, 0, 0, 6, nrow * ncol * nlayer){
  m_nrow   = nrow;
  m_ncol   = ncol;
  m_nlayer = nlayer;
  
  for(int k(0); k < m_numComp; k++){
    m_comp->at(k) = new ProstCell(this);
    m_numOut += (m_comp->at(k))->getNumOut();
  }  
  m_tumEdge   = new vector<ProstCell *>((unsigned int)0,0);
  m_deadCells = new vector<ProstCell *>((unsigned int)0,0);
  m_treatment = 0;
}


Gen3DProstTissue::Gen3DProstTissue(const int nrow, const int ncol,
				   const int nlayer,
				   const string nFInPO2,
				   const string nFInTum,
				   const string nFInVes,
				   Treatment *const treatment) :
  Model(DESS, 0, 0, 0, 6, nrow * ncol * nlayer){

  double inputPO2, inputTum, inputVes;
  vector<Model **> map2D;
  ifstream fInPO2(nFInPO2.c_str());
  ifstream fInTum(nFInTum.c_str());
  ifstream fInVes(nFInVes.c_str());
  
  m_nrow = nrow;
  m_ncol = ncol;
  m_nlayer = nlayer;

  m_treatment=treatment;
    
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
  
  m_tumEdge = new vector<ProstCell *>((unsigned int)0,0);
  m_deadCells = new vector<ProstCell *>((unsigned int)0,0);
  
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
    ProstCell *cell;
    for(int k(0); k < m_numComp; k++){
      cell = ((ProstCell *)m_comp->at(k));
      if(fInPO2 >> inputPO2){
	cell->setInPO2(inputPO2);
      }
      else{
	cout << "Insufficient data in PO2 file"<<endl;
	break;
      }
      if(fInTum >> inputTum){
	setInTum(cell, inputTum);
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
	  setInVes(cell, inputVes);
	}
      }
      else{
	cout << "Insufficient data in vessel file" << endl;
	break;
      }
      cell->updateModel();
      setInTum(cell, 0.0);
      setInVes(cell, 0.0);
    }
    fInPO2.close();
    fInTum.close();
    fInVes.close();
  }
}


Gen3DProstTissue::~Gen3DProstTissue(){
  delete m_tumEdge;
  delete m_deadCells;
}


int Gen3DProstTissue::calcModelOut(){
  for(int k(0); k < m_numComp; k++){
    (m_comp->at(k))->calcModelOut();
  }
  return 0;
}


int Gen3DProstTissue::initModel(const double DT){
  int deadTime, doubTime;

  doubTime = ((ProstCell *)m_comp->at(0))->getDoubTime();
  deadTime = ((ProstCell *)m_comp->at(0))->getDeadTime();
  PAR_PF = pow(2.0, DT / doubTime);
  PAR_INIT_NUM_TUM = getNumTum();
  PAR_NUM_TUM = PAR_INIT_NUM_TUM;
  PAR_RF = 1.0 - pow(2.0, -DT / deadTime);
  PAR_NUM_DEAD = getNumDead();
  PAR_ACC_DOSE = 0.0;
  m_flag = 0;
  srand(time(NULL));
  
  cout <<"Total number of cells = " << m_numComp << endl;
  cout << "Initial number of living cells = " <<
    getNumAlive() << endl;
  cout << "Initial number of tumor cells = " << getNumTum() << endl;
  cout << "Initial tumor density: " <<
    (double)getNumTum() / (double)m_numComp * 100.0 << endl;
  cout << "Initial number of vessels = " << getNumVes() << endl;
  cout << "Initial vascular density: " <<
    (double)getNumVes() / (double)m_numComp * 100.0 << endl;
  cout << "Initial number of dead cells = " << getNumDead() << endl;
  cout << "Initial tumor edge size = " <<
    m_tumEdge->size() << endl;
  cout << "---------------------------------------------" << endl;
  
  return 0;
}


//It does nothing for the moment
int Gen3DProstTissue::startModel(){
  for (int k(0); k < m_numComp; k++){
    (m_comp->at(k))->startModel();
  }
  return 0;
}


int Gen3DProstTissue::terminateModel(){
  for(int k(0); k < m_numComp; k++){
    (m_comp->at(k))->terminateModel();
  }
  cout << "---------------------------------------------" << endl;
  cout << "Final number of living cells = " <<
    getNumAlive() << endl;
  cout << "Final number of tumor cells = " << getNumTum() << endl;
  cout << "Final number of vessels = " << getNumVes() << endl;
  cout << "Final number of dead cells = " << getNumDead() << endl;
  cout << "Final tumor edge size = " << m_tumEdge->size() << endl;
  
  return 0;
}


int Gen3DProstTissue::updateModel(const double currentTime,
				  const double DT){
  calcTumGrowth();
  if(m_treatment){
    if(fmod(currentTime, m_treatment->getInterval()) == 0){
      int i(currentTime / m_treatment->getInterval());
      if((m_treatment->getSchedule()).at(i)){
	calcRespToIrr();
      }
    }
    
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
  calcCellsResor();
  return 0;
}


void Gen3DProstTissue::addToDeadCells(ProstCell *const cell){
  m_deadCells->push_back(cell);
}


void Gen3DProstTissue::addToTumEdge(ProstCell *const cell){
  bool alive(false), notInEdge(true);
 
  for(int i(0); i < m_tumEdge->size(); i++){
    if(cell == m_tumEdge->at(i)){
      notInEdge = false;
      break;
    }
  }
    
  alive = cell->getAlive() || cell->getInAlive();
    
  if(notInEdge && alive){
    m_tumEdge->push_back(cell);
  }
}


void Gen3DProstTissue::calcCellsResor(){
  int numDead;
  double m;
  ProstCell * cell;
  
  numDead = getNumDead();
  PAR_NUM_DEAD -= PAR_NUM_DEAD * PAR_RF;
  
  for(int i((int)PAR_NUM_DEAD + 1); i < numDead; i++){
    m = rand() % m_deadCells->size();
    cell = m_deadCells->at(m);
    setInAlive(cell, 1.0);
    cell->updateModel();
    setInAlive(cell, 0.0);
  }
}


void Gen3DProstTissue::calcRespToIrr(){
  int preTumSize;
  double p;
  ProstCell *cell;

  preTumSize = getNumTum();
  for(int k(0); k < m_numComp; k++){
    cell = ((ProstCell*)m_comp->at(k));
    p = (double)rand() / (double)(RAND_MAX);
    if(cell->calcSF(m_treatment->getFraction()) < p){
      setInDead(cell, 1.0);
      cell->updateModel();
      setInDead(cell, 0.0);
    }
  }
	
  PAR_NUM_TUM -= preTumSize - getNumTum();
  PAR_NUM_DEAD += preTumSize - getNumTum();
  PAR_ACC_DOSE += m_treatment->getFraction();
}


void Gen3DProstTissue::calcTumGrowth(){
  int m;
  int numTum;
  ProstCell *cell;
  
  PAR_NUM_TUM *= PAR_PF;
  numTum = getNumTum();
  for(int i(numTum); i < (int)PAR_NUM_TUM; i++){
    if(m_tumEdge->size() > 0){
      m = rand() % m_tumEdge->size();
      cell = m_tumEdge->at(m);
      setInTum(cell, 1.0);
      cell->updateModel();
      setInTum(cell, 0.0);
    }
  }
}


double Gen3DProstTissue::getAlpha() const{
  return PAR_ALPHA;
}


double Gen3DProstTissue::getBeta() const{
  return PAR_BETA;
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


int Gen3DProstTissue::getNumTum() const{
  int count(0);
  for(int k(0); k<m_numComp; k++){
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


Treatment * Gen3DProstTissue::getTreatment() const{
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
  case 3:
    perc = "95";
    break;
  case 4:
    perc = "99";
  case 5:
    perc = "99.9";
    break;
  }
 
  cout << "Total dose needed to kill " << perc <<
    "% of tumor cells = " << PAR_ACC_DOSE << endl;
}

  
void Gen3DProstTissue::rmFromDeadCells(ProstCell *const cell){
  vector<ProstCell *>::iterator it;

  it = find(m_deadCells->begin(), m_deadCells->end(), cell);
  if(it != m_deadCells->end()){
    m_deadCells->erase(it);
  }
}


void Gen3DProstTissue::rmFromTumEdge(ProstCell *const cell){
  vector<ProstCell *>::iterator it;
  
  it = find(m_tumEdge->begin(), m_tumEdge->end(), cell);
  if(it != m_tumEdge->end()){
    m_tumEdge->erase(it);
  }   
}


void Gen3DProstTissue::setInAlive(ProstCell *const cell,
				  const double input){
  cell->setInAlive(input);
  
  if(input){
    rmFromDeadCells(cell);
    vector<ProstCell *> *cellEdge;
    cellEdge = cell->getEdge();
    for(int i(0); i < cellEdge->size(); i++){
      if(cellEdge->at(i)->getTum()){
	addToTumEdge(cell);
	break;
      }
    }
  }
}


void Gen3DProstTissue::setInDead(ProstCell *const cell,
				 const double input){
  cell->setInDead(input);
  
  if(input){
    addToDeadCells(cell);
    vector<ProstCell *> *cellEdge;
    vector<ProstCell *> *neiCellEdge; //neighborCellEdge
    bool nTiN; //not tumor in neighborhood
    cellEdge = cell->getEdge();
    for(int i(0); i < cellEdge->size(); i++){
      nTiN = true;
      neiCellEdge = cellEdge->at(i)->getEdge();
      for(int ii(0); ii < neiCellEdge->size(); ii++){
	if(neiCellEdge->at(ii)->getTum() &&
	   neiCellEdge->at(ii) != cell){
	  nTiN = false;
	  break;
	}
      }
      if(nTiN){
	rmFromTumEdge(cellEdge->at(i));
      }
    }
  }
}
      


void Gen3DProstTissue::setInTum(ProstCell *cell,
				const double input){
  cell->setInTum(input);
  
  if(input){
    rmFromTumEdge(cell);
    vector<ProstCell *> *cellEdge;
    cellEdge = cell->getEdge();
    for(int i(0); i < cellEdge->size(); i++){
      addToTumEdge(cellEdge->at(i));
    }
  }
}


void Gen3DProstTissue::setInVes(ProstCell *cell,
				const double input){
  cell->setInVes(input);
  
  if(input){
    rmFromTumEdge(cell);
  }
}








