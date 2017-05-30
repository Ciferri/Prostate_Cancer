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
  m_tumEdge   = new vector<int>(0,0);
  m_deadCells = new vector<int>(0,0);
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
  ifstream fInPO2(nFInPO2.c_str());
  ifstream fInTum(nFInTum.c_str());
  ifstream fInVes(nFInVes.c_str());
  
  m_nrow = nrow;
  m_ncol = ncol;
  m_nlayer = nlayer;
  
  for(int k(0); k < m_numComp; k++){
    m_comp->at(k) = new ProstCell(this);
    m_numOut += (m_comp->at(k))->getNumOut();
  }
  
  m_tumEdge = new vector<int>(0,0);
  m_deadCells = new vector<int>(0,0);

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
    for(int k(0); k < m_numComp; k++){
      if(fInPO2 >> inputPO2){
	((ProstCell *)m_comp->at(k))->setInPO2(inputPO2);
      }
      else{
	cout << "Insufficient data in PO2 file"<<endl;
	break;
      }
      if(fInTum >> inputTum){
	setInTum(k, inputTum);
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
	  setInVes(k, inputVes);
	}
      }
      else{
	cout << "Insufficient data in vessel file" << endl;
	break;
      }
      (m_comp->at(k))->updateModel();
      setInTum(k, 0.0);
      setInVes(k, 0.0);
    }
    fInPO2.close();
    fInTum.close();
    fInVes.close();
  }

  m_treatment=treatment;
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
  PAR_NUM_SESSION = 0.0;
  m_flag = 0;
  srand(time(NULL));
  
  cout <<"Total number of cells = " << m_numComp << endl;
  cout << "Initial number of living cells = " <<
    getNumAlive() << endl;
  cout << "Initial number of tumor cells = " << getNumTum() << endl;
  cout << "Initial number of vessels = " << getNumVes() << endl;
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
  int k, m;
  int numTum, numDead;
  
  PAR_NUM_TUM *= PAR_PF;
  numTum = getNumTum();
  for(int i(numTum); i < (int)PAR_NUM_TUM; i++){
    if(m_tumEdge->size() > 0){
      m = rand() % m_tumEdge->size();
      k = m_tumEdge->at(m);
      setInTum(k, 1.0);
      (m_comp->at(k))->updateModel();
      setInTum(k, 0.0);
    }
  }
  
  if(m_treatment){
    if(fmod(currentTime, m_treatment->getInterval()) == 0){
      int i(currentTime / m_treatment->getInterval());

      if((m_treatment->getSchedule()).at(i)){
	int preTumSize;
	double n;

	preTumSize = getNumTum();
	for(int k(0); k < m_numComp; k++){
	  n = (double)rand() / (double)(RAND_MAX);
	  if(((ProstCell *)m_comp->at(k))->calcSF() < n){
	    setInDead(k, 1.0);
	    (m_comp->at(k))->updateModel();
	    setInDead(k, 0.0);
	  }
	}
	
	PAR_NUM_TUM -= preTumSize - getNumTum();
	PAR_NUM_DEAD += preTumSize - getNumTum();
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
    else if((getNumTum()) / PAR_INIT_NUM_TUM < 0.05 && m_flag == 3){
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

    numDead = getNumDead();
    PAR_NUM_DEAD -= PAR_NUM_DEAD * PAR_RF;
  
    for(int i((int)PAR_NUM_DEAD + 1); i < numDead; i++){
      m = rand() % m_deadCells->size();
      k = m_deadCells->at(m);
      setInAlive(k, 1.0);
      (m_comp->at(k))->updateModel();
      setInAlive(k, 0.0);
    }
  }
  return 0;
}


void Gen3DProstTissue::addToDeadCells(const int k){
  m_deadCells->push_back(k);
}


void Gen3DProstTissue::addToEdge(const int x, const int y,
				 const int z){
  int k;
  bool alive(false), notInEdge(true);
 
  k = xyzTok(x, y, z);
  if(k != -1){
    for(int i(0); i < m_tumEdge->size(); i++){
      if(k == m_tumEdge->at(i)){
	notInEdge = false;
	i = m_tumEdge->size();
      }
    }
    
    alive = (((ProstCell *)m_comp->at(k))->getAlive() ||
	     ((ProstCell *)m_comp->at(k))->getInAlive());
    
    if(notInEdge && alive){
      m_tumEdge->push_back(k);
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


void Gen3DProstTissue::removeFromDeadCells(const int k){
  vector<int>::iterator it;

  it = find(m_deadCells->begin(), m_deadCells->end(), k);
  if(it != m_deadCells->end()){
    m_deadCells->erase(it);
  }
}


void Gen3DProstTissue::removeFromEdge(const int k){
  vector<int>::iterator it;
  
  it = find(m_tumEdge->begin(), m_tumEdge->end(), k);
  if(it != m_tumEdge->end()){
    m_tumEdge->erase(it);
  }   
}


void Gen3DProstTissue::setInAlive(const int k, const double input){
  int x, y, z;
  int m;
  int *coordxyz;
  
  ((ProstCell *)m_comp->at(k))->setInAlive(input);
  
  if(input){
    removeFromDeadCells(k);
    coordxyz = kToXyz(k);
    x = coordxyz[0];
    y = coordxyz[1];
    z = coordxyz[2];
    for(int i(-1); i < 2; i++){
      for(int j(-1); j < 2; j++){
	for(int l(-1); l < 2; l++){
	  if(i != 0 || j != 0 || l != 0){
	    m = xyzTok(x + i, y + j, z + l);
	    if(m != -1){
	      if(((ProstCell *)m_comp->at(m))->getTum()){
		addToEdge(x, y, z);
		i = 2;
		j = 2;
		l = 2;
	      }
	    }
	  }
	}
      }
    }
  }
}


void Gen3DProstTissue::setInDead(const int k, const double input){
  int x, xx, y, yy, z, zz;
  int m, n;
  int *coordxyz;
  bool nTiNCSDC; //not Tumor in Neighborhood of the Cells
  //Surrounding the Dead Cell
  
  ((ProstCell *)m_comp->at(k))->setInDead(input);
  
  if(input){
    addToDeadCells(k);
    coordxyz = kToXyz(k);
    x = coordxyz[0];
    y = coordxyz[1];
    z = coordxyz[2];
    for(int i(-1); i < 2; i++){
      for(int j(-1); j < 2; j++){
	for(int l(-1); l < 2; l++){
	  if(i != 0 || j != 0 || l != 0){
	    m = xyzTok(x + i, y + j, z + l);
	    if(m != -1){
	      nTiNCSDC = true;
	      xx = x + i;
	      yy = y + j;
	      zz = z + l;
	      for(int ii(-1); ii < 2; ii++){
		for(int jj(-1); jj < 2; jj++){
		  for(int ll(-1); ll < 2; ll++){
		    if((ii != 0 || jj != 0 || ll != 0) &&
		       (ii != -i || jj !=-j || ll != -l)){
		      n = xyzTok(xx + ii, yy + jj, zz + ll);
		      if(n != -1){
			nTiNCSDC = nTiNCSDC &&
			  !(((ProstCell *)m_comp->at(n))
			    ->getTum());		       
		      }
		    }
		  }
		}  
	      }
	      if(nTiNCSDC){
		removeFromEdge(m);
	      }
	    }
	  }
	}
      }
    }
  }
}


void Gen3DProstTissue::setInTum(const int k, const double input){
  int x, y, z;
  int *coordxyz;
  
  ((ProstCell *)m_comp->at(k))->setInTum(input);
  
  if(input){
    removeFromEdge(k);
    coordxyz = kToXyz(k);
    x = coordxyz[0];
    y = coordxyz[1];
    z = coordxyz[2];
    for(int i(-1); i < 2; i++){
      for(int j(-1); j < 2; j++){
	for(int l(-1); l < 2; l++){
	  if(i != 0 || j != 0 || l != 0){
	    addToEdge(x+i, y+j, z+l);
	  }
	}
      }
    }
  }
}


void Gen3DProstTissue::setInVes(const int k, const double input){
  ((ProstCell *)m_comp->at(k))->setInVes(input);
  
  if(input){
    removeFromEdge(k);
  }
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






