/**
 * @file ProstCell.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.19.17 
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "ProstCell.hpp"

using namespace std;

ProstCell::ProstCell() : Model(DESS, 6, 8, 1, 12, 0){
  ST_ALIVE = 1.0;		
  ST_DEAD  = 0.0;
  ST_TUM   = 0.0;		
  ST_VES   = 0.0;

  PAR_TIMER    = 0; //h
  PAR_DOUBTIME = 1008; //h 
  PAR_DEADTIME = 234; //h 
  PAR_LIMG1S   = 0.55 * PAR_DOUBTIME;
  PAR_LIMSG2   = 0.75 * PAR_DOUBTIME;
  PAR_LIMG2M   = 0.9 * PAR_DOUBTIME;
  
  ST_G1 = 0 <= PAR_TIMER && PAR_TIMER < PAR_LIMG1S;
  ST_S  = PAR_LIMG1S <= PAR_TIMER && PAR_TIMER < PAR_LIMSG2;
  ST_G2 = PAR_LIMSG2 <= PAR_TIMER && PAR_TIMER < PAR_LIMG2M;
  ST_M  = PAR_LIMG2M <= PAR_TIMER && PAR_TIMER <= PAR_DOUBTIME;

  PAR_M     = 3.0; //adim.
  PAR_K     = 3.0; //mmHg
  PAR_PO2   = 3.5; //mmHg
  PAR_ALPHA = 0; //Gy^-1
  PAR_BETA  = 0; //Gy^-2

  PAR_ACC_DOSE = 0; //Gy
  
  m_parent = 0;
  m_edge = new vector<ProstCell *>((unsigned int)0, 0);
  m_treatment = 0;
  
}


ProstCell::ProstCell(Model *const parent) :
  Model(DESS, 6, 8, 1, 12, 0){
  ST_ALIVE = 1.0;		
  ST_DEAD  = 0.0;
  ST_TUM   = 0.0;		
  ST_VES   = 0.0;

  PAR_TIMER    = 0; //h
  PAR_DOUBTIME = 1008; //h 
  PAR_DEADTIME = 234; //h
  PAR_LIMG1S   = 0.55 * PAR_DOUBTIME;
  PAR_LIMSG2   = 0.75 * PAR_DOUBTIME;
  PAR_LIMG2M   = 0.9 * PAR_DOUBTIME;
  
  ST_G1 = 0 <= PAR_TIMER && PAR_TIMER < PAR_LIMG1S;
  ST_S  = PAR_LIMG1S <= PAR_TIMER && PAR_TIMER < PAR_LIMSG2;
  ST_G2 = PAR_LIMSG2 <= PAR_TIMER && PAR_TIMER < PAR_LIMG2M;
  ST_M  = PAR_LIMG2M <= PAR_TIMER && PAR_TIMER <= PAR_DOUBTIME;
  
  PAR_M     = 3.0; //adim.
  PAR_K     = 3.0; //mmHg
  PAR_PO2   = 3.5; //mmHg
  PAR_ALPHA = 0; //Gy^-1
  PAR_BETA  = 0; //Gy^-2

  PAR_ACC_DOSE = 0; //Gy
  
  m_parent = parent;
  m_edge = new vector<ProstCell *>((unsigned int)0, 0);
  m_treatment = ((Gen3DProstTissue *)m_parent)->getTreatment();
}


ProstCell::~ProstCell(){
}


int ProstCell::calcModelOut(){
  OUT_STATE = ST_ALIVE + 2*ST_TUM + 3*ST_VES + 4*ST_DEAD;	    
  return 0;
}


int ProstCell::initModel(const double DT){
  ST_ALIVE = (ST_ALIVE || IN_ALIVE)  && !IN_TUM && !IN_DEAD &&
    !IN_VES;   
  ST_TUM   = (ST_TUM || IN_TUM) && !ST_ALIVE && !IN_DEAD;
  ST_DEAD  = (ST_DEAD || IN_DEAD) && !ST_ALIVE && !ST_TUM;
  ST_VES   = (ST_VES || IN_VES) && !IN_DEAD;
  setInTum(0.0);
  setInDead(0.0);
  setInVes(0.0);
  
  PAR_TIMER = IN_TIMER;
  ST_G1 = !ST_DEAD && 0 <= PAR_TIMER && PAR_TIMER < PAR_LIMG1S;
  ST_S  = !ST_DEAD && PAR_LIMG1S <= PAR_TIMER &&
    PAR_TIMER < PAR_LIMSG2;
  ST_G2 = !ST_DEAD && PAR_LIMSG2 <= PAR_TIMER &&
    PAR_TIMER < PAR_LIMG2M;
  ST_M  = !ST_DEAD && PAR_LIMG2M <= PAR_TIMER &&
    PAR_TIMER <= PAR_DOUBTIME; 

  if(ST_ALIVE){
    PAR_ALPHA = 0;
    PAR_BETA  = 0;
  }
  if(ST_DEAD){
    PAR_ALPHA = 0;
    PAR_BETA  = 0;
  }
  if(ST_VES){
    PAR_ALPHA = 0;
    PAR_BETA  = 0;
  }
  if(ST_TUM){
    if(ST_G1){
      PAR_ALPHA = 0.158;
      PAR_BETA  = 0.051;
    }
    if(ST_S){
      PAR_ALPHA = 0.113;
      PAR_BETA  = 0.037;
    }
    if(ST_G2){
      PAR_ALPHA = 0.169;
      PAR_BETA  = 0.055;
    }
    if(ST_M){
      PAR_ALPHA = 0.189;
      PAR_BETA  = 0.061;
    }
  }

  PAR_PO2 = IN_PO2;
  
  return 0;
}


int ProstCell::startModel(){
  return 0;
}


int ProstCell::terminateModel(){
  return 0;
}


int ProstCell::updateModel(const double currentTime,
			   const double DT){
  int i;
  double p;
  ProstCell *newTumCell(0);
  
  ST_ALIVE = (ST_ALIVE || IN_ALIVE)  && !IN_TUM && !IN_DEAD &&
    !IN_VES;   
  ST_TUM   = (ST_TUM || IN_TUM) && !ST_ALIVE && !IN_DEAD;
  ST_DEAD  = (ST_DEAD || IN_DEAD) && !ST_ALIVE && !ST_TUM;
  ST_VES   = (ST_VES || IN_VES) && !IN_DEAD;
  setInAlive(0.0);
  setInTum(0.0);
  setInDead(0.0);
  setInVes(0.0);

  PAR_TIMER += DT;
  ST_G1 = !ST_DEAD && 0 <= PAR_TIMER && PAR_TIMER < PAR_LIMG1S;
  ST_S  = !ST_DEAD && PAR_LIMG1S <= PAR_TIMER &&
    PAR_TIMER < PAR_LIMSG2;
  ST_G2 = !ST_DEAD && PAR_LIMSG2 <= PAR_TIMER &&
    PAR_TIMER < PAR_LIMG2M;
  ST_M  = !ST_DEAD && PAR_LIMG2M <= PAR_TIMER &&
    PAR_TIMER <= PAR_DOUBTIME;
    
  if(ST_ALIVE){
    PAR_ALPHA = 0;
    PAR_BETA  = 0;
  }
  if(ST_DEAD){
    PAR_ALPHA = 0;
    PAR_BETA  = 0;
  }
  if(ST_VES){
    PAR_ALPHA = 0;
    PAR_BETA  = 0;
  }
  if(ST_TUM){
    if(ST_G1){
      PAR_ALPHA = 0.158;
      PAR_BETA  = 0.051;
    }
    if(ST_S){
      PAR_ALPHA = 0.113;
      PAR_BETA  = 0.037;
    }
    if(ST_G2){
      PAR_ALPHA = 0.169;
      PAR_BETA  = 0.055;
    }
    if(ST_M){
      PAR_ALPHA = 0.189;
      PAR_BETA  = 0.061;
    }
  }

  PAR_PO2 = IN_PO2;
  
  //Tumor growth
  if(PAR_TIMER >= PAR_DOUBTIME){
    PAR_TIMER = 0;
    if(ST_TUM && ((Gen3DProstTissue *)m_parent)->getNumAlive()){
      newTumCell = searchSpace();
      newTumCell->setInTum(1.0);
      newTumCell->PAR_TIMER = 0;
    }     
  }

  if(m_treatment){
    //Response to radiation
    if(fmod(currentTime, m_treatment->getInterval()) == 0){
      i = currentTime / m_treatment->getInterval();
      if((m_treatment->getSchedule()).at(i)){
	p = (double)rand() / (double)(RAND_MAX);
	if(calcSF() < p){
	  setInDead (1.0);
	  p = (double)rand() / (double)(RAND_MAX);
	  if(p < 0.8){
	    PAR_DEADTIME = 234; //apoptotic
	  }
	  else{
	    PAR_DEADTIME = 468; //necrotic
	  }
	  PAR_TIMER = 0;
	}
	PAR_ACC_DOSE += m_treatment->getFraction();
      }
    }
    //Resorption
    if(ST_DEAD){
      p = (double)rand() / (double)(RAND_MAX);
      if(calcRF(DT) > p){
	PAR_TIMER = 0;
	setInAlive(1.0);
      }
    }
  }
  return 0;
}


void ProstCell::addToEdge(ProstCell *const cell){
  m_edge->push_back(cell);
}


double ProstCell::calcOER() const{
  double OER;
  
  OER = (PAR_M * PAR_PO2 + PAR_K) / (PAR_PO2 + PAR_K); //mmHg
  return OER;
}


double ProstCell::calcRF(const double DT) const{
  double RF;

  RF = 1.0 - pow(2.0, -DT / PAR_DEADTIME);
  return RF;
}


double ProstCell::calcSF() const{
  double fraction, SF;
  
  fraction = m_treatment->getFraction();
  SF = exp(-PAR_ALPHA / PAR_M * fraction * calcOER() -
	   PAR_BETA / (PAR_M * PAR_M) * fraction * fraction *
	   calcOER() * calcOER());
  return SF;
}


double ProstCell::getAlive() const{
  return ST_ALIVE;
}


double ProstCell::getAccDose() const{
  return PAR_ACC_DOSE;
}


double ProstCell::getDead() const{
  return ST_DEAD;
}


double ProstCell::getDeadTime() const{
  return PAR_DEADTIME;
}


double ProstCell::getDoubTime() const{
  return PAR_DOUBTIME;
}


vector<ProstCell *> *ProstCell::getEdge() const{
  return m_edge;
}


double ProstCell::getG1() const{
  return ST_G1;
}


double ProstCell::getG2() const{
  return ST_G2;
}


double ProstCell::getM() const{
  return ST_M;
}


int ProstCell::getOutState() const{
  return (int)OUT_STATE;
}


double ProstCell::getS() const{
  return ST_S;
}


double ProstCell::getTum() const{
  return ST_TUM;
}


double ProstCell::getVes() const{
  return ST_VES;
}


ProstCell *ProstCell::searchSpace() const{
  int edgeSize, m;
  
  edgeSize = m_edge->size();
  m = rand() % edgeSize;
  for(int n(0); n < edgeSize; n++){
    if(((ProstCell *)m_edge->at(m))->ST_ALIVE){
      return ((ProstCell *)m_edge->at(m));
    }
    m++;
    if(m == edgeSize){
      m = 0;
    }
  }
  return ((ProstCell *)m_edge->at(m))->searchSpace();
}


void ProstCell::setInAlive(const double input){
  IN_ALIVE = input;
}


void ProstCell::setInDead(const double input){
  IN_DEAD = input;
}


void ProstCell::setInPO2(const double input){
  IN_PO2 = input;
  
}


void ProstCell::setInTimer(const double input){
  IN_TIMER = input;
}


void ProstCell::setInTum(const double input){
  IN_TUM = input;
}


void ProstCell::setInVes(const double input){
  IN_VES = input;
}









