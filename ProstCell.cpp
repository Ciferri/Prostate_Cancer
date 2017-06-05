/**
 * @file ProstCell.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.19.17 
 */

#include "ProstCell.hpp"
#include "Treatment.hpp"

#include <iostream>
#include <math.h>

using namespace std;

ProstCell::ProstCell() : Model(DESS, 5, 4, 1, 8, 0){
  ST_ALIVE = 1.0;		
  ST_DEAD  = 0.0;
  ST_TUM   = 0.0;		
  ST_VES   = 0.0;

  PAR_DOUBTIME = 1008; //h 
  PAR_DEADTIME = 234; //h

  PAR_M     = 3.0; //adim.
  PAR_K     = 3.0; //mmHg
  PAR_PO2   = 3.5; //mmHg
  PAR_ALPHA = 0; //Gy^-1
  PAR_BETA  = 0; //Gy^-2
  
  m_parent = 0;
  m_edge = new vector<ProstCell *>((unsigned int)0, 0);
}

ProstCell::ProstCell(Model *const parent) :
  Model(DESS, 5, 4, 1, 8, 0){
  ST_ALIVE = 1.0;		
  ST_DEAD  = 0.0;
  ST_TUM   = 0.0;		
  ST_VES   = 0.0;
  
  PAR_DOUBTIME = 1008; //h
  PAR_DEADTIME = 234; //h

  PAR_M     = 3.0; //adim.
  PAR_K     = 3.0; //mmHg
  PAR_PO2   = 3.5; //mmHg
  PAR_ALPHA = 0; //Gy^-1
  PAR_BETA  = 0; //Gy^-2
  
  m_parent  = parent;
  m_edge = new vector<ProstCell *>((unsigned int)0, 0);
}


ProstCell::~ProstCell(){
}


int ProstCell::calcModelOut(){
  OUT_STATE = ST_ALIVE + 2*ST_TUM + 3*ST_VES + 4*ST_DEAD;	    
  return 0;
}


int ProstCell::initModel(const double DT){	
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
  ST_ALIVE = (ST_ALIVE || IN_ALIVE)  && !IN_TUM && !IN_DEAD &&
    !IN_VES;   
  ST_TUM   = (ST_TUM || IN_TUM) && !ST_ALIVE && !IN_DEAD;
  ST_DEAD  = (ST_DEAD || IN_DEAD) && !ST_ALIVE && !ST_TUM;
  ST_VES   = (ST_VES || IN_VES) && !IN_DEAD;

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
    PAR_ALPHA = 0.15;
    PAR_BETA  = 0.048;
  }
  
  PAR_PO2 = IN_PO2;
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


double ProstCell::calcSF(const double fraction) const{
  double  SF;
  
  SF = exp(-PAR_ALPHA / PAR_M * fraction * calcOER() -
	   PAR_BETA / (PAR_M * PAR_M) * fraction * fraction *
	   calcOER() * calcOER());
}


double ProstCell::getAlive() const{
  return ST_ALIVE;
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


double ProstCell::getInAlive() const{
  return IN_ALIVE;
}


double ProstCell::getTum() const{
  return ST_TUM;
}


double ProstCell::getVes() const{
  return ST_VES;
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


void ProstCell::setInTum(const double input){
  IN_TUM = input;
}


void ProstCell::setInVes(const double input){
  IN_VES = input;
}









