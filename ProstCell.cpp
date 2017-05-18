/*
 *  prostCell.cpp
 *
 *  Copyright (c) 2001 INSERM. All rights reserved.
 *	03.23.2017
 */

#include "ProstCell.hpp"
#include "Treatment.hpp"

#include <iostream>
#include <math.h>

using namespace std;

ProstCell::ProstCell() : Model(DESS, 5, 4, 1, 8, 0){
  ST_ALIVE = 1.0;		
  ST_DEAD = 0.0;
  ST_TUMOR = 0.0;		
  ST_VES = 0.0;
  PAR_DOUBTIME = 1008; //h 
  PAR_DEADTIME = 168; //h
  PAR_M = 3.0; //adim.
  PAR_K = 3.0; //mmHg
  PAR_PO2 = 3.5; //mmHg
  PAR_ALPHA = 0.15; //Gy^-1
  PAR_BETA = 0.048; //Gy^-2
  m_parent = 0;
}

ProstCell::ProstCell(Model *parent) : Model(DESS, 5, 4, 1, 8, 0){
  ST_ALIVE = 1.0;		
  ST_DEAD = 0.0;
  ST_TUMOR = 0.0;		
  ST_VES = 0.0;
  PAR_DOUBTIME = 1008; //h
  PAR_DEADTIME = 168; //h
  PAR_M = 3.0; //adim.
  PAR_K = 3.0; //mmHg
  PAR_PO2 = 3.5; //mmHg
  PAR_ALPHA = 0.15; //Gy^-1
  PAR_BETA = 0.048; //Gy^-2
  m_parent = parent;
}


ProstCell::~ProstCell(){
}


int ProstCell::ModelInitSim(double DT){	
  return 0;
}


int ProstCell::ModelOut(){
  OUT_STATE = ST_ALIVE + 2*ST_TUMOR + 3*ST_VES + 4*ST_DEAD;	    
  return 0;
}


int ProstCell::ModelStart(){
  return 0;
}


int ProstCell::ModelTerminate(){
  return 0;
}


int ProstCell::ModelUpdate(double currentTime, double DT){
  ST_ALIVE = (ST_ALIVE || IN_ALIVE)  && !IN_TUMOR && !IN_VES;   
  ST_TUMOR = (ST_TUMOR || IN_TUMOR) && !IN_DEAD;
  ST_DEAD  = (ST_DEAD || IN_DEAD) && !ST_ALIVE && !ST_VES;
  ST_VES = ST_VES || IN_VES;

  if(ST_ALIVE){
    PAR_ALPHA = 0;
    PAR_BETA = 0;
  }
  if(ST_DEAD){
    PAR_ALPHA = 0;
    PAR_BETA = 0;
  }
  if(ST_VES){
    PAR_ALPHA = 0;
    PAR_BETA = 0;
  }
  if(ST_TUMOR){
    PAR_ALPHA = 0.15;
    PAR_BETA = 0.048;
  }
  
  PAR_PO2 = IN_PO2;
  return 0;
}


double ProstCell::CalcOER() const{
  double OER;
  
  OER = (PAR_M * PAR_PO2 + PAR_K)/(PAR_PO2 + PAR_K); //mmHg
  return OER;
}


double ProstCell::CalcSF() const{
  double fraction, SF;
  fraction = ((Gen3DProstTissue *)m_parent)->getTreatment()->
    getFraction();
  SF = exp(-PAR_ALPHA/PAR_M * fraction * CalcOER() -
	   PAR_BETA/(PAR_M*PAR_M) * fraction*fraction *
	   CalcOER()*CalcOER());
  return SF;
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


double ProstCell::getInAlive() const{
  return IN_ALIVE;
}


double ProstCell::getTumor() const{
  return ST_TUMOR;
}


double ProstCell::getVes() const{
  return ST_VES;
}


void ProstCell::setInAlive(double input){
  IN_ALIVE = input;
}


void ProstCell::setInDead(double input){
  IN_DEAD = input;
}


void ProstCell::setInPO2(double input){
  IN_PO2 = input;
  
}


void ProstCell::setInTumor(double input){
  IN_TUMOR = input;
}


void ProstCell::setInVes(double input){
  IN_VES = input;
}









