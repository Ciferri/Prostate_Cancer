/**
 * @file SimMeth.cpp
 * @brief
 * @author Alfredo Hernandez
 */

#include "SimMeth.hpp"

#include <stdlib.h>

SimMeth::SimMeth(methName meth, Model* model){
  m_meth = meth;
  m_model = model;  
}


SimMeth::SimMeth(methName meth){
  m_meth = meth;
}

SimMeth::~SimMeth(){
  delete m_model;
}

void SimMeth::setModel(Model* model){
  m_model = model;
}


methName SimMeth::getMeth() const{
  return m_meth;
}


int SimMeth::setMeth(methName meth){
  m_meth = meth;
  return 0;
}


void SimMeth::initMeth(){
}


void SimMeth::simulate(double globalTime, double DT){
  switch(m_meth){
    /*case rk4:
    break;
  case ac:
    break;
  case euler:
    euler(globalTime, DT);
    break;*/
  default:
    euler(globalTime, DT);
    break;
  }
}


/*void SimMeth::simRk4(double x, double h){
  int i;
  double xh,hh,h6;
  ParamVect *dyma,*dyt,*yt;
  int n;
  ParamVect *y,*dydx,*yout;
  

  n = modelSim->getNumStates();
  y = yout = modelSim->getStates();
  dydx = modelSim->getDerivStates();

  dyma = new ParamVect (n,0);
  dyt = new ParamVect (n,0);
  yt = new ParamVect (n,0);

  hh = h*0.5;
  h6 = h/6.0;
  xh = x+hh;
                                                                                                                
  for (i=0;i<n;i++)
    yt->at(i) = y->at(i)+hh*dydx->at(i);

  modelSim->ModelDerivatives (xh,yt,dyt);

  for (i=0;i<n;i++)
    yt->at(i) = y->at(i)+hh*dyt->at(i);

  modelSim->ModelDerivatives(xh,yt,dyma);

  for (i=0;i<n;i++) {
    yt->at(i) = y->at(i)+h*dyma->at(i);
    dyma->at(i) += dyt->at(i);
  }

  modelSim->ModelDerivatives(x+h,yt,dyt);

  for (i=0;i<n;i++)
    yout->at(i) = y->at(i)+h6*(dydx->at(i)+dyt->at(i)+2.0*dyma->at(i));

  delete(yt);
  delete(dyt);
  delete(dyma);
  }*/

void SimMeth::euler(double x, double h){
  /* int i;
  int n;
  ParamVect *y, *dydx, *yout;

  n = m_model->getNumSt();

  y = yout = m_model->getSt();
  dydx = m_model->getDerivSt();

  //Evaluacion de las derivadas en el punto actual AIHR
  m_model->ModelDeriv(x, y, dydx);

  for(i=0;i<n;i++){
    yout->at(i)=y->at(i)+ h*(dydx->at(i));
    }*/
}




