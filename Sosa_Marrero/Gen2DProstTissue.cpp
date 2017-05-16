/*
 *  generic2dtissueProstate.cpp
 *
 *  Copyright (c) 2001 INSERM. All rights reserved.
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "Gen2DProstTissue.hpp"

using namespace std;

Gen2DProstTissue::Gen2DProstTissue() :
  Model(DESS, 0, 0, 0, 4, TISSUEROW*TISSUECOL){
  int k(0);
  double input;
  
  //Creation of the cells composing the tissue model
  for(int k=0;k<m_numComp;k++){
    m_comp->at(k) = new ProstCell(this);
    m_comp->at(k)->ModelOut();
    m_numOut += (m_comp->at(k))->getNumOut();
  }
 
  for(int i=0;i<TISSUEROW;i++){
    for(int j=0;j<TISSUECOL;j++){
      m_tissue[i][j] = m_comp->at(k);
      k++;
    }
  }
  m_tumorEdge = new vector<int>(0,0);
  m_deadCells = new vector<int>(0,0);
  srand(time(NULL));

  m_treatment = 0;
}


Gen2DProstTissue::Gen2DProstTissue(string nFInPO2, string nFInTum,
				   string nFInVes,
				   Treatment *treatment) :
  Model(DESS, 0, 0, 0, 4, TISSUEROW*TISSUECOL){
  int k(0);
  double input;
  
  //Creation of the cells composing the tissue model
  for(int k=0;k<m_numComp;k++){
    m_comp->at(k) = new ProstCell(this);
    m_comp->at(k)->ModelOut();
    m_numOut += (m_comp->at(k))->getNumOut();
  }
 
  for(int i=0;i<TISSUEROW;i++){
    for(int j=0;j<TISSUECOL;j++){
      m_tissue[i][j] = m_comp->at(k);
      k++;
    }
  }
  m_tumorEdge = new vector<int>(0,0);
  m_deadCells = new vector<int>(0,0);
  srand(time(NULL));

  //Initialization of the PO2
  k = 0;
  ifstream fInPO2 (nFInPO2.c_str());
  if(fInPO2.is_open()){
    while(fInPO2>>input){
      ((ProstCell *)m_comp->at(k))->setInPO2(input);
      (m_comp->at(k))->ModelUpdate();
      //((ProstCell *)m_comp->at(k))->setInPO2(0.0);
      (m_comp->at(k))->ModelOut();
      k++;
    }
    fInPO2.close();
  }
  else{
    cout<<"An error occurred while opening initial PO2 data file"
	<<endl;
  }

  //Initialization of the tumor cells
  k = 0;
  ifstream fInTum (nFInTum.c_str());
  if(fInTum.is_open()){
    while(fInTum>>input){
      setInTumor(k, input);
      (m_comp->at(k))->ModelUpdate();
      setInTumor(k, 0.0);
      (m_comp->at(k))->ModelOut();
      k++;
    }
    fInTum.close();
  }
  else{
    cout<<"An error occurred while opening initial tumor data file"
	<<endl;
  }

  //Initialization of vessels
  k = 0;
  ifstream fInVes (nFInVes.c_str());
  if(fInVes.is_open()){
    while(fInVes>>input){
      setInVes(k, input);
      (m_comp->at(k))->ModelUpdate();
      setInVes(k, 0.0);
      (m_comp->at(k))->ModelOut();
      k++;
    }
    fInVes.close();
  }
  else{
    cout<<"An error occurred while opening initial vessel data file"
	<<endl;
  }
  m_treatment=treatment;
}


Gen2DProstTissue::~Gen2DProstTissue(){
  delete m_tumorEdge;
  delete m_deadCells;
}


int Gen2DProstTissue::ModelInitSim(double DT){
  PAR_PF = pow(2.0, DT/(((ProstCell *)m_comp->at(0))->
			getDoubTime()));
  PAR_NUM_TUMOR = getNumTumor();
  PAR_RF = 1 - pow(2.0, -DT/(((ProstCell *)m_comp->at(0))->
			     getDeadTime()));
  PAR_NUM_DEAD = getNumDead();
  
  cout<<"Initial number of living cells = "<<getNumAlive()<<endl;
  cout<<"Initial number of tumor cells = "<<getNumTumor()<<endl;
  cout<<"Initial number of vessels = "<<getNumVes()<<endl;
  cout<<"Initial number of dead cells = "<<getNumDead()<<endl;
  cout<<"Initial tumor edge size = "<<m_tumorEdge->size()<<endl;
  
  return 1;
}


int Gen2DProstTissue::ModelOut(){
  for(int k=0;k<m_numComp;k++){
    (m_comp->at(k))->ModelOut();
  }
  return 0;
}


//It does nothing for the moment
int Gen2DProstTissue::ModelStart(){
  for (int k=0;k<m_numComp;k++){
    (m_comp->at(k))->ModelStart();
  }
  return 0;
}


int Gen2DProstTissue::ModelTerminate(){
  for(int k=0;k<m_numComp;k++){
    (m_comp->at(k))->ModelTerminate();
  }
  
  cout<<"Final number of living cells = "<<getNumAlive()<<endl;
  cout<<"Final number of tumor cells = "<<getNumTumor()<<endl;
  cout<<"Final number of vessels = "<<getNumVes()<<endl;
  cout<<"Final number of dead cells = "<<getNumDead()<<endl;
  cout<<"Final tumor edge size = "<<m_tumorEdge->size()<<endl;
  
  return 0;
}


int Gen2DProstTissue::ModelUpdate(double currentTime, double DT){
  int k, l;
  int numTumor, numDead;
  
  PAR_NUM_TUMOR *= PAR_PF;
  numTumor = getNumTumor();

  for(int i=numTumor;i<(int)PAR_NUM_TUMOR;i++){
    if(m_tumorEdge->size()>0){
      l=rand()%m_tumorEdge->size();
      k=m_tumorEdge->at(l);

      setInTumor(k, 1.0);
      (m_comp->at(k))->ModelUpdate(currentTime, DT);
      setInTumor(k, 0.0);
      (m_comp->at(k))->ModelOut();
    }
  }
  
  if(m_treatment!=0){
    if(fmod(currentTime,m_treatment->getInterval())==0){
      int i(currentTime/m_treatment->getInterval());

      if((m_treatment->getSchedule()).at(i)){
	int preTumorSize;
	double m;

	preTumorSize = getNumTumor();

	for(int k=0;k<m_numComp;k++){
	  m=(double)rand()/(double)(RAND_MAX);
	  if(((ProstCell *)m_comp->at(k))->CalcSF()<m){
	    setInDead(k, 1.0);
	    (m_comp->at(k))->ModelUpdate(currentTime, DT);
	    setInDead(k, 0.0);
	    (m_comp->at(k))->ModelOut();
	  }
	}
	PAR_NUM_TUMOR -= preTumorSize - getNumTumor();
	PAR_NUM_DEAD += preTumorSize - getNumTumor();
      }
    }
  
    numDead = getNumDead();
    PAR_NUM_DEAD -= PAR_NUM_DEAD * PAR_RF;
  
    for(int i=(int)PAR_NUM_DEAD+1;i<numDead;i++){
      l=rand()%m_deadCells->size();
      k=m_deadCells->at(l);
      
      setInAlive(k, 1.0);
      (m_comp->at(k))->ModelUpdate(currentTime, DT);
      setInAlive(k, 0.0);
      (m_comp->at(k))->ModelOut();
    }
  }

return 0;
}


void Gen2DProstTissue::AddToDeadCells(int k){
  m_deadCells->push_back(k);
}


void Gen2DProstTissue::AddToEdge(int x, int y){
  int k;
  bool alive(false), notInEdge(true);
 
  k = XYTok(x, y);
  if(k !=-1){
    for(int i=0;i<m_tumorEdge->size();i++){
      if(k==m_tumorEdge->at(i)){
	notInEdge = false;
	i = m_tumorEdge->size();
      }
    }
    
    alive = (((ProstCell *)m_comp->at(k))->getOutState()==1 ||
	     ((ProstCell *)m_comp->at(k))->getInAlive()==1);
    
    if(notInEdge && alive){
      m_tumorEdge->push_back(k);
    }
  }
}


double Gen2DProstTissue::getAlpha() const{
  return PAR_ALPHA;
}


double Gen2DProstTissue::getBeta() const{
  return PAR_BETA;
}


int Gen2DProstTissue::getNumAlive() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getOutState()==1){
      count++;
    }
  }
  return count;
}


int Gen2DProstTissue::getNumDead() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getOutState()==4){
      count++;
    }
  }
  return count;
}


int Gen2DProstTissue::getNumTumor() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getOutState()==2){
      count++;
    }
  }
  return count;
}


int Gen2DProstTissue::getNumVes() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getOutState()==3){
      count++;
    }
  }
  return count;
}


Treatment * Gen2DProstTissue::getTreatment() const{
  return m_treatment;
}


int *Gen2DProstTissue::kToXY(int k) const{
  int *tab= new int[2];
  if(k>-1 && k<TISSUEROW*TISSUECOL){
    tab[0] = k/TISSUECOL;
    tab[1] = k%TISSUECOL;
  }
  else{
    tab[0] = -1;
    tab[1] = -1;
  }
  return tab;
}


void Gen2DProstTissue::setInAlive(int k, double input){
  int x, y;
  int k8[8];
  int *coordxy;
  
  ((ProstCell *)m_comp->at(k))->setInAlive(input);
  if(input){
    RemoveFromDeadCells(k);
    coordxy = kToXY(k);
    x = coordxy[0];
    y = coordxy[1];  
    k8[0] = XYTok(x, y+1);
    k8[1] = XYTok(x, y-1);
    k8[2] = XYTok(x+1, y);
    k8[3] = XYTok(x+1, y+1);
    k8[4] = XYTok(x+1, y-1);
    k8[5] = XYTok(x-1, y);
    k8[6] = XYTok(x-1, y+1);
    k8[7] = XYTok(x-1, y-1);
    
    for(int i=0;i<8;i++){
      if(k8[i]!=-1){
	if(((ProstCell *)m_comp->at(k8[i]))->getOutState()==2){
	  AddToEdge(x,y);
	  i=8;
	}
      }
    }
  }
}


void Gen2DProstTissue::setInDead(int k, double input){
  int x, y;
  int k24[24];
  int *coordxy;
  bool notTumor24[24];
  
  ((ProstCell *)m_comp->at(k))->setInDead(input);
  if (input && ((ProstCell *)m_comp->at(k))->getOutState()==2){
    AddToDeadCells(k);
    coordxy = kToXY(k);
    x = coordxy[0];
    y = coordxy[1];  
    k24[0] = XYTok(x, y+1);
    k24[1] = XYTok(x, y-1);
    k24[2] = XYTok(x+1, y);
    k24[3] = XYTok(x+1, y+1);
    k24[4] = XYTok(x+1, y-1);
    k24[5] = XYTok(x-1, y);
    k24[6] = XYTok(x-1, y+1);
    k24[7] = XYTok(x-1, y-1);

    k24[8] = XYTok(x, y+2);
    k24[9] = XYTok(x, y-2);
    k24[10] = XYTok(x+1, y+2);
    k24[11] = XYTok(x+1, y-2);
    k24[12] = XYTok(x+2, y);
    k24[13] = XYTok(x+2, y+1);
    k24[14] = XYTok(x+2, y+2);
    k24[15] = XYTok(x+2, y-1);
    k24[16] = XYTok(x+2, y-2);
    k24[17] = XYTok(x-1, y+2);
    k24[18] = XYTok(x-1, y-2);
    k24[19] = XYTok(x-2, y);
    k24[20] = XYTok(x-2, y+1);
    k24[21] = XYTok(x-2, y+2);
    k24[22] = XYTok(x-2, y-1);
    k24[23] = XYTok(x-2, y-2);
 
    for(int i=0;i<24;i++){
      if(k24[i]==-1){
	notTumor24[i] = true;
      }
      else{
	notTumor24[i] = ((ProstCell *)m_comp->at(k24[i]))
	  ->getOutState()!=2;
      }
    }

    if(notTumor24[2] && notTumor24[3] && notTumor24[10]
       && notTumor24[8]  && notTumor24[17] && notTumor24[6]
       && notTumor24[5]){
      RemoveFromEdge(k24[0]);
    }

    if(notTumor24[0] && notTumor24[8] && notTumor24[17]
       && notTumor24[21]  && notTumor24[20] && notTumor24[19]
       && notTumor24[5]){
      RemoveFromEdge(k24[6]);
    }

    if(notTumor24[0] && notTumor24[6] && notTumor24[20]
       && notTumor24[19]  && notTumor24[22] && notTumor24[7]
       && notTumor24[1]){
      RemoveFromEdge(k24[5]);
    }

    if(notTumor24[5] && notTumor24[19] && notTumor24[22]
       && notTumor24[23]  && notTumor24[18] && notTumor24[9]
       && notTumor24[1]){
      RemoveFromEdge(k24[7]);
    }

    if(notTumor24[5] && notTumor24[7] && notTumor24[18]
       && notTumor24[9]  && notTumor24[11] && notTumor24[4]
       && notTumor24[2]){
      RemoveFromEdge(k24[1]);
    }

    if(notTumor24[1] && notTumor24[9] && notTumor24[11]
       && notTumor24[16]  && notTumor24[15] && notTumor24[12]
       && notTumor24[2]){
      RemoveFromEdge(k24[4]);
    }
    if(notTumor24[1] && notTumor24[4] && notTumor24[15]
       && notTumor24[12]  && notTumor24[13] && notTumor24[3]
       && notTumor24[0]){
      RemoveFromEdge(k24[2]);
    }

    if(notTumor24[2] && notTumor24[12] && notTumor24[13]
       && notTumor24[14]  && notTumor24[10] && notTumor24[8]
       && notTumor24[0]){
      RemoveFromEdge(k24[3]);
    }
  }
}


void Gen2DProstTissue::setInTumor(int k, double input){
  int x, y;
  int *coordxy;
  
  ((ProstCell *)m_comp->at(k))->setInTumor(input);
  if(input){
    RemoveFromEdge(k);
    coordxy = kToXY(k);
    x = coordxy[0];
    y = coordxy[1];
    AddToEdge(x, y+1);
    AddToEdge(x, y-1);
    AddToEdge(x+1, y);
    AddToEdge(x+1, y+1);
    AddToEdge(x+1, y-1);
    AddToEdge(x-1, y);
    AddToEdge(x-1, y+1);
    AddToEdge(x-1, y-1);
  }
}


void Gen2DProstTissue::setInVes(int k, double input){
  ((ProstCell *)m_comp->at(k))->setInVes(input);
  if(input){
    RemoveFromEdge(k);
  }
}


void Gen2DProstTissue::RemoveFromDeadCells(int k){
  for(int i=0;i<m_deadCells->size();i++){
    if(m_deadCells->at(i)==k){
      for(int ii=i;ii<m_deadCells->size()-1;ii++){
	m_deadCells->at(ii)=m_deadCells->at(ii+1);
      }
      m_deadCells->pop_back();
      i=m_deadCells->size();
    }
  }
}


void Gen2DProstTissue::RemoveFromEdge(int k){
  for(int i=0;i<m_tumorEdge->size();i++){
    if(m_tumorEdge->at(i)==k){
      for(int ii=i;ii<m_tumorEdge->size()-1;ii++){
	m_tumorEdge->at(ii)=m_tumorEdge->at(ii+1);
      }
      m_tumorEdge->pop_back();
      i=m_tumorEdge->size();
    }
  }
}


int Gen2DProstTissue::XYTok(int x, int y) const{
  if(x>-1 && x<TISSUEROW && y>-1 && y<TISSUECOL){
    return x*TISSUECOL + y;
  }
  else{
    return -1;
  }
}






