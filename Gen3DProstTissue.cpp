/*
 *  generic3dtissueProstate.cpp
 *
 *  Copyright (c) 2001 INSERM. All rights reserved.
 *
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
  Model(DESS, 0, 0, 0, 4, TISSUEROW*TISSUECOL*TISSUELAYER){
  int k(0);
  double input;
  
  //Creation of the cells composing the tissue model
  for(int i=0;i<TISSUEROW;i++){
    for(int j=0;j<TISSUECOL;j++){
      for(int l=0;l<TISSUELAYER;l++){
	m_comp->at(k) = new ProstCell(this);
	m_numOut += (m_comp->at(k))->getNumOut();
	m_tissue[i][j][l] = m_comp->at(k);
	k++;
      }  
    }
  }
  m_tumorEdge = new vector<int>(0,0);
  m_deadCells = new vector<int>(0,0);
  m_treatment = 0;
}


Gen3DProstTissue::Gen3DProstTissue(string nFInPO2, string nFInTum,
				   string nFInVes,
				   Treatment *treatment) :
  Model(DESS, 0, 0, 0, 4, TISSUEROW*TISSUECOL*TISSUELAYER){
  int k(0);
  double input;
  
  //Creation of the cells composing the tissue model
  for(int i=0;i<TISSUEROW;i++){
    for(int j=0;j<TISSUECOL;j++){
      for(int l=0;l<TISSUELAYER;l++){
	m_comp->at(k) = new ProstCell(this);
	m_numOut += (m_comp->at(k))->getNumOut();
	m_tissue[i][j][l] = m_comp->at(k);
	k++;
      }
    }
  }
  
  m_tumorEdge = new vector<int>(0,0);
  m_deadCells = new vector<int>(0,0);

  //Initialization of the PO2
  k = 0;
  ifstream fInPO2 (nFInPO2.c_str());
  if(fInPO2.is_open()){
    while(fInPO2>>input && k<m_numComp){
      ((ProstCell *)m_comp->at(k))->setInPO2(input);
      (m_comp->at(k))->ModelUpdate();
      ((ProstCell *)m_comp->at(k))->setInPO2(0.0);
      k++;
    }
    fInPO2.close();
  }
  else{
    cout<<"An error occurred while opening initial PO2 data file"
	<<endl;
  }

  //Initialization of the tumor cells
  k = 0*TISSUEROW*TISSUECOL;
  ifstream fInTum (nFInTum.c_str());
  if(fInTum.is_open()){
    while(fInTum>>input){
      setInTumor(k, input);
      (m_comp->at(k))->ModelUpdate();
      setInTumor(k, 0.0);
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
  ifstream fInVes;
  for(int l=0;l<TISSUELAYER;l++){
    fInVes.open(nFInVes.c_str());
    if(fInVes.is_open()){
      while(fInVes>>input){
	setInVes(k, input);
	(m_comp->at(k))->ModelUpdate();
	setInVes(k, 0.0);
	k++;
      }
      fInVes.close();
    }
    else{
      cout<<"An error occurred while opening initial"<<
	"vessel data file"<<endl;
    }
  } 
  m_treatment=treatment;
}


Gen3DProstTissue::~Gen3DProstTissue(){
  delete m_tumorEdge;
  delete m_deadCells;
}


int Gen3DProstTissue::ModelInitSim(double DT){
  int deadTime, doubTime;

  doubTime = ((ProstCell *)m_comp->at(0))->getDoubTime();
  deadTime = ((ProstCell *)m_comp->at(0))->getDeadTime();
  PAR_PF = pow(2.0, DT/doubTime);
  PAR_NUM_TUMOR = getNumTumor();
  PAR_RF = 1 - pow(2.0, -DT/deadTime);
  PAR_NUM_DEAD = getNumDead();
  srand(time(NULL));
  
  cout<<"Total number of cells = "<<m_numComp<<endl;
  cout<<"Initial number of living cells = "<<getNumAlive()<<endl;
  cout<<"Initial number of tumor cells = "<<getNumTumor()<<endl;
  cout<<"Initial number of vessels = "<<getNumVes()<<endl;
  cout<<"Initial number of dead cells = "<<getNumDead()<<endl;
  cout<<"Initial tumor edge size = "<<m_tumorEdge->size()<<endl;
  cout<<"---------------------------------------------"<<endl;
  
  return 0;
}


int Gen3DProstTissue::ModelOut(){
  for(int k=0;k<m_numComp;k++){
    (m_comp->at(k))->ModelOut();
  }
  return 0;
}


//It does nothing for the moment
int Gen3DProstTissue::ModelStart(){
  for (int k=0;k<m_numComp;k++){
    (m_comp->at(k))->ModelStart();
  }
  return 0;
}


int Gen3DProstTissue::ModelTerminate(){
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


int Gen3DProstTissue::ModelUpdate(double currentTime, double DT){
  int k, m;
  int numTumor, numDead;
  
  PAR_NUM_TUMOR *= PAR_PF;
  numTumor = getNumTumor();

  for(int i=numTumor;i<(int)PAR_NUM_TUMOR;i++){
    if(m_tumorEdge->size()>0){
      m=rand()%m_tumorEdge->size();
      k=m_tumorEdge->at(m);
      setInTumor(k, 1.0);
      (m_comp->at(k))->ModelUpdate();
      setInTumor(k, 0.0);
    }
  }
  
  if(m_treatment!=0){
    if(fmod(currentTime,m_treatment->getInterval())==0){
      int i(currentTime/m_treatment->getInterval());

      if((m_treatment->getSchedule()).at(i)){
	int preTumorSize;
	double n;

	preTumorSize = getNumTumor();
	for(int k=0;k<m_numComp;k++){
	  n=(double)rand()/(double)(RAND_MAX);
	  if(((ProstCell *)m_comp->at(k))->CalcSF()<n){
	    setInDead(k, 1.0);
	    (m_comp->at(k))->ModelUpdate();
	    setInDead(k, 0.0);
	  }
	}
	
	PAR_NUM_TUMOR -= preTumorSize - getNumTumor();
	PAR_NUM_DEAD += preTumorSize - getNumTumor();
      }
    }
  
    numDead = getNumDead();
    PAR_NUM_DEAD -= PAR_NUM_DEAD * PAR_RF;
  
    for(int i=(int)PAR_NUM_DEAD+1;i<numDead;i++){
      m=rand()%m_deadCells->size();
      k=m_deadCells->at(m);
      setInAlive(k, 1.0);
      (m_comp->at(k))->ModelUpdate();
      setInAlive(k, 0.0);
    }
  }
  return 0;
}


void Gen3DProstTissue::AddToDeadCells(int k){
  m_deadCells->push_back(k);
}


void Gen3DProstTissue::AddToEdge(int x, int y, int z){
  int k;
  bool alive(false), notInEdge(true);
 
  k = XYZTok(x, y, z);
  if(k !=-1){
    for(int i=0;i<m_tumorEdge->size();i++){
      if(k==m_tumorEdge->at(i)){
	notInEdge = false;
	i = m_tumorEdge->size();
      }
    }
    
    alive = (((ProstCell *)m_comp->at(k))->getAlive() ||
	     ((ProstCell *)m_comp->at(k))->getInAlive()==1);
    
    if(notInEdge && alive){
      m_tumorEdge->push_back(k);
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
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getAlive()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumDead() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getDead()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumTumor() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getTumor()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumVes() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getVes()){
      count++;
    }
  }
  return count;
}


Treatment * Gen3DProstTissue::getTreatment() const{
  return m_treatment;
}


int *Gen3DProstTissue::kToXYZ(int k) const{
  int *tab= new int[3];
  if(k>-1 && k<TISSUEROW*TISSUECOL*TISSUELAYER){
    tab[0] = (k%(TISSUEROW*TISSUECOL))/TISSUECOL;
    tab[1] = (k%(TISSUEROW*TISSUECOL))%TISSUECOL;
    tab[2] = k/(TISSUEROW*TISSUECOL);
  }
  else{
    tab[0] = -1;
    tab[1] = -1;
    tab[2] = -1;
  }
  return tab;
}


void Gen3DProstTissue::setInAlive(int k, double input){
  int x, y, z;
  int m;
  int *coordxyz;
  
  ((ProstCell *)m_comp->at(k))->setInAlive(input);
  
  if(input){
    RemoveFromDeadCells(k);
    coordxyz = kToXYZ(k);
    x = coordxyz[0];
    y = coordxyz[1];
    z = coordxyz[2];
    for(int i=-1;i<2;i++){
      for(int j=-1;j<2;j++){
	for(int l=-1;l<2;l++){
	  if(i!=0 || j!=0 || l!=0){
	    m = XYZTok(x+i, y+j, z+l);
	    if(m!=-1){
	      if(((ProstCell *)m_comp->at(m))->getTumor()){
		AddToEdge(x, y, z);
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


void Gen3DProstTissue::setInDead(int k, double input){
  int x, xx, y, yy, z, zz;
  int m, n;
  int *coordxyz;
  bool nTiNCSDC; //not Tumor in Neighborhood of the Cells
  //Surrounding the Dead Cell
  
  ((ProstCell *)m_comp->at(k))->setInDead(input);
  
  if(input){
    AddToDeadCells(k);
    coordxyz = kToXYZ(k);
    x = coordxyz[0];
    y = coordxyz[1];
    z = coordxyz[2];
    for(int i=-1;i<2;i++){
      for(int j=-1;j<2;j++){
	for(int l=-1;l<2;l++){
	  if(i!=0 || j!=0 || l!=0){
	    m = XYZTok(x+i, y+j, z+l);
	    if(m!=-1){
	      nTiNCSDC = true;
	      xx = x+i;
	      yy = y+j;
	      zz = z+l;
	      for(int ii=-1;ii<2;ii++){
		for(int jj=-1;jj<2;jj++){
		  for(int ll=-1;ll<2;ll++){
		    if((ii!=0 || jj!=0 || ll!=0) &&
		       (ii != -i || jj !=-j || ll != -l)){
		      n = XYZTok(xx+ii, yy+jj, zz+ll);
		      if(n!=-1){
			nTiNCSDC = nTiNCSDC &&
			  !(((ProstCell *)m_comp->at(n))
			  ->getTumor());		       
		      }
		    }
		  }
		}  
	      }
	      if(nTiNCSDC){
		RemoveFromEdge(m);
	      }
	    }
	  }
	}
      }
    }
  }
}


void Gen3DProstTissue::setInTumor(int k, double input){
  int x, y, z;
  int *coordxyz;
  
  ((ProstCell *)m_comp->at(k))->setInTumor(input);
  
  if(input){
    RemoveFromEdge(k);
    coordxyz = kToXYZ(k);
    x = coordxyz[0];
    y = coordxyz[1];
    z = coordxyz[2];
    for(int i=-1;i<2;i++){
      for(int j=-1;j<2;j++){
	for(int l=-1;l<2;l++){
	  if(i != 0 || j != 0 || l != 0){
	    AddToEdge(x+i, y+j, z+l);
	  }
	}
      }
    }
  }
}


void Gen3DProstTissue::setInVes(int k, double input){
  ((ProstCell *)m_comp->at(k))->setInVes(input);
  
  if(input){
    RemoveFromEdge(k);
  }
}


void Gen3DProstTissue::RemoveFromDeadCells(int k){
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


void Gen3DProstTissue::RemoveFromEdge(int k){
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


int Gen3DProstTissue::XYZTok(int x, int y, int z) const{
  if(x>-1 && x<TISSUEROW && y>-1 && y<TISSUECOL && z>-1 &&
     z<TISSUELAYER){
    return x*TISSUECOL + y + z*TISSUEROW*TISSUECOL;
  }
  else{
    return -1;
  }
}






