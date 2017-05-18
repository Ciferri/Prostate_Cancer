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
  Model(DESS, 0, 0, 0, 0, TISSUEROW*TISSUECOL*TISSUELAYER){
  int k;
  double input;
  
  //Creation of the cells composing the tissue model
  for(int k=0;k<m_numComp;k++){
    m_comp->at(k) = new ProstCell(this);
    m_comp->at(k)->ModelOut();
    m_numOut += (m_comp->at(k))->getNumOut();
  }

  k = 0;
  for(int i=0;i<TISSUEROW;i++){
    for(int j=0;j<TISSUECOL;j++){
      for(int l=0;l<TISSUELAYER;l++){
	m_tissue[i][j][l] = m_comp->at(k);
	k++;
	
      }
    }
  }
  m_treatment = 0;
}


Gen3DProstTissue::Gen3DProstTissue(string nFInPO2, string nFInTum,
				   string nFInVes,
				   Treatment *treatment) :
  Model(DESS, 0, 0, 0, 0, TISSUEROW*TISSUECOL*TISSUELAYER){
  int k;
  int selInitPhase;
  double doubTime, input;

  m_treatment=treatment;
  //Creation of the cells composing the tissue model
  for(k=0;k<m_numComp;k++){
    m_comp->at(k) = new ProstCell(this);
    m_comp->at(k)->ModelOut();
    m_numOut += (m_comp->at(k))->getNumOut();
  }
  k = 0;
  for(int i=0;i<TISSUEROW;i++){
    for(int j=0;j<TISSUECOL;j++){
      for(int l=0;l<TISSUELAYER;l++){
	m_tissue[i][j][l] = m_comp->at(k);
	k++;
      }
    }
  }
  k = 0;
  for(int i=0;i<TISSUEROW;i++){
    for(int j=0;j<TISSUECOL;j++){
      for(int l=0;l<TISSUELAYER;l++){
	for(int ii=-1;ii<=1;ii++){
	  for(int jj=-1;jj<=1;jj++){
	    for(int ll=-1;ll<=1;ll++){
	      if(ii!=0 || jj!=0 || ll!=0){
		if(i+ii>=0 && j+jj>=0 && l+ll>=0 && i+ii<TISSUEROW
		   && j+jj<TISSUECOL && l+ll<TISSUELAYER){
		  ((ProstCell *)m_tissue[i][j][l])
		    ->AddToEdge(((ProstCell *)m_tissue[i+ii]
				 [j+jj][l+ll]));
		}
	      }
	    }
	  }
	}
	k++;
      }
    }
  }
  
  //Initialization of the PO2
  k = 0;
  ifstream fInPO2 (nFInPO2.c_str());
  if(fInPO2.is_open()){
    while(fInPO2>>input){
      ((ProstCell *)m_comp->at(k))->setInPO2(input);
      (m_comp->at(k))->ModelInit();
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
      ((ProstCell *)m_comp->at(k))->setInTumor(input);
      (m_comp->at(k))->ModelInit();
      ((ProstCell *)m_comp->at(k))->setInTumor(0.0);
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
	((ProstCell *)m_comp->at(k))->setInVes(input);
	(m_comp->at(k))->ModelInit();
	((ProstCell *)m_comp->at(k))->setInVes(0.0);
	k++;
      }
      fInVes.close();
    }
    else{
      cout<<"An error occurred while opening initial"<<
	"vessel data file"<<endl;
    }
  }
  
  //Initialization of timer
  srand(time(NULL));
  for(k=0;k<m_numComp;k++){
    doubTime = ((ProstCell *)m_comp->at(k))->getDoubTime();
    selInitPhase = rand()%200;
    if(selInitPhase<120){
      input = rand()%(int)(0.55*doubTime);
    }
    else if(selInitPhase<170){
      input = (int)(0.55*doubTime) + rand()%(int)(0.2*doubTime);
    }
    else if(selInitPhase<185){
      input = (int)(0.75*doubTime) + rand()%(int)(0.15*doubTime);
    }
    else{
      input = (int)(0.9*doubTime) + rand()%(int)(0.1*doubTime);
    }
    ((ProstCell *)m_comp->at(k))->setInTimer(input);
    (m_comp->at(k))->ModelInit();
  }
}


Gen3DProstTissue::~Gen3DProstTissue(){
}


int Gen3DProstTissue::ModelInit(double DT){
  cout<<"Total number of cells = "<<m_numComp<<endl;
  cout<<"Initial number of cells at G1 = "<<getNumG1()<<endl;
  cout<<"Initial number of cells at S = "<<getNumS()<<endl;
  cout<<"Initial number of cells at G2 = "<<getNumG2()<<endl;
  cout<<"Initial number of cells at M = "<<getNumM()<<endl;
  cout<<"Initial number of living cells = "<<getNumAlive()<<endl;
  cout<<"Initial number of tumor cells = "<<getNumTumor()<<endl;
  cout<<"Initial number of vessels = "<<getNumVes()<<endl;
  cout<<"Initial number of dead cells = "<<getNumDead()<<endl;
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
 
  cout<<"Final number of cells at G1 = "<<getNumG1()<<endl;
  cout<<"Final number of cells at S = "<<getNumS()<<endl;
  cout<<"Final number of cells at G2 = "<<getNumG2()<<endl;
  cout<<"Final number of cells at M = "<<getNumM()<<endl;
  cout<<"Final number of living cells = "<<getNumAlive()<<endl;
  cout<<"Final number of tumor cells = "<<getNumTumor()<<endl;
  cout<<"Final number of vessels = "<<getNumVes()<<endl;
  cout<<"Final number of dead cells = "<<getNumDead()<<endl;
  
  return 0;
}


int Gen3DProstTissue::ModelUpdate(double currentTime, double DT){
  for(int k=0;k<m_numComp;k++){
    (m_comp->at(k))->ModelUpdate(currentTime, DT);
    (m_comp->at(k))->ModelOut();
  }
  return 0;
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


int Gen3DProstTissue::getNumG1() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getG1()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumG2() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getG2()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumM() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getM()){
      count++;
    }
  }
  return count;
}


int Gen3DProstTissue::getNumS() const{
  int count(0);
  for(int k=0;k<m_numComp;k++){
    if(((ProstCell *)m_comp->at(k))->getS()){
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


Treatment *Gen3DProstTissue::getTreatment() const{
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


int Gen3DProstTissue::XYZTok(int x, int y, int z) const{
  if(x>-1 && x<TISSUEROW && y>-1 && y<TISSUECOL && z>-1 &&
     z<TISSUELAYER){
    return x*TISSUECOL + y + z*TISSUEROW*TISSUECOL;
  }
  else{
    return -1;
  }
}






