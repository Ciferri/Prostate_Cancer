/*
 *  main.cpp
 *
 *  Created by Alfredo Hernández on Fri Jul 26 2002.
 *  Copyright (c) 2001 INSERM. All rights reserved.
 *
 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include <string>

#include "Gen3DProstTissue.hpp"
#include "Model.hpp"
#include "ProstCell.hpp"
#include "Simulator.hpp"
#include "Treatment.hpp"

using namespace std;

int main(int argc, char *argv[]){
  int nrow, ncol, nlayer;
  double simTime;
  string nFInPO2, nFInTum, nFInVes, nFTissueDim;
  Model *model;
  Simulator *sim;
  Treatment *treatment;

  nFTissueDim = "tissueDim.dat";
  nFInTum = "inTum.dat";
  nFInVes = "inVes.dat";
  nFInPO2 = "inPO2.dat";

  treatment = new Treatment();
  cout<<treatment<<endl;
  //treatment = 0;
  ifstream fTissueDim(nFTissueDim.c_str());
  if(fTissueDim.is_open()){
    fTissueDim >> nrow;
    fTissueDim >> ncol;
    fTissueDim >> nlayer;
    fTissueDim.close();
    cout << "Tissue dimensions: "  << endl;
    cout << nrow << " rows" << endl;
    cout << ncol << " columns" << endl;
    cout << nlayer << " layers" << endl;
    cout << "---------------------------------------------" << endl;
  }
  else{
    cout << "An error occurred while opening tissue dimensions file"
	 << endl;
  }
  model = new Gen3DProstTissue(nrow, ncol, nlayer, nFInPO2, nFInTum,
			       nFInVes, treatment);
  sim = new Simulator();
  simTime = treatment->getDuration();
  //simTime = 2016;
  sim->setModel(model);
  sim->start(simTime);
  return EXIT_SUCCESS;
}
