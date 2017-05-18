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
  Model *model;
  Simulator *sim;
  Treatment *treatment;
  double simTime;
  string nFInPO2, nFInTum, nFInVes;
 
  nFInTum = "inTum.dat";
  nFInVes = "inVes.dat";
  nFInPO2 = "inPO2.dat";
  treatment = new Treatment();
  cout<<treatment<<endl;
  //treatment = 0;
  model = new Gen3DProstTissue(nFInPO2, nFInTum, nFInVes,
			       treatment);
  sim = new Simulator();
  simTime = treatment->getDuration();
  //simTime = 2016;
  sim->setModel(model);
  sim->start(simTime);
  return EXIT_SUCCESS;
}
