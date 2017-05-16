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
#include "Model.hpp"
#include "Simulator.hpp"
#include "Gen2DProstTissue.hpp"
#include "ProstCell.hpp"
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
  //treatment = 0;
  model = new Gen2DProstTissue(nFInPO2, nFInTum, nFInVes,
			       treatment);
  sim = new Simulator();
  simTime = treatment->getDuration();
  //simTime = 3456000; //s (40 days)
  sim->setModel(model);
  sim->start(simTime);
  return EXIT_SUCCESS;
}
