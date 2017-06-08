/**
 * @file main.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.19.17 
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
  int nrow, ncol, nlayer;
  double apopDeadTime, apopProb, doubTime, necDeadTime;
  string nFInPO2, nFInTum, nFInVes, nFTissueDim;
  vector<double> alpha(7, 0.0);
  vector<double> beta(7, 0.0);
  vector<double> cycDistrib(4, 0.0);
  vector<double> cycDur(4, 0.0);
  Treatment *treatment;

  Simulator *sim;
  double simTime;


  nFTissueDim = "tissueDim.dat";
  nFInTum = "inTum.dat";
  nFInVes = "inVes.dat";
  nFInPO2 = "inPO2.dat";

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
  doubTime = 1008;

  cycDur.at(0) = 0.55;
  cycDur.at(1) = 0.2;
  cycDur.at(2) = 0.15;
  cycDur.at(3) = 0.1;

  cycDistrib.at(0) = 0.6;
  cycDistrib.at(1) = 0.25;
  cycDistrib.at(2) = 0.075;
  cycDistrib.at(3) = 0.075;

  apopProb     = 0.8;
  apopDeadTime = 234;
  necDeadTime  = 468;
  
  alpha.at(1) = 0.158;
  alpha.at(2) = 0.113;
  alpha.at(3) = 0.169;
  alpha.at(4) = 0.189;

  beta.at(1) = 0.051;
  beta.at(2) = 0.037;
  beta.at(3) = 0.055;
  beta.at(4) = 0.061;
  
  treatment = new Treatment();
  cout<<treatment<<endl;
  //treatment = 0;
  model = new Gen3DProstTissue(nrow, ncol, nlayer, nFInPO2, nFInTum,
			       nFInVes, doubTime, cycDur,
			       cycDistrib, apopDeadTime,
			       necDeadTime, apopProb, alpha, beta,
			       treatment);
  sim = new Simulator();
  simTime = treatment->getDuration();
  //simTime = 2016;
  sim->setModel(model);
  sim->start(simTime);
  return EXIT_SUCCESS;
}
