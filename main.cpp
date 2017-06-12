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

#include "coupler.hpp"
#include "diffusion3D.hpp"
#include "Gen3DProstTissue.hpp"
#include "Model.hpp"
#include "prostateCell.hpp"
#include "ProstCell.hpp"
#include "rootSimulator.hpp"
#include "Simulator.hpp"
#include "Treatment.hpp"

using namespace std;

int main(int argc, char *argv[]){
  Model *coupler, *model1, *model2;
  int nrow, ncol, nlayer;
  double doubTime;
  double apopDeadTime, apopProb, necDeadTime;
  double KConso, Vmax;
  string nFInPO2, nFInTum, nFInVes, nFTissueDim;
  vector<double> alpha(7, 0.0);
  vector<double> beta(7, 0.0);
  vector<double> cycDistrib(4, 0.0);
  vector<double> cycDur(4, 0.0);
  Treatment *treatment;
  RootSimulator *sim;
  double DT1, DT2, simTime;

  if(argc != 29){
    cout << "Incorrect number of parameters. 28 parameters expected"
	 << endl;
    cout << "doubTime G1Dur SDur G2Dur MDur G1Distrib SDistrib" <<
      "G2Distrib MDistrib alphaAlive alphaTumG1 alphaTumS" <<
      "alphaTumG2 alphaTumM alphaDead alphaVes betaAlive" <<
      "betaTumG1 betaTumS betaTumG2 betaTumM betaDead betaVes" <<
      "apopProb apopDeadTime necDeadTime Vmax KConso" << endl;
    return EXIT_FAILURE;
    }
  
  
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
  doubTime = atof(argv[1]);

  cycDur[0] = atof(argv[2]);
  cycDur[1] = atof(argv[3]);
  cycDur[2] = atof(argv[4]);
  cycDur[3] = atof(argv[5]);

  cycDistrib[0] = atof(argv[6]);
  cycDistrib[1] = atof(argv[7]);
  cycDistrib[2] = atof(argv[8]);
  cycDistrib[3] = atof(argv[9]);
  
  alpha[0] = atof(argv[10]);
  alpha[1] = atof(argv[11]);
  alpha[2] = atof(argv[12]);
  alpha[3] = atof(argv[13]);
  alpha[4] = atof(argv[14]);
  alpha[5] = atof(argv[15]);
  alpha[6] = atof(argv[16]);

  beta[0] = atof(argv[17]);
  beta[1] = atof(argv[18]);
  beta[2] = atof(argv[19]);
  beta[3] = atof(argv[20]);
  beta[4] = atof(argv[21]);
  beta[5] = atof(argv[22]);
  beta[6] = atof(argv[23]);

  apopProb     = atof(argv[24]);
  apopDeadTime = atof(argv[25]);
  necDeadTime  = atof(argv[26]);

  Vmax   = atof(argv[27]);
  KConso = atof(argv[28]);
  
  treatment = new Treatment();
  cout<<treatment<<endl;
  //treatment = 0;
  model1 = new Gen3DProstTissue(nrow, ncol, nlayer, nFInPO2,
				nFInTum, nFInVes, doubTime, cycDur,
				cycDistrib, apopDeadTime,
				necDeadTime, apopProb, alpha, beta,
				treatment);
  model2 = new diffusion3D(nrow, ncol, nlayer, nFInVes, Vmax,
			   KConso);
  coupler = new Coupler(model1, model2);
  DT1 = 1; //h;
  DT2 = 1; //s;
  sim = new RootSimulator(coupler, DT1, DT2);
  simTime = treatment->getDuration();
  //simTime = 2016;
  sim->initSim();
  cout << "test2" << endl;
  sim->simulate(0.0, simTime);
  return EXIT_SUCCESS;
}
