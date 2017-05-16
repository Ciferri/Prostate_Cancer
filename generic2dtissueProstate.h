/***************************************************************************
                          generic2dtissueProstate.h  -  description
                             -------------------
    begin                : mer jan 28 2004
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef generic2dtissueProstate_H
#define generic2dtissueProstate_H

#include "Model.h"
#include "prostateCell.h"
#include <iostream>
#include "math.h"
#include <cmath>
#include <fstream>
#include <sstream>
//~ #include <itkImage.h>
//~ #include <itkTimeProbe.h>
//~ #include <itkRGBAPixel.h>
//~ #include <itkImageFileReader.h>
//~ #include <itkImageRegionIterator.h>

#define TISSUESIZE_COL 98 
#define TISSUESIZE_ROW 63 
#define FILTERSIZE 3
#define INIT_VASCULAR_PO2 42.0

using namespace std;

typedef vector< vector<Model *> > TwoDTissue;

class generic2dtissueProstate : public Model {
public:

//variables
//double d, Rc, Rm, Cm;
  double Tm, Lm;
  double coupParam;
  //double count1;
  double Vmax;		//
  double K_conso;		//
  double coeff_diff;
  double size;
  double time_step_s;
  double mat_diff;
  double pO2;
  double conso_pO2;	//Consomation of pO2
  double delta;

  double filterIn[FILTERSIZE][FILTERSIZE];  
  double filterOut[FILTERSIZE][FILTERSIZE];

  Model *tissue[TISSUESIZE_COL][TISSUESIZE_ROW];
  double *gradTissue[TISSUESIZE_COL][TISSUESIZE_ROW];
  double gradLin[TISSUESIZE_COL*TISSUESIZE_ROW][1];
  double resDiffusionIn[TISSUESIZE_COL][TISSUESIZE_ROW];  
  double resDiffusionOut[TISSUESIZE_COL][TISSUESIZE_ROW];
  double initImageTissue[TISSUESIZE_COL][TISSUESIZE_ROW];
  double Conso_PO2[TISSUESIZE_COL][TISSUESIZE_ROW];
  double vectTissue[TISSUESIZE_COL*TISSUESIZE_ROW];

// constructor / destructor
	generic2dtissueProstate();
	~generic2dtissueProstate();

// methods
  virtual int ModelInitSim();
  virtual int ModelInitSim(int x,int y,double pO2);
  //virtual int ModelInitSim(double *tab,int size_x, int size_y);
  virtual int NextSampleHit();
  virtual int ModelOutputs();
  virtual int ModelUpdate(double time);
  virtual int ModelDerivatives(double time, ParamVect *variables, ParamVect *derivadas);
  virtual int ModelTerminate();
  virtual int ModelStart();
  virtual int Coord_XY_to_K(int x, int y);
  //virtual int *Coord_K_to_XY(int k);

};


#endif
