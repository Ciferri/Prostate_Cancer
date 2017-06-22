/**
 * @file diffusion3D.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#ifndef diffusion3D_H
#define diffusion3D_H

//#include "Model.hpp"
//#include "prostateCell.hpp"
#include "prostateCell.hpp"
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


#define FILTERSIZE 3
#define INIT_VASCULAR_PO2 42.0

using namespace std;

//typedef vector<string> TwoDTissue;
//typedef vector< vector< vector<Model *> > > Tissue;

class diffusion3D : public Model {

	public:
	//static const int nrow =37, ncol=55, nlayer=3;
	// constructor / destructor
	diffusion3D(const int nrow, const int ncol, const int nlayer, const string nFInVes
							,double Vmax, double K_conso);
	~diffusion3D();

	//variables
	//double d, Rc, Rm, Cm;
  double Tm, Lm;
  double coupParam;
  //double count1;
  double m_Vmax;		//
  double m_K_conso;		//
  double coeff_diff;
  double size;
  double time_step_s;
  double mat_diff;
  double pO2;
  double conso_pO2;	//Consomation of pO2
  double delta;

  double filterIn[FILTERSIZE][FILTERSIZE][FILTERSIZE];
  double filterOut[FILTERSIZE][FILTERSIZE][FILTERSIZE];



	// methods
  virtual int initModel(const double DT);
  virtual int initModel(int x,int y,int z,double pO2);
  //virtual int initModel(double *tab,int size_x, int size_y);
  virtual int calcModelOut();
  virtual int updateModel(double time, const double DT);
  virtual int terminateModel();
  virtual int Coord_XY_to_K(int x, int y);
  //virtual int *Coord_K_to_XY(int k);

	int m_ncol, m_nlayer, m_nrow;
	vector<double> vectorTissue;
	vector<Model *> Tissue;

};


#endif
