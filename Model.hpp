/*
 *  Model.hpp
 *  AlfSimLib
 *
 *  Created by Alfredo Hern‡ndez on Fri Jul 26 2002.
 *  Copyright (c) 2001 INSERM. All rights reserved.
 *
 */

#ifndef DEF_MODEL
#define DEF_MODEL

#include <vector>

//define paramters storage
typedef std::vector<double> DoubleVect;

//model types
enum modelType {DEVS,DESS,DTSS};

class Model {
public :
  Model(const modelType type, const int numIn, const int numSt,
	const int numOut, const int numParam, const int numComp);
  virtual ~Model(); 
  virtual int ModelInit(const double DT=0)=0;
  virtual int ModelOut()=0;
  virtual int ModelUpdate(const double currentTime=0,
			  const double DT=0)=0;
  virtual int ModelTerminate();
  virtual int ModelStart();
  DoubleVect *getIn() const;
  DoubleVect *getSt() const;
  DoubleVect *getDerivSt() const;
  DoubleVect *getOut() const;
  DoubleVect *getParam() const;
  modelType getModelType() const;
  std::vector<Model*> *getComp() const;
  int getNumIn() const;
  int getNumSt() const;
  int getNumOut() const;
  int getNumParam() const;
  int getNumComp() const;

protected:
  modelType m_typeModel;
  int m_numIn;
  int m_numSt;
  int m_numOut;
  int m_numParam;
  int m_numComp;
  DoubleVect *m_in;	
  DoubleVect *m_st;   
  DoubleVect *m_derivSt; //cache the value of states during integ
  DoubleVect *m_out; 
  DoubleVect *m_param; 
  std::vector<Model *>  *m_comp;
  Model	*m_parent; 
};

#endif
