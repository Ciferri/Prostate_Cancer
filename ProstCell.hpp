/*
 *  prostCell.hpp
 *
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 * 	03.23.2017
 *
 */

#ifndef DEF_PROSTCELL
#define DEF_PROSTCELL

#include "Gen3DProstTissue.hpp"
#include "Model.hpp"
#include "Treatment.hpp"

//Inputs
#define IN_TUMOR      m_in->at(0)
#define IN_VES        m_in->at(1)
#define IN_DEAD       m_in->at(2)
#define IN_ALIVE      m_in->at(3)
#define IN_PO2        m_in->at(4)

//State variables
#define ST_ALIVE      m_st->at(0)
#define ST_DEAD       m_st->at(1)
#define ST_TUMOR      m_st->at(2)
#define ST_VES        m_st->at(3)

//Outputs
#define OUT_STATE     m_out->at(0)

//Internal Parameters
#define PAR_DOUBTIME  m_param->at(0)
#define PAR_DEADTIME  m_param->at(1)
#define PAR_M         m_param->at(2)
#define PAR_K         m_param->at(3)
#define PAR_PO2       m_param->at(4)
#define PAR_ALPHA     m_param->at(5)
#define PAR_BETA      m_param->at(6)

class ProstCell: public Model{
public :
  ProstCell();
  ProstCell(Model *parent);
  virtual ~ProstCell();
  virtual int ModelInitSim(double DT);
  virtual int ModelOut();
  virtual int ModelStart();
  virtual int ModelTerminate();
  virtual int ModelUpdate(double currentTime, double DT);
  double CalcOER() const;
  double CalcSF() const;
  double getAlive() const;
  double getDead() const;
  double getDeadTime() const;
  double getDoubTime() const;
  double getInAlive() const;
  int getOutState() const;
  double getTumor() const;
  double getVes() const;
  void setInAlive(double input);
  void setInDead(double input);
  void setInPO2(double input);
  void setInTumor(double input);
  void setInVes(double input);
};
#endif
