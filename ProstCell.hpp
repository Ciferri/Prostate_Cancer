/**
 * @file ProstCell.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.19.17 
 */

#ifndef DEF_PROSTCELL
#define DEF_PROSTCELL

#include "Gen3DProstTissue.hpp"
#include "Model.hpp"
#include "Treatment.hpp"

//Inputs
#define IN_ALIVE     m_in->at(0)
#define IN_TUM       m_in->at(1)
#define IN_DEAD      m_in->at(2)
#define IN_VES       m_in->at(3)
#define IN_PO2       m_in->at(4)
#define IN_TIMER     m_in->at(5)

//State variables
#define ST_ALIVE     m_st->at(0)
#define ST_TUM       m_st->at(1)
#define ST_DEAD      m_st->at(2)
#define ST_VES       m_st->at(3)
#define ST_G1        m_st->at(4)
#define ST_S         m_st->at(5)
#define ST_G2        m_st->at(6)
#define ST_M         m_st->at(7)

//Outputs
#define OUT_STATE    m_out->at(0)

//Internal Parameters
#define PAR_TIMER    m_param->at(0)
#define PAR_DOUBTIME m_param->at(1)
#define PAR_LIMG1S   m_param->at(2)
#define PAR_LIMSG2   m_param->at(3)
#define PAR_LIMG2M   m_param->at(4)
#define PAR_DEADTIME m_param->at(5)
#define PAR_M        m_param->at(6)
#define PAR_K        m_param->at(7)
#define PAR_PO2      m_param->at(8)
#define PAR_ALPHA    m_param->at(9)
#define PAR_BETA     m_param->at(10)
#define PAR_ACC_DOSE m_param->at(11)

class ProstCell: public Model{
public :
  ProstCell();
  ProstCell(Model *const parent);
  virtual ~ProstCell();
  virtual int calcModelOut();
  virtual int initModel(const double DT);
  virtual int startModel();
  virtual int terminateModel();
  virtual int updateModel(const double currentTime,
			  const double DT);
  void addToEdge(ProstCell *const cell);
  void calcDeadCellsResor(const double DT);
  double calcOER() const;
  void calcRespToIrr();
  double calcRF (const double DT) const;
  double calcSF() const;
  void calcTumGrowth();
  double getAlive() const;
  double getAccDose() const;
  double getDead() const;
  double getDeadTime() const;
  double getDoubTime() const;
  std::vector<ProstCell *> *getEdge() const;
  double getG1() const;
  double getG2() const;
  double getM() const;
  int getOutState() const;
  double getS() const;
  double getTum() const;
  double getVes() const;
  ProstCell *searchSpace() const;
  void setInAlive(const double input);
  void setInDead(const double input);
  void setInPO2(const double input);
  void setInTimer(const double input);
  void setInTum(const double input);
  void setInVes(const double input);

protected:
  std::vector<ProstCell *> *m_edge;
  Treatment *m_treatment;
};
#endif
