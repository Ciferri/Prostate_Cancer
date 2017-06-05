/**
 * @file Gen3DProstTissue.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Alfredo Hernandez
 * @date 05.19.17 
 */

#ifndef DEF_Gen3DProstTissue
#define DEF_Gen3DProstTissue

#include <vector>
#include <string>

#include "Model.hpp"
#include "ProstCell.hpp"
#include "Treatment.hpp"

//Internal parameters
#define PAR_INIT_NUM_TUM m_param->at(0)

class Gen3DProstTissue : public Model {
public:
  Gen3DProstTissue(const int nrow, const int ncol,
		   const int nlayer);
  Gen3DProstTissue(const int nrow, const int ncol, const int nlayer,
		   const std::string nFInPO2,
		   const std::string nFInTum,
		   const std::string nFInVes,
		   Treatment *const treatment);
  ~Gen3DProstTissue();
  virtual int calcModelOut();
  virtual int initModel(const double DT);
  virtual int startModel();
  virtual int terminateModel();
  virtual int updateModel(const double currentTime,
			  const double DT);
  int getNumAlive() const;
  int getNumDead() const;
  int getNumG1() const;
  int getNumG2() const;
  int getNumM() const;
  int getNumS() const;
  int getNumTum() const;
  int getNumVes() const;
  Treatment *getTreatment() const;
  void printNeededDose() const;
  
protected:
  int m_ncol, m_nlayer, m_nrow;
  std::vector<Model ***> m_map;
  int m_flag;
  Treatment *m_treatment;
};

#endif
