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

#define TISSUEROW 57
#define TISSUECOL 97
#define TISSUELAYER 1

//Internal Parameters
#define PAR_INIT_NUM_TUMOR m_param->at(0)
#define PAR_NUM_DEAD       m_param->at(1)
#define PAR_NUM_SESSION    m_param->at(2)
#define PAR_NUM_TUMOR      m_param->at(3)
#define PAR_PF             m_param->at(4)
#define PAR_RF             m_param->at(5)

class Gen3DProstTissue : public Model {
 public:
  Gen3DProstTissue();
  Gen3DProstTissue(const std::string nFInPO2,
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
  void addToDeadCells(const int k);
  void addToEdge(const int x, const int y, const int z);
  double getAlpha() const;
  double getBeta() const;
  int getNumAlive() const;
  int getNumDead() const;
  int getNumTumor() const;
  int getNumVes() const;
  Treatment *getTreatment() const;
  int *kToXyz(const int k) const;
  void setInAlive(const int k, const double input);
  void setInDead(const int k, const double input);
  void setInTumor(const int k, const double input);
  void setInVes(const int k, const double input);
  void removeFromDeadCells(const int k);
  void removeFromEdge(const int k);
  int xyzTok(const int x, const int y, const int z) const;

protected:
  int m_flag;
  Model *m_tissue[TISSUEROW][TISSUECOL][TISSUELAYER];
  Treatment *m_treatment;
  std::vector<int> *m_deadCells;
  std::vector<int> *m_tumorEdge;
};

#endif
