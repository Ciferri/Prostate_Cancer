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
#define TISSUECOL 94
#define TISSUELAYER 1

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
  int getNumAlive() const;
  int getNumDead() const;
  int getNumG1() const;
  int getNumG2() const;
  int getNumM() const;
  int getNumS() const;
  int getNumTumor() const;
  int getNumVes() const;
  Treatment *getTreatment() const;

protected:
  Model *m_tissue[TISSUEROW][TISSUECOL][TISSUELAYER];
  Treatment *m_treatment;
};

#endif
