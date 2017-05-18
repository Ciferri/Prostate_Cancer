/***************************************************************************
                          gen2dtissueProst.hpp  -  description
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

//Internal Parameters
#define PAR_PF        m_param->at(0)
#define PAR_NUM_TUMOR m_param->at(1)
#define PAR_RF        m_param->at(2)
#define PAR_NUM_DEAD  m_param->at(3)

class Gen3DProstTissue : public Model {
 public:
  Gen3DProstTissue();
  Gen3DProstTissue(std::string nFInPO2, std::string nFInTum,
		   std::string nFInVes, Treatment *treatment);
  ~Gen3DProstTissue();
  virtual int ModelInitSim(double DT);
  virtual int ModelOut();
  virtual int ModelStart();
  virtual int ModelTerminate();
  virtual int ModelUpdate(double currentTime, double DT);
  void AddToDeadCells(int k);
  void AddToEdge(int x, int y, int z);
  double getAlpha() const;
  double getBeta() const;
  int getNumAlive() const;
  int getNumDead() const;
  int getNumTumor() const;
  int getNumVes() const;
  Treatment *getTreatment() const;
  int *kToXYZ(int k) const;
  void setInAlive(int k, double input);
  void setInDead(int k, double input);
  void setInTumor(int k, double input);
  void setInVes(int k, double input);
  void RemoveFromDeadCells(int k);
  void RemoveFromEdge(int k);
  int XYZTok(int x, int y, int z) const;

protected:
  Model *m_tissue[TISSUEROW][TISSUECOL][TISSUELAYER];
  Treatment *m_treatment;
  std::vector<int> *m_deadCells;
  std::vector<int> *m_tumorEdge;
};

#endif
