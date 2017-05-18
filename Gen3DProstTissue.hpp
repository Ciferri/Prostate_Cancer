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

class Gen3DProstTissue : public Model {
 public:
  Gen3DProstTissue();
  Gen3DProstTissue(const std::string nFInPO2,
		   const std::string nFInTum,
		   const std::string nFInVes,
		   Treatment *const treatment);
  ~Gen3DProstTissue();
  virtual int ModelInit(const double DT);
  virtual int ModelOut();
  virtual int ModelStart();
  virtual int ModelTerminate();
  virtual int ModelUpdate(const double currentTime,
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
