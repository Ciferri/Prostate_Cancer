/**
 * @file coupler.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17 
 */

#ifndef DEF_COUPLER
#define DEF_COUPLER

#include "prostTissue.hpp"
#include "model.hpp"

//#define MODEL1 m_comp->at(0)
//#define MODEL2 m_comp->at(1)

class Coupler : public Model{
public :
  Coupler(Model *model1, Model *model2);
  virtual ~Coupler();
  virtual int calcModelOut();
  virtual int initModel(const double DT);
  virtual int startModel();
  virtual int terminateModel();
  virtual int updateModel(const double currentTime,
			  const double DT);
  Model *getModel1() const;
  Model *getModel2() const;
};
#endif
