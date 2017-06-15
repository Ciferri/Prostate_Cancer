/**
 * @file treatment.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.19.17 
 */

#include <iostream>
#include <vector>

#include "treatment.hpp"

using namespace std;

Treatment::Treatment(){
  int accSession(0), i(0), numSession;
  
  m_fraction  = 2.0; //Gy
  m_interval  = 24.0; //h
  m_totalDose = 80.0; //Gy

  //Sessions scheduled Mon-Fri
  numSession=(int)(m_totalDose/m_fraction);
  while(accSession<=numSession){
    if((i+1)%7==6||(i+1)%7==0){
      m_schedule.push_back(false);
    }
    else{
      m_schedule.push_back(true);
      accSession++;
    }
    i++;
  }
}


Treatment::Treatment(const double fraction, const double totalDose,
		     const double interval,
		     const vector<bool> schedule){
  m_fraction  = fraction;
  m_interval  = interval;
  m_totalDose = totalDose;
  m_schedule  = schedule;
}


double Treatment::getDuration() const{
  return m_schedule.size() * m_interval;
}


double Treatment::getFraction() const{
  return m_fraction;
}


double Treatment::getInterval() const{
  return m_interval;
}


vector<bool> Treatment::getSchedule() const{
  return m_schedule;
}


double Treatment::getTotalDose() const{
  return m_totalDose;
}


ostream &operator<<(ostream &stream, Treatment *treatment){
  stream << "Total dose = " << treatment->getTotalDose() <<
    " Gy" << endl;
  stream << "Fraction = " << treatment->getFraction() << " Gy" <<
    endl;
  stream << "---------------------------------------------";
  return stream;
}
