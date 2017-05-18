#include <iostream>
#include <vector>
#include "Treatment.hpp"

using namespace std;

Treatment::Treatment(){
  m_fraction = 2.0; //Gy
  m_interval = 24; //s
  m_totalDose = 60.0; //Gy

  //Sessions scheduled Mon-Fri
  int accSession(0), i(0), numSession;
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


Treatment::Treatment(double fraction, double totalDose,
		     double interval, vector<bool> schedule){
  m_fraction = fraction;
  m_interval = interval;
  m_totalDose = totalDose;
  m_schedule = schedule;
}


double Treatment::getDuration() const{
  return m_schedule.size()*m_interval;
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
  stream<<"Total dose = "<<treatment->getTotalDose()<<" Gy"<<endl;
  stream<<"Fraction = "<<treatment->getFraction()<<" Gy"<<endl;
  stream<<"---------------------------------------------";
  return stream;
}
