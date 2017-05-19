/**
 * @file Treatment.hpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.19.17 
 */

#ifndef DEF_Treatment
#define DEF_Treatment

#include <iostream>
#include <vector>

class Treatment{
public:
  Treatment();
  Treatment(const double fraction, const double totalDose,
	    const double interval,
	    const std::vector<bool> schedule);
  double getDuration() const;
  double getFraction() const;
  double getInterval() const;
  std::vector<bool> getSchedule() const;
  double getTotalDose() const;

private:
  double m_fraction;
  double m_interval;
  double m_totalDose;
  std::vector<bool> m_schedule;
};

std::ostream &operator<<(std::ostream &stream,
			 Treatment *treatment);
#endif
