#ifndef DEF_Treatment
#define DEF_Treatment

#include <vector>

class Treatment{
public:
  Treatment();
  Treatment(double fraction, double totalDose, double interval,
	    std::vector<bool> schedule);
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
