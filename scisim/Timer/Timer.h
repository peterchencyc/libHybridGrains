#ifndef TIMER_H
#define TIMER_H

#include <string>
#include <unordered_map>

namespace TimerCollection {

const std::unordered_map<std::string, double> &timerValues();

void startTimer(const std::string &name);
void stopTimer(const std::string &name);

void enableTimers();
void disabletimers();

void printTotalTimerStats();

} // namespace TimerCollection

class RAIITimer final {

public:
  RAIITimer(const std::string &name);

  RAIITimer(const RAIITimer &) = delete;
  RAIITimer(RAIITimer &&) = delete;

  RAIITimer &operator=(const RAIITimer &) = delete;
  RAIITimer &operator=(RAIITimer &&) = delete;

  ~RAIITimer();

private:
  const std::string m_name;
};

#endif
