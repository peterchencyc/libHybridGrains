// TimeUtils.h
//
// Breannan Smith
// Last updated: 09/03/2015

#ifndef TIME_UTILS_H
#define TIME_UTILS_H

#include <string>
#include <sys/time.h>
#include <time.h>

namespace TimeUtils {

// Returns current time as: Day Month DayOfMonth Hour:Minute:Second Year.
//   For example: Fri Jun 21 14:55:24 2013.
std::string currentTime();

} // namespace TimeUtils

class TimingTools final {

public:
  TimingTools();

  void start();

  void stop(const std::string &caption, const bool print = true);

  void serialize(std::ostream &output_stream) const;
  void deserialize(std::istream &input_stream);

  const double &elapsedTime() const;

private:
  timeval s;
  timeval e;
  double elapsed_sec;
};

#endif
