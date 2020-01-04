// TimeUtils.cpp
//
// Breannan Smith
// Last updated: 09/03/2015

#include "TimeUtils.h"

#include "scisim/StringUtilities.h"
#include "scisim/Utilities.h"
#include <ctime>
#include <iostream>

std::string TimeUtils::currentTime() {
  std::time_t result = std::time(nullptr);
  return StringUtilities::trim(std::ctime(&result));
}

static timeval getTimeVal() {
  timeval tod;
  gettimeofday(&tod, nullptr);
  return tod;
}

TimingTools::TimingTools() : s(getTimeVal()), e(s), elapsed_sec(0.0) {}

void TimingTools::start() { s = getTimeVal(); }

void TimingTools::stop(const std::string &caption, const bool print) {
  e = getTimeVal();
  elapsed_sec = (double(e.tv_sec) + double(e.tv_usec) * 0.000001) -
                (double(s.tv_sec) + double(s.tv_usec) * 0.000001);
  if (print) {
    std::cout << caption << ": " << elapsed_sec << "(sec)" << std::endl;
  }
}

const double &TimingTools::elapsedTime() const { return elapsed_sec; }

void TimingTools::serialize(std::ostream &output_stream) const {
  Utilities::serializeBuiltInType(s.tv_sec, output_stream);
  Utilities::serializeBuiltInType(s.tv_usec, output_stream);
  Utilities::serializeBuiltInType(e.tv_sec, output_stream);
  Utilities::serializeBuiltInType(e.tv_usec, output_stream);
}

void TimingTools::deserialize(std::istream &input_stream) {
  s.tv_sec = Utilities::deserialize<time_t>(input_stream);
  s.tv_usec = Utilities::deserialize<suseconds_t>(input_stream);
  e.tv_sec = Utilities::deserialize<time_t>(input_stream);
  e.tv_usec = Utilities::deserialize<suseconds_t>(input_stream);
}
