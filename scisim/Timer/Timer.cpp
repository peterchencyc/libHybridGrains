#include "Timer.h"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <stack>
#include <vector>

class Timer final {

public:
  Timer(const std::string &name)
      : m_name(name), m_running(false), m_elapsed(0.0) {}

  void start() {
    if (m_running) {
      return;
    }
    m_running = true;
    m_start = std::chrono::high_resolution_clock::now();
  }

  void stop() {
    if (!m_running) {
      return;
    }
    m_running = false;
    const std::chrono::time_point<std::chrono::high_resolution_clock> end{
        std::chrono::high_resolution_clock::now()};
    const std::chrono::duration<double> duration{end - m_start};
    m_elapsed += duration.count();
  }

  const std::string &name() const { return m_name; }

  const double &elapsedTime() const { return m_elapsed; }

private:
  const std::string m_name;
  bool m_running;
  double m_elapsed;
  std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
};

static std::unordered_map<std::string, double> g_timer_values;
static std::stack<Timer> g_timers;
static bool g_enabled{false};

void TimerCollection::startTimer(const std::string &name) {
  if (!g_enabled) {
    return;
  }

  // Stop the previous timer
  if (!g_timers.empty()) {
    g_timers.top().stop();
  }
  // Start the new timer
  g_timers.emplace(name);
  g_timers.top().start();
  auto timer_value = g_timer_values.find(name);
  if (timer_value == g_timer_values.end()) {
    g_timer_values.insert(std::make_pair(name, 0.0));
  }
}

void TimerCollection::stopTimer(const std::string &name) {
  if (!g_enabled) {
    return;
  }

  // Stop the current timer
  if (g_timers.empty()) {
    std::cerr << "Error, timer start/stop mismatch for timer: " << name
              << std::endl;
    std::cerr << "Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  g_timers.top().stop();
  if (g_timers.top().name() != name) {
    std::cerr << "Error, timer start/stop mismatch for timer: " << name
              << std::endl;
    std::cerr << "Exiting." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  g_timer_values[name] += g_timers.top().elapsedTime();
  g_timers.pop();
  // Start the previous timer
  if (!g_timers.empty()) {
    g_timers.top().start();
  }
}

const std::unordered_map<std::string, double> &TimerCollection::timerValues() {
  return g_timer_values;
}

void TimerCollection::enableTimers() { g_enabled = true; }

void TimerCollection::disabletimers() { g_enabled = false; }

void TimerCollection::printTotalTimerStats() {
  // Sum all timers
  double total_time = 0.0;
  for (const auto &timer : TimerCollection::timerValues()) {
    total_time += timer.second;
  }

  // Sort in decreasing order
  std::vector<std::pair<double, std::string>> timing_stats;
  timing_stats.reserve(TimerCollection::timerValues().size());
  for (const auto &timer : TimerCollection::timerValues()) {
    timing_stats.emplace_back(timer.second, timer.first);
  }
  // std::sort( timing_stats.begin(), timing_stats.end(), []( const
  // std::pair<std::string,double>& a, const std::pair<std::string,double>& b ){
  // return a.second >= b.second; } );
  std::sort(timing_stats.begin(), timing_stats.end());
  std::reverse(timing_stats.begin(), timing_stats.end());

  // Print timing stats
  std::cout << "Total time: " << total_time << std::endl;
  std::cout << std::setw(80) << std::left << "Timer" << std::setw(20)
            << std::left << "Seconds" << std::setw(20) << "Percent"
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "---------------------------------------------"
            << std::endl;
  for (const auto &timer : timing_stats) {
    std::cout << std::setw(80) << std::left << timer.second << std::setw(20)
              << std::left << timer.first << std::setw(20)
              << 100.0 * (timer.first / total_time) << std::endl;
  }
}

RAIITimer::RAIITimer(const std::string &name) : m_name(name) {
  TimerCollection::startTimer(m_name);
}

RAIITimer::~RAIITimer() { TimerCollection::stopTimer(m_name); }
