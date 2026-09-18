// Standalone concurrency regression (run from the repository root):
// c++ -std=c++17 -O1 -g -pthread -I dtcc_core/cpp/include \
//   tests/cpp/test_timer.cpp -o /tmp/dtcc-test-timer && /tmp/dtcc-test-timer
// Add -fsanitize=thread to check concurrent updates and reports for data races.

#include <atomic>
#include <cassert>
#include <sstream>
#include <thread>
#include <vector>

#include "Timer.h"

int main()
{
  constexpr int worker_count = 4;
  constexpr int iterations = 2000;
  std::atomic<int> ready{0};
  std::atomic<bool> start{false};
  std::vector<std::thread> workers;
  std::ostringstream output;
  auto *original_output = std::cout.rdbuf(output.rdbuf());

  for (int i = 0; i < worker_count; ++i)
  {
    workers.emplace_back([&, i]() {
      const auto name = "worker-" + std::to_string(i);
      ++ready;
      while (!start.load())
        std::this_thread::yield();
      for (int j = 0; j < iterations; ++j)
      {
        DTCC_BUILDER::Timer shared("shared");
        DTCC_BUILDER::Timer individual(name);
      }
    });
  }
  while (ready.load() != worker_count)
    std::this_thread::yield();
  start = true;
  for (int i = 0; i < 100; ++i)
    DTCC_BUILDER::Timer::report("Concurrent snapshot");
  for (auto &worker : workers)
    worker.join();

  output.str("");
  DTCC_BUILDER::Timer::report("Final counts");
  std::cout.rdbuf(original_output);

  // Inspect the public report, without exposing timer internals for testing.
  std::map<std::string, size_t> counts;
  std::istringstream lines(output.str());
  std::string line;
  while (std::getline(lines, line))
  {
    std::istringstream fields(line);
    std::string timestamp, component, level, name;
    double mean, total;
    size_t count;
    if (fields >> timestamp >> component >> level >> name >> mean >> total >> count)
      counts[name] = count;
  }
  assert(counts.size() == worker_count + 1);
  assert(counts.at("shared") == worker_count * iterations);
  for (int i = 0; i < worker_count; ++i)
    assert(counts.at("worker-" + std::to_string(i)) == iterations);
  std::cout << "Concurrent timer updates and reports passed\n";
}
