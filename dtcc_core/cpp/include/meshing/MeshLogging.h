#ifndef MESH_LOGGER_H
#define MESH_LOGGER_H

#include "Table.h"
#include "Logging.h"
#include "Timer.h"

#include <iostream>


namespace DTCC_BUILDER
{

class MeshLogger
{
public:

  size_t step;

  MeshLogger()
    : _table("Volume Mesh Generation Metrics"),
      _timer("Volume Mesh Generation Metrics", false),
      _running(false),
      step(0)
  {
    // add header row
    _table.rows.push_back({ "Step", "Name", "Time (s)",
                             "Min AR", "Median AR", "Max AR" });
  }

  inline void step_begin(const std::string &step_name){
    
    if (_running)
    {
      _timer.stop();
      warning("step_begin() called before previous step_stop(); auto-stopping timer.");
    }

    _current_step_name = step_name;
    info(_current_step_name);
    step ++;
    _running = true;
    _timer.start();
  }

  inline void step_stop(){
    
    if (!_running)
    {
      warning("step_stop() called but timer is not running.");
      return;
    }

    _timer.stop();
    _running = false;
  }

  /// Call this after each step, passing in a vector of all element ARs
  void log_step(const double min_ar,
               const double median_ar,
               const double max_ar)
  { 
    if (_running) {
      // auto-stop
      _timer.stop();
      _running = false;
    }
    auto elapsed = _timer.time(); 

    _table.rows.push_back({
      std::to_string(step),
      _current_step_name,
      format(elapsed),
      format(min_ar),
      format(median_ar),
      format(max_ar)
    });
  }

  /// Print the whole table to stdout
  void summary() const
  {
    std::cout << _table.__str__() << std::endl;
  }

private:
  Table _table;

  Timer _timer;

  bool _running;

  std::string _current_step_name;

  inline static std::string format(double x)
  {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(4) << x;
    return ss.str();
  }

  
};

}
#endif // MESH_LOGGER_H