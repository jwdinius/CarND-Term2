#include <uWS/uWS.h>
#include <iostream>
#include <fstream>
#include "json.hpp"
#include "PID.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <chrono>
#include <list>

// Uncomment this define to start to train PID controller.
//#define enable_tuning_mode


// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
double deg2rad(double x) { return x * M_PI / 180; }
double rad2deg(double x) { return x * 180 / M_PI; }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
std::string hasData(std::string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_last_of("]");
  if (found_null != std::string::npos) {
    return "";
  }
  else if (b1 != std::string::npos && b2 != std::string::npos) {
    return s.substr(b1, b2 - b1 + 1);
  }
  return "";
}

namespace
{
  enum class Mode { driving, tuning };

  class Timer
  {
  public:
    using high_resolution_clock = std::chrono::high_resolution_clock;
    using ms = std::chrono::milliseconds;

    Timer(bool run = true) {
      if (run)
        start();
    }

    void start() {
      start_ = high_resolution_clock::now();
    }

    ms duration() {
      high_resolution_clock::time_point now = high_resolution_clock::now();
      ms dur = std::chrono::duration_cast<ms>(now - start_);
      start_ = now;
      return dur;
    }

  private:
    high_resolution_clock::time_point start_;
  };

  class PIDTuning
  {
    // best params: iteration 181
    // .206488, 0, 1.63698

    static constexpr double track_distance_ = 2900.;
    static constexpr double inc_step = 1.1;
    static constexpr double dec_step = 0.9;
  public:
    PIDTuning(PID& pid, double maxThrottle = 0.3, Mode mode = Mode::driving);
    ~PIDTuning();

    bool isTuning() const { return mode_ == Mode::tuning; }
    double getBestError() const { return best_err_; }
    double getCurrThrottle() const { return curThrottle_; }

    double integrateVelocity(Timer::ms duration);
    bool trackPassedOrErrorIncreasing();
    void startTwiddle();
    void accumError(double error);
    bool testErrorIncreasing();
    void setSpeed(double speed);
    size_t getIters() { return iters_; }

    void printState();

  private:
    void updateParams();
    void reset();
    void nextParamToTune();
    bool trackWasPassed() const { return total_dist_ > track_distance_; }
    double calc_err() const;

  private:
    double total_time_;
    double maxThrottle_;
    double curThrottle_;
    Mode mode_;
    double total_dist_;
    double best_err_;
    double err_;
    double err_ticks_;
    PID& pid_;
    size_t iters_; // # of iterations
    int tuned_param_ind_;
    size_t sim_count_;
    double speed_;
    std::list<double> err_list_;
    std::ofstream outFile;


    enum ParamsState { inc_param, calc_best_err_after_inc, calc_best_err_after_dec };

    int state_; // 0-run after inc pd, 1-run after decrease dp

    // P.I.D. sequence
    std::vector<double> params_;
    std::vector<double> d_params_;
    std::vector<double> best_params_;

    static constexpr double tol_ = 0.2;
  };


  void reset_simulator(uWS::WebSocket<uWS::SERVER>& ws)
  {
    // reset
    std::string msg("42[\"reset\", {}]");
    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
  }
}


// Implementation of PIDTuning
namespace
{
  PIDTuning::PIDTuning(PID& pid, double maxThrottle, Mode mode)
    : total_time_(0)
    , maxThrottle_(maxThrottle)
    , curThrottle_(0.8)
    , mode_(mode)
    , total_dist_(0)
    , best_err_(0)
    , err_(0)
    , err_ticks_(1)
    , pid_(pid)
    , iters_(0)
    , tuned_param_ind_(0)
    , sim_count_(0)
    , speed_(0.)
    , state_(inc_param)
    , params_{ 0., 0., 0. }
//    , d_params_{ 1., 1., 1. } // is better to use if curThrottle_=0.2
    , d_params_{ 0.5, 0.5, 0.5 }  // is better to use if curThrottle_=0.3
    , best_params_{ 0., 0., 0. }
  {
    curThrottle_ = std::min(curThrottle_, maxThrottle_);
    outFile.open("Successes.txt");
  }
  PIDTuning::~PIDTuning(){
    outFile.close();
  }

  void PIDTuning::setSpeed(double speed)
  {
    speed_ = speed;
    return;
  }
  double PIDTuning::integrateVelocity(Timer::ms duration)
  {
    const auto dur_sec = (double)duration.count() / 1000.;
    total_time_ += dur_sec;
    total_dist_ += speed_ * dur_sec;
    return total_dist_;
  }

  void PIDTuning::reset()
  {
    total_time_ = 0;
    total_dist_ = 0;
    speed_ = 0.;
    sim_count_  = 0;
    err_list_.resize(0);
    iters_++;
  }

  void PIDTuning::startTwiddle()
  {
    if (isTuning()) {
      pid_.Init(params_[0], params_[1], params_[2]);
      printState();
    }
  }

  // add logic to end run after several frames where error doesn't decrease
  bool PIDTuning::trackPassedOrErrorIncreasing()
  {
    bool errorIncreasing = testErrorIncreasing();
    bool trackPassed = false;

    if (!errorIncreasing) {
      trackPassed = trackWasPassed();
    }

    if (trackPassed){
      // write output to file
      outFile << "Iter: " << iters_ << " best_err: " << best_err_ << " err: " << calc_err() << std::endl;
      outFile << "params: " << params_[0] << ", " << params_[1] << ", " << params_[2] << ", learning:" << isTuning() << std::endl;
    }

    if (errorIncreasing || trackPassed) {
      updateParams();
      reset();
    }

    return errorIncreasing || trackPassed;
  }

  bool PIDTuning::testErrorIncreasing()
  {
    double max_err;
    if (sim_count_ < 50) {
      err_list_.push_front(calc_err());
      sim_count_++;
      return false;
    }
    else {
      err_list_.push_front(calc_err());
      err_list_.pop_back();

    }
    sim_count_++;

    max_err = std::max(err_list_.front(), err_list_.back());
    bool output = (max_err > 5.) && ((err_list_.front() > err_list_.back()) ? true : false);
    output = output && (total_time_ > 5.);
    output = output && (speed_ > .1);
    // not only check for error not increasing, but we may have the case where we're stuck 
    // (i.e. error isn't too large but it's not decreasing after awhile)
    return output;
  }

  double PIDTuning::calc_err() const
  {
    double err = err_ / err_ticks_;
    err *= track_distance_ / total_dist_; // adds penalty if car cannot drive along the track
    return err;
  }

  void PIDTuning::accumError(double error)
  {
    err_ticks_++;
    err_ += error * error;
  }

  void PIDTuning::updateParams()
  {
    double err = calc_err();
    assert(err > 0);

    // Check if tuning is finished
    double p_sum = 0.;
    for (auto par : d_params_) {
      p_sum += par;
    }
    if (p_sum < tol_) {
      // tuning is finished for current speed!
      double throttle = std::min(curThrottle_ + 0.05, maxThrottle_);
      if (curThrottle_ == throttle) {
        mode_ = Mode::driving;
        return;
      }
      else { // continue tuning
        curThrottle_ = throttle;
        if (curThrottle_ == maxThrottle_) {

        }
      }
    }

    if (!iters_) {
      best_err_ = err;
    }

    {
      if (state_ == inc_param) {
        params_[tuned_param_ind_] += d_params_[tuned_param_ind_];
        pid_.Init(params_[0], params_[1], params_[2]);
        state_ = calc_best_err_after_inc; // run next error calculation on the same track
      }
      else if (state_ == calc_best_err_after_inc) {
        if (err < best_err_) {
          best_err_ = err;
          d_params_[tuned_param_ind_] *= inc_step;
          best_params_ = params_;

          nextParamToTune();

          params_[tuned_param_ind_] += d_params_[tuned_param_ind_];
          pid_.Init(params_[0], params_[1], params_[2]);
          state_ = calc_best_err_after_inc; // run next error calculation on the same track
        }
        else {
          params_[tuned_param_ind_] -= 2* d_params_[tuned_param_ind_];
          pid_.Init(params_[0], params_[1], params_[2]);
          state_ = calc_best_err_after_dec; // run next error calculation on the same track
        }
      }
      else if (state_ == calc_best_err_after_dec) {
        if (err < best_err_) {
          best_err_ = err;
          d_params_[tuned_param_ind_] *= inc_step;
          best_params_ = params_;
        }
        else {
          params_[tuned_param_ind_] += d_params_[tuned_param_ind_];
          d_params_[tuned_param_ind_] *= dec_step;
        }

        nextParamToTune();

        params_[tuned_param_ind_] += d_params_[tuned_param_ind_];
        pid_.Init(params_[0], params_[1], params_[2]);
        state_ = calc_best_err_after_inc; // run next error calculation on the same track
      }
    }

    err_ = 0;
    err_ticks_ = 1;
  }

  void PIDTuning::nextParamToTune()
  {
    ++tuned_param_ind_;

    if (tuned_param_ind_ == 1) {
      ++tuned_param_ind_; // skip I-term;
    }

    if (tuned_param_ind_ >= 3) {
      iters_++;
      tuned_param_ind_ = 0;
    }
  }

  void PIDTuning::printState()
  {
    if (mode_ == Mode::tuning) {
      std::cout << "Iter: " << iters_ << " best_err: " << best_err_ << " err: " << calc_err() << std::endl;
      std::cout << "params: " << params_[0] << ", " << params_[1] << ", " << params_[2] << ", learning:" << isTuning() << std::endl;
    }
  }

}


int main()
{
    uWS::Hub h;

  PID pid;
  // TODO: Initialize the pid variable.  Optimized by #defining enable_tuning_mode above
  pid.Init(0.206488, 0., 1.63698);
  
  Mode mode = Mode::driving;
  double maxThrottle = 0.5;

#ifdef enable_tuning_mode
  mode = Mode::tuning;
  maxThrottle = 0.5;
#endif

  Timer timer;
  PIDTuning tuning(pid, maxThrottle, mode);

  if (tuning.isTuning()) {
    tuning.startTwiddle();
  }

  h.onMessage([&pid, &tuning, &timer, maxThrottle](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (length && length > 2 && data[0] == '4' && data[1] == '2')
    {
      auto s = hasData(std::string(data));
      if (s != "") {
        auto j = json::parse(s);
        std::string event = j[0].get<std::string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          double cte = std::stod(j[1]["cte"].get<std::string>());
          tuning.setSpeed(std::stod(j[1]["speed"].get<std::string>()));
          /*
          * TODO: Calcuate steering value here, remember the steering value is
          * [-1, 1].
          * NOTE: Feel free to play around with the throttle and speed. Maybe use
          * another PID controller to control the speed!
          */

      pid.UpdateError(cte);
      const double steer_value = pid.TotalError();

      double dist = 0;
      Timer::ms duration = timer.duration();
      if (tuning.isTuning()) {
        dist = tuning.integrateVelocity(duration);
        tuning.accumError(cte);
        if (tuning.trackPassedOrErrorIncreasing()) {
          reset_simulator(ws);
          timer.start();
        }
        if (tuning.getIters() > 300) {
          ws.terminate();
          std::cout << "done" << std::endl;
        }
      }

      tuning.printState();

      double throttle = !tuning.isTuning() ? maxThrottle : tuning.getCurrThrottle();

      if (!tuning.isTuning()) {
        if (std::abs(cte) > 1) {
          throttle = -0.001;
        }
      }

          json msgJson;
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle;// 0.3;
          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
    }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1)
    {
      res->end(s.data(), s.length());
    }
    else
    {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port))
  {
    std::cout << "Listening to port " << port << std::endl;
  }
  else
  {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}