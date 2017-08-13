#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

MPC tmp;

// TODO: Set the timestep length and duration
size_t N = 30;
double dt = 0.1;

// helpful for indexing
size_t x_start = 0;
size_t y_start = N + x_start;
size_t psi_start = N + y_start;
size_t v_start = N + psi_start;
size_t cte_start = N + v_start;
size_t epsi_start = N + cte_start;
size_t delta_start = N + epsi_start;
size_t a_start = N-1 + delta_start;

// reference speed (tunable)
double ref_v = 45.*tmp.mph2mps(); // MPH -> m/s

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    // the cost is in element 0
    fg[0] = 0.;

    // reference state info
    for (unsigned int t = 0; t < N; t++) {
      fg[0] += CppAD::pow(vars[cte_start + t], 2);
      fg[0] += CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += CppAD::pow(vars[v_start + t] - ref_v, 2);
    }
    // Minimize control effort.
    for (unsigned int t = 0; t < N - 1; t++) {
      fg[0] += 6500.*CppAD::pow(vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the difference between control iterations.
    for (unsigned int t = 0; t < N - 2; t++) {
      fg[0] += 20000.*CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    // Add constraints
    fg[1 + x_start]    = vars[x_start];
    fg[1 + y_start]    = vars[y_start];
    fg[1 + psi_start]  = vars[psi_start];
    fg[1 + v_start]    = vars[v_start];
    fg[1 + cte_start]  = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    //  Rest of the constraints
    for (unsigned t = 1; t < N; t++) {
      // The state at next time (t+1) .
      AD<double> x_tp1    = vars[x_start + t];
      AD<double> y_tp1    = vars[y_start + t];
      AD<double> psi_tp1  = vars[psi_start + t];
      AD<double> v_tp1    = vars[v_start + t];
      AD<double> cte_tp1  = vars[cte_start + t];
      AD<double> epsi_tp1 = vars[epsi_start + t];

      // The state at time t.
      AD<double> x_t    = vars[x_start + t - 1];
      AD<double> y_t    = vars[y_start + t - 1];
      AD<double> psi_t  = vars[psi_start + t - 1];
      AD<double> v_t    = vars[v_start + t - 1];
      AD<double> cte_t  = vars[cte_start + t - 1];
      AD<double> epsi_t = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta_t = vars[delta_start + t - 1];
      AD<double> a_t     = vars[a_start + t - 1];

      AD<double> x2 = x_t * x_t;
      AD<double> x3 = x_t * x2;
      
      AD<double> f_t       = coeffs[0] + coeffs[1] * x_t + coeffs[2] * x2 + coeffs[3] * x3;
      AD<double> psi_des_t = CppAD::atan(coeffs[1] + 2. * coeffs[2] * x_t + 3. * coeffs[3] * x2);

      // Recall the equations for the model:
      // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
      // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
      // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
      // v_[t] = v[t-1] + a[t-1] * dt
      // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
      // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
      //
      // remember: the goal here is enforce the constraint that the updated state
      // equals the results of the update equations
      fg[1 + x_start + t]    = x_tp1    - (x_t + v_t * CppAD::cos(psi_t) * dt);
      fg[1 + y_start + t]    = y_tp1    - (y_t + v_t * CppAD::sin(psi_t) * dt);
      fg[1 + psi_start + t]  = psi_tp1  - (psi_t + v_t * (delta_t / Lf) * dt);
      fg[1 + v_start + t]    = v_tp1    - (v_t + a_t * dt);
      fg[1 + cte_start + t]  = cte_tp1  - ((f_t - y_t) + (v_t * CppAD::sin(epsi_t) * dt));
      fg[1 + epsi_start + t] = epsi_tp1 - ((psi_t - psi_des_t) + v_t * (delta_t / Lf) * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = N * 6 + (N - 1) * 2;
  //Set the number of constraints
  size_t n_constraints = N*6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set the initial variable values
  vars[x_start]    = state(0);
  vars[y_start]    = state(1);
  vars[psi_start]  = state(2);
  vars[v_start]    = state(3);
  vars[cte_start]  = state(4);
  vars[epsi_start] = state(5);

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  // [-25deg,+25deg] for steering
  for (i = delta_start; i < a_start; i++) {
      vars_lowerbound[i] = deg2rad(-25.);
      vars_upperbound[i] = deg2rad( 25.);
  }
  // Set lower and upper limits for a.
  for (i = a_start; i < n_vars; i++) {
      vars_lowerbound[i] = -1.;
      vars_upperbound[i] =  1.;
  }

  // Set all state upper and lower limits
  // to something big
  for (i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.e10;
    vars_upperbound[i] =  1.e10;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // Set lower and upper limits for the initial state constraints
  // (we are bound to start from where we're at, and there is nothing
  // to be done about it!)
  constraints_lowerbound[x_start]    = state(0);
  constraints_upperbound[x_start]    = state(0);
  constraints_lowerbound[y_start]    = state(1);
  constraints_upperbound[y_start]    = state(1);
  constraints_lowerbound[psi_start]  = state(2);
  constraints_upperbound[psi_start]  = state(2);
  constraints_lowerbound[v_start]    = state(3);
  constraints_upperbound[v_start]    = state(3);
  constraints_lowerbound[cte_start]  = state(4);
  constraints_upperbound[cte_start]  = state(4);
  constraints_lowerbound[epsi_start] = state(5);
  constraints_upperbound[epsi_start] = state(5);

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // if the solution's not cool, just get out of here
  if (!ok) return {};

  // otherwise, get the solution
  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> sol;
  
  sol.push_back(solution.x[delta_start]);
  sol.push_back(solution.x[a_start]);

  // Add x and y values for optimal chosen trajectory
  for (i = 0; i < N; ++i) {
    sol.push_back(solution.x[x_start + i]);
    sol.push_back(solution.x[y_start + i]);
  }

  return sol;
}
