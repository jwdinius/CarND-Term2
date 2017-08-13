#ifndef MPC_H
#define MPC_H

#include <vector>
#include <math.h>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();
  double mph2mps() { return 0.44704; }
  // For converting back and forth between radians and degrees.
  double pi() { return M_PI; }
  double deg2rad(double x) { return x * pi() / 180; }
  double rad2deg(double x) { return x * 180 / pi(); }

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
