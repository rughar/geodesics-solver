#include <geodesics/solver.hpp>
#include <iostream>
#include <iomanip>   // std::scientific, std::setprecision
#include <fstream>
#include <cmath>

// Indices convention for Schwarzchild coordinates within plane (time, radius, angle)
constexpr auto T = 0;
constexpr auto R = 1;
constexpr auto P = 2;

template<class U>
class SchwarzschildPlane : public StormerVerletCore<U>
{
private:
  U h = 1;                      // Schwarzchild parameter  h = horizon radius

  using base = StormerVerletCore<U>;
  
  // Optional override. Default is Minkowski (-1,1,1,...)
  // Define metric coefficents from position coordinates.
  // In our case it is needed for functions get_norm, get_L, get_U, set_from_turning_points
  // and metric coefficents are also used in Christoffel_symbols.
  void metric()
  {
    base::g[T][T] = -1 + h / base::x[R];
    base::g[R][R] = -1 / base::g[T][T];
    base::g[P][P] = base::x[R] * base::x[R];
  }

  // Must have override. Define all non-zero Christoffel symbols. 
  void Christoffel_symbols() 
  {
    base::check_metric();       // Call this only if you want to use metric coefficents

    U g_tt_r = -h / (base::x[R] * base::x[R]);

    base::Gamma[R][T][T] = base::g[T][T] * g_tt_r / 2;
    base::Gamma[R][R][R] = base::g[R][R] * g_tt_r / 2;
    base::Gamma[R][P][P] = base::g[T][T] * base::x[R];

    base::Gamma[T][T][R] = base::Gamma[T][R][T] = -base::Gamma[R][R][R];
    base::Gamma[P][P][R] = base::Gamma[P][R][P] = 1 / base::x[R];
  }

public:
  // Set number of dimensions into default constructor as following.
  SchwarzschildPlane()
  {
    base::set_dimension(3);
  }

  // Following functions are optional for this specific case
  U get_L()
  {
    base::check_metric();
    return base::g[P][P] * base::u[P];
  }

  U get_E()
  {
    base::check_metric();
    return -base::g[T][T] * base::u[T];
  }

  // Function to set initial position and veloctity from radius and turning points
  void set_from_turning_points(const U r_init, const U r_min, const U r_max)
  {
    base::x[T] = 0.0;
    base::x[R] = r_init;
    base::x[P] = 0.0;
    base::check_metric();

    const U t1 = -1 + h / r_min;
    const U t2 = -1 + h / r_max;
    const U tmp1 = -1 + h / (r_max + r_min) - t1 - t2;
    const U tmp2 = h * r_min * r_max / (r_min + r_max);

    base::u[T] = -sqrt(t1 * t2 / tmp1) / base::g[T][T];    
    base::u[P] = sqrt(tmp2 / tmp1) / base::g[P][P];
    base::set_ui_from_metric(R, -1.0);
  }
};


int main(void) {
  
  std::ofstream filep("particle_data.txt");
  std::ofstream filei("simulation_data.txt");
  std::ofstream files("stepsize_curve_data.txt");

  // Set scientific format with enough precision
  filep.imbue(std::locale::classic());
  filei.imbue(std::locale::classic());
  files.imbue(std::locale::classic());
  filep << std::scientific << std::setprecision(16);
  filei << std::scientific << std::setprecision(16);
  files << std::scientific << std::setprecision(16);

  SchwarzschildPlane <double> core;

  // Set initial values computed from turning points and initial r
  // r_start = 15.0 , r_min = 8.0 , r_max = 50.0
  core.set_from_turning_points(15.0, 8.0, 50.0);

  auto L0 = core.get_L();
  auto E0 = core.get_E();
  auto G0 = core.get_norm();

  double stepsize_curve_param = 0.45;

  std::vector<double> x_ref = { 1.0, 50.0, 6.28 };
  std::vector<double> u_ref = { 1.0, 0.1, 0.05 };

  double dtcore = stepsize_curve_param * core.x[R];
  double dt = dtcore;
  core.step_1(dt);
  dtcore = (dtcore + stepsize_curve_param * core.x[R]);
  core.step_1(-dt);
  dt = dtcore / 2;

  for (size_t i = 0; i < 250; i++) {
    core.step_3(dt, true);
    //dt = core.suggest_stepsize_3(x_ref, u_ref, 10e-8);
    dt = 2 * stepsize_curve_param * core.x[R] - dt;

    // save X,Y coordinates in Schwarzchild plane
    filep << core.x[R] * cos(core.x[P]) << " " << core.x[R] * sin(core.x[P]) << '\n';

    // save invariants  norm, E, L
    filei << (core.get_norm() - G0) / G0 << " " << (core.get_E() - E0) / E0 << " " << (core.get_L() - L0) / L0 << '\n';

    // save stepsize curve
    files << core.x[R] << " " << dt << " " << dt / core.x[R] << '\n';
  }

  return 0;
}