#pragma once
#include "ling.hpp"
#include <vector>
#include <cmath>

// =============================================================================
//  FILE: solver.h - symplectic (Stormer-Verlet) core for geodesic integration
// =============================================================================
//  High level summary
//  ------------------
//    * Template parameter  U  - floating point type (e.g. double, long double).
//    * N dimensional phase space: position  x[i]  and four velocity  u[i].
//    * Physics specific details (metric  g  and Christoffel symbols Gamma) are
//      provided by a derived class via two virtual functions:
//          - metric()               -> fills   g   from current   x.
//          - Christoffel_symbols()  -> fills   Gamma   (may access current g).
//    * Time integration:
//          -  step_1  = 2nd order velocity Verlet (one kick drift kick).
//          -  step_2  = 4th order Yoshida composition (quartic coefficients).
//          -  step_3  = 6th order Yoshida composition (sextic coefficients).
//
//  The file keeps the interface minimal: user code sets the dimension, initial
//  state, then repeatedly calls  step_*() .
// =============================================================================
template<class U>
class StormerVerletCore {

private:

  // ---- Yoshida composition coefficients (publications 1990) ---------------
  // https://aiichironakano.github.io/phys516/Yoshida-symplectic-PLA00.pdf
  const U quartic[2] = {
    1.351207191959657634047687808971460826922,
    -1.702414383919315268095375617942921653844
  };
  const U sectic[4] = {
    0.7845136104775572638194976338663498757768,
    0.2355732133593581336847931829785346016865,
    -1.177679984178871006946415680964315734639,
    1.315186320683911218884249728238862514352
  };

  std::vector<std::vector<U>> mat;                  // n*n temporary matrix for implicit boost
  bool metric_updated;                              // lazy evaluation flag

  void push(const U dt);                            // Drift  (x := x + u * dt)
  void boost(const U dt, const bool control=false); // Kick   (u := inv[I + dt * Gamma.u].u)
  void update_mat(const U dt);                      // Assemble (I + dt * Gamma.u) into  mat

  std::vector<U> x_e, u_e;                          // Used for evaluation of highest eigenvalue


protected:

  size_t n;                                         // dimension
  std::vector<std::vector<U>> g;                    // metric  g_ij(x)
  std::vector<std::vector<std::vector<U>>> Gamma;   // Christoffels Gamma^i_jk

  // Hooks for the derived concrete spacetime --------------------------------
  void check_metric();
  virtual void Christoffel_symbols() = 0;           // Christoffel symbols definition. Must be overriden.
  virtual void metric();                            // Optional routine to compute metric if necessary 
                                                    // (see below which functions require metric).

public:
  U omega;
  std::vector<U> x, u;                              // position, velocity

  // ---- integrator entry‑points --------------------------------------------
  void step_1(const U dt, const bool control = false); // 2nd‑order
  void step_2(const U dt, const bool control = false); // 4th‑order via quartic composition
  void step_3(const U dt, const bool control = false); // 6th‑order via sextic composition
  
  // ---- configuration ------------------------------------------------------
  void set_dimension(const size_t dim);             // allocate containers
  void set_ui_from_metric(size_t i, const U norm);  // adjust u[i] s.t. g.u.u = norm, require metric

  // ---- diagnostics --------------------------------------------------------
  U get_norm();                                     // get norm = g.u.u, require metric

  // ---- error analysos -----------------------------------------------------
  template<class V> U get_top_eigen(const V& x_ref, const V& u_ref) const;                    // estimate of highest eigenvalue                
  template<class V> U suggest_stepsize_1(const V& x_ref, const V& u_ref, const U tol) const;  // suggest stepsize for step_1
  template<class V> U suggest_stepsize_2(const V& x_ref, const V& u_ref, const U tol) const;  // suggest stepsize for step_2
  template<class V> U suggest_stepsize_3(const V& x_ref, const V& u_ref, const U tol) const;  // suggest stepsize for step_3
};

// ---------------------------------------------------------------------------
//  IMPLEMENTATION (header‑only template)
// ---------------------------------------------------------------------------

// push  — drift:  x := x + u * dt  ------------------------------------------
template<class U>
inline void StormerVerletCore<U>::push(const U dt)
{
  for (size_t i = n; i--;)
    x[i] += u[i] * dt;
  metric_updated = false;                           // geometry depends on x
}

// boost — implicit velocity update via linear solve -------------------------
//   (I + dt * Gamma.u_{old}).u_{new} = u_{old}
//
//   1. assemble  mat = I + dt * Gamma.u
//   2. inplace LU  (O(n^3))
//   3. solve  mat.u_{new} = u_{old}
template<class U>
inline void StormerVerletCore<U>::boost(const U dt, const bool control)
{
  update_mat(dt);

  if (control) {
    for (size_t i = n; i--;) {
      U tmp = 2 * u[i];
      for (size_t j = n; j--;)
        tmp -= mat[i][j] * u[j];
      u_e[i] = tmp;
      x_e[i] = x[i] + u[i] * dt / 2;
    }
  }

  lu_naive(n, mat);
  fb_naive(n, mat, u);
}

// update_mat — fill I + dt * Gamma.u into  mat  -----------------------------
template<class U>
inline void StormerVerletCore<U>::update_mat(const U dt)
{
  Christoffel_symbols();
  omega = 0.0;
  for (size_t i = n; i--;)
    for (size_t j = n; j--;) {
      U tmp = U(0);
      for (size_t k = n; k--;)
        tmp += Gamma[i][j][k] * u[k];
      mat[i][j] = tmp * dt + (i == j ? 1 : 0);
      omega = std::max(omega, fabs(tmp));
    }
}

// ensure  g  is up‑to‑date --------------------------------------------------
template<class U>
inline void StormerVerletCore<U>::check_metric()
{
  if (!metric_updated)
    metric();
  metric_updated = true;
}

// Default metric, set as Minkowski
template<class U>
inline void StormerVerletCore<U>::metric()
{
  g[0][0] = -1;
  for (size_t i = 1; i < n; i++)
    g[i][i] = 1;
}

// 2nd‑order velocity‑Verlet -------------------------------------------------
//   kick‑half -> drift -> kick‑half
template<class U>
inline void StormerVerletCore<U>::step_1(const U dt, const bool control)
{
  push(dt / 2);
  boost(dt, control);
  push(dt / 2);

  if (control) {
    for (size_t i = n; i--;) {
      x_e[i] = (x_e[i] - x[i]) / (dt * dt / 2);
      u_e[i] = (u_e[i] - u[i]) / (dt * dt / 2);
    }
  }
}

// Yoshida 4th‑order composition ---------------------------------------------
template<class U>
inline void StormerVerletCore<U>::step_2(const U dt, const bool control)
{
  for (int i : {0, 1})
    step_1(quartic[i] * dt, false);
  step_1(quartic[0] * dt, control);
}

// Yoshida 6th‑order composition ---------------------------------------------
template<class U>
inline void StormerVerletCore<U>::step_3(const U dt, const bool control)
{
  for (int i : {0, 1, 2, 3, 2, 1})
    step_1(sectic[i] * dt, false);
  step_1(sectic[0] * dt, control);
}

// set_dimension — allocate & zero the containers ----------------------------
template<class U>
inline void StormerVerletCore<U>::set_dimension(const size_t dim)
{
  n = dim;
  metric_updated = false;
  Gamma.assign(n, std::vector<std::vector<U>>(n, std::vector<U>(n, U{})));
  mat.assign(n, std::vector<U>(n, U{}));
  g.assign(n, std::vector<U>(n, U{}));
  x.assign(n, U{});
  u.assign(n, U{});
  x_e.assign(n, U{});
  u_e.assign(n, U{});
}

// set u[i] so that g_{ij} u^i u^j = norm ------------------------------------
//  Solves a quadratic equation in u[i].  If no real solution exists, u[i]
//  remains unchanged.
template<class U>
inline void StormerVerletCore<U>::set_ui_from_metric(const size_t i, const U norm)
{
  u[i] = 0;

  U gtot = get_norm();
  U ltot = 0;

  check_metric();
  for (size_t k = n; k--;)
    ltot -= g[i][k] * u[k];

  if (gtot == norm)                                 // already satisfied (rare)
    return;

  U A = gtot - norm;
  U D = ltot * ltot - g[i][i] * A;
  
  if (D < 0)                                        // no real root
    return;

  D = sqrt(D);

  U x1 = (norm - gtot) / (ltot + D);
  U x2 = (norm - gtot) / (ltot - D);
  
  u[i] = std::max(x1, x2);                          // choose positive root by convention
}

// returns   g_{ij} u^i u^j   (should be −1 for timelike) --------------------
template<class U>
inline U StormerVerletCore<U>::get_norm()
{
  check_metric();
 
  U gtot = 0;
  for (size_t k = n; k--;)
    for (size_t l = n; l--;)
      gtot += g[k][l] * u[k] * u[l];
  return gtot;
}

// Subtract rough estimate of highest eigen value (assume stiff problem]
template<class U>
template<class V>
inline U StormerVerletCore<U>::get_top_eigen(const V& x_ref, const V& u_ref) const
{
  U tmp = 0.0;
  for (size_t k = 0; k < 3; k++) {
    U tmpx = abs(x_e[k]) / x_ref[k];
    U tmpu = abs(u_e[k]) / u_ref[k];
    tmp = std::max(tmp, std::max(tmpx, tmpu));
  };
  return tmp;
}

template<class U>
template<class V>
inline U StormerVerletCore<U>::suggest_stepsize_1(const V& x_ref, const V& u_ref, const U tol) const
{
  return pow(tol * 12 * pow(get_top_eigen(x_ref, u_ref), -4.0 / 3.0), 1.0 / 3.0);
}

template<class U>
template<class V>
inline U StormerVerletCore<U>::suggest_stepsize_2(const V& x_ref, const V& u_ref, const U tol) const
{
  return pow(tol * 360 * pow(get_top_eigen(x_ref, u_ref), -6.0 / 3.0), 1.0 / 5.0);
}

template<class U>
template<class V>
inline U StormerVerletCore<U>::suggest_stepsize_3(const V& x_ref, const V& u_ref, const U tol) const
{
  return pow(tol * 20160 * pow(get_top_eigen(x_ref, u_ref), -8.0 / 3.0), 1.0 / 7.0);
}
