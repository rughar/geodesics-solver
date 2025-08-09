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
  void boost(const U dt);                           // Kick   (u := inv[I + dt * Gamma.u].u)
  void get_vel_jac();                               // Assemble Gamma.u into mat

  // ---- error analysis -----------------------------------------------------
  U get_top_eigen(const U dt);                      // estimate of highest eigenvalue   

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
  std::vector<U> x, u;                              // position, velocity

  // ---- integrator entry‑points --------------------------------------------
  void step_1(const U dt);                          // 2nd‑order
  void step_2(const U dt);                          // 4th‑order via quartic composition
  void step_3(const U dt);                          // 6th‑order via sextic composition
  
  // ---- configuration ------------------------------------------------------
  void set_dimension(const size_t dim);             // allocate containers
  void set_ui_from_metric(size_t i, const U norm);  // adjust u[i] s.t. g.u.u = norm, require metric

  // ---- diagnostics --------------------------------------------------------
  template<class V1, class V2>
  U dot_product(const V1 &a, const V2 &b);          // Get scalar product = g.a.b, require metric    

  U suggest_stepsize(const U omegadt, const U dtold);
  U suggest_stepsize_init(const U omegadt, const U dtref);
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
inline void StormerVerletCore<U>::boost(const U dt)
{
  get_vel_jac();

  for(size_t i = n; i--;)
    for (size_t j = n; j--;)
      mat[i][j] = (i == j ? 1 : 0) - dt * mat[i][j];
 
  lu_naive(n, mat);
  fb_naive(n, mat, u);
}

// ensure  g  is up‑to‑date --------------------------------------------------
template<class U>
inline void StormerVerletCore<U>::check_metric()
{
  if (!metric_updated)
    metric();
  metric_updated = true;
}

// Default metric
template<class U>
inline void StormerVerletCore<U>::metric()
{
  for (size_t i = 0; i < n; i++)
    g[i][i] = 1;
}

template<class U>
inline void StormerVerletCore<U>::get_vel_jac()
{
  Christoffel_symbols();
  for (size_t i = n; i--;)
    for (size_t j = n; j--;) {
      U tmp = U(0);
      for (size_t k = n; k--;)
        tmp -= Gamma[i][j][k] * u[k];
      mat[i][j] = tmp;
    }
}

// 2nd‑order velocity‑Verlet -------------------------------------------------
//   kick‑half -> drift -> kick‑half
template<class U>
inline void StormerVerletCore<U>::step_1(const U dt)
{
  push(dt / 2);
  boost(dt);
  push(dt / 2);
}

// Yoshida 4th‑order composition ---------------------------------------------
template<class U>
inline void StormerVerletCore<U>::step_2(const U dt)
{
  for (int i : {0, 1, 0})
    step_1(quartic[i] * dt);
}

// Yoshida 6th‑order composition ---------------------------------------------
template<class U>
inline void StormerVerletCore<U>::step_3(const U dt)
{
  for (int i : {0, 1, 2, 3, 2, 1, 0})
    step_1(sectic[i] * dt);
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
}

// set u[i] so that g_{ij} u^i u^j = norm ------------------------------------
//  Solves a quadratic equation in u[i].  If no real solution exists, u[i]
//  remains unchanged.
template<class U>
inline void StormerVerletCore<U>::set_ui_from_metric(const size_t i, const U norm)
{
  u[i] = 0;
  check_metric();
  U gtot = dot_product(u, u);

  if (gtot == norm)                                 // already satisfied
    return;

  U ltot = 0;
  for (size_t k = n; k--;)
    ltot -= g[i][k] * u[k];

  U A = gtot - norm;
  U D = ltot * ltot - g[i][i] * A;
  
  if (D < 0)                                        // no real root
    return;

  D = sqrt(D);

  U x1 = (norm - gtot) / (ltot + D);
  U x2 = (norm - gtot) / (ltot - D);
  
  u[i] = std::max(x1, x2);                          // choose positive root by convention
}

template<class U>
inline U StormerVerletCore<U>::suggest_stepsize(const U omegadt, const U dtold)
{
  auto tmp = omegadt / get_top_eigen(0.1 * dtold);
  return tmp * tmp / dtold;
}

template<class U>
inline U StormerVerletCore<U>::suggest_stepsize_init(const U omegadt, const U dtref)
{
  auto omega = get_top_eigen(dtref);
  auto dt = omegadt / omega;
  step_1(dt);
  omega = get_top_eigen(dtref);
  step_1(-dt);
  return sqrt((omegadt / omega) * dt);
}

// Rough estimate of highest eigen value
template<class U>
inline U StormerVerletCore<U>::get_top_eigen(const U dt)
{
  U t1 = 0;
  U t2 = 0;
  std::vector<U> a(n);

  push(dt);
  get_vel_jac();
  for (size_t i = n; i--;) {
    U tmp = 0;
    t1 += mat[i][i];
    for (size_t j = n; j--;) {
      t2 += 2 * mat[i][j] * mat[j][i];
      tmp -= mat[i][j] * u[j];
    };
    a[i] = tmp;
  };

  push(-2 * dt);
  get_vel_jac();
  for (size_t i = n; i--;) {
    U tmp = 0;
    t1 += mat[i][i];
    for (size_t j = n; j--;) {
      t2 += 2 * mat[i][j] * mat[j][i];
      tmp -= mat[i][j] * u[j];
    };
    a[i] -= tmp;
  };

  push(dt);
 
  t2 += dot_product(a, u) / (2 * dt * dot_product(u, u));
  U det = 2 * t2 - t1 * t1;

  if (det > 0) {
    det = sqrt(det);
    return std::max(std::abs(t1 + det), std::abs(t1 - det)) / 2;
  }

  return sqrt(std::abs(t2 / 2));
}


template<class U>
template<class V1, class V2>
inline U StormerVerletCore<U>::dot_product(const V1& a, const V2& b)
{
  check_metric();

  U tmp = 0;
  for (size_t k = n; k--;)
    for (size_t l = n; l--;)
      tmp += g[k][l] * a[k] * b[l];
  return tmp;
}
