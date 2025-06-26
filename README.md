# Geodesic‑type ODE Integrator

> Symplectic, time‑reversible integration of systems of the form\
> \(\ddot x^i + \Gamma^{i}_{\;jk}\,\dot x^j\dot x^k = 0\,,\qquad i=0,\dots,n-1\)\
> Implemented in C++ as a header‑only library that couples a **Störmer–Verlet** core with higher‑order **Yoshida** compositions [[Yoshida 1990](https://doi.org/10.1016/0375-9601\(90\)90092-3)].

---

## 1  Mathematical background

### 1.1  Geodesic equation and its cousins

The second‑order autonomous system above is best known from **general relativity**, but the *algebraic structure* also covers many mechanical or field‑theory problems whose acceleration is quadratic in the velocities:

- **Free rigid‑body rotation** (body angular momentum expressed in principal axes)
- **Motion in an electromagnetic field** — see § 1.2
- **Planar Kepler motion** in *spherical* coordinates

Whenever the RHS can be written as a bilinear map \(\dot u^i = -\Gamma^{i}_{\;jk} u^j u^k,\qquad u^i = \dot x^i,\) the library can integrate it without modification.

> **Split formulation.**  With phase vector \$(x^i,u^i)\$:\
> \(\dot x^i = u^i,\qquad \dot u^i = -\Gamma^{i}_{\;jk} u^j u^k.\)\
> This split is **symplectic** and **time‑reversible**, hence the choice of Verlet.

### 1.2  Electromagnetic field as a \$(n!+!1)\$‑D geodesic

Add one auxiliary coordinate \$x^{c}\$ with \$u^{c}=1\$ and set

- \$\Gamma^{c}\_{;ij}=0\$  (extra dimension is cyclic)
- \$\Gamma^{i}\_{;cc}\$ reproduces the *electric‑field* term
- \$\Gamma^{i}\_{;ck}\$ produces the magnetic term \$q,\mathbf u\times\mathbf B\$

Lorentz force = **quadratic geodesic part** + **linear correction**.  For *purely* linear Lorentz motion use the simpler, volume‑preserving **Boris leap‑frog** [[Boris 1970](https://ntrs.nasa.gov/citations/19710026052)].

### 1.3  Norms and invariants

Monitor e.g.\
\(g_{ij}u^i u^j = \text{const}\quad(=-1 \text{ for time‑like})\)\
and, where applicable, specific energy \$E\$ and angular momentum \$L\$.  Helpers `get_norm()`, `get_E()`, `get_L()` are provided. fileciteturn0file0turn0file2

---

## 2  Integrator hierarchy

| Order | Method                      | Entry‑point  |
| ----- | --------------------------- | ------------ |
| 2     | velocity Störmer–Verlet     | `step_1(dt)` |
| 4     | Yoshida quartic composition | `step_2(dt)` |
| 6     | Yoshida sextic composition  | `step_3(dt)` |

`step_2` and `step_3` are thin wrappers around `step_1` with Yoshida coefficients.

---

## 3  Advanced Störmer–Verlet formulation

### 3.1  Classic velocity‑Verlet

$$
\begin{aligned}
\mathbf x_{1/2}&=\mathbf x_0+\tfrac{\Delta t}{2}\,\mathbf u_0,\\[4pt]
\mathbf u_{1}&=\mathbf u_0+\Delta t\,\mathbf f(\mathbf x_{1/2},\mathbf u_0,\mathbf u_1),\\[4pt]
\mathbf x_{1}&=\mathbf x_{1/2}+\tfrac{\Delta t}{2}\,\mathbf u_1.
\end{aligned}
$$

Symmetric kicks

$$
\mathbf f_A=\ddot{\mathbf x}(\mathbf x_{1/2},(\mathbf u_0+\mathbf u_1)/2),\qquad
\mathbf f_B=\tfrac12[\ddot{\mathbf x}(\mathbf x_{1/2},\mathbf u_0)+\ddot{\mathbf x}(\mathbf x_{1/2},\mathbf u_1)].
$$

Take\
\(\boxed{\mathbf f=2\mathbf f_A-\mathbf f_B}\) so that for geodesics\
\$ f^i=-\Gamma^{i}\_{;jk}u\_0^{,j}u\_1^{,k}\$.

Resulting step

$$
\boxed{\begin{aligned}
 x_{1/2}^i &= x_0^i+\frac{\Delta t}{2}u_0^i,\\[6pt]
 u_1^i &= u_0^i-\Delta t\,\Gamma^{i}_{\;jk}u_0^{\,j}u_1^{\,k},\\[6pt]
 x_1^i &= x_{1/2}^i+\frac{\Delta t}{2}u_1^i
\end{aligned}}
$$

— implicit but **linear** in \$u\_1\$.

### 3.2  LU solve per step

Linear form\
\(\mathbf u_1=(\mathbf I+\Delta t\,\Gamma(\mathbf u_0))^{-1}\,\mathbf u_0,\;\Gamma(\mathbf u_0)^{i}{}_{k}=\Gamma^{i}_{\;jk}u_0^{\,j}.\)\
`StormerVerletCore` LU‑factorises this \$n\times n\$ matrix once per step (\$\mathcal O(n^3)\$).

---

## 4  Adaptive step‑size strategy

### 4.1  Why universal controllers fail

Reversibility demands a **palindromic** step list; classic PID or embedded‑RK controllers violate this and leak invariants.

### 4.2  Eigenvalue‑based heuristic

For Jacobian \$J\$:\
\(f\approx1/\lambda_{\max},\qquad \Delta t_{\text{new}}=2f-\Delta t_{\text{old}}.\) `StormerVerletCore` provides quick asymmetric probes:

```cpp
suggest_stepsize_k(x_ref,u_ref,tol); // k = 1,2,3
```

### 4.3  Recommended workflow

1. **Probe run** — call `step_3(dt,true)` and log the suggested `dt`s.
2. **Fit** an empirical map \$(\mathbf x,\mathbf u)\mapsto f\$.
3. **Production run** — feed that \$f\$ into the symmetric update formula; long‑term drifts disappear while stiff zones still adapt.

---

## 5  Quick start — Euler’s rigid‑body equations

```cpp
#include <geodesics/solver.hpp>
#include <array>
#include <fstream>

struct RigidBody : public StormerVerletCore<double> {
  std::array<double,3> I{2.0,1.0,0.5};
  std::array<double,3> M{0.0,0.0,0.1};

  RigidBody(){ set_dimension(4); }

  void Christoffel_symbols() override {
    // quadratic ω×Iω terms
    Gamma[idx(1,2,3)]=Gamma[idx(1,3,2)]=-(I[1]-I[2])/(2*I[0]);
    Gamma[idx(2,3,1)]=Gamma[idx(2,1,3)]=-(I[2]-I[0])/(2*I[1]);
    Gamma[idx(3,1,2)]=Gamma[idx(3,2,1)]=-(I[0]-I[1])/(2*I[2]);
    // constant torque via u^0
    Gamma[idx(1,0,0)]=-M[0]/I[0];
    Gamma[idx(2,0,0)]=-M[1]/I[1];
    Gamma[idx(3,0,0)]=-M[2]/I[2];
  }
  static constexpr size_t idx(int i,int j,int k){ return 16*i+4*j+k; }
};

int main(){
  RigidBody rb;
  rb.u={1.0,0.3,0.4,9.0}; // u^0=1 keeps torque active
  rb.x.assign(4,0.0);

  std::ofstream log("omega.txt");
  double dt=1e-2;
  for(size_t n=0;n<500;++n){
    log<<rb.u[1]<<' '<<rb.u[2]<<' '<<rb.u[3]<<'\n';
    rb.step_3(dt);
  }
}
```

---

## 6  Limitations & future work

- No event detection (periapsis, horizon crossing, …)
- LU without pivoting — safe only if \$(\mathbf I+\Delta t,\Gamma u)\$ is diagonally dominant. fileciteturn0file1
- Currently single‑threaded; integrations are embarrassingly parallel.

---

© 2025 — MIT Licence.  Contributions welcome!

