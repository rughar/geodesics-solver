# Geodesic‑type ODE Integrator

A header‑only C++ library for **symplectic, time‑reversible** integration of ordinary differential equations that can be cast as a *quadratic‑in‑velocity* geodesic equation.

$$
\ddot x^i + {\Gamma^{i}}_{jk}\\dot x^j \dot x^k = 0,\qquad i = 0,\dots,n-1 .
$$

The core algorithm is the velocity **Störmer–Verlet** step, lifted to higher orders by the 4th and 6th‑order **Yoshida** compositions [[Yoshida 1990](https://doi.org/10.1016/0375-9601\(90\)90092-3)].

---

## 1 Mathematical background

### 1.1 Geodesic equation & related systems

Many seemingly unrelated problems share the geodesic algebraic structure, e.g.

- **Free rigid‑body rotation** (Euler equations in principal axes),
- **Planar Kepler motion** in *spherical* coordinates,
- **Charged particle in an electromagnetic field** (§ 1.2).

Writing \$u^i \equiv \dot x^i\$, the first‑order form is

$$
\dot x^i = u^i, \qquad \dot u^i = -{\Gamma^{i}}_{jk}u^j u^k.
$$

The split \$(x,u)\$ system is **symplectic** and **time‑reversible** – perfect for drift–kick–drift (Verlet) integrators.

### 1.2 Electromagnetic field as an \$(n+1)\$‑dimensional geodesic

Add a cyclic coordinate \$x^c\$ with fixed velocity \$u^c=1\$ and set

| Symbol                  | Role                                         |
| ----------------------- | -------------------------------------------- |
| \${\Gamma^{c}}_{ij}=0\$ | extra dimension remains unchanged            |
| \${\Gamma^{i}}_{cc}\$   | electric‑field part                          |
| \${\Gamma^{i}}_{ck}\$   | magnetic term \$\mathbf u\times\mathbf B\$   |

For purely linear Lorentz motion the classic **Boris** leap‑frog is simpler [[Boris](https://www.sciencedirect.com/science/article/abs/pii/0375960190900923)]. This section demonstrates how to include Lorentz force into the system with possibly additional **quadratic** term (geodesic‑like) correction, eg. geodesic equation for particle in additional electromagnetic field.

### 1.3 Norms and invariants

Typical invariants to monitor

$$
 g_{ij}u^i u^j = \text{const}\quad(=-1 \text{ for time‑like}),
$$

plus specific energy \$E\$ and angular momentum \$L\$ in stationary or axisymmetric metrics.

---

## 2 Integrator hierarchy

| Order | Method                      | Entry‑point  |
| ----- | --------------------------- | ------------ |
| 2     | velocity Störmer–Verlet     | `step_1(dt)` |
| 4     | Yoshida quartic composition | `step_2(dt)` |
| 6     | Yoshida sextic composition  | `step_3(dt)` |

`step_2` and `step_3` are wrappers that repeatedly call `step_1` with Yoshida coefficients.

---

## 3 Generalised Störmer–Verlet scheme

### 3.1 Classic velocity‑Verlet

$$
\begin{aligned}
\mathbf x_{1/2} &= \mathbf x_0 + \tfrac{\Delta t}{2}\\mathbf u_0,\\
\mathbf u_1 &= \mathbf u_0 + \Delta t \mathbf f(\mathbf x_{1/2},\mathbf u_0,\mathbf u_1),\\
\mathbf x_1 &= \mathbf x_{1/2} + \tfrac{\Delta t}{2}\\mathbf u_1.
\end{aligned}
$$

Time symmetry requires the kick \$\mathbf f\$ to use the same midpoint position \$\mathbf x\_{1/2}\$. The choise of velocity dependence between \$\mathbf u_0\$ and \$\mathbf u_1\$ may vary. For symmetric scheme we have two common choices

$$
\begin{aligned}
\mathbf f_A &= \ddot{\mathbf x}\bigl(\mathbf x_{1/2},\tfrac12(\mathbf u_0+\mathbf u_1)\bigr),\\
\mathbf f_B &= \tfrac12\Bigl[\ddot{\mathbf x}(\mathbf x_{1/2},\mathbf u_0)+\ddot{\mathbf x}(\mathbf x_{1/2},\mathbf u_1)\Bigr].
\end{aligned}
$$

Choose a linear combination that annihilates the quadratic dependence on \$\mathbf u\_1\$

$$
 \mathbf f = 2\mathbf f_A - \mathbf f_B.
$$

For geodesic acceleration this becomes

$$
f^i = -{\Gamma^{i}}_{jk} u_0^{j} u_1^{k}
$$

leading to

$$
\boxed{%
\begin{aligned}
\mathbf x_{1/2} &= \mathbf x_0 + \tfrac{\Delta t}{2}\\mathbf u_0,\\
\mathbf u_1 &= \mathbf u_0 - \Delta t \Gamma^i_{jk} u_0^{j} u_1^{k} ,\\
\mathbf x_1 &= \mathbf x_{1/2} + \tfrac{\Delta t}{2}\\mathbf u_1.
\end{aligned}
}
$$

Only **linear** in the unknown \$u\_1\$ → one LU solve per step. Complexity per step: \$\mathcal O(n^3)\$.

---

## 4 Adaptive step‑size workflow

Estimate stiffness via the largest eigenvalue \$\lambda\_{\max}\$ of the Jacobian \$J\$:

$$
h \approx \frac{1}{\lambda_{\max}}, \qquad \Delta t_{\text{new}} = 2h - \Delta t_{\text{old}}.
$$

Because the same \$h\$ is reused in both halves of the drift‑kick‑drift palindrome, reversibility is preserved. Our solver itselef contains fast (however non symmetry perserving) method for suggesting stepsize:

```cpp
suggest_stepsize_*(x_ref, u_ref, tol);  // * = 1,2,3
```

**Typical workflow**

1. **Probe run** – `step_3(dt,true)`; record suggested \$\Delta t\$ via `suggest_stepsize_*(...)` subroutine.
2. **Fit** – build \$(\mathbf x,\mathbf u) \mapsto \Delta t\$.
3. **Production run** – use that map in the symmetric update above.

---

## 5 Quick start – Euler rigid‑body with constant torque

```cpp
#include <geodesics/solver.hpp>
#include <array>
#include <fstream>

struct RigidBody : public StormerVerletCore<double> {
  using base = StormerVerletCore<double>;

  std::array<double, 3> I{ 2.0, 1.0, 0.5 };  // principal moments
  std::array<double, 3> M{ 0.0, 0.0, 0.1 };  // constant body torque

  RigidBody() { set_dimension(4); }

  void Christoffel_symbols() override {
    // quadratic ω×Iω terms
    base::Gamma[1][2][3] = base::Gamma[1][3][2] = -(I[1] - I[2]) / (2 * I[0]);
    base::Gamma[2][3][1] = base::Gamma[2][1][3] = -(I[2] - I[0]) / (2 * I[1]);
    base::Gamma[3][1][2] = base::Gamma[3][2][1] = -(I[0] - I[1]) / (2 * I[2]);
    // constant torque via u^0
    base::Gamma[1][0][0] = -M[0] / I[0];
    base::Gamma[2][0][0] = -M[1] / I[1];
    base::Gamma[3][0][0] = -M[2] / I[2];
  }
};

int main() {
  RigidBody rb;
  rb.u = { 1.0, 0.3, 0.4, 9.0 };   // u^0 = 1 keeps torque active
  rb.x = { 0.0, 0.0, 0.0, 0.0 };   // Initial positions (not used in this example)

  std::ofstream out("omega.txt");
  const double dt = 1.0e-2;
  for (size_t n = 0; n < 500; ++n) {
    out << rb.u[1] << ' ' << rb.u[2] << ' ' << rb.u[3] << '\n';
    rb.step_3(dt);
  }
  return 0;
}
```

