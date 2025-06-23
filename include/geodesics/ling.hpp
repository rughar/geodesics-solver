#pragma once

// =============================================================================
//  FILE: ling.h  -  lightweight naive LU solver
// =============================================================================
//  Purpose :
//    * lu_naive  - factorises a diagonally dominant matrix A into L*U.
//    * fb_naive  - solves A*x = b in place on v by forward + backward
//      substitution using the factors produced by lu_naive.
//
//  Limitations :
//    - No pivoting (numeric stability relies on diagonal dominance).
//    - Works on user supplied containers that support A[i][j] style access
//      and are mutable for the LU stage.
// =============================================================================
 
// ---------------------------------------------------------------------------
//  lu_naive
// ---------------------------------------------------------------------------
//  In‑place LU factorisation of an  n*n  matrix  A  without pivoting.
//  After the call:
//      • U  is stored on and above the main diagonal.
//      • L  (without the unit diagonal) is stored below the diagonal.
//  PRECONDITION:  A  must be diagonally dominant so that no element on the
//  diagonal becomes zero.
// ---------------------------------------------------------------------------
template<class Um>
void lu_naive(const size_t n, Um& A)
{
  for (size_t i = 0; i < n; i++) {
    A[i][i] = 1 / A[i][i];
    for (size_t j = i + 1; j < n; j++) {
      A[j][i] *= A[i][i];
      for (size_t k = i + 1; k < n; k++)
        A[j][k] -= A[j][i] * A[i][k];
    }
  }
}

// ---------------------------------------------------------------------------
//  fb_naive
// ---------------------------------------------------------------------------
//  Forward + backward substitution for  A.x = v .
//  The matrix  A  must already be factorised by  lu_naive.
//  The solution overwrites  v.
// ---------------------------------------------------------------------------
template<class Um, class Uv>
void fb_naive(const size_t n, const Um& A, Uv& v)
{
  for (size_t i = 1; i < n; i++)
    for (size_t j = 0; j < i; j++)
      v[i] -= A[i][j] * v[j];
  for (size_t i = n; i--;) {
    for (size_t j = i + 1; j < n; j++)
      v[i] -= A[i][j] * v[j];
    v[i] *= A[i][i];
  }
}