/* $Id: shoot_regularisers.h 6799 2016-05-20 16:50:25Z john $ */
/* (c) John Ashburner (2011) */

/** \function trapprox
 * Pointwise approximation of the trace of (L + H), where L is the penalty 
 * operator and H is a symmetric tensor. This trace is approximated by 
 * Tr(L + H) = Tr(diag(L) + H).
 *
 * \param[in]  dm  Dimension of the lattice
 * \param[in]  a   Symmetric tensor field
 * \param[out] s   Allocated array of dimension dm where to store the 
 *                 output traces.
 *
 * \note As usual, all arrays ar F-ordered with 0-indexing, meaning that
 * a(i, j, k, r, c) = a + i + dm[0]*(j + dm[1]*(k + dm[2]*(r + dm[3]*c)))
 */
extern double trapprox(mwSize dm[], float a[], double s[]);

/** \function vel2mom
 * Computes the momentum based on the velocity field and the regularization 
 * matrix (m = L * v)
 *
 * \param[in]  dm     Dimensions of the velocity lattice (voxels)
 * \param[in]  f      Velocities computed on the above lattice
 * \param[in]  s[0:2] Voxel size of the lattice (voxels/mm)
 * \param[in]  s[3]   Parameter of the absolute displacement penalty
 * \param[in]  s[4]   Parameter of the membrane energy (penalizes 
 *                    elements of the Jacobian matrix -> 1st order 
 *                    smoothness)
 * \param[in]  s[5]   Parameter of the bending energy (penalizes 
 *                    elements of the Hessian matrix -> 2nd order 
 *                    smoothness)
 * \param[in]  s[6]   Parameter of the linear elastic energy (penalizes 
 *                    elements of the symmetric part of the Jacobian matrix 
 *                    -> penalizes scaling and shearing)
 * \param[in]  s[7]   Parameter of the linear elastic energy (penalizes 
 *                    the divergence of the Jacobian matrix -> preserves 
 *                    volumes)
 * \param[out] g      Allocated array in which to store the output momentum 
 *
 * \note As usual, all arrays ar F-ordered with 0-indexing, meaning that
 * a(i, j, k, r, c) = a + i + dm[0]*(j + dm[1]*(k + dm[2]*(r + dm[3]*c)))
 */
extern void vel2mom(mwSize dm[], float f[], double s[], float g[]);

/** \function relax
 * Relaxation iterations when the matrix to invert is known to decompose 
 * locally into L + H, where L is the penalty operator (which can be 
 * implemented as a convolution) and H is a symmetric tensor field.
 * . We solve locally for (L + H) * u = b
 * . F = nondiag(L) (i.e., L without its diagonal elements)
 * . E = H + diag(L) + sI (to insure diagonal dominance)
 * . u = E^{-1} * ( b - F * u )
 * Updates are performed in place (Gauss-Siediel) rather than after each 
 * iteration (Jacobi).
 *
 * \param[in]    dm     Dimension of the lattice.
 * \param[in]    a      (Optional) Symmetric tensor field H (i.e. a 3x3   
 *                      symmetric matrix at each point of the lattice.).
 * \param[in]    b      Point at which to solve the system.
 * \param[in]    s[0:2] Voxel size of the lattice (voxels/mm)
 * \param[in]    s[3]   Parameter of the absolute displacement penalty
 * \param[in]    s[4]   Parameter of the membrane energy (penalizes 
 *                      elements of the Jacobian matrix -> 1st order 
 *                      smoothness)
 * \param[in]    s[5]   Parameter of the bending energy (penalizes 
 *                      elements of the Hessian matrix -> 2nd order 
 *                      smoothness)
 * \param[in]    s[6]   Parameter of the linear elastic energy (penalizes 
 *                      elements of the symmetric part of the Jacobian  
 *                      matrix -> penalizes scaling and shearing)
 * \param[in]    s[7]   Parameter of the linear elastic energy (penalizes 
 *                      the divergence of the Jacobian matrix -> preserves 
 *                      volumes)
 * \param[in]    nit    Number of relaxation iterations.
 * \param[inout] u      [in] Initial guess for the solution.
 *                      [out] Output relaxation solution.
 *
 * \note As usual, all arrays ar F-ordered with 0-indexing, meaning that
 * a(i, j, k, r, c) = a + i + dm[0]*(j + dm[1]*(k + dm[2]*(r + dm[3]*c)))
 */
extern void relax(mwSize dm[], float a[], float b[], double s[], int nit, float u[]);

/** \function Atimesp
 * Pointwise matrix multiplication in the momentum space.
 * Returns A[i] * vel2mom(p[i]), where A contains M symmetric 3*3 matrix 
 * (tensor field) and p contains M 3*1 vectors (velocity field).
 * 
 * \param[in]  dm  Dimension of the lattice (M = prod(dm))
 * \param[in]  A   Symmetric tensor field
 * \param[in]  p   Vector field
 * \param[out] Ap  Result of the pointwise matrix multiplication (A * p)
 *
 * \note As usual, all arrays ar F-ordered with 0-indexing, meaning that
 * a(i, j, k, r, c) = a + i + dm[0]*(j + dm[1]*(k + dm[2]*(r + dm[3]*c)))
 */
extern void Atimesp(mwSize dm[], float A[], double param[], float p[], float Ap[]);

/** \function sumsq
 * Returns the sum of squares of (b - (L + H) * u)
 *
 * \param[in]    a      (Optional) Symmetric tensor field H (i.e. a 3x3   
 *                      symmetric matrix at each point of the lattice.).
 * \param[in]    b      Point at which to solve the system.
 * \param[in]    s[0:2] Voxel size of the lattice (voxels/mm)
 * \param[in]    s[3]   Parameter of the absolute displacement penalty
 * \param[in]    s[4]   Parameter of the membrane energy (penalizes 
 *                      elements of the Jacobian matrix -> 1st order 
 *                      smoothness)
 * \param[in]    s[5]   Parameter of the bending energy (penalizes 
 *                      elements of the Hessian matrix -> 2nd order 
 *                      smoothness)
 * \param[in]    s[6]   Parameter of the linear elastic energy (penalizes 
 *                      elements of the symmetric part of the Jacobian  
 *                      matrix -> penalizes scaling and shearing)
 * \param[in]    s[7]   Parameter of the linear elastic energy (penalizes 
 *                      the divergence of the Jacobian matrix -> preserves 
 *                      volumes)
 * \param[in]    u      Approximate solution.
 *
 * \note As usual, all arrays ar F-ordered with 0-indexing, meaning that
 * a(i, j, k, r, c) = a + i + dm[0]*(j + dm[1]*(k + dm[2]*(r + dm[3]*c)))
 */
extern double sumsq(mwSize dm[], float a[], float b[], double s[], float u[]);

/** \function kernel
 * Returns the velocity-to-momentum convolution kernel.
 *
 * \param[in]  dm     Dimensions of the velocity lattice (voxels)
 * \param[in]  s[0:2] Voxel size of the lattice (voxels/mm)
 * \param[in]  s[3]   Parameter of the absolute displacement penalty
 * \param[in]  s[4]   Parameter of the membrane energy (penalizes 
 *                    elements of the Jacobian matrix -> 1st order 
 *                    smoothness)
 * \param[in]  s[5]   Parameter of the bending energy (penalizes 
 *                    elements of the Hessian matrix -> 2nd order 
 *                    smoothness)
 * \param[in]  s[6]   Parameter of the linear elastic energy (penalizes 
 *                    elements of the symmetric part of the Jacobian matrix 
 *                    -> penalizes scaling and shearing)
 * \param[in]  s[7]   Parameter of the linear elastic energy (penalizes 
 *                    the divergence of the Jacobian matrix -> preserves 
 *                    volumes)
 * \param[out] f      Allocated array of size dm in which to store the 
 *                    kernel
 *
 * \note As usual, all arrays ar F-ordered with 0-indexing, meaning that
 * a(i, j, k, r, c) = a + i + dm[0]*(j + dm[1]*(k + dm[2]*(r + dm[3]*c)))
 */
extern void kernel(mwSize dm[], double s[], float f[]);

