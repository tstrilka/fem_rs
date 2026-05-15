# 1D Poisson with linear FEM — derivations

Reference for the TODOs in `src/element.rs`, `src/assembly.rs`, `src/bc.rs`.
Notation follows Larson & Bengzon, *The Finite Element Method*, Ch. 2.

## 1. Strong form

Find `u(x)` on `(0, 1)` such that

    -u''(x) = f(x)        on (0, 1)
     u(0)   = u(1) = 0

## 2. Weak form

Multiply by a test function `v(x)` with `v(0) = v(1) = 0`, integrate over `(0, 1)`,
and integrate by parts:

    ∫₀¹ -u'' v dx = ∫₀¹ f v dx
    [-u' v]₀¹ + ∫₀¹ u' v' dx = ∫₀¹ f v dx

The boundary term vanishes because `v` is zero at both endpoints. So:

    a(u, v) := ∫₀¹ u' v' dx  =  ∫₀¹ f v dx  =: L(v)        for all admissible v

Galerkin's idea: pick a finite-dimensional space `V_h ⊂ V` and look for `u_h ∈ V_h`
with `a(u_h, v_h) = L(v_h)` for all `v_h ∈ V_h`.

## 3. P1 (linear) finite element space

Partition `[0, 1]` into `N` elements with nodes `0 = x₀ < x₁ < … < x_N = 1`.
Define hat functions `φ_i(x)` — piecewise linear, equal to 1 at node `i` and 0 at
every other node. Then `V_h = span{φ₀, …, φ_N}`, and any `u_h` can be written as

    u_h(x) = Σᵢ Uᵢ φᵢ(x)

where `Uᵢ` are the **nodal values** (= the unknowns we solve for).

Plug `u_h = Σⱼ Uⱼ φⱼ` and `v = φᵢ` into the weak form:

    Σⱼ Uⱼ ∫₀¹ φⱼ' φᵢ' dx = ∫₀¹ f φᵢ dx
            └────┬─────┘     └────┬────┘
              K_{ij}              F_i

That gives the linear system `K U = F`.

## 4. Element-local view

Both `K` and `F` are sums of element contributions. On element `e = [xₗ, xᵣ]`
with length `h = xᵣ - xₗ`, only the two hat functions `φ_ℓ` and `φ_r` are nonzero,
so the element only contributes a 2×2 block to `K` and a 2-vector to `F`.

Map to the reference element `ξ ∈ [0, 1]` via `x = xₗ + ξ h`. The two shape
functions are

    N₀(ξ) = 1 - ξ          N₁(ξ) = ξ

with derivatives w.r.t. `x`:

    dN₀/dx = -1/h          dN₁/dx = 1/h

### 4a. Element stiffness `K_e`

    K_e[a, b] = ∫_{xₗ}^{xᵣ} (dN_a/dx)(dN_b/dx) dx

For a P1 line element the integrands are constants, so each integral is just
that constant times `h`:

    K_e[0,0] = (-1/h)(-1/h) · h = 1/h
    K_e[0,1] = (-1/h)( 1/h) · h = -1/h
    K_e[1,0] = -1/h
    K_e[1,1] = 1/h

→  **K_e = (1/h) · [[1, -1], [-1, 1]]**

This is what `element_stiffness(h)` should return.

### 4b. Element load `F_e`

    F_e[a] = ∫_{xₗ}^{xᵣ} N_a(x) f(x) dx

For general `f` we use quadrature. The simplest accurate-enough choice for P1
is the **midpoint rule** (1-point Gauss on `[0, 1]` in reference coords):

    ∫_{xₗ}^{xᵣ} g(x) dx  ≈  h · g(x_mid),    x_mid = (xₗ + xᵣ) / 2

At the midpoint, `N₀ = N₁ = 1/2`, so

    F_e[0] ≈ h · (1/2) · f(x_mid) = (h/2) · f(x_mid)
    F_e[1] ≈ (h/2) · f(x_mid)

→  **F_e ≈ (h/2) · f(x_mid) · [1, 1]ᵀ**

Note: midpoint rule integrates linear `f` exactly. For higher-order `f` or P2
shape functions you'd move to 2-point Gauss.

## 5. Assembly (scatter)

Loop over elements `e = 0 .. N-1`. Element `e` has global node indices
`(i, j) = (e, e+1)`. Add the 2×2 block into the global matrix:

    K_global[i, i] += K_e[0, 0]
    K_global[i, j] += K_e[0, 1]
    K_global[j, i] += K_e[1, 0]
    K_global[j, j] += K_e[1, 1]
    F_global[i]    += F_e[0]
    F_global[j]    += F_e[1]

Same pattern in 2D: a triangle has 3 nodes, you scatter a 3×3 block.

## 6. Dirichlet BC by row/column elimination

We have `u(node) = value` as a constraint. Naive approach: drop that row and
column, modify the RHS. Easier in code: keep matrix size fixed, do this in place.

Given `K U = F` with constraint `U[k] = g`:

1. **Move the constraint into the RHS** so the other equations remain consistent:

       for each row i ≠ k:   F[i] -= K[i, k] · g

2. **Zero out row k and column k** of `K`.
3. **Set `K[k, k] = 1` and `F[k] = g`.**

Order matters — step 1 must happen before step 2, otherwise you've lost the
column you need to subtract.

Why this works: the new system says `1 · U[k] = g` for row `k`, and every other
row `i` has had the `K[i, k] · g` term moved across, so it now solves the
original equation with `U[k]` fixed.

## 7. Convergence (sanity check)

For a smooth solution and P1 elements:

    ‖u - u_h‖_{L∞}  =  O(h²)     (also true at nodes, and often "superconvergent")

So when you double `N` (halve `h`), the max nodal error should drop by ~4×.
That's exactly what `tests/poisson_1d.rs` checks: ratio between 3.5 and 4.5.

For the test problem `f = π² sin(πx)`, the exact solution is `u = sin(πx)`.
Expected max nodal error:

    N = 16  →  ~1.6e-3
    N = 32  →  ~4.0e-4
    N = 64  →  ~1.0e-4

## 8. Phase 2 — 2D Poisson on a triangular mesh

Strong form:

    -Δu(x, y) = f(x, y)   on Ω ⊂ R²
     u = 0  on ∂Ω

Weak form (same derivation as 1D, integration by parts on `∫∫ -Δu v dA`):

    a(u, v) := ∫∫_Ω ∇u · ∇v dA  =  ∫∫_Ω f v dA  =: L(v)

The boundary term vanishes for `v|_∂Ω = 0`.

## 9. P1 triangle: barycentric shape functions

For a triangle with vertices `p_0, p_1, p_2`, define barycentric coordinates
`λ_0, λ_1, λ_2` by

    λ_i(p_j) = δ_{ij}        λ_0 + λ_1 + λ_2 = 1

`λ_i` is affine in `(x, y)` → its gradient is **constant** on the triangle.

**Closed-form gradients.** Let `A` = signed area of the triangle:

    A = (1/2) · ((x_1 - x_0)(y_2 - y_0) - (x_2 - x_0)(y_1 - y_0))

Then

    ∇λ_0 = (1/(2A)) · (y_1 - y_2,  x_2 - x_1)
    ∇λ_1 = (1/(2A)) · (y_2 - y_0,  x_0 - x_2)
    ∇λ_2 = (1/(2A)) · (y_0 - y_1,  x_1 - x_0)

Geometric reading: `∇λ_i` is perpendicular to the edge opposite to vertex `i`,
pointing inward, scaled so `λ_i` rises from 0 (on that edge) to 1 (at vertex i).

This is implemented in `shape_gradients(verts)` — no TODO.

### 9a. Element stiffness `K_e`

    K_e[i, j] = ∫∫_T (∇λ_i · ∇λ_j) dA

Both gradients are constant on the triangle, so the integral is trivial:

    K_e[i, j] = A · (∇λ_i · ∇λ_j)

→ a 3×3 matrix per triangle. This is what `element_stiffness(verts)` should
return.

### 9b. Element load `F_e`

    F_e[i] = ∫∫_T λ_i · f(x, y) dA

For P1 the cheapest correct quadrature is **1-point Gauss at the centroid**:

    ∫∫_T g dA  ≈  A · g(centroid),    centroid = (p_0 + p_1 + p_2) / 3

At the centroid, `λ_0 = λ_1 = λ_2 = 1/3`. So

    F_e[i] ≈ A · (1/3) · f(centroid) = (A / 3) · f(centroid)

Exact for linear `f`, second-order accurate for general smooth `f` — enough
for P1. Bump to 3-point Gauss when you go to P2.

## 10. Assembly (2D scatter)

Same shape as 1D, just 3 local DOFs per element. For triangle `t` with node
indices `[n_0, n_1, n_2]`:

    for a in 0..3:
        F_global[n_a] += F_e[a]
        for b in 0..3:
            K_global[n_a, n_b] += K_e[a, b]

The pattern generalizes to any element type — just iterate over the local
DOFs and the connectivity table.

## 11. Dirichlet BCs (2D)

Identical to 1D. The `apply_dirichlet` function from Phase 1 already handles
the general case — just loop over `mesh.boundary_nodes` and call it for each.

## 12. Expected convergence

For P1 on a quasi-uniform mesh with smooth `u`:

    ‖u - u_h‖_{L∞}  =  O(h²)

For our manufactured solution `u = sin(πx)sin(πy)`:

    nx = ny = 16  →  ~5e-3
    nx = ny = 32  →  ~1.3e-3
    nx = ny = 64  →  ~3.2e-4

Doubling resolution should drop the error by ~4×.

## 13. Sparse storage and the conjugate gradient solver

For 2D Poisson on a regular triangulation, the global stiffness matrix `K`
has ~7 nonzeros per row (the node plus its ~6 neighbors). On a 64×64 mesh,
that's 29 057 nonzeros in a 4225×4225 matrix — only 0.16% dense. Dense
storage is wasteful; dense LU is O(n³) and quickly unusable.

### Assembly into a sparse matrix

Use a **COO (coordinate) triplet list** during assembly:

    for each element:
        for each (a, b) in 3×3:
            coo.push(local[a], local[b], K_e[a, b])

`push` allows duplicate `(row, col)` pairs. When you convert COO → CSR, the
duplicates are **summed**, which is exactly the scatter semantics we want
— `nalgebra_sparse::CsrMatrix::from(&coo)` does this in one call.

### Penalty method for Dirichlet BCs (sparse-friendly)

The row/col elimination from §6 needs random-access updates to a sparse
matrix — expensive. The **penalty method** is simpler:

For the constraint `u[k] = g`, add a huge value `P` (e.g. `1e30`) to the
diagonal and update the RHS:

    K[k, k]  += P
    F[k]      = P · g

The equation for row `k` is now approximately

    P · u[k] + (small terms) = P · g    →    u[k] ≈ g

The "small terms" don't get zeroed out, but they're dominated by `P`.

Pros: O(1) per BC, no structural change to the sparsity pattern.
Cons: condition number gets worse (κ scales with P). Fine for direct or
preconditioned iterative solvers on well-conditioned problems; not great
for stiff/ill-conditioned ones.

The cleaner alternative is **row/col elimination** that reduces the system
to free DOFs only. Slightly more code, no conditioning cost. Worth doing
if you hit numerical issues.

### Conjugate Gradient (CG)

For symmetric positive-definite systems — exactly the case for our Poisson
matrix — **CG** is the canonical iterative solver. The algorithm minimizes
the energy `½ x^T A x - b^T x` over a sequence of A-conjugate directions:

    r₀ = b - A x₀                       # initial residual
    p₀ = r₀                             # search direction
    for k = 0, 1, …:
        α_k = (r_k · r_k) / (p_k · A p_k)   # step length (line search exact)
        x_{k+1} = x_k + α_k p_k
        r_{k+1} = r_k - α_k A p_k
        if ‖r_{k+1}‖ < tol: done
        β_k = (r_{k+1} · r_{k+1}) / (r_k · r_k)
        p_{k+1} = r_{k+1} + β_k p_k     # A-conjugate to all previous p_j

Each iteration costs **one sparse matvec** `A p` plus a few dot products
and axpy operations — all O(nnz). Convergence: ‖e_k‖_A ≤ 2·((√κ - 1)/
(√κ + 1))^k · ‖e_0‖_A, where κ is the condition number. For our 64×64
Poisson, κ is mild and CG converges in a few hundred iterations.

For ill-conditioned problems you'd add a **preconditioner** M ≈ A⁻¹
(diagonal/Jacobi, incomplete Cholesky, multigrid). Save for later phases.

### Why CG over sparse LU here?

- LU/Cholesky: O(nnz · √n) work in 2D, O(n^(4/3)) in 3D — fast but factor
  storage gets fat (fill-in).
- CG: O(nnz · √κ) work, O(nnz) memory — leaner, scales better in 3D.

For this learning project CG is the better choice because it's
algorithmically transparent (~30 lines of code) and works on plain sparse
matrices without symbolic factorization machinery.

## 14. What still needs doing in Phase 2

- **gmsh `.msh` v4 ASCII parser** — replace the structured `unit_square`
  mesh with arbitrary meshes (e.g. L-shape).
- **Interactive viewer** — `egui` window with pan/zoom instead of static SVG.

Both are infrastructure — no new FEM theory required.

## 15. Phase 3 — 2D linear elasticity (plane stress)

### Setup

Unknown: displacement field **u**(x, y) = [u_x(x, y), u_y(x, y)]ᵀ.
2 DOFs per node → global system has size 2N.

Strain (small-strain, Voigt notation, **engineering** shear γ_xy = 2 ε_xy):

    ε_xx = ∂u_x/∂x
    ε_yy = ∂u_y/∂y
    γ_xy = ∂u_x/∂y + ∂u_y/∂x

    ε = [ε_xx, ε_yy, γ_xy]ᵀ

Constitutive law for **plane stress**, isotropic material (E, ν):

    σ = D · ε,    D = (E / (1 - ν²)) · [[1, ν,        0       ],
                                         [ν, 1,        0       ],
                                         [0, 0, (1-ν)/2 ]]

Equilibrium (strong form):    ∇·σ + b = 0   on Ω
plus Dirichlet BC `u = ū` on Γ_D and traction BC `σ·n = t̄` on Γ_N.

### Weak form

Multiply by a test displacement **v** that vanishes on Γ_D, integrate, and
use the divergence theorem:

    ∫_Ω σ(u)ᵀ ε(v) dA  =  ∫_Ω b·v dA  +  ∫_{Γ_N} t̄·v dS

In FEM-matrix form using the strain–displacement operator **B** that
maps element DOFs to strain:

    K_e u_e = f_e,    K_e = ∫_T Bᵀ D B dA,    f_e = ∫_T Nᵀ b dA + ∫_{∂T∩Γ_N} Nᵀ t̄ dS

### P1 triangle, the B matrix

DOF order: u_e = [u_x⁰, u_y⁰, u_x¹, u_y¹, u_x², u_y²]ᵀ (interleaved per node).

B is **3 × 6** (3 strain components × 6 element DOFs). Recall from §9 that
∇λ_i = (∂λ_i/∂x, ∂λ_i/∂y) is constant on a P1 triangle. For node i:

    column 2i   (u_x_i):    Bᵢ_x = [∂λ_i/∂x,      0,    ∂λ_i/∂y]ᵀ
    column 2i+1 (u_y_i):    Bᵢ_y = [     0, ∂λ_i/∂y,    ∂λ_i/∂x]ᵀ

Why those entries? `ε_xx = ∂u_x/∂x` picks up `∂λ_i/∂x` for the u_x DOF only;
`γ_xy = ∂u_x/∂y + ∂u_y/∂x` picks up `∂λ_i/∂y` for u_x and `∂λ_i/∂x` for u_y.

### Element stiffness, the simple way

Since both B and D are **constant** on a P1 triangle:

    K_e = ∫_T Bᵀ D B dA  =  area · Bᵀ D B           # 6×6 matrix

This is what `element_stiffness(verts, mat)` returns.

### Element body force

1-point Gauss at the centroid (`λ_i = 1/3` there):

    F_e[2i]   = (area/3) · b_x(centroid)
    F_e[2i+1] = (area/3) · b_y(centroid)

### Assembly and DOF mapping

For triangle `t` with node indices [n₀, n₁, n₂], the 6 global DOFs are

    [2n₀, 2n₀+1, 2n₁, 2n₁+1, 2n₂, 2n₂+1]

Scatter the 6×6 K_e into K_global using that index list. Same pattern as
Poisson, just doubled.

### Boundary conditions

- **Clamped node** (u_x = u_y = 0): penalty on both DOFs.
- **Roller** (u_x = 0, u_y free, or vice versa): penalty on a single DOF.
- **Point load** P at a node: add to RHS directly,
      F[2n]   += P_x
      F[2n+1] += P_y
- **Distributed traction** on Γ_N: edge-quadrature on each boundary segment.
  Not implemented here — for point loads only.

## 16. Stress recovery and von Mises

For each triangle, compute strain from element displacements:

    ε = B u_e

(B is constant on a P1 triangle, so strain is per-element, not per-vertex.)

Stress:    σ = D ε    →    [σ_xx, σ_yy, σ_xy].

**von Mises** scalar stress (for plane stress with σ_zz = 0):

    σ_vm = sqrt( σ_xx² - σ_xx σ_yy + σ_yy² + 3 σ_xy² )

Used as a single scalar "is the material in trouble here?" indicator.
Standard yield criterion: σ_vm > σ_yield → plastic deformation.

## 17. Cantilever beam verification

Test problem: clamped at x=0, point load P (downward) at the tip x=L.
**Euler-Bernoulli** (slender beam, ignores shear deformation):

    δ_EB = -P L³ / (3 E I),    I = h³/12 (unit thickness)

For our setup L=8, h=1, E=10⁴, ν=0.3, P=1:
    I = 1/12,    δ_EB = -64/(3·10⁴·1/12) = -64/2500 = -0.2048

### P1 bending locking — observed convergence

P1 triangles overestimate bending stiffness because their constant-strain
assumption can't represent linear strain gradients across the beam depth.
Measured FEM/EB ratios on this problem:

    nx × ny    FEM/EB ratio   notes
    32 ×  4    0.83           severe locking
    64 ×  8    0.96           mild locking
    128× 16    0.996          essentially converged
    256× 32    1.006          Timoshenko shear deformation kicks in,
                              FEM > EB by ~0.6 %

The 1.006 ratio at the finest mesh is the *correct* physics — plane stress
elasticity includes shear deformation that EB ignores. EB is the
approximation; FEM is closer to truth.

### Remedies for locking (not implemented here)

- **Quadratic (P2) triangles** — captures linear strain.
- **Selective reduced integration** — under-integrate the volumetric
  (incompressible) part. More effective for plane strain near ν → 0.5.
- **Mixed formulations** — separate displacement + stress unknowns.

For learning purposes, just refine the mesh. P1 + enough elements through
thickness is fine.

## 18. What still needs doing in Phase 3

- Distributed traction (Neumann) BCs via edge quadrature.
- Body force support (gravity demos).
- Heat equation (time-dependent) or P2 triangles, as alternative paths.

Phase-3 core (math + assembly + viewer + verified test) is done.
