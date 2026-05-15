# fem_rs

A 2D finite element method (FEM) solver written in Rust from scratch as a
learning project. Covers 1D Poisson, 2D Poisson on triangular meshes, and
2D linear elasticity (plane stress), with a gmsh `.msh` parser, sparse
assembly + conjugate-gradient solver, static SVG plotting, and an
interactive `egui` viewer.

See [`PLAN.md`](PLAN.md) for the phased roadmap and [`NOTES.md`](NOTES.md)
for the math derivations behind every implementation step.

## What's implemented

- **1D Poisson**, P1 elements, dense LU, manufactured-solution test
  with O(h²) convergence ratio.
- **2D Poisson** on triangles, P1 elements with closed-form barycentric
  gradients. Sparse assembly via COO → CSR, conjugate-gradient solve,
  penalty-method Dirichlet BCs.
- **gmsh `.msh` v4 ASCII parser** with topological boundary detection
  (edges appearing in exactly one triangle).
- **2D linear elasticity (plane stress)**, P1 triangles, Voigt-notation
  B-matrix, isotropic D matrix, sparse assembly, clamped / roller / point
  load / distributed-edge-traction BCs, per-triangle stress recovery and
  von Mises.
- **Interactive viewer** (`egui` + `eframe`): pan, cursor-anchored zoom,
  wireframe toggle, reset, optional deformation overlay and triangle
  field colouring.
- **Static SVG plots** via `plotters`.

## Quickstart

```sh
cargo test                                              # all unit + integration tests
cargo run --release                                     # 2D Poisson, sin(πx)sin(πy) on a 64×64 unit square; writes out.svg
cargo run --release --example poisson_1d                # 1D Poisson convergence demo
cargo run --release --example poisson_msh               # 2D Poisson on assets/lshape.msh
cargo run --release --bin viewer                        # interactive viewer, built-in mesh
cargo run --release --bin viewer assets/lshape.msh      # interactive viewer, loaded mesh
cargo run --release --bin cantilever                    # plane-stress cantilever vs Euler-Bernoulli
cargo run --release --bin disc                          # Brazilian disc compression
cargo run --release --bin plate_hole                    # plate with central hole under tension
```

## Layout

```
src/
  lib.rs                  module root
  mesh.rs, element.rs, assembly.rs, solve.rs, bc.rs    1D Poisson
  mesh_2d.rs, element_2d.rs, assembly_2d.rs            2D Poisson
  gmsh.rs                 .msh v4 ASCII parser
  elasticity_2d.rs        plane-stress elasticity (D, B, K_e, stress recovery)
  viz_2d.rs               static SVG field plot
  viewer.rs               egui interactive viewer
  bin/
    viewer.rs, cantilever.rs, disc.rs, plate_hole.rs
  main.rs                 default 2D Poisson demo
examples/
  poisson_1d.rs, poisson_msh.rs
tests/
  poisson_1d.rs, poisson_2d.rs, gmsh_parser.rs, elasticity.rs
assets/
  lshape.{geo,msh}, disc.{geo,msh}, plate_hole.{geo,msh}
```

## Verification results

- 1D Poisson, P1: max nodal error ratio ≈ 4.0 when halving `h` (matches
  the O(h²) prediction).
- 2D Poisson on 32×32 unit square, manufactured `u = sin(πx)sin(πy)`:
  max nodal error 1.34e-3.
- L-shape gmsh demo (1489 nodes, 2816 triangles): parse + assemble + CG
  solve in under 5 ms (release).
- Cantilever (L=8, h=1, E=10⁴, ν=0.3, P=1) versus Euler-Bernoulli
  `δ = -PL³/(3EI)`: FEM/EB = 0.996 at 128×16, asymptoting to 1.006 at
  256×32 (Timoshenko shear correction is the right physics; EB is the
  approximation).
- Brazilian disc, FEM stress at the centre within 0.5 % of the
  closed-form `σ_xx = +2P/(πD)`, `σ_yy = -6P/(πD)`.
- Plate with central hole under tension, FEM stress concentration
  K_t = 3.37 (Kirsch infinite-plate value 3.0, Howland finite-width
  ≈ 3.05).

## Dependencies

`nalgebra`, `nalgebra-sparse`, `plotters` (SVG backend only),
`eframe` / `egui`, `anyhow`. Dev: `approx`. Edition 2024.
