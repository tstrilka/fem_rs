# fem_rs — learning roadmap

Build a 2D FEM solver in Rust from scratch, with mesh I/O, custom mesher,
and interactive visualization. Pure learning project.

Start: 2026-05-15

## How to use this file

Check off boxes as you finish. Each phase ends in a runnable, working
checkpoint — feel free to stop after any phase and take the win.

## Sub-agents available

Project-scoped agents live in `.claude/agents/`. Invoke by mentioning the
name or by asking the matching question:

- **fem-math-tutor** — "why does X work?", "derive Y"
- **numerics-reviewer** — "review my element_stiffness", run after each TODO
- **fem-debugger** — "my test fails", "convergence ratio is wrong"

## Reading list

- **Larson & Bengzon**, *The Finite Element Method: Theory, Implementation,
  and Applications* — free PDF, primary reference. Ch. 1–3 cover Phase 1–2.
- **Hughes**, *The Finite Element Method* — denser, optional, great on
  isoparametric elements (Phase 3).
- **Shewchuk's "Triangle" paper** — Delaunay refinement (Phase 4).

---

## Phase 1 — 1D Poisson, no mesh, no GUI ✅ DONE

**Goal:** prove the FEM loop end-to-end on the simplest problem.

- [x] Project scaffold (`cargo new`, modules, deps, plotting)
- [x] `Mesh1D::uniform`
- [x] `solve_dense` via `nalgebra` LU
- [x] 1D demo (moved to `examples/poisson_1d.rs`)
- [x] `tests/poisson_1d.rs` O(h²) convergence test
- [x] `NOTES.md` derivations §1–§7
- [x] `element_stiffness(h)`
- [x] `element_load(x_left, x_right, f)`
- [x] `assemble`
- [x] `apply_dirichlet`
- [x] `cargo run --example poisson_1d` → max nodal error 4.019e-4
- [x] `cargo test` green (convergence ratio ≈ 4)

---

## Phase 2 — 2D Poisson on triangles, gmsh mesh 🟡 IN PROGRESS

**Goal:** real FEM on a real (loaded) mesh, with real visualization.

Milestone A — math working on a built-in mesh: ✅ DONE
- [x] `Mesh2D` struct: nodes, triangles, boundary_nodes (`mesh_2d.rs`)
- [x] Built-in `unit_square(nx, ny)` mesh generator
- [x] `shape_gradients(verts)` — closed-form `∇λ_i` and area
- [x] `viz_2d::render_field` — flat-shaded SVG of scalar field
- [x] `main.rs` 2D demo, sin(πx)sin(πy) manufactured solution
- [x] NOTES.md §8–§13
- [x] `element_stiffness(verts)`
- [x] `element_load(verts, f)`
- [x] `assemble` scatter
- [x] `cargo run` → max nodal error 1.339e-3 at 32×32 (matches §12)
- [x] `tests/poisson_2d.rs` — O(h²) convergence test green

Milestone B — gmsh + sparse + interactive:
- [x] Switch global K to sparse (`nalgebra-sparse` CSR via COO triplets)
- [x] Conjugate Gradient sparse solver (`solve_cg`)
- [x] Penalty-method Dirichlet BC for sparse (`apply_dirichlet_penalty`)
- [x] NOTES.md §13 — sparse + CG derivations
- [x] 64×64 mesh solves in 17ms (release): assemble 3.2ms, CG 13.3ms
- [x] gmsh `.msh` v4 ASCII parser (`src/gmsh.rs`)
- [x] Boundary nodes computed topologically (edge appearing once)
- [x] `tests/gmsh_parser.rs` — parses inline mesh, solves Poisson on it
- [x] `assets/lshape.geo` + `examples/poisson_msh.rs` ready for real meshes
- [x] L-shape demo: 1489 nodes, 2816 triangles, total < 5ms (parse+assemble+solve)
- [x] `egui` + `eframe` interactive viewer (`src/viewer.rs`, `src/bin/viewer.rs`)
- [x] Pan (drag), zoom (scroll, cursor-anchored), wireframe toggle, reset view, stats sidebar

**Done when:** can load `lshape.msh`, solve Poisson, view in an interactive
window.

**Time estimate:** 2–3 weekends.

---

## Phase 3 — physics that looks like something

### Option A: linear elasticity (plane stress) ✅ DONE
- [x] Vector-valued field (2 DOF per node)
- [x] Plane-stress constitutive matrix D (`src/elasticity_2d.rs`)
- [x] Strain–displacement matrix B (3×6)
- [x] Element K via `K_e = area · Bᵀ D B`
- [x] Sparse assembly + clamp/point-load BC helpers
- [x] Stress recovery + von Mises (per triangle)
- [x] Viewer extension: builder API, optional displacement + triangle field
- [x] Deformed-shape rendering with undeformed wireframe overlay
- [x] Cantilever demo (`src/bin/cantilever.rs`)
- [x] Convergence test against EB analytic — 128×16 hits 0.996 ratio
- [x] NOTES.md §15–§18 — elasticity derivations + locking discussion

Optional follow-ons (skipped for now):
- [ ] Distributed traction BCs via edge quadrature
- [ ] Body force demos (gravity)
- [ ] P2 elements (removes bending locking)

Additional cases built on Phase 3:
- [x] Brazilian disc compression (`src/bin/disc.rs`, `assets/disc.geo`)
       FEM σ at center within 0.5% of closed-form (σ_xx = +2P/(πD), σ_yy = -6P/(πD))
- [x] Plate with hole under tension (`src/bin/plate_hole.rs`, `assets/plate_hole.geo`)
       Added `Mesh2D::boundary_edges` and `apply_edge_traction` for distributed loads.
       FEM K_t = 3.37 (Kirsch infinite-plate K_t = 3.0, Howland finite-width ≈ 3.05)

### Option B: heat equation (time-dependent) — not started
### Option C: P2 elements — not started

**Time estimate:** 2–3 weekends. (Took ~3 hours conversation time.)

---

## Phase 4 — write your own mesher (optional, big chunk)

**Goal:** stop depending on gmsh; produce meshes from polygons.

- [ ] Half-edge / DCEL data structure
- [ ] Robust geometric predicates (use `robust` crate, or accept fragility)
- [ ] Bowyer–Watson incremental Delaunay
- [ ] Constrained Delaunay (boundary edge insertion)
- [ ] Ruppert's algorithm for quality refinement
- [ ] Replace gmsh in Phase 2 demo with your own mesh

**Time estimate:** ~1 month.

This is genuinely hard. The trap is non-robust predicates — when you do
floating-point orientation tests and miss edge cases, the algorithm
silently produces degenerate triangles. Read Shewchuk before starting.

---

## Phase 5 — stretch goals

- [ ] 3D tetrahedral meshes
- [ ] `wgpu`-based 3D renderer with proper camera, lighting
- [ ] Nonlinear material (Newton iteration, line search)
- [ ] Parallel assembly via `rayon`
- [ ] GPU-side assembly + solve for fun
- [ ] CLI: `fem-rs solve config.toml` with problem definition in TOML

---

## Crate choices (lock in by Phase 2)

- **Linear algebra/sparse solver:** `faer` (modern, fast, ergonomic). Fall
  back to `nalgebra-sparse` if you hit limitations.
- **2D GUI:** `egui` + `eframe`. Quick to set up, immediate-mode.
- **3D GUI (Phase 5):** `wgpu` direct, or `bevy` if you want batteries.
- **Static plots:** `plotters` with `svg_backend` (no system deps).
- **Mesh format:** roll your own gmsh parser (`mshio` exists if needed).
- **Robust predicates (Phase 4):** `robust` crate.

## Working agreements with myself

- One phase at a time. Don't skip ahead.
- After every TODO: run `cargo test`, then ask numerics-reviewer to look.
- Compare numerical results to an analytic solution at every phase.
- Commit at each green checkpoint so it's easy to roll back.
- It's a learning project — favor clarity over performance until Phase 5.
