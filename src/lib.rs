// 1D (Phase 1)
pub mod assembly;
pub mod bc;
pub mod element;
pub mod mesh;
pub mod solve;

// 2D (Phase 2)
pub mod assembly_2d;
pub mod element_2d;
pub mod gmsh;
pub mod mesh_2d;
pub mod viewer;
pub mod viz_2d;

// Phase 3 — linear elasticity
pub mod elasticity_2d;
