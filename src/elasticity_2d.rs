//! 2D linear elasticity (plane stress) on P1 triangles.
//!
//! DOF layout: 2*k = u_x(node k), 2*k+1 = u_y(node k). Total DOFs = 2N.
//! Strain in Voigt notation: ε = [εxx, εyy, γxy]ᵀ (engineering shear).
//! See NOTES.md §15–§19 for derivations.

use crate::element_2d::shape_gradients;
use crate::mesh_2d::Mesh2D;
use nalgebra::{DVector, Matrix3, SMatrix};
use nalgebra_sparse::{CooMatrix, CsrMatrix, SparseEntryMut};

// ---------------------------------------------------------------------------
// Material
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
pub struct PlaneStress {
    pub e: f64,  // Young's modulus
    pub nu: f64, // Poisson's ratio
}

impl PlaneStress {
    /// Plane-stress constitutive matrix D (3×3) mapping ε → σ.
    pub fn d_matrix(&self) -> Matrix3<f64> {
        let c = self.e / (1.0 - self.nu * self.nu);
        Matrix3::new(
            c,
            c * self.nu,
            0.0,
            c * self.nu,
            c,
            0.0,
            0.0,
            0.0,
            c * (1.0 - self.nu) / 2.0,
        )
    }
}

// ---------------------------------------------------------------------------
// Element-level math
// ---------------------------------------------------------------------------

/// Strain–displacement matrix B (3×6) for a P1 triangle.
/// Maps element DOFs [u_x0, u_y0, u_x1, u_y1, u_x2, u_y2]ᵀ to strain.
pub fn strain_displacement(verts: [[f64; 2]; 3]) -> SMatrix<f64, 3, 6> {
    let (grads, _area) = shape_gradients(verts);
    let mut b = SMatrix::<f64, 3, 6>::zeros();
    for i in 0..3 {
        let dx = grads[i][0];
        let dy = grads[i][1];
        b[(0, 2 * i)] = dx;
        b[(1, 2 * i + 1)] = dy;
        b[(2, 2 * i)] = dy;
        b[(2, 2 * i + 1)] = dx;
    }
    b
}

/// Element stiffness: K_e = ∫_T Bᵀ D B dA. Constant integrand on a P1
/// triangle, so K_e = area · Bᵀ D B.
pub fn element_stiffness(verts: [[f64; 2]; 3], mat: &PlaneStress) -> SMatrix<f64, 6, 6> {
    let (_, area) = shape_gradients(verts);
    let b = strain_displacement(verts);
    let d = mat.d_matrix();
    area * b.transpose() * d * b
}

/// Element body-force vector: F_e[2i]   = ∫_T λ_i b_x dA,
///                            F_e[2i+1] = ∫_T λ_i b_y dA.
/// 1-point centroid quadrature.
pub fn element_body_force<F>(verts: [[f64; 2]; 3], body: &F) -> SMatrix<f64, 6, 1>
where
    F: Fn(f64, f64) -> [f64; 2],
{
    let (_, area) = shape_gradients(verts);
    let [p0, p1, p2] = verts;
    let cx = (p0[0] + p1[0] + p2[0]) / 3.0;
    let cy = (p0[1] + p1[1] + p2[1]) / 3.0;
    let [bx, by] = body(cx, cy);
    let scale = area / 3.0;
    let mut f = SMatrix::<f64, 6, 1>::zeros();
    for i in 0..3 {
        f[(2 * i, 0)] = scale * bx;
        f[(2 * i + 1, 0)] = scale * by;
    }
    f
}

// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------

/// Assemble global K (sparse 2N × 2N) and F (2N × 1) for plane-stress
/// elasticity with given material and body-force function.
pub fn assemble<B>(mesh: &Mesh2D, mat: &PlaneStress, body: &B) -> (CsrMatrix<f64>, DVector<f64>)
where
    B: Fn(f64, f64) -> [f64; 2],
{
    let n = mesh.num_nodes();
    let dof = 2 * n;
    let mut coo = CooMatrix::<f64>::new(dof, dof);
    let mut f_global = DVector::<f64>::zeros(dof);

    for t in 0..mesh.num_triangles() {
        let verts = mesh.triangle_coords(t);
        let local = mesh.triangles[t];

        let k_e = element_stiffness(verts, mat);
        let f_e = element_body_force(verts, body);

        let gdof = [
            2 * local[0],
            2 * local[0] + 1,
            2 * local[1],
            2 * local[1] + 1,
            2 * local[2],
            2 * local[2] + 1,
        ];

        for a in 0..6 {
            f_global[gdof[a]] += f_e[(a, 0)];
            for b in 0..6 {
                coo.push(gdof[a], gdof[b], k_e[(a, b)]);
            }
        }
    }

    (CsrMatrix::from(&coo), f_global)
}

// ---------------------------------------------------------------------------
// Boundary conditions
// ---------------------------------------------------------------------------

const PENALTY: f64 = 1e30;

/// Clamp a node: both u_x = u_y = 0. Penalty method on both DOFs.
pub fn clamp_node(k: &mut CsrMatrix<f64>, f: &mut DVector<f64>, node: usize) {
    pin_dof(k, f, 2 * node, 0.0);
    pin_dof(k, f, 2 * node + 1, 0.0);
}

/// Pin a single DOF (e.g. u_x at one node but free u_y) to a value.
pub fn pin_dof(k: &mut CsrMatrix<f64>, f: &mut DVector<f64>, dof: usize, value: f64) {
    match k.get_entry_mut(dof, dof) {
        Some(SparseEntryMut::NonZero(v)) => *v += PENALTY,
        _ => panic!("expected nonzero diagonal at DOF {dof}"),
    }
    f[dof] = PENALTY * value;
}

/// Apply a concentrated point load at a node.
pub fn apply_point_load(f: &mut DVector<f64>, node: usize, fx: f64, fy: f64) {
    f[2 * node] += fx;
    f[2 * node + 1] += fy;
}

/// Apply a constant traction `t = (tx, ty)` to every boundary edge passing
/// `edge_filter`. Linear (P1) edge lumping: each endpoint of an edge of
/// length `L` receives `t · L / 2`. Equivalent to ∫_edge λ_i · t dS for
/// linear shape functions.
pub fn apply_edge_traction<F>(
    mesh: &Mesh2D,
    f: &mut DVector<f64>,
    traction: [f64; 2],
    edge_filter: F,
) where
    F: Fn([f64; 2], [f64; 2]) -> bool,
{
    for (a, b) in mesh.boundary_edges() {
        let pa = mesh.nodes[a];
        let pb = mesh.nodes[b];
        if !edge_filter(pa, pb) {
            continue;
        }
        let len = ((pb[0] - pa[0]).powi(2) + (pb[1] - pa[1]).powi(2)).sqrt();
        let half = 0.5 * len;
        f[2 * a] += traction[0] * half;
        f[2 * a + 1] += traction[1] * half;
        f[2 * b] += traction[0] * half;
        f[2 * b + 1] += traction[1] * half;
    }
}

// ---------------------------------------------------------------------------
// Post-processing: strain, stress, von Mises
// ---------------------------------------------------------------------------

/// Strain [εxx, εyy, γxy] at triangle `t` (constant on P1 element).
pub fn strain_at_triangle(mesh: &Mesh2D, u: &[f64], t: usize) -> [f64; 3] {
    let verts = mesh.triangle_coords(t);
    let local = mesh.triangles[t];
    let b = strain_displacement(verts);
    let u_e = SMatrix::<f64, 6, 1>::from_column_slice(&[
        u[2 * local[0]],
        u[2 * local[0] + 1],
        u[2 * local[1]],
        u[2 * local[1] + 1],
        u[2 * local[2]],
        u[2 * local[2] + 1],
    ]);
    let eps = b * u_e;
    [eps[(0, 0)], eps[(1, 0)], eps[(2, 0)]]
}

/// Stress σ = D ε (plane stress).
pub fn stress_at_triangle(mat: &PlaneStress, strain: [f64; 3]) -> [f64; 3] {
    let d = mat.d_matrix();
    let eps = SMatrix::<f64, 3, 1>::from_column_slice(&strain);
    let s = d * eps;
    [s[(0, 0)], s[(1, 0)], s[(2, 0)]]
}

/// Plane-stress von Mises stress:
///   σ_vm = sqrt(σxx² - σxx·σyy + σyy² + 3·σxy²).
pub fn von_mises_2d(sigma: [f64; 3]) -> f64 {
    let [sxx, syy, sxy] = sigma;
    (sxx * sxx - sxx * syy + syy * syy + 3.0 * sxy * sxy).sqrt()
}

/// von Mises stress per triangle.
pub fn stress_field(mesh: &Mesh2D, u: &[f64], mat: &PlaneStress) -> Vec<f64> {
    (0..mesh.num_triangles())
        .map(|t| {
            let eps = strain_at_triangle(mesh, u, t);
            let sig = stress_at_triangle(mat, eps);
            von_mises_2d(sig)
        })
        .collect()
}

/// Reshape a flat 2N displacement vector into per-node [u_x, u_y] pairs.
pub fn displacement_to_pairs(u: &[f64]) -> Vec<[f64; 2]> {
    u.chunks_exact(2).map(|c| [c[0], c[1]]).collect()
}
