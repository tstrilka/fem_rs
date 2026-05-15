use crate::element::{element_load, element_stiffness};
use crate::mesh::Mesh1D;
use nalgebra::{DMatrix, DVector};

/// Assemble global stiffness K and load F by scattering element contributions.
///
/// TODO:
///   Loop over each element, compute element_stiffness(h) and element_load(x0, x1, f),
///   then add into K[i, j] and F[i] using the connectivity from mesh.element_nodes(e).
///
/// For 1D with sequential numbering, element e has nodes (e, e+1).
/// In 2D this becomes a triangle with 3 nodes — same scatter pattern.
pub fn assemble<F: Fn(f64) -> f64>(mesh: &Mesh1D, f: &F) -> (DMatrix<f64>, DVector<f64>) {
    let n = mesh.num_nodes();
    let mut k_global = DMatrix::<f64>::zeros(n, n);
    let mut f_global = DVector::<f64>::zeros(n);

    for e in 0..mesh.num_elements() {
        let (i, j) = mesh.element_nodes(e);
        let h = mesh.element_length(e);
        let x_left = mesh.nodes[i];
        let x_right = mesh.nodes[j];

        let k_e = element_stiffness(h);
        let f_e = element_load(x_left, x_right, f);

        let local = [i, j];
        for a in 0..2 {
            f_global[local[a]] += f_e[a];
            for b in 0..2 {
                k_global[(local[a], local[b])] += k_e[(a, b)];
            }
        }
    }

    (k_global, f_global)
}
