use crate::element_2d::{element_load, element_stiffness};
use crate::mesh_2d::Mesh2D;
use nalgebra::DVector;
use nalgebra_sparse::{CooMatrix, CsrMatrix};

/// Assemble global stiffness K (sparse CSR) and load F (dense vector) for a
/// 2D Poisson problem on a triangular mesh.
///
/// Build a COO triplet list during assembly — `CooMatrix::push` allows
/// duplicates, and the COO→CSR conversion sums them automatically, which is
/// exactly the scatter semantics we want.
pub fn assemble<F: Fn(f64, f64) -> f64>(mesh: &Mesh2D, f: &F) -> (CsrMatrix<f64>, DVector<f64>) {
    let n = mesh.num_nodes();
    let mut coo = CooMatrix::<f64>::new(n, n);
    let mut f_global = DVector::<f64>::zeros(n);

    for t in 0..mesh.num_triangles() {
        let verts = mesh.triangle_coords(t);
        let local = mesh.triangles[t];

        let k_e = element_stiffness(verts);
        let f_e = element_load(verts, f);

        for a in 0..3 {
            f_global[local[a]] += f_e[a];
            for b in 0..3 {
                coo.push(local[a], local[b], k_e[(a, b)]);
            }
        }
    }

    (CsrMatrix::from(&coo), f_global)
}
