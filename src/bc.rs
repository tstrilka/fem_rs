use nalgebra::{DMatrix, DVector};
use nalgebra_sparse::{CsrMatrix, SparseEntryMut};

/// Apply a Dirichlet BC u(node) = value by row/column elimination on a dense
/// matrix.
///
/// Standard technique:
///   1. Subtract value * column(node) from F  (so other equations stay consistent)
///   2. Zero out row(node) and column(node)
///   3. Set K[node, node] = 1, F[node] = value
pub fn apply_dirichlet(k: &mut DMatrix<f64>, f: &mut DVector<f64>, node: usize, value: f64) {
    let n = k.nrows();

    for i in 0..n {
        if i != node {
            f[i] -= k[(i, node)] * value;
        }
    }

    for i in 0..n {
        k[(node, i)] = 0.0;
        k[(i, node)] = 0.0;
    }

    k[(node, node)] = 1.0;
    f[node] = value;
}

/// Apply a Dirichlet BC u(node) = value via the **penalty method** on a
/// sparse matrix.
///
/// Adds a large penalty `p` to the diagonal entry K[node, node] and sets
/// F[node] = p * value. The equation for row `node` becomes approximately
/// `p · u[node] = p · value`, so u[node] ≈ value.
///
/// Pros: one-line BC application, no structural changes to the sparsity
/// pattern (diagonal is already nonzero from assembly).
/// Cons: worse conditioning than row/col elimination. Acceptable for direct
/// or well-preconditioned solvers on Poisson-like problems.
pub fn apply_dirichlet_penalty(
    k: &mut CsrMatrix<f64>,
    f: &mut DVector<f64>,
    node: usize,
    value: f64,
) {
    const PENALTY: f64 = 1e30;
    match k.get_entry_mut(node, node) {
        Some(SparseEntryMut::NonZero(v)) => *v += PENALTY,
        _ => panic!(
            "expected nonzero diagonal at node {node} \
             (every assembled node should touch its own row)"
        ),
    }
    f[node] = PENALTY * value;
}
