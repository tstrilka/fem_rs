use anyhow::{Result, anyhow};
use nalgebra::{DMatrix, DVector};
use nalgebra_sparse::CsrMatrix;

/// Solve K u = F using a dense LU. Used by the 1D demo where N is tiny.
pub fn solve_dense(k: &DMatrix<f64>, f: &DVector<f64>) -> Result<DVector<f64>> {
    k.clone()
        .lu()
        .solve(f)
        .ok_or_else(|| anyhow!("linear system has no solution (singular matrix?)"))
}

/// Conjugate Gradient solver for symmetric positive-definite sparse systems.
///
/// The canonical iterative solver for FEM on elliptic PDEs. For our Poisson
/// matrix (SPD after Dirichlet BCs), CG converges in O(√κ) iterations where
/// κ is the condition number — typically a few hundred for 32×32 meshes.
///
/// Reference: Shewchuk, *An Introduction to the Conjugate Gradient Method
/// Without the Agonizing Pain*.
///
/// Algorithm (no preconditioner):
///     r₀ = b - A x₀
///     p₀ = r₀
///     for k = 0, 1, …:
///         α_k = (r_k · r_k) / (p_k · A p_k)
///         x_{k+1} = x_k + α_k p_k
///         r_{k+1} = r_k - α_k A p_k
///         if ‖r_{k+1}‖ < tol: break
///         β_k = (r_{k+1} · r_{k+1}) / (r_k · r_k)
///         p_{k+1} = r_{k+1} + β_k p_k
pub fn solve_cg(
    a: &CsrMatrix<f64>,
    b: &DVector<f64>,
    tol: f64,
    max_iter: usize,
) -> Result<DVector<f64>> {
    assert_eq!(a.nrows(), a.ncols(), "CG requires a square matrix");
    assert_eq!(a.nrows(), b.len(), "matrix/RHS dimension mismatch");

    let n = b.len();
    let mut x = DVector::<f64>::zeros(n);
    let mut r = b - a * &x;
    let mut p = r.clone();
    let mut rs_old = r.dot(&r);

    let tol2 = tol * tol;
    if rs_old < tol2 {
        return Ok(x);
    }

    for _ in 0..max_iter {
        let ap: DVector<f64> = a * &p;
        let denom = p.dot(&ap);
        if denom.abs() < 1e-300 {
            return Err(anyhow!("CG breakdown: p·Ap ≈ 0"));
        }
        let alpha = rs_old / denom;
        x.axpy(alpha, &p, 1.0);
        r.axpy(-alpha, &ap, 1.0);
        let rs_new = r.dot(&r);
        if rs_new < tol2 {
            return Ok(x);
        }
        let beta = rs_new / rs_old;
        p = &r + beta * &p;
        rs_old = rs_new;
    }
    Err(anyhow!(
        "CG did not converge in {max_iter} iterations (final residual {})",
        rs_old.sqrt()
    ))
}
