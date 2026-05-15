use approx::assert_abs_diff_eq;
use fem_rs::{assembly::assemble, bc::apply_dirichlet, mesh::Mesh1D, solve::solve_dense};
use std::f64::consts::PI;

/// Convergence check: doubling N should cut the L∞ error by ~4 for P1 elements
/// (second-order convergence in mesh size h).
#[test]
fn poisson_converges_second_order() {
    let f = |x: f64| PI * PI * (PI * x).sin();
    let exact = |x: f64| (PI * x).sin();

    let solve_for = |n: usize| -> f64 {
        let mesh = Mesh1D::uniform(0.0, 1.0, n);
        let (mut k, mut rhs) = assemble(&mesh, &f);
        apply_dirichlet(&mut k, &mut rhs, 0, 0.0);
        apply_dirichlet(&mut k, &mut rhs, mesh.num_nodes() - 1, 0.0);
        let u = solve_dense(&k, &rhs).unwrap();
        mesh.nodes
            .iter()
            .zip(u.iter())
            .map(|(&x, &uh)| (uh - exact(x)).abs())
            .fold(0.0_f64, f64::max)
    };

    let e_coarse = solve_for(16);
    let e_fine = solve_for(32);
    let ratio = e_coarse / e_fine;
    assert!(
        ratio > 3.5 && ratio < 4.5,
        "expected ~4x error reduction (P1, h -> h/2), got {ratio:.2}"
    );
    assert_abs_diff_eq!(e_fine, 0.0, epsilon = 1e-2);
}
