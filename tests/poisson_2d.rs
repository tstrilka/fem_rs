use approx::assert_abs_diff_eq;
use fem_rs::{
    assembly_2d::assemble, bc::apply_dirichlet_penalty, mesh_2d::unit_square, solve::solve_cg,
};
use std::f64::consts::PI;

/// 2D convergence: P1 on quasi-uniform mesh → O(h²) in L∞.
/// Doubling resolution should drop max nodal error by ~4×.
#[test]
fn poisson_2d_converges_second_order() {
    let f = |x: f64, y: f64| 2.0 * PI * PI * (PI * x).sin() * (PI * y).sin();
    let exact = |x: f64, y: f64| (PI * x).sin() * (PI * y).sin();

    let solve_for = |n: usize| -> f64 {
        let mesh = unit_square(n, n);
        let (mut k, mut rhs) = assemble(&mesh, &f);
        for &node in &mesh.boundary_nodes {
            apply_dirichlet_penalty(&mut k, &mut rhs, node, 0.0);
        }
        let u = solve_cg(&k, &rhs, 1e-10, 100_000).unwrap();
        mesh.nodes
            .iter()
            .zip(u.iter())
            .map(|(&[x, y], &uh)| (uh - exact(x, y)).abs())
            .fold(0.0_f64, f64::max)
    };

    let e_coarse = solve_for(32);
    let e_fine = solve_for(64);
    let ratio = e_coarse / e_fine;
    assert!(
        ratio > 3.5 && ratio < 4.5,
        "expected ~4x error reduction (P1, h -> h/2), got {ratio:.2} (coarse={e_coarse:.3e}, fine={e_fine:.3e})"
    );
    assert_abs_diff_eq!(e_fine, 0.0, epsilon = 5e-3);
}
