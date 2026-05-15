use anyhow::Result;
use fem_rs::{
    assembly_2d::assemble, bc::apply_dirichlet_penalty, mesh_2d::unit_square, solve::solve_cg,
    viz_2d::render_field,
};
use std::f64::consts::PI;
use std::time::Instant;

/// Problem:
///   -Δu(x, y) = 2π² sin(πx) sin(πy)   on (0, 1)²
///   u = 0 on ∂Ω
/// Exact solution: u(x, y) = sin(πx) sin(πy).
fn main() -> Result<()> {
    let nx = 64;
    let ny = 64;

    let t_start = Instant::now();
    let mesh = unit_square(nx, ny);

    let f = |x: f64, y: f64| 2.0 * PI * PI * (PI * x).sin() * (PI * y).sin();
    let (mut k, mut rhs) = assemble(&mesh, &f);
    let t_assemble = t_start.elapsed();

    for &node in &mesh.boundary_nodes {
        apply_dirichlet_penalty(&mut k, &mut rhs, node, 0.0);
    }

    let t_solve_start = Instant::now();
    let u = solve_cg(&k, &rhs, 1e-10, 10_000)?;
    let t_solve = t_solve_start.elapsed();

    let max_err = mesh
        .nodes
        .iter()
        .zip(u.iter())
        .map(|(&[x, y], &u_h)| (u_h - (PI * x).sin() * (PI * y).sin()).abs())
        .fold(0.0_f64, f64::max);
    println!(
        "nx = {nx}, ny = {ny}, #nodes = {}, nnz = {}, max nodal error = {max_err:.3e}",
        mesh.num_nodes(),
        k.nnz()
    );
    println!(
        "timings: assemble = {:?}, solve (CG) = {:?}",
        t_assemble, t_solve
    );

    render_field(&mesh, u.as_slice(), "out.svg", "2D Poisson, u = sin(πx)sin(πy)")?;
    println!("Wrote out.svg");

    Ok(())
}
