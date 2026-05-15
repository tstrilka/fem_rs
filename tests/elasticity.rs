use fem_rs::{
    elasticity_2d::{PlaneStress, apply_point_load, assemble, clamp_node},
    mesh_2d::{Mesh2D, unit_square},
    solve::solve_cg,
};

fn rectangle_mesh(length: f64, height: f64, nx: usize, ny: usize) -> Mesh2D {
    let mut mesh = unit_square(nx, ny);
    for node in mesh.nodes.iter_mut() {
        node[0] *= length;
        node[1] = (node[1] - 0.5) * height;
    }
    mesh
}

/// Cantilever beam tip deflection vs Euler-Bernoulli analytic.
///
/// Mesh study (L=8, h=1, ν=0.3):
///   32×4:    FEM/EB = 0.83  (severe P1 bending locking)
///   64×8:    FEM/EB = 0.96  (mild locking)
///   128×16:  FEM/EB = 0.996 (essentially converged)
///   256×32:  FEM/EB = 1.006 (Timoshenko shear deformation kicks in,
///                            FEM correctly exceeds the simplified EB result)
///
/// We use 128×16 for the test — converged, runs fast in release.
#[test]
fn cantilever_tip_deflection() {
    let length = 8.0;
    let height = 1.0;
    let nx = 128;
    let ny = 16;

    let mesh = rectangle_mesh(length, height, nx, ny);
    let material = PlaneStress { e: 1.0e4, nu: 0.3 };
    let p_total = 1.0;

    let body = |_x: f64, _y: f64| [0.0_f64, 0.0_f64];
    let (mut k, mut rhs) = assemble(&mesh, &material, &body);

    for (i, &[x, _]) in mesh.nodes.iter().enumerate() {
        if x.abs() < 1e-9 {
            clamp_node(&mut k, &mut rhs, i);
        }
    }

    let tip_nodes: Vec<usize> = mesh
        .nodes
        .iter()
        .enumerate()
        .filter(|&(_, &[x, _])| (x - length).abs() < 1e-9)
        .map(|(i, _)| i)
        .collect();
    let load_per_node = -p_total / tip_nodes.len() as f64;
    for &n in &tip_nodes {
        apply_point_load(&mut rhs, n, 0.0, load_per_node);
    }

    let u = solve_cg(&k, &rhs, 1e-12, 200_000).unwrap();
    let delta_fem: f64 =
        tip_nodes.iter().map(|&n| u[2 * n + 1]).sum::<f64>() / tip_nodes.len() as f64;
    let i_moment = height.powi(3) / 12.0;
    let delta_eb: f64 = -p_total * length.powi(3) / (3.0 * material.e * i_moment);

    let ratio = delta_fem / delta_eb;
    assert!(
        ratio > 0.99 && ratio < 1.05,
        "tip deflection FEM/EB ratio {ratio:.4} outside [0.99, 1.05] \
         (FEM = {delta_fem:+.4e}, EB = {delta_eb:+.4e})"
    );
}
