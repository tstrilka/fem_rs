use anyhow::Result;
use fem_rs::{
    elasticity_2d::{
        PlaneStress, apply_point_load, assemble, clamp_node, displacement_to_pairs, stress_field,
    },
    mesh_2d::{Mesh2D, unit_square},
    solve::solve_cg,
    viewer::{FemViewer, ViewerApp},
};
use std::time::Instant;

/// Cantilever beam: clamped at the left edge, point load at the tip.
///
/// Domain: rectangle [0, L] × [-h/2, h/2].
/// BCs: u_x = u_y = 0 on x = 0.
/// Load: total P in -y direction, distributed evenly across tip nodes (x = L).
/// Analytic (Euler-Bernoulli, unit thickness):  δ = P L³ / (3 E I),  I = h³/12.
fn rectangle_mesh(length: f64, height: f64, nx: usize, ny: usize) -> Mesh2D {
    let mut mesh = unit_square(nx, ny);
    for node in mesh.nodes.iter_mut() {
        node[0] *= length;
        node[1] = (node[1] - 0.5) * height;
    }
    mesh
}

fn main() -> Result<()> {
    let length: f64 = 8.0;
    let height: f64 = 1.0;
    let nx = 64;
    let ny = 8;

    let mesh = rectangle_mesh(length, height, nx, ny);
    let material = PlaneStress { e: 1.0e4, nu: 0.3 };
    let p_total: f64 = 1.0;

    let body = |_x: f64, _y: f64| [0.0_f64, 0.0_f64];

    let t0 = Instant::now();
    let (mut k, mut rhs) = assemble(&mesh, &material, &body);

    let mut clamped = 0;
    for (i, &[x, _]) in mesh.nodes.iter().enumerate() {
        if x.abs() < 1e-9 {
            clamp_node(&mut k, &mut rhs, i);
            clamped += 1;
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
    let t_assemble = t0.elapsed();

    let t1 = Instant::now();
    let u: Vec<f64> = solve_cg(&k, &rhs, 1e-12, 200_000)?
        .iter()
        .copied()
        .collect();
    let t_solve = t1.elapsed();

    let delta_fem: f64 =
        tip_nodes.iter().map(|&n| u[2 * n + 1]).sum::<f64>() / tip_nodes.len() as f64;
    let i_moment = height.powi(3) / 12.0;
    let delta_analytic = -p_total * length.powi(3) / (3.0 * material.e * i_moment);

    println!("cantilever: L={length}, h={height}, mesh = {nx}×{ny}");
    println!(
        "  #nodes={}, #tri={}, clamped={clamped}, tip nodes={}",
        mesh.num_nodes(),
        mesh.num_triangles(),
        tip_nodes.len()
    );
    println!(
        "  tip deflection:  FEM = {delta_fem:+.5e},  EB analytic = {delta_analytic:+.5e},  \
         ratio FEM/EB = {:.4}",
        delta_fem / delta_analytic
    );
    println!("  timings: assemble+BC = {t_assemble:?}, CG = {t_solve:?}");

    let stress = stress_field(&mesh, &u, &material);
    let disp = displacement_to_pairs(&u);

    let max_disp: f64 = disp
        .iter()
        .map(|d| (d[0] * d[0] + d[1] * d[1]).sqrt())
        .fold(0.0_f64, f64::max);
    let target_visual = 0.2 * height;
    let viz_scale = if max_disp > 0.0 {
        target_visual / max_disp
    } else {
        1.0
    };

    let viewer = FemViewer::new(mesh)
        .with_displacement(disp, viz_scale)
        .with_triangle_field(stress);

    let app = ViewerApp {
        label: format!(
            "Cantilever L={length}, h={height}\nE={:.0e}, ν={}\nP={p_total} at tip\nviz scale ≈ {viz_scale:.0}",
            material.e, material.nu
        ),
        viewer,
    };

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 600.0])
            .with_title("fem_rs cantilever"),
        ..Default::default()
    };

    eframe::run_native("cantilever", options, Box::new(|_cc| Ok(Box::new(app))))
        .map_err(|e| anyhow::anyhow!("eframe error: {e}"))?;

    Ok(())
}
