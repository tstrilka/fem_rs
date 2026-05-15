use anyhow::Result;
use fem_rs::{
    elasticity_2d::{
        PlaneStress, apply_edge_traction, assemble, displacement_to_pairs, pin_dof,
        strain_at_triangle, stress_at_triangle,
    },
    gmsh::parse_msh,
    mesh_2d::Mesh2D,
    solve::solve_cg,
    viewer::{FemViewer, ViewerApp},
};
use std::time::Instant;

/// Rectangular plate with a circular hole under uniform tension.
///
/// Geometry: rectangle [-W/2, W/2] × [-H/2, H/2], hole of radius r at origin.
/// Loading: uniform traction σ_∞ in +x on the right edge.
/// BCs: u_x = 0 on the left edge (mirrors the loading); a single corner
/// node is also pinned in u_y to remove the remaining y-translation.
///
/// Kirsch's solution for an infinite plate under remote tension σ_∞:
///     σ_xx at (r, ±π/2)  =  +3 σ_∞      ← stress concentration K_t = 3
///     σ_xx at (±r, 0)    =  -σ_∞        ← compressive lobe along loading axis
///     σ_xx → σ_∞         as r/R → 0     ← far field
///
/// For our finite plate (d/H = 0.3), the true K_t is slightly above 3
/// (Howland correction). Expect FEM peak σ_xx/σ_∞ ≈ 3.0–3.2.
fn main() -> Result<()> {
    let half_w = 2.0_f64;
    let half_h = 0.5_f64;
    let sigma_inf = 1.0_f64;

    let msh_path = "assets/plate_hole.msh";
    if !std::path::Path::new(msh_path).exists() {
        eprintln!("Missing {msh_path}. Generate with:");
        eprintln!("  gmsh -2 -format msh4 -o assets/plate_hole.msh assets/plate_hole.geo");
        std::process::exit(1);
    }
    let mesh = parse_msh(msh_path)?;
    let material = PlaneStress { e: 1.0e4, nu: 0.3 };

    let t0 = Instant::now();
    let (mut k, mut rhs) = assemble(&mesh, &material, &|_, _| [0.0, 0.0]);

    // BC: u_x = 0 on the entire left edge.
    let mut left_count = 0;
    let mut left_nodes: Vec<usize> = Vec::new();
    for (i, &[x, _]) in mesh.nodes.iter().enumerate() {
        if (x + half_w).abs() < 1e-6 {
            pin_dof(&mut k, &mut rhs, 2 * i, 0.0);
            left_count += 1;
            left_nodes.push(i);
        }
    }
    // Additionally pin u_y at the bottom-left corner to remove y-translation.
    let corner = *left_nodes
        .iter()
        .min_by(|&&a, &&b| {
            mesh.nodes[a][1]
                .partial_cmp(&mesh.nodes[b][1])
                .unwrap()
        })
        .expect("at least one left-edge node");
    pin_dof(&mut k, &mut rhs, 2 * corner + 1, 0.0);

    // Loading: uniform σ_inf traction on the right edge.
    apply_edge_traction(&mesh, &mut rhs, [sigma_inf, 0.0], |pa, pb| {
        (pa[0] - half_w).abs() < 1e-6 && (pb[0] - half_w).abs() < 1e-6
    });
    let t_asm = t0.elapsed();

    let t1 = Instant::now();
    let u: Vec<f64> = solve_cg(&k, &rhs, 1e-12, 200_000)?
        .iter()
        .copied()
        .collect();
    let t_solve = t1.elapsed();

    // Per-triangle σ_xx.
    let s_xx_field: Vec<f64> = (0..mesh.num_triangles())
        .map(|t| {
            let eps = strain_at_triangle(&mesh, &u, t);
            stress_at_triangle(&material, eps)[0]
        })
        .collect();

    // Peak σ_xx (we expect it near the hole top/bottom, y = ±r).
    let (peak_t, &peak_val) = s_xx_field
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap();
    let peak_centroid = centroid(&mesh, peak_t);
    let k_t = peak_val / sigma_inf;

    // Stress at the top of the hole, picked as the triangle whose centroid
    // is closest to (0, r) where r is the hole radius — computed below.
    println!("plate with hole: half-W={half_w}, half-H={half_h}, σ_∞={sigma_inf}");
    println!(
        "  #nodes={}, #tri={}, left edge nodes={left_count}",
        mesh.num_nodes(),
        mesh.num_triangles()
    );
    println!("  timings: assemble+BC = {t_asm:?}, CG = {t_solve:?}");
    println!(
        "  peak σ_xx = {peak_val:+.4} at centroid {:?}",
        peak_centroid
    );
    println!("  stress concentration K_t = peak σ_xx / σ_∞ = {k_t:.3}");
    println!("  (Kirsch infinite-plate prediction: K_t = 3.000)");

    // Clamp the color range for visualization: -σ_∞ to 3.5 σ_∞ covers the
    // Kirsch field nicely.
    let s_xx_clamped: Vec<f64> = s_xx_field
        .iter()
        .map(|&v| v.clamp(-sigma_inf, 3.5 * sigma_inf))
        .collect();

    let disp = displacement_to_pairs(&u);
    let max_disp: f64 = disp
        .iter()
        .map(|d| (d[0] * d[0] + d[1] * d[1]).sqrt())
        .fold(0.0_f64, f64::max);
    let viz_scale = if max_disp > 0.0 {
        0.1 * half_h / max_disp
    } else {
        1.0
    };

    let viewer = FemViewer::new(mesh)
        .with_displacement(disp, viz_scale)
        .with_triangle_field(s_xx_clamped);

    let app = ViewerApp {
        label: format!(
            "Plate with hole, uniaxial tension\nσ_∞ = {sigma_inf}, E = {:.0e}, ν = {}\nK_t (FEM) ≈ {k_t:.2}\nfield: σ_xx",
            material.e, material.nu
        ),
        viewer,
    };

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 500.0])
            .with_title("fem_rs plate_hole"),
        ..Default::default()
    };

    eframe::run_native(
        "plate_hole",
        options,
        Box::new(|_cc| Ok(Box::new(app))),
    )
    .map_err(|e| anyhow::anyhow!("eframe error: {e}"))?;

    Ok(())
}

fn centroid(mesh: &Mesh2D, t: usize) -> [f64; 2] {
    let [p0, p1, p2] = mesh.triangle_coords(t);
    [
        (p0[0] + p1[0] + p2[0]) / 3.0,
        (p0[1] + p1[1] + p2[1]) / 3.0,
    ]
}
