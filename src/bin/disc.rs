use anyhow::Result;
use fem_rs::{
    elasticity_2d::{
        PlaneStress, apply_point_load, assemble, clamp_node, displacement_to_pairs, pin_dof,
        strain_at_triangle, stress_at_triangle,
    },
    gmsh::parse_msh,
    mesh_2d::Mesh2D,
    solve::solve_cg,
    viewer::{FemViewer, ViewerApp},
};
use std::f64::consts::PI;
use std::time::Instant;

/// Brazilian (diametral compression) test on a disc.
///
/// Domain: disc of radius R centered at origin.
/// Loading: downward point load P at top (0, +R). Bottom (0, -R) is pinned.
/// To prevent rotation: u_x = 0 at the top node as well.
///
/// Analytical stress at the center (closed form for an infinite-thickness disc):
///     σ_xx(0, 0) = +2P / (π D)   (tensile)
///     σ_yy(0, 0) = -6P / (π D)   (compressive)
///     σ_xy(0, 0) = 0
/// where D = 2R is the diameter.
///
/// The tensile σ_xx at the center is the whole point: brittle materials
/// fail in tension first, so a compression test indirectly measures
/// tensile strength.
fn main() -> Result<()> {
    let radius = 0.5_f64;
    let diameter = 2.0 * radius;
    let p_total = 1.0_f64;

    let msh_path = "assets/disc.msh";
    if !std::path::Path::new(msh_path).exists() {
        eprintln!("Missing {msh_path}. Generate with:");
        eprintln!("  gmsh -2 -format msh4 -o assets/disc.msh assets/disc.geo");
        std::process::exit(1);
    }
    let mesh = parse_msh(msh_path)?;
    let material = PlaneStress { e: 1.0e4, nu: 0.3 };

    let top = closest_node(&mesh, [0.0, radius]);
    let bottom = closest_node(&mesh, [0.0, -radius]);
    let center = closest_node(&mesh, [0.0, 0.0]);

    let t0 = Instant::now();
    let (mut k, mut rhs) = assemble(&mesh, &material, &|_, _| [0.0, 0.0]);

    clamp_node(&mut k, &mut rhs, bottom);
    pin_dof(&mut k, &mut rhs, 2 * top, 0.0); // u_x of top = 0
    apply_point_load(&mut rhs, top, 0.0, -p_total);
    let t_asm = t0.elapsed();

    let t1 = Instant::now();
    let u: Vec<f64> = solve_cg(&k, &rhs, 1e-12, 200_000)?
        .iter()
        .copied()
        .collect();
    let t_solve = t1.elapsed();

    // Stress at the (closest-to-origin) center node: average over the
    // triangles touching it for a cleaner estimate.
    let (sxx_c, syy_c, sxy_c) = stress_at_node(&mesh, &u, &material, center);
    let sxx_analytic = 2.0 * p_total / (PI * diameter);
    let syy_analytic = -6.0 * p_total / (PI * diameter);

    println!("disc Brazilian test: R={radius}, P={p_total}");
    println!(
        "  #nodes={}, #tri={}, top={top}, bottom={bottom}, center={center}",
        mesh.num_nodes(),
        mesh.num_triangles()
    );
    println!("  timings: assemble+BC = {t_asm:?}, CG = {t_solve:?}");
    println!("  σ at center (FEM vs analytic):");
    println!(
        "    σ_xx = {sxx_c:+.4}  (analytic {sxx_analytic:+.4}, ratio {:.3})",
        sxx_c / sxx_analytic
    );
    println!(
        "    σ_yy = {syy_c:+.4}  (analytic {syy_analytic:+.4}, ratio {:.3})",
        syy_c / syy_analytic
    );
    println!("    σ_xy = {sxy_c:+.4}  (analytic 0)");

    // Per-triangle σ_xx for visualization. Clamp the color range so the
    // singularities at the load points don't wash out the interesting band.
    let s_xx_field: Vec<f64> = (0..mesh.num_triangles())
        .map(|t| {
            let eps = strain_at_triangle(&mesh, &u, t);
            stress_at_triangle(&material, eps)[0]
        })
        .collect();
    let limit = 1.5 * sxx_analytic;
    let s_xx_clamped: Vec<f64> = s_xx_field.iter().map(|&v| v.clamp(-limit, limit)).collect();

    let disp = displacement_to_pairs(&u);
    let max_disp: f64 = disp
        .iter()
        .map(|d| (d[0] * d[0] + d[1] * d[1]).sqrt())
        .fold(0.0_f64, f64::max);
    let viz_scale = if max_disp > 0.0 {
        0.05 * radius / max_disp
    } else {
        1.0
    };

    let viewer = FemViewer::new(mesh)
        .with_displacement(disp, viz_scale)
        .with_triangle_field(s_xx_clamped);

    let app = ViewerApp {
        label: format!(
            "Brazilian disc compression\nR={radius}, P={p_total}, E={:.0e}, ν={}\nfield: σ_xx (red = tensile)",
            material.e, material.nu
        ),
        viewer,
    };

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([900.0, 900.0])
            .with_title("fem_rs disc"),
        ..Default::default()
    };

    eframe::run_native("disc", options, Box::new(|_cc| Ok(Box::new(app))))
        .map_err(|e| anyhow::anyhow!("eframe error: {e}"))?;

    Ok(())
}

/// Closest node to `target`, restricted to nodes that are part of at least
/// one triangle. (gmsh sometimes emits construction points that aren't
/// triangulated — closest_node would pick those and break BCs/stress.)
fn closest_node(mesh: &Mesh2D, target: [f64; 2]) -> usize {
    let mut in_triangle = vec![false; mesh.num_nodes()];
    for tri in &mesh.triangles {
        for &n in tri {
            in_triangle[n] = true;
        }
    }
    (0..mesh.num_nodes())
        .filter(|&i| in_triangle[i])
        .min_by(|&a, &b| {
            let da = (mesh.nodes[a][0] - target[0]).powi(2)
                + (mesh.nodes[a][1] - target[1]).powi(2);
            let db = (mesh.nodes[b][0] - target[0]).powi(2)
                + (mesh.nodes[b][1] - target[1]).powi(2);
            da.partial_cmp(&db).unwrap()
        })
        .unwrap()
}

/// Area-weighted average of stress over all triangles incident to a node.
fn stress_at_node(
    mesh: &Mesh2D,
    u: &[f64],
    mat: &PlaneStress,
    node: usize,
) -> (f64, f64, f64) {
    let mut sum_a = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    let mut sxy = 0.0;
    for t in 0..mesh.num_triangles() {
        let tri = mesh.triangles[t];
        if tri.contains(&node) {
            let a = mesh.signed_area(t).abs();
            let eps = strain_at_triangle(mesh, u, t);
            let s = stress_at_triangle(mat, eps);
            sxx += a * s[0];
            syy += a * s[1];
            sxy += a * s[2];
            sum_a += a;
        }
    }
    (sxx / sum_a, syy / sum_a, sxy / sum_a)
}
