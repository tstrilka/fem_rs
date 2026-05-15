use anyhow::Result;
use fem_rs::{
    assembly_2d::assemble,
    bc::apply_dirichlet_penalty,
    gmsh::parse_msh,
    mesh_2d::{Mesh2D, unit_square},
    solve::solve_cg,
    viewer::{FemViewer, ViewerApp},
};
use std::f64::consts::PI;
use std::path::Path;

/// Interactive FEM viewer.
///
/// Usage:
///   cargo run --release --bin viewer                       # built-in unit square demo
///   cargo run --release --bin viewer assets/lshape.msh     # load a .msh file
fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let (label, mesh, u) = match args.get(1) {
        Some(path) => {
            let mesh = parse_msh(path)?;
            let f = |_x: f64, _y: f64| 1.0;
            let u = solve_poisson(&mesh, &f)?;
            let name = Path::new(path)
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("mesh")
                .to_string();
            (format!("loaded: {name}\nf = 1"), mesh, u)
        }
        None => {
            let mesh = unit_square(32, 32);
            let f = |x: f64, y: f64| 2.0 * PI * PI * (PI * x).sin() * (PI * y).sin();
            let u = solve_poisson(&mesh, &f)?;
            (
                "unit_square(32×32)\nf = 2π² sin(πx) sin(πy)".to_string(),
                mesh,
                u,
            )
        }
    };

    println!(
        "solved: {} nodes, {} triangles, u range = [{:.3e}, {:.3e}]",
        mesh.num_nodes(),
        mesh.num_triangles(),
        u.iter().copied().fold(f64::INFINITY, f64::min),
        u.iter().copied().fold(f64::NEG_INFINITY, f64::max),
    );

    let app = ViewerApp {
        label,
        viewer: FemViewer::new(mesh).with_node_scalar(u),
    };

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1000.0, 750.0])
            .with_title("fem_rs viewer"),
        ..Default::default()
    };

    eframe::run_native(
        "fem_rs viewer",
        options,
        Box::new(|_cc| Ok(Box::new(app))),
    )
    .map_err(|e| anyhow::anyhow!("eframe error: {e}"))?;

    Ok(())
}

fn solve_poisson<F: Fn(f64, f64) -> f64>(mesh: &Mesh2D, f: &F) -> Result<Vec<f64>> {
    let (mut k, mut rhs) = assemble(mesh, f);
    for &node in &mesh.boundary_nodes {
        apply_dirichlet_penalty(&mut k, &mut rhs, node, 0.0);
    }
    let u = solve_cg(&k, &rhs, 1e-10, 100_000)?;
    Ok(u.iter().copied().collect())
}
