use anyhow::{Context, Result};
use fem_rs::{
    assembly_2d::assemble, bc::apply_dirichlet_penalty, gmsh::parse_msh, solve::solve_cg,
    viz_2d::render_field,
};
use std::path::Path;
use std::time::Instant;

/// Solve `-Δu = 1` with `u = 0` on the boundary, on an arbitrary mesh
/// loaded from a gmsh `.msh` v4 ASCII file.
///
/// Usage:
///   cargo run --release --example poisson_msh -- path/to/mesh.msh
///
/// Generate an L-shape mesh first:
///   gmsh -2 -format msh4 -o assets/lshape.msh assets/lshape.geo
fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let msh_path = args
        .get(1)
        .map(String::as_str)
        .unwrap_or("assets/lshape.msh");

    if !Path::new(msh_path).exists() {
        eprintln!("File not found: {msh_path}");
        eprintln!();
        eprintln!("Install gmsh, then generate the demo mesh:");
        eprintln!("  sudo apt install gmsh   # Ubuntu/Debian");
        eprintln!("  gmsh -2 -format msh4 -o assets/lshape.msh assets/lshape.geo");
        eprintln!();
        eprintln!("Then re-run:");
        eprintln!("  cargo run --release --example poisson_msh");
        std::process::exit(1);
    }

    let t0 = Instant::now();
    let mesh = parse_msh(msh_path).with_context(|| format!("parsing {msh_path}"))?;
    let t_parse = t0.elapsed();
    println!(
        "loaded {msh_path}: {} nodes, {} triangles, {} boundary nodes  ({:?})",
        mesh.num_nodes(),
        mesh.num_triangles(),
        mesh.boundary_nodes.len(),
        t_parse,
    );

    let f = |_x: f64, _y: f64| 1.0;

    let t1 = Instant::now();
    let (mut k, mut rhs) = assemble(&mesh, &f);
    let t_asm = t1.elapsed();

    for &node in &mesh.boundary_nodes {
        apply_dirichlet_penalty(&mut k, &mut rhs, node, 0.0);
    }

    let t2 = Instant::now();
    let u = solve_cg(&k, &rhs, 1e-10, 100_000)?;
    let t_solve = t2.elapsed();

    let umin = u.iter().cloned().fold(f64::INFINITY, f64::min);
    let umax = u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!(
        "u range: [{umin:.4e}, {umax:.4e}]   assemble: {t_asm:?}, CG: {t_solve:?}"
    );

    let stem = Path::new(msh_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("out");
    let out = format!("{stem}.svg");
    render_field(&mesh, u.as_slice(), &out, &format!("Poisson on {stem}"))?;
    println!("Wrote {out}");
    Ok(())
}
