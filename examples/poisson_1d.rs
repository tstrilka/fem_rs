use anyhow::Result;
use fem_rs::{assembly::assemble, bc::apply_dirichlet, mesh::Mesh1D, solve::solve_dense};
use plotters::prelude::*;
use std::f64::consts::PI;

/// Problem:
///   -u''(x) = pi^2 * sin(pi x)   on (0, 1)
///   u(0) = u(1) = 0
/// Exact solution: u(x) = sin(pi x). Use this to sanity-check convergence.
fn main() -> Result<()> {
    let num_elements = 32;
    let mesh = Mesh1D::uniform(0.0, 1.0, num_elements);

    let f = |x: f64| PI * PI * (PI * x).sin();
    let (mut k, mut rhs) = assemble(&mesh, &f);

    apply_dirichlet(&mut k, &mut rhs, 0, 0.0);
    apply_dirichlet(&mut k, &mut rhs, mesh.num_nodes() - 1, 0.0);

    let u = solve_dense(&k, &rhs)?;

    let max_err = mesh
        .nodes
        .iter()
        .zip(u.iter())
        .map(|(&x, &u_h)| (u_h - (PI * x).sin()).abs())
        .fold(0.0_f64, f64::max);
    println!("N = {num_elements}, max nodal error = {max_err:.3e}");

    plot(&mesh.nodes, u.as_slice(), "out.svg")?;
    println!("Wrote out.svg");

    Ok(())
}

fn plot(xs: &[f64], u_h: &[f64], path: &str) -> Result<()> {
    let root = SVGBackend::new(path, (800, 500)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .caption("1D Poisson: FEM vs exact", ("sans-serif", 24))
        .build_cartesian_2d(0f64..1f64, -0.1f64..1.2f64)?;
    chart.configure_mesh().draw()?;

    let n = 200;
    let exact = (0..=n).map(|i| {
        let x = i as f64 / n as f64;
        (x, (PI * x).sin())
    });
    chart
        .draw_series(LineSeries::new(exact, BLACK.stroke_width(1)))?
        .label("exact")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], BLACK));

    let approx: Vec<(f64, f64)> = xs.iter().zip(u_h.iter()).map(|(&x, &u)| (x, u)).collect();
    chart
        .draw_series(LineSeries::new(approx.iter().cloned(), RED.stroke_width(2)))?
        .label("FEM")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], RED));
    chart.draw_series(approx.iter().map(|&p| Circle::new(p, 3, RED.filled())))?;

    chart.configure_series_labels().border_style(BLACK).draw()?;
    root.present()?;
    Ok(())
}
