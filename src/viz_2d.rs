use crate::mesh_2d::Mesh2D;
use anyhow::Result;
use plotters::prelude::*;

/// Render a scalar field `u` (one value per node) on a triangular mesh as an
/// SVG. Flat shading: each triangle gets the color of its centroid value.
pub fn render_field(mesh: &Mesh2D, u: &[f64], path: &str, title: &str) -> Result<()> {
    let (xmin, xmax, ymin, ymax) = mesh.nodes.iter().fold(
        (f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY),
        |(xa, xb, ya, yb), &[x, y]| (xa.min(x), xb.max(x), ya.min(y), yb.max(y)),
    );
    let umin = u.iter().cloned().fold(f64::INFINITY, f64::min);
    let umax = u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let root = SVGBackend::new(path, (700, 700)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .caption(title, ("sans-serif", 20))
        .build_cartesian_2d(xmin..xmax, ymin..ymax)?;
    chart.configure_mesh().draw()?;

    for t in 0..mesh.num_triangles() {
        let [a, b, c] = mesh.triangles[t];
        let avg = (u[a] + u[b] + u[c]) / 3.0;
        let color = blue_to_red(normalize(avg, umin, umax));
        let pts = [mesh.nodes[a], mesh.nodes[b], mesh.nodes[c]];
        chart.draw_series(std::iter::once(Polygon::new(
            pts.iter().map(|&[x, y]| (x, y)).collect::<Vec<_>>(),
            color.filled(),
        )))?;
        chart.draw_series(std::iter::once(PathElement::new(
            [pts[0], pts[1], pts[2], pts[0]]
                .iter()
                .map(|&[x, y]| (x, y))
                .collect::<Vec<_>>(),
            BLACK.mix(0.2),
        )))?;
    }

    root.present()?;
    Ok(())
}

fn normalize(v: f64, lo: f64, hi: f64) -> f64 {
    if hi - lo < 1e-15 {
        0.5
    } else {
        ((v - lo) / (hi - lo)).clamp(0.0, 1.0)
    }
}

/// Simple linear blue → white → red colormap.
fn blue_to_red(t: f64) -> RGBColor {
    if t < 0.5 {
        let s = t * 2.0;
        RGBColor((255.0 * s) as u8, (255.0 * s) as u8, 255)
    } else {
        let s = (t - 0.5) * 2.0;
        RGBColor(255, (255.0 * (1.0 - s)) as u8, (255.0 * (1.0 - s)) as u8)
    }
}
