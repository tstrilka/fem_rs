use nalgebra::{Matrix3, Vector3};

/// P1 (linear) triangular element with 3 nodes.
///
/// Shape functions are barycentric coordinates λ_0, λ_1, λ_2.
/// On a triangle with vertices p_0, p_1, p_2:
///     λ_i(p_j) = δ_{ij}
///     λ_0 + λ_1 + λ_2 = 1 everywhere
///
/// Because λ_i is affine in (x, y), its gradient is constant on the triangle.
///
/// See NOTES.md §8–§10 for the derivation.

/// Returns ([∇λ_0, ∇λ_1, ∇λ_2], signed_area).
///
/// Closed-form gradients (rotated edge normals divided by 2A):
///     ∇λ_0 = (1/(2A)) · (y_1 - y_2,  x_2 - x_1)
///     ∇λ_1 = (1/(2A)) · (y_2 - y_0,  x_0 - x_2)
///     ∇λ_2 = (1/(2A)) · (y_0 - y_1,  x_1 - x_0)
pub fn shape_gradients(verts: [[f64; 2]; 3]) -> ([[f64; 2]; 3], f64) {
    let [p0, p1, p2] = verts;
    let area = 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
    let two_a = 2.0 * area;
    let grads = [
        [(p1[1] - p2[1]) / two_a, (p2[0] - p1[0]) / two_a],
        [(p2[1] - p0[1]) / two_a, (p0[0] - p2[0]) / two_a],
        [(p0[1] - p1[1]) / two_a, (p1[0] - p0[0]) / two_a],
    ];
    (grads, area)
}

/// Element stiffness matrix for `-Δu = f` on a P1 triangle.
///
/// TODO:
///   K_e[i][j] = ∫_T (∇λ_i · ∇λ_j) dA
///   Gradients are constant on the triangle, so the integral is just
///   the dot product times the area:
///       K_e[i][j] = area * (∇λ_i · ∇λ_j)
///   Use `shape_gradients(verts)` to get the gradients and area.
pub fn element_stiffness(verts: [[f64; 2]; 3]) -> Matrix3<f64> {
    let (grads, area) = shape_gradients(verts);
    let mut k_e = Matrix3::<f64>::zeros();
    for i in 0..3 {
        for j in 0..3 {
            let dot = grads[i][0] * grads[j][0] + grads[i][1] * grads[j][1];
            k_e[(i, j)] = area * dot;
        }
    }
    k_e
}

/// Element load vector ∫_T λ_i · f(x, y) dA on a P1 triangle.
///
/// TODO:
///   Use 1-point Gauss quadrature at the centroid (exact for linear f,
///   sufficient for P1 demo problems).
///   At the centroid, λ_0 = λ_1 = λ_2 = 1/3, so:
///       F_e[i] ≈ (area / 3) · f(centroid)
///   Centroid = (p_0 + p_1 + p_2) / 3.
///   Later: bump to 3-point Gauss for higher accuracy.
pub fn element_load<F: Fn(f64, f64) -> f64>(verts: [[f64; 2]; 3], f: &F) -> Vector3<f64> {
    let [p0, p1, p2] = verts;
    let (_, area) = shape_gradients(verts);
    let cx = (p0[0] + p1[0] + p2[0]) / 3.0;
    let cy = (p0[1] + p1[1] + p2[1]) / 3.0;
    let c = (area / 3.0) * f(cx, cy);
    Vector3::new(c, c, c)
}
