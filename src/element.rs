use nalgebra::{Matrix2, Vector2};

/// Linear (P1) line element with two nodes.
///
/// Shape functions on reference element xi in [0, 1]:
///     N_0(xi) = 1 - xi
///     N_1(xi) = xi
///
/// You have everything you need below to fill these in by hand.
/// Reference: Larson & Bengzon, Chapter 2.

/// Element stiffness matrix for `-u''(x) = f(x)` on an element of length `h`.
///
/// TODO:
///   Derive K_e = ∫ B^T B dx where B = dN/dx.
///   For a P1 line element this is closed-form — no quadrature needed.
///   Expected result: K_e = (1/h) * [[ 1, -1],
///                                   [-1,  1]]
pub fn element_stiffness(h: f64) -> Matrix2<f64> {
    let c = 1.0 / h;
    Matrix2::new(c, -c, -c, c)
}

/// Element load vector ∫ N_i * f(x) dx over the element.
///
/// TODO:
///   Use midpoint rule for f (1-point Gauss on [0, 1] in reference coords).
///   For midpoint rule: F_e ≈ (h/2) * f(x_mid) * [1, 1].
///   Later: swap to 2-point Gauss for P2 elements / higher-order f.
pub fn element_load<F: Fn(f64) -> f64>(x_left: f64, x_right: f64, f: &F) -> Vector2<f64> {
    let h = x_right - x_left;
    let x_mid = 0.5 * (x_left + x_right);
    let c = 0.5 * h * f(x_mid);
    Vector2::new(c, c)
}
