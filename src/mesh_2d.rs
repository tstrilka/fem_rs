/// 2D triangular mesh.
///
/// `nodes[i]` is the (x, y) coordinate of node i.
/// `triangles[t]` is a triple of node indices (CCW orientation expected).
/// `boundary_nodes` is the list of node indices on the domain boundary —
/// where Dirichlet BCs will be applied.
#[derive(Debug, Clone)]
pub struct Mesh2D {
    pub nodes: Vec<[f64; 2]>,
    pub triangles: Vec<[usize; 3]>,
    pub boundary_nodes: Vec<usize>,
}

impl Mesh2D {
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn num_triangles(&self) -> usize {
        self.triangles.len()
    }

    /// Vertices of triangle `t` in world coordinates, in the order stored.
    pub fn triangle_coords(&self, t: usize) -> [[f64; 2]; 3] {
        let [a, b, c] = self.triangles[t];
        [self.nodes[a], self.nodes[b], self.nodes[c]]
    }

    /// Signed area of triangle `t` (positive if vertices are CCW).
    pub fn signed_area(&self, t: usize) -> f64 {
        let [p0, p1, p2] = self.triangle_coords(t);
        0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
    }

    /// All edges that appear in exactly one triangle — the boundary edges.
    /// Each returned tuple `(a, b)` is normalised with `a < b`.
    pub fn boundary_edges(&self) -> Vec<(usize, usize)> {
        use std::collections::HashMap;
        let mut count: HashMap<(usize, usize), u32> = HashMap::new();
        for tri in &self.triangles {
            for k in 0..3 {
                let a = tri[k];
                let b = tri[(k + 1) % 3];
                let key = if a < b { (a, b) } else { (b, a) };
                *count.entry(key).or_insert(0) += 1;
            }
        }
        count
            .into_iter()
            .filter(|&(_, c)| c == 1)
            .map(|(k, _)| k)
            .collect()
    }
}

/// Structured triangular mesh on [0, 1] × [0, 1].
/// Each grid cell is split into 2 triangles (lower-left + upper-right).
/// `nx`, `ny` = number of cells along each axis (so #nodes = (nx+1)*(ny+1)).
pub fn unit_square(nx: usize, ny: usize) -> Mesh2D {
    assert!(nx >= 1 && ny >= 1);

    let idx = |i: usize, j: usize| j * (nx + 1) + i;

    let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
    for j in 0..=ny {
        for i in 0..=nx {
            nodes.push([i as f64 / nx as f64, j as f64 / ny as f64]);
        }
    }

    let mut triangles = Vec::with_capacity(2 * nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let n00 = idx(i, j);
            let n10 = idx(i + 1, j);
            let n11 = idx(i + 1, j + 1);
            let n01 = idx(i, j + 1);
            // CCW orientation
            triangles.push([n00, n10, n11]);
            triangles.push([n00, n11, n01]);
        }
    }

    let mut boundary_nodes = Vec::new();
    for (k, &[x, y]) in nodes.iter().enumerate() {
        if x == 0.0 || x == 1.0 || y == 0.0 || y == 1.0 {
            boundary_nodes.push(k);
        }
    }

    Mesh2D {
        nodes,
        triangles,
        boundary_nodes,
    }
}
