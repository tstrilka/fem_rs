/// 1D mesh: uniformly spaced nodes on [x_min, x_max].
/// In Phase 2 this becomes 2D — keep the names generic-ish.
#[derive(Debug, Clone)]
pub struct Mesh1D {
    pub nodes: Vec<f64>,
}

impl Mesh1D {
    pub fn uniform(x_min: f64, x_max: f64, num_elements: usize) -> Self {
        assert!(num_elements >= 1);
        assert!(x_max > x_min);
        let h = (x_max - x_min) / num_elements as f64;
        let nodes = (0..=num_elements).map(|i| x_min + i as f64 * h).collect();
        Self { nodes }
    }

    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub fn num_elements(&self) -> usize {
        self.nodes.len() - 1
    }

    /// Returns (node_i, node_j) — the two endpoints of element `e`.
    pub fn element_nodes(&self, e: usize) -> (usize, usize) {
        (e, e + 1)
    }

    pub fn element_length(&self, e: usize) -> f64 {
        self.nodes[e + 1] - self.nodes[e]
    }
}
