use fem_rs::gmsh::parse_msh_str;

/// Minimal gmsh v4.1 ASCII mesh: unit square with 4 nodes and 2 triangles.
const UNIT_SQUARE_MSH: &str = "\
$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 1 0
1 0 0 0 1 1 0 0
$EndEntities
$Nodes
1 4 1 4
2 1 0 4
1
2
3
4
0 0 0
1 0 0
1 1 0
0 1 0
$EndNodes
$Elements
1 2 1 2
2 1 2 2
1 1 2 3
2 1 3 4
$EndElements
";

#[test]
fn parse_unit_square() {
    let mesh = parse_msh_str(UNIT_SQUARE_MSH).expect("parse");

    assert_eq!(mesh.num_nodes(), 4, "expected 4 nodes");
    assert_eq!(mesh.num_triangles(), 2, "expected 2 triangles");

    assert_eq!(mesh.nodes[0], [0.0, 0.0]);
    assert_eq!(mesh.nodes[1], [1.0, 0.0]);
    assert_eq!(mesh.nodes[2], [1.0, 1.0]);
    assert_eq!(mesh.nodes[3], [0.0, 1.0]);

    // All 4 corners are on the boundary.
    assert_eq!(mesh.boundary_nodes, vec![0, 1, 2, 3]);

    // Triangles should have positive (CCW) area.
    for t in 0..mesh.num_triangles() {
        let area = mesh.signed_area(t);
        assert!(area > 0.0, "triangle {t} has non-positive area {area}");
    }
}

/// Sanity check: solve Poisson on the parsed unit square, compare to the
/// answer from the structured `unit_square` mesh generator at the same
/// resolution. The two should agree exactly (same mesh topology).
#[test]
fn parsed_mesh_solves_poisson() {
    use fem_rs::{
        assembly_2d::assemble, bc::apply_dirichlet_penalty, mesh_2d::unit_square, solve::solve_cg,
    };
    use std::f64::consts::PI;

    let f = |x: f64, y: f64| 2.0 * PI * PI * (PI * x).sin() * (PI * y).sin();

    let parsed = parse_msh_str(UNIT_SQUARE_MSH).unwrap();
    let generated = unit_square(1, 1);

    assert_eq!(parsed.num_nodes(), generated.num_nodes());
    assert_eq!(parsed.num_triangles(), generated.num_triangles());

    for mesh in [&parsed, &generated] {
        let (mut k, mut rhs) = assemble(mesh, &f);
        for &n in &mesh.boundary_nodes {
            apply_dirichlet_penalty(&mut k, &mut rhs, n, 0.0);
        }
        let u = solve_cg(&k, &rhs, 1e-10, 100).unwrap();
        // Every node is on the boundary in a 1×1 mesh — solution should be ≈ 0.
        let max = u.iter().cloned().fold(0.0_f64, |a, b| a.max(b.abs()));
        assert!(max < 1e-20, "expected zero solution, got {max:e}");
    }
}
