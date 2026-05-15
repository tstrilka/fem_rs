// L-shape domain for Poisson demo.
//
// Generate the mesh with:
//   gmsh -2 -format msh4 -o assets/lshape.msh assets/lshape.geo
//
// Then run:
//   cargo run --release --example poisson_msh -- assets/lshape.msh
//
// `h` controls target element size. Smaller h → finer mesh.

h = 0.05;

Point(1) = {0, 0, 0, h};
Point(2) = {2, 0, 0, h};
Point(3) = {2, 1, 0, h};
Point(4) = {1, 1, 0, h};
Point(5) = {1, 2, 0, h};
Point(6) = {0, 2, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};
