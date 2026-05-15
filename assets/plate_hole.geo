// Rectangular plate with a circular hole, loaded by uniform tension.
//
// Generate with:
//   gmsh -2 -format msh4 -o assets/plate_hole.msh assets/plate_hole.geo

W = 4.0;     // plate length (x direction, loading axis)
H = 1.0;     // plate height (y direction)
r = 0.15;    // hole radius

h_far  = 0.08;   // mesh size away from hole
h_near = 0.012;  // mesh size at the hole (resolves stress concentration)

// Outer corners
Point(1) = {-W/2, -H/2, 0, h_far};
Point(2) = { W/2, -H/2, 0, h_far};
Point(3) = { W/2,  H/2, 0, h_far};
Point(4) = {-W/2,  H/2, 0, h_far};

// Hole (construction point + 4 arc endpoints)
Point(5) = { 0,   0, 0, h_near};
Point(6) = { r,   0, 0, h_near};
Point(7) = { 0,   r, 0, h_near};
Point(8) = {-r,   0, 0, h_near};
Point(9) = { 0,  -r, 0, h_near};

Line(1) = {1, 2};  // bottom
Line(2) = {2, 3};  // right (tension applied here)
Line(3) = {3, 4};  // top
Line(4) = {4, 1};  // left  (u_x = 0)

Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Curve Loop(1) = {1, 2, 3, 4};     // outer
Curve Loop(2) = {5, 6, 7, 8};     // hole
Plane Surface(1) = {1, 2};
