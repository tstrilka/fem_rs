// Disc for Brazilian (diametral compression) test.
//
// Generate with:
//   gmsh -2 -format msh4 -o assets/disc.msh assets/disc.geo

h = 0.025;
R = 0.5;

Point(1) = {0, 0, 0, h};      // center (not used in arcs)
Point(2) = {R, 0, 0, h};      // right (3 o'clock)
Point(3) = {0, R, 0, h};      // top (12 o'clock)
Point(4) = {-R, 0, 0, h};     // left (9 o'clock)
Point(5) = {0, -R, 0, h};     // bottom (6 o'clock)

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
