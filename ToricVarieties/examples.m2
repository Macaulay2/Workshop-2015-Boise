--input: Nothing
--output: Created examples
needsPackage("Polyhedra");
needsPackage("NormalToricVarieties");

P2 = projectiveSpace(2);
A1 = affineSpace(1);
X = P2**A1;
N = matrix{{1,0},{0,1},{0,0}};

H4 = hirzebruchSurface(4);
P1 = projectiveSpace(1);
M = matrix{{0,1}};

H4;
P2;
L = matrix{{1,0},{0,1}};

