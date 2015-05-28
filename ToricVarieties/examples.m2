--input: Nothing
--output: Created examples
needsPackage("Polyhedra");
needsPackage("NormalToricVarieties");

P2 = projectiveSpace(2);
A1 = affineSpace(1);
X = P2**A1;
N = matrix{{1,0},{0,1},{0,0}};
