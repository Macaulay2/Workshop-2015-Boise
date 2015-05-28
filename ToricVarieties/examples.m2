--input: Nothing
--output: Created examples
needsPackcage("Polyhedra");

P2 = projectiveSpace(2);
X = normalToricVariety({{1, 0, 0},{0, 1, 0}, {-1, -1, 0}, {0, 0, 1}}, {{0,1,3},{0,2,3},{1,2,3}});
N = matrix{{1,0},{0,1},{0,0}};
