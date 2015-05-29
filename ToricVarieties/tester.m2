needsPackage("ToricMaps");
needsPackage("Polyhedra");
needsPackage("NormalToricVarieties");
P2 = projectiveSpace(2);
N = normalToricVariety({{1,0},{1,1},{0,1},{-1,-1}},{{0,1},{1,2},{2,3},{3,0}});
f = toricMap(P2,N,matrix({{1,0},{0,1}}))
D = P2_0 + 2*P2_1 - P2_2;
p = pullback(f,D)
print p;
