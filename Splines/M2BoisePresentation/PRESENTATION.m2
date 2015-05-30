------ Examples for SplineMatrix
------ SmallExample:Star of Vertex 2D--
R = QQ[x,y];
V={{-1,-1},{0,1},{1,-1},{-2,-2},{0,2},{2,-2}};
F={{0,1,2},{0,1,3,4},{1,2,4,5},{0,2,3,5}};
E={{0,1},{0,2},{1,2},{0,3},{1,4},{2,5},{3,4},{4,5},{3,5}};
L={V,F,E}
splineMatrix(V,F,E,1)
M=splineModule(V,F,E,2)
splineDimTable(0,8,M)
posNum M

adjFacets = {{0,1},{1,2},{2,3},{3,4},{0,4}};
forms = {y,x-y,x+y,x-2*y,x};
ideals = apply(apply(forms,f -> f^3),ideal);
N=generalizedSplines(adjFacets,ideals)
