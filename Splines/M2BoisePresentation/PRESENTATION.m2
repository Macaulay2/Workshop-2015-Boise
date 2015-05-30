------ Examples for cube
R = QQ[x,y];
V={{-1,-1},{-1,1},{1,1},{1,-1},{-2,-2},{-2,2},{2,2},{2,-2}};
F={{0,1,2,3},{0,1,4,5},{1,2,5,6},{2,3,6,7},{0,3,4,7}};
E={{0,1},{1,2},{2,3},{0,3},{0,4},{1,5},{2,6},{3,7}};

splineMatrix(V,F,E,2)
M=splineModule(V,F,E,2)
splineDimTable(0,8,M)
posNum M

adjFacets = {{0,1},{1,2},{2,3},{3,4},{0,4}};
forms = {y,x-y,x+y,x-2*y,x};
ideals = apply(apply(forms,f -> f^3),ideal);
N=generalizedSplines(adjFacets,ideals)
