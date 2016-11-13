restart
installPackage("AlgebraicSplines",FileName=>"/Users/whieldon/Workshop-2015-Boise/Splines/AlgebraicSplines.m2")

V = {{0,0},{1,0},{1,1},{-1,1},{-2,-1},{0,-1}};
F = {{0,2,1},{0,2,3},{0,3,4},{0,4,5},{0,1,5}};
E = {{0,1},{0,2},{0,3},{0,4},{0,5}};
S = splines(V,F,E,1) -- splines in R^2 with smoothness 1

splines(V,F,E,0,BaseRing=>QQ[x,y,z])
splineModule(V,F,E,0,BaseRing=>QQ[x,y,z])
splineMatrix(V,F,E,0,BaseRing=>QQ[x,y,z])

splineComplex(V,F,1,BaseRing=>QQ[x,y,z])
splineComplex(V,F,1)



R = QQ[x,y]

V = {{-1,0},{1,0},{2,1},{-2,1},{-2,-1},{2,-1}}
F = {{0,1,2,3},{1,2,5},{0,3,4},{0,1,4,5}}
E = {{0,1},{1,2},{2,3},{0,3},{3,4},{0,4},{4,5},{1,5},{2,5}}
