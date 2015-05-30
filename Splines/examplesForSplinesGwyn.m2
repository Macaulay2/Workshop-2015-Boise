uninstallPackage("Splines")
restart
installPackage("Splines",FileName=>"/Users/whieldon/Macaulay2/Workshop-2015-Boise/Splines/Splines.m2")
viewHelp Splines

V = {{0,0},{1,0},{1,1},{-1,1},{-2,-1},{0,-1}};
F = {{0,2,1},{0,2,3},{0,3,4},{0,4,5},{0,1,5}};
E = {{0,1},{0,2},{0,3},{0,4},{0,5}};

S = splineSet(V,F,E,0)
peek S
gens S.SplineModule
