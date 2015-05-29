V = {{0,0},{1,-1},{1,1},{2,0}};
F = {{0,1,2},{1,2,3}};
adjFaces = {{0,1}};
R = QQ[x,y];
forms = {x-1};
M = image splineModule(V,F,1)
N = image splineModule(adjFaces,forms,1,InputType=>"ByLinearForms")
hilbertSeries M
hilbertSeries N

E = {{0,1}};
S = QQ[x,y,w];
ideals = {ideal (x-w)^2}
hilbertSeries generalizedSplines(E,ideals,S)
