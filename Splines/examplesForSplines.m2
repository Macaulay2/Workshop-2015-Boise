restart
--Put your own Splines address here before running.
installPackage("Splines",FileName=>"/Users/whieldon/Macaulay2/Workshop-2015-Boise/Splines/Splines.m2")


V = {{0,0},{1,0},{1,1},{-1,1},{-2,-1},{0,-1}}
F = {{0,2,1},{0,2,3},{0,3,4},{0,4,5},{0,1,5}}
E = {{0,1},{0,2},{0,3},{0,4},{0,5}}

splineMatrix(V,F,E,1,CoefficientRing=>ZZ)
splineMatrix(V,F,E,1,Homogenize=>false)


V={{-5,0},{-3,0},{-1,-4},{-1,4},{-1,-2},{-1,2},{0,-1},{0,1},{1,-2},{1,2},{1,-4},{1,4},{3,0},{5,0}}
F={{0, 1, 4, 2}, {0, 1, 5, 3}, {8, 10, 13, 12}, {9, 11, 13, 12}, {1, 4, 6, 7, 5}, {2, 4, 6, 8, 10}, {3, 5, 7, 9, 11}, {6, 7, 9, 12, 8}}
E={{0, 1}, {0, 2}, {0, 3}, {1, 4}, {1, 5}, {2, 4}, {2, 10}, {3, 5}, {3, 11}, {4, 6}, {5, 7}, {6, 7}, {6, 8}, {7, 9}, {8, 10}, {8, 12}, {9, 11}, {9, 12}, {10, 13}, {11, 13}, {12, 13}}

splineMatrix(V,F,E,1,CoefficientRing=>QQ)
splineMatrix(V,F,E,1,Homogenize=>false)



d = # (first V)
facetEdgeH = apply(#E, e-> positions(F, f-> all(E_e,v-> member(v,f))))
indx = positions(facetEdgeH, i-> #i === 2)
E = E_indx;
facetEdgeH = facetEdgeH_indx
BM = matrix apply(facetEdgeH, i-> apply(#F, j-> if (j === first i) then 1 else if (j===last i) then -1 else 0))
S = QQ[t_1..t_d]
varlist = vars S
varCol = transpose varlist
M = (transpose(matrix(S,V)))
mM = numrows M
minorList = apply(E, e-> gens gb minors(mM,M_e|varCol))
if any(minorList, I-> ideal I === sub(ideal 1,S)) then (
    error "Some vertices on entered face are not in codimension 1 face."
