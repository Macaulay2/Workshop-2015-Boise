------ Examples for SplineMatrix
------ SmallExample:Two Dimensional Star of Vertex--
--Input Coordinate List
V={{-1,-1},{0,1},{1,-1},{-2,-2},{0,2},{2,-2}};
--Input List of Facets
F={{0,1,2},{0,1,3,4},{1,2,4,5},{0,2,3,5}};
--Input List of Edges
E={{0,1},{0,2},{1,2},{0,3},{1,4},{2,5},{3,4},{4,5},{3,5}};
--Compute the Matrix whose kernel is the algebra of splines with smoothness 1
--By inputting list of vertices, facets, and edges
splineMatrix(V,F,E,1)
--Can just take vertex and facet data--
splineMatrix(V,F,1)
--Alternatively, user can input list of adjacent facets and equation of linear form between them
--List of adjacent facets
adjFacets = {{0,1},{1,2},{2,3},{3,4},{0,4}};
--List of forms between adjacent facets
R = QQ[x,y];
forms = {y,x-y,x+y,x-2*y,x};
splineMatrix(adjFacets,forms,1)
--Note the user needs to define an ambient ring for this method
--We get generators for the module of splines by calling
splineModule(V,F,1)

--Let's do the example from the slides
--Schlegel Diagram of Cube--
V={{-1,-1},{-1,1},{1,1},{1,-1},{-2,-2},{-2,2},{2,2},{2,-2}};
F={{0,1,2,3},{0,1,4,5},{1,2,5,6},{2,3,6,7},{0,3,4,7}};
E={{0,1},{1,2},{2,3},{0,3},{0,4},{1,5},{2,6},{3,7}};
M=splineModule(V,F,E,0);
--We output the dimensions of the graded pieces of M in a nice way
splineDimTable(0,8,M)
--and we compare these dimensions to the values of the hilbert polynomial
hilbertCTable(0,8,M)
--the values start agreeing when we pass the posulation number of M--
posNum M

--Generalized Splines--
ideals = apply(apply(forms,f -> f^3),ideal);
N=generalizedSplines(adjFacets,ideals)
