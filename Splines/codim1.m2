INTgetCodim1Intersections = method();
-- Input: facets of a pure simplicial complex (as lists of vertices)
-- Output: the codimension-1 intersections
INTgetCodim1Intersections(List) := List => facets ->(
    G = {};
    for i from 0 to #facets-2 do(
    	f = facets_i;
    	for j from 0 to #f-1 do(
            g = drop(f,{j,j});
            if not instance(position(drop(facets,i+1),B ->
		    isSubset(g,B)),Nothing) then G = append(G,g);
    	    ) -- end for
	); -- end for
    G
)

INTgetSize = method();
INTgetSize(List) := ZZ => vectors ->(
    n := #(vectors_0);
    if instance(position(vectors,v->#v != n),Nothing) then return n
    else(
	return null
    )
)

INTisSimplicial = method();
-- Assumes that the inputted complex is pure
INTisSimplicial(List,List) := Boolean => (vertices, facets) ->(
    n := getSize(vertices);
    f := getSize(facets);
    if not instance(n, Nothing) and not instance(f,Nothing) and n + 1 == f then true
    else(
	if instance(n, Nothing) then print "Vertices have inconsistent dimension."
	else false
    )
)


