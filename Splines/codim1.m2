-- Input: facets of a pure simplicial complex (as lists of vertices)
-- Output: the codimension-1 intersections
getCodim1Intersections = Facets ->(
    G = {};
    for i from 0 to #Facets-2 do(
    	f = Facets_i;
    	for j from 0 to #f-1 do(
            g = drop(f,{j,j});
            if not instance(position(drop(Facets,i+1),B ->
		    isSubset(g,B)),Nothing) then G = append(G,g);
    	    ) -- end for
	); -- end for
    G
)
