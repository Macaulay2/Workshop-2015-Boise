
getSize = method();
getSize(List) := ZZ => vectors ->(
    if all(vectors, v-> #v == #(vectors_0)) then #vectors_0 else null
)

V={{-1, -1, -1}, {3, -1, -1}, {-1, 3, -1}, {-1, -1, 3}, {8, 8, 8}, {-24, 8, 8}, {8, -24, 8}, {8, 8, -24}};
F ={{0, 1, 2, 3}, {0, 5, 6, 7}, {1, 4, 6, 7}, {2, 4, 5, 7}, {3, 4, 5, 6}, {1, 2, 3, 4}, {0, 2, 3, 5}, {0, 1, 3, 6}, {0, 1, 2, 7}, {0, 1, 6, 7}, {0, 2, 5, 7}, {0, 3, 5, 6}, {1, 2, 4, 7}, {1, 3, 4, 6}, {2, 3, 4, 5}};
E=unique flatten apply(F,s->subsets(s,3));


getCodim1IntersectionsPolytope(List,List) := List => (V,F) ->(
    n := #F;
    d := getSize(V);
    --For each non-final facet, construct all codimension 1 subsets.
    codim1faces := apply(n-1, i -> subsets(F_i,d-1));
    --Check if a codimension 1 subset is contained in another facet,
    --store it as a codim 1 intersection.
    sort flatten apply(#codim1faces, i -> 
	select(codim1faces_i, 
	    s -> any(F_{i+1..n-1}, 
		f-> all(s, v-> member(v,f)))))
)

getCodimIFaces = method()
getCodimIFaces(List,ZZ) := List => (F,i) -> (
    d := getSize(F);
    unique flatten apply(F, f-> subsets(f,d-i))
    )

