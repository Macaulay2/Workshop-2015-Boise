
getSize = method();
getSize(List) := ZZ => vectors ->(
    if all(vectors, v-> #v == #(vectors_0)) then #vectors_0 else null
)

V={{-1, -1, -1}, {3, -1, -1}, {-1, 3, -1}, {-1, -1, 3}, {8, 8, 8}, {-24, 8, 8}, {8, -24, 8}, {8, 8, -24}};
F ={{0, 1, 2, 3}, {0, 5, 6, 7}, {1, 4, 6, 7}, {2, 4, 5, 7}, {3, 4, 5, 6}, {1, 2, 3, 4}, {0, 2, 3, 5}, {0, 1, 3, 6}, {0, 1, 2, 7}, {0, 1, 6, 7}, {0, 2, 5, 7}, {0, 3, 5, 6}, {1, 2, 4, 7}, {1, 3, 4, 6}, {2, 3, 4, 5}};
E=unique flatten apply(F,s->subsets(s,3));

restart
installPackage("Splines")

V={{0,0,-1},{-1,-1,-2},{1,-1,-2},{1,1,-2},{-1,1,-2},{0,0,1},{-1,-1,2},{1,-1,2},{1,1,2},{-1,1,2}};
F={{0,1,2,3,4},{0,1,4,5,6,9},{0,1,2,5,6,7},{0,2,3,5,7,8},{0,3,4,5,8,9},{5,6,7,8,9}};
E={{0,1,2},{0,2,3},{0,3,4},{0,1,4},{0,1,5,6},{0,2,5,7},{0,3,5,8},{0,4,5,9},{5,6,7},{5,7,8},{5,6,9},{5,8,9}};

sort interiorFaces(F,E)
time sort getCodimIFacesPolytope(F,1)
time sort interiorFaces(F,getCodimIFacesSimplicial(F,1))

(sort interiorFaces(F,E)) == (sort getCodimIFacesPolytope(F,1))


V={{-1, -1, -1}, {3, -1, -1}, {-1, 3, -1}, {-1, -1, 3}, {8, 8, 8}, {-24, 8, 8}, {8, -24, 8}, {8, 8, -24}};
F ={{0, 1, 2, 3}, {0, 5, 6, 7}, {1, 4, 6, 7}, {2, 4, 5, 7}, {3, 4, 5, 6}, {1, 2, 3, 4}, {0, 2, 3, 5}, {0, 1, 3, 6}, {0, 1, 2, 7}, {0, 1, 6, 7}, {0, 2, 5, 7}, {0, 3, 5, 6}, {1, 2, 4, 7}, {1, 3, 4, 6}, {2, 3, 4, 5}};
E=unique flatten apply(F,s->subsets(s,3));


getCodim1FacesPolytope = method()
getCodim1FacesPolytope(List) := List => F ->(
    --This function ASSUMES that the polytopal complex 
    --considered is hereditary.
    n := #F;
    --For each pair of facets, take their intersection:
    intersectFacets := unique flatten apply(#F-1, i-> apply(toList(i+1..#F-1), j-> sort select(F_i, v-> member(v,F_j))));
    --Remove any non-maximal faces in this intersections:
    select(intersectFacets, f -> (
    	(number(intersectFacets, g-> all(f, j-> member(j,g)))) === 1
    ))
)

getCodimIFacesPolytope = method()
getCodimIFacesPolytope(List,ZZ) := List => (F,d) ->(
    Fcodim := F;
    apply(d, i-> Fcodim = getCodim1FacesPolytope(Fcodim));
    Fcodim
    )

getCodimIFacesSimplicial = method()
getCodimIFacesSimplicial(List,ZZ) := List => (F,i) -> (
    d := getSize(F);
    unique flatten apply(F, f-> subsets(f,d-i))
    )

