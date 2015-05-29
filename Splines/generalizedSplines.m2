generalizedSplines = method()
--assume vertices are 0,...,n-1
generalizedSplines(List,List,Ring) := Module => (edges,ideals,S) ->(
    edges = apply(edges,sort);
    vertices := unique flatten edges;
    n := #vertices;
    T := directSum(apply(ideals,I->coker gens I))
    M := matrix apply(edges,
	e->apply(n,
	    v->if(v===first e) then 1
	    else if(v===last e) then -1
	    else 0))
   ker(map(T,S^n,sub(M,S)))
);
