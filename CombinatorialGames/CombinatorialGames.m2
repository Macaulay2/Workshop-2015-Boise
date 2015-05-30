-- Copyright 2015: Svenja Huntemann, Gwyn Whieldon
-- You may redistribute this file under the terms of the GNU General Public
-- License as published by the Free Software Foundation, either version 2
-- of the License, or any later version.

------------------------------------------
------------------------------------------
-- Header
------------------------------------------
------------------------------------------

if version#"VERSION" <= "1.4" then (
    needsPackage "SimplicialComplexes"
    )

newPackage select((
    "CombinatorialGames",
        Version => "0.0.1", 
        Date => "29. May 2015",
        Authors => {
            {Name => "Svenja Huntemann", Email => "svenja.huntemann@dal.ca", HomePage => "http://mathstat.dal.ca/~svenjah/"},
            {Name => "Gwyn Whieldon", Email => "whieldon@hood.edu", HomePage => "http://cs.hood.edu/~whieldon"}
	},
        Headline => "Package for computing combinatorial game representations of bipartitioned simplicial complexes.",
        Configuration => {},
        DebuggingMode => true,
        if version#"VERSION" > "1.4" then PackageExports => {
	    "SimplicialComplexes"
	    }
        ), x -> x =!= null)

if version#"VERSION" <= "1.4" then (
    needsPackage "SimplicialComplexes"
    )

export {
    "gameRepresentation",
    "facetIdeal",
    "facetComplex",
    "coverDual"
    }

------------------------------------------
------------------------------------------
-- Methods
------------------------------------------
------------------------------------------

------------------------------------------
--Internal method to concatenate strings
--from left/right move sets.
------------------------------------------
cgtJoin = method()
------------------------------------------
cgtJoin(List,HashTable) := String => (L,dummyH) -> (
    if #L === 0 then (
	"") else (
    	fold(concatenate,apply(#L, i-> if i=!=#L-1 then concatenate(dummyH#(L_i),",") else dummyH#(L_i)))
    	)
    )
------------------------------------------


gameRepresentation = method()
gameRepresentation(SimplicialComplex,List,List) := String => (Delta,L,R) -> (
    S := ring Delta;
    L = apply(L,v-> sub(v,S));
    R = apply(R,v-> sub(v,S));
    V := S_*;
    if ((#L+#R) === #V) and (all(join(L,R),v->member(v,V))) then (
	S = QQ[L,R,Degrees=>join(apply(L,i->{1,0}),apply(R,i->{0,1}))];
	L = (S_*)_(toList(0..#L-1));
	R = (S_*)_(toList(#L..#V-1));
	Delta = sub(Delta,S);
	d := dim(Delta);
	Fvec := apply(toList(0..d+1), i->flatten entries faces(d-i,Delta));
	H := hashTable apply(Fvec_0, f-> f=>"{|}");
	E := hashTable {{}=>""};
	dummyH := merge(H,E,join);
	tempStringL := {};
	tempStringR := {};
	Hnew :={};
	for F in drop(Fvec,1) do (
	    tempStringL = apply(apply(F, m-> select(keys H, k-> (degree(k)-degree(m) == {1,0}) and (k % m == 0))), i-> if i=={} then {} else i);
	    tempStringR = apply(apply(F, m-> select(keys H, k-> (degree(k)-degree(m) == {0,1}) and (k % m == 0))), i-> if i=={} then {} else i);
	    Hnew = hashTable apply(#F, i-> F_i => concatenate("{",cgtJoin(tempStringL_i,dummyH),"|",cgtJoin(tempStringR_i,dummyH),"}"));
    	    H = merge(H,Hnew,join);
	    --print H;
    	    dummyH = merge(dummyH,Hnew,join);
	    );
	--print H;
	H#(sub(1,ring Delta))
	)
    else "Variables of Delta not bi-partitioned."
    )

facetIdeal = method()
facetIdeal(SimplicialComplex) := Ideal => Delta -> (
    ideal flatten entries facets(Delta)    
    )

facetComplex = method()
facetComplex(MonomialIdeal) := SimplicialComplex => facetIdeal ->(
    simplicialComplex facetIdeal_*	
    )

coverDual = method()
coverDual(SimplicialComplex) := SimplicialComplex => Delta ->(
    I:=facetIdeal(Delta);
    L:=primaryDecomposition I;
    K:=for i from 0 to #L-1 list product((L_i)_*);
    simplicialComplex K
    )
------------------------------------------
------------------------------------------
-- Documentation
------------------------------------------
------------------------------------------

beginDocumentation()

-- Front Page
doc ///
    Key
        CombinatorialGames
    Headline
        a package for outputting combinatorial game suite formatted games
    Description
        Text
            @SUBSECTION "Definitions"@
	    Let $\Delta$ be a simplicial complex with vertices labeled by
	    a partition of variables $L = \{x_1,x_2,...,x_n\}$ and
	    $R = \{y_1,y_2,...,y_m\}$ where $L\cup R = V$.
	Text
	
	Text
	    Put description of how this constructor works in here.
        Text
            @SUBSECTION "Other acknowledgements"@
            --
            This package started at Macaulay2 Workshop in Boise, supported by
	    NSF grant and organized by Zach Teitler, Hirotachi Abo and 
	    Frank Moore.
///

------------------------------------------
-- Data type & constructor
------------------------------------------

-- gameRepresentation
doc ///
    Key
        gameRepresentation
	(gameRepresentation,SimplicialComplex,List,List)
    Headline
        compute the game representation to export to cgsuite
    Usage
    	G = gameRepresentation(Delta,L,R)
    Inputs
    	Delta:SimplicialComplex
	    simplicial complex in bipartitioned variables
	L:List
	    list of left-labeled vertices in Delta
	R:List
	    list of right-labeled vertices in Delta
    Outputs
    	G:String
	  string to export to cgsuite
    Description
        Text
	    This takes a simplicial complex and a partition of its
	    vertices into two subsets and outputs the game representation
	    of Delta in the format understood by open-source software 
	    @HREF("http://cgsuite.sourceforge.net/", "cgsuite")@.
	Example
	    S = QQ[x_1,x_2,x_3,y_1];-- create ring for simplicial complex Delta
            Delta = simplicialComplex(apply({{x_1,x_2,x_3},{x_2,x_3,y_1}},product));  -- create Delta, here input by list of facet vertices
            L = {x_1,x_2,x_3} -- list of left player (L) vertex labels
	    R = {y_1} -- list of right player (R) vertex labels
    	    gameRepresentation(Delta,L,R)
        Text
            Here calling the variables $x_1,...,x_n$ and $y_1,...,y_m$ is
	    for convenience, and is not necessary in inputting game.
	Example
	    S = QQ[a,b,c,d,e,f]
	    F = {{a,b,c},{a,c,d},{a,d,f},{c,d,e},{b,c,e},{d,e,f},{b,e,f}}
	    Delta = simplicialComplex(apply(F,product))
	    L = {a,b,d,f}
	    R = {c,e}
	    gameRepresentation(Delta,L,R)

///


end

