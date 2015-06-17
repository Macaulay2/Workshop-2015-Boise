-- Copyright 2014-2015: Mike Dipasquale
-- You may redistribute this file under the terms of the GNU General Public
-- License as published by the Free Software Foundation, either version 2
-- of the License, or any later version.

------------------------------------------
------------------------------------------
-- Header
------------------------------------------
------------------------------------------

if version#"VERSION" <= "1.4" then (
    needsPackage "Graphs",
    needsPackage "Polyhedra"
    )

newPackage select((
    "AlgebraicSplines",
        Version => "0.1.0", 
        Date => "27. May 2015",
        Authors => {
            {Name => "Mike DiPasquale", Email => "midipasq@gmail.com", HomePage => "http://illinois.edu/~dipasqu1"},
            {Name => "Gwyn Whieldon", Email => "whieldon@hood.edu", HomePage => "http://cs.hood.edu/~whieldon"},
	    {Name => "Eliana Duarte", Email => "emduart2@illinois.edu", HomePage => "http://illinois.edu/~emduart2"},
	    {Name => "Daniel Irving Bernstein", Email=> "dibernst@ncsu.edu", HomePage =>"http://www4.ncsu.edu/~dibernst"}
        },
        Headline => "Package for computing topological boundary maps and piecewise continuous splines on polyhedral complexes.",
        Configuration => {},
        DebuggingMode => true,
        if version#"VERSION" > "1.4" then PackageExports => {
	    "Polyhedra",
	    "Graphs"
	    }
        ), x -> x =!= null)

if version#"VERSION" <= "1.4" then (
    needsPackage "Polyhedra",
    needsPackage "Graphs"
    )

export {
   "Splines",
   "VertexCoordinates",
   "Regions",
   "SplineModule",
   "splines",
   "Spline",
   "spline",
   "isTPure",
   "getDim",
   "formsList",
   "splineMatrix",
   "splineModule",
   "InputType",
   "RingType",
   "ByFacets",
   "ByLinearForms",
   "isHereditary",
   "CheckHereditary",
   "Homogenize",
   "VariableName",
   "interiorFaces",
   "splineDimTable",
   "posNum",
   "hilbertCTable",
   "hilbertPolyEval",
   "generalizedSplines",
   "cellularComplex",
   "idealsComplex",
   "splineComplex",
   --get rid of these once testing is done
   "getCodim1Intersections",
   "getCodimDIntersections",
   "getCodimDFacesSimplicial",
   "issimplicial",
   "simpBoundary",
   "boundaryComplex",
   "subsetL",
   "codim1Cont",
   "orient",
   "polyBoundaryPair",
   "polyBoundary"
    }

------------------------------------------
------------------------------------------
-- Data Types and Constructors
------------------------------------------
------------------------------------------

--Create an object that gives ALL splines
--on a given subdivision.
Splines = new Type of HashTable
splines = method(Options => {
	symbol InputType => "ByFacets", 
	symbol CheckHereditary => false, 
	symbol Homogenize => true, 
	symbol VariableName => getSymbol "t",
	symbol CoefficientRing => QQ})

splines(List,List,List,ZZ) := Matrix => opts -> (V,F,E,r) -> (
    	AD := splineMatrix(V,F,E,r,opts);
	K := ker AD;
	b := #F;
    	new Splines from {
	    symbol cache => new CacheTable from {"name" => "Unnamed Spline"},
	    symbol VertexCoordinates => V,
	    symbol Regions => F,
	    symbol SplineModule => image submatrix(gens K, toList(0..b-1),)
	}
)


net Splines := S -> S.SplineModule

Spline = new Type of HashTable
spline = method()

spline(Splines,List) := (S,L) -> (
    M := S.SplineModule;
    )
   


------------------------------------------
------------------------------------------
-- Methods
------------------------------------------
------------------------------------------


------------------------------------------
subsetL=method()
------------------------------------------
--Containment function for lists--

subsetL(List,List):=List=>(L1,L2)->(
    all(L1,f->member(f,L2))
    )

-----------------------------------------
isTPure=method()
-----------------------------------------
--Inputs:
--V = list of vertices
--F = list of facets
-----------------------------------------
--Outputs:
--True, if every facet has same dimension
--AS THE AMBIENT SPACE
--False, if some facets have dimension different
--from the ambient space
-----------------------------------------

isTPure(List,List):= Boolean =>(V,F)->(
    d := #(first V);
    V = apply(V,v->prepend(1,v));
    if all(F,f->((rank matrix(V_f))==d+1)) then (true) else (false)
    )

-----------------------------------------
getDim=method()
-----------------------------------------
--Input: V= list of vertices
-----------------------------------------
--Output: Dimension of affine span of V
-----------------------------------------
getDim(List):=ZZ=> V ->(
    V=apply(V,v->prepend(1,v));
    (rank matrix V)-1
    )


-----------------------------------------
-----------------------------------------
interiorFaces = method()
-----------------------------------------
-----------------------------------------
--Inputs: 
-----------------------------------------
--F = list of facets
--E = list of codimension 1 faces
-- (possibly including non-interior)
-----------------------------------------
-----------------------------------------
--Outputs:
-----------------------------------------
--E' = list of interior codimension 1 faces
-----------------------------------------
interiorFaces(List,List) := List => (F,E) -> (
    --Compute which facets are adjacent to each edge:
    facetEdgeH := apply(#E, e-> positions(F, f-> all(E_e,v-> member(v,f))));
    --Compute indices of interior edges, and replace edge list and 
    --facet adjacencies to only include these interior edges:
    indx := positions(facetEdgeH, i-> #i === 2);
    E_indx
    )

------------------------------------------
------------------------------------------
getCodim1Intersections = method(Options=>{
	symbol InputType => "Polyhedral"
	}
    )
------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--F = list of facets of a polytopal 
--or a simplicial complex
------------------------------------------
--Outputs:
-----------------------------------------
--E = largest (under containment) intersections
--of facets.  If input is hereditary, this is
--the list of interior codim 1 intersections
-----------------------------------------

getCodim1Intersections(List) := List => opts -> F ->(
     n := #F;
    if opts.InputType==="Polyhedral" then(
    	--For each pair of facets, take their intersection:
     	intersectFacets := unique flatten apply(#F-1, i-> 
    	    apply(toList(i+1..#F-1), 
	    	j-> sort select(F_i,
		    v-> member(v,F_j))));
	--Remove any non-maximal faces in this intersections:
    	codim1int:=select(intersectFacets, f -> (
    		(number(intersectFacets, g-> all(f, j-> member(j,g)))) === 1
    		))
    	) else if opts.InputType==="Simplicial" then(
    	d := #(F_0);
    	--For each non-final facet, construct all codimension 1 subsets.
    	codim1faces := apply(n-1, i -> subsets(F_i,d-1));
    	--Check if a codimension 1 subset is contained in another facet,
    	--store it as a codim 1 intersection.
    	codim1int=sort flatten apply(#codim1faces, i -> 
	    select(codim1faces_i, 
	    	s -> any(F_{i+1..n-1}, 
		    f-> subsetL(s,f))));
	);
    codim1int
)

------------------------------------------
------------------------------------------
getCodimDIntersections = method(Options=>{
	symbol InputType => "Polyhedral"
	}
    )
------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--F = list of faces of a polytope
--d = desired codimesion
------------------------------------------
--Outputs:
-----------------------------------------
--E = list of (interior) codim d faces
-----------------------------------------
getCodimDIntersections(List,ZZ) := List => opts->(F,d) ->(
    if opts.InputType === "Polyhedral" then(
    	Fcodim := F;
    	--Iteratively compute intersections up to codim d --
    	apply(d, i-> Fcodim = getCodim1Intersections(Fcodim))
    ) else if opts.InputType === "Simplicial" then(
	--Get all faces of codimension d--
	Fcodim = getCodimDFacesSimplicial(F,d);
	--Get boundary faces of codimension d--
	boundaryF := boundaryComplex(F);
	boundaryCodim := getCodimDFacesSimplicial(boundaryF,d);
	--Select nonBoundary faces of codimension d--
	Fcodim = select(Fcodim, f-> not member(f, boundaryCodim)
	    )
    	);
    Fcodim
    )

------------------------------------------
------------------------------------------
getCodimDFacesSimplicial = method()
------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--F = list of facets of a simplicial complex
--d = desired codimension
------------------------------------------
--Outputs:
-----------------------------------------
--E = list of (all) codim d faces
-----------------------------------------
getCodimDFacesSimplicial(List,ZZ) := List => (F,D) -> (
    d := #first(F);
    unique flatten apply(F, f-> subsets(f,d-D))
    )

------------------------------------------
verifyComplex= method()
------------------------------------------
-- This method should check everything.
------------------------------------------
--Inputs:
--V = list of vertices
--F = list of facets
------------------------------------------
--Outputs:
--True or false
--Also should pinpoint what the problem is.
--This won't detect everything... the user
--can't screw up too much...
------------------------------------------

verifyComplex(List,List):=Boolean=>(V,F)->(
    --Check that (V,F) is pure of correct dim
    --Check that intersections of facets have correct dimension
    --Limited Check that facets don't overlap (check along codim 1 intersections)
    --Check hereditary
    )

------------------------------------------
------------------------------------------
isHereditary= method()
------------------------------------------
------------------------------------------
-- This method checks if the PURE polyhedral
-- complex with vertices V and facets F
-- is hereditary.  Dimension of polyhedral
-- complex must be same as ambient dimension.
-- Needs correction for 3 dim and higher!
------------------------------------------
--Inputs: 
------------------------------------------
--F = ordered lists of facets
------------------------------------------
--Outputs:
--Boolean, if complex is hereditary
------------------------------------------

isHereditary(List,List) := Boolean => (V,F) -> (
    d := #(first V);
    -- Checks that all maximal intersections of facets have dimension d-1
    E := getCodim1Intersections(F);
    if any(E,e->(getDim(V_e)!= d-1)) then (
	bool:=false
	) else (
        -- Checks non-branching condition: all codim 1 faces are in two facets
	dualE := apply(E,e->select(#F,f->subsetL(e,F_f)));
	if any(dualE,f->(#f>2)) then (
	    bool = false
	    )else (
	    dualG := graph(dualE,EntryMode=>"edges");
	    linkH := hashTable apply(#V, v-> v=>select(#F, f -> member(v,F_f)));
	    -- Checks if the dual graph of the star of each each vertex is connected. Note: need to extend
	    -- this to link of every face... my bad (Mike D.)
	    bool = all(keys linkH, k-> isConnected inducedSubgraph(dualG,linkH#k))
	    ));
    bool
)

isHereditary(List) := Boolean => F -> (
    V := unique flatten join F;
    E := getCodimDFacesSimplicial(F,1);
    dualV := toList(0..#F-1);
    dualE := apply(#E, e-> positions(F, f-> all(E_e,v-> member(v,f))));
    if not all(dualE,e-> #e <= 2) then (
	false -- Checks pseudo manifold condition
      ) else (
      dualG := graph(dualE,EntryMode=>"edges");
      linkH := hashTable apply(V, v-> v=>select(#F, f -> member(v,F_f)));
      -- Checks if the link of each vertex is connected.
      all(keys linkH, k-> isConnected inducedSubgraph(dualG,linkH#k))
      )
)

-----------------------------------------

formsList=method(Options=>{
	symbol InputType => "ByFacets", 
	symbol CheckHereditary => false, 
	symbol Homogenize => true, 
	symbol VariableName => getSymbol "t",
	symbol CoefficientRing => QQ}
    )
----------------------------------------------------
--This method returns a list of forms corresponding to codimension one faces
--------
--Input:
--V= vertex list
--E= codim 1 face list
--r= smoothness parameter
--------
--Output:
--List of forms defining input list of codimension one faces 
--raised to (r+1) power
------------------------------------------------------------

formsList(List,List,ZZ):=List=>opts->(V,E,r)->(
    --To homogenize, we append a 1 as the final coordinate of each vertex coord in list V.
    --If not homogenizing, still need 1s for computing equations
    d := #(first V);
    t := opts.VariableName;
    V = apply(V, v-> append(v,1));
    if opts.Homogenize then (
	    S := (opts.CoefficientRing)[t_0..t_d];
	    varlist := (vars S)_(append(toList(1..d),0));
	    ) else (
	    S = (opts.CoefficientRing)[t_1..t_d];
	    varlist = (vars S)|(matrix {{sub(1,S)}});
	    );
    varCol := transpose varlist;
    M := (transpose(matrix(S,V)));
    mM := numrows M;
    minorList := apply(E, e-> gens gb minors(mM,matrix(M_e)|varCol));
    if any(minorList, I-> ideal I === ideal 1) then (
    	error "Some vertices on entered face are not in codimension 1 face."
	    );
    flatten apply(minorList, m -> (m_(0,0))^(r+1))
)


-----------------------------------------
-----------------------------------------
splineMatrix = method(Options => {
	symbol InputType => "ByFacets", 
	symbol CheckHereditary => false, 
	symbol Homogenize => true, 
	symbol VariableName => getSymbol "t",
	symbol CoefficientRing => QQ})
------------------------------------------
------------------------------------------

------------------------------------------
------------------------------------------
-- splineMatrix "ByFacets"
------------------------------------------
--Inputs: 
------------------------------------------
--("ByFacets")
--L = {L_0,L_1,L_2} (i.e. {V,F,E})
--r = degree of desired continuity 
--
-- OR!!!!
--
--("ByLinearForms")
--L = {L_0,L_1} (i.e. {B,L})
--r = degree of desired continuity
------------------------------------------
--Outputs:
-- BM = matrix with columns corresponding
-- to facets and linear forms separating facets.
------------------------------------------
splineMatrix(List,ZZ) := Matrix => opts -> (L,r) -> (
    --Use this if your list L = {V,F,E} contains
    --The inputs as a single list L.
    if opts.InputType === "ByFacets" then (
	splineMatrix(L_0,L_1,L_2,r)
	);
    if opts.InputType == "ByLinearForms" then (
	splineMatrix(L_0,L_1,r,InputType=>"ByLinearForms")
	)
    )

------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--("ByFacets")
--V = list of coordinates of vertices
--F = ordered lists of facets
--E = list of edges
--r = degree of desired continuity
------------------------------------------
--Outputs:
-- BM = matrix with columns corresponding
-- to facets and linear forms separating facets.
------------------------------------------
splineMatrix(List,List,List,ZZ) := Matrix => opts -> (V,F,E,r) -> (
    if opts.InputType === "ByFacets" then (
		if opts.CheckHereditary then (
	    	    if not isHereditary(F,E) then (
			error "Not hereditary."
			);
	    	    );
	d := # (first V);
	--Compute which facets are adjacent to each edge:
	facetEdgeH := apply(#E, e-> positions(F, f-> all(E_e,v-> member(v,f))));
	--Compute indices of interior edges, and replace edge list and 
	--facet adjacencies to only include these interior edges:
	indx := positions(facetEdgeH, i-> #i === 2);
	E = E_indx;
	facetEdgeH = facetEdgeH_indx;
	--Compute top boundary map for complex:
	BM := matrix apply(
	    facetEdgeH, i-> apply(
		#F, j-> if (
		    j === first i) then 1 else if (
		    j===last i) then -1 else 0));
	--List of forms definining interior codim one faces (raised to (r+1) power)
	flist := formsList(V,E,r,opts);
	T := diagonalMatrix(flist);
	splineM := BM|T;
	) else if opts.InputType === "ByLinearForms" then (
	 print "Wrong inputs, put in lists of adjacent facets and linear forms and continuity r."
    	 );
    splineM
)

------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
------------------------------------------
-- ("ByFacets")
--V = list of vertex coordinates
--F = list of facets
--r = degree of desired continuity
--
--    OR!!!!
--
--("ByLinearForms")
--B = list of adjacent facets
--L = list of ordered linear forms
--defining codim 1 faces, ordered as in B
--r = degree of desired continuity
------------------------------------------
--Outputs:
-- BM = matrix with columns corresponding
-- to facets and linear forms separating facets.
------------------------------------------
splineMatrix(List,List,ZZ) := Matrix => opts -> (V,F,r) ->(
    --Warn user if they are accidentally using ByFacets method with too few inputs.
    --This code assumes that the polytopal complex is hereditary.
    if opts.InputType === "ByFacets" then (
	if issimplicial(V,F) then(
	    E := getCodim1Intersections(F,InputType=>"Simplicial");
	    SM := splineMatrix(V,F,E,r,opts)  
	    )
	else(
	    E = getCodim1Intersections(F);
	    SM = splineMatrix(V,F,E,r,opts)
	    );
	);
    --If user DOES want to define complex by regions and dual graph.
    if opts.InputType === "ByLinearForms" then (
	B := V;
	L := F;
	m := max flatten B;
	A := matrix apply(B, i-> apply(toList(0..m), j-> 
		if (j=== first i) then 1 
		else if (j===last i) then -1 
		else 0));
	D := matrix apply(#L, i-> apply(#L, j-> if i===j then L_i^(r+1) else 0));
	SM = A|D;
    );
    SM
)


------------------------------------------
------------------------------------------
splineModule = method(Options => {
	symbol InputType => "ByFacets",
	symbol CheckHereditary => false, 
	symbol Homogenize => true, 
	symbol VariableName => getSymbol "t",
	symbol CoefficientRing => QQ}
    )
------------------------------------------
------------------------------------------
-- This method computes the splineModule
-- of a complex Delta, given by either
-- facets, codim 1 faces, and vertex coors,
-- or by pairs of adjacent faces and
-- linear forms.
------------------------------------------
--Inputs: 
------------------------------------------
--V = list of vertices
--F = list of facets
--E = list of edges
--r = desired continuity of splines
------------------------------------------
--Outputs:
--Spline module S^r(Delta)
------------------------------------------
splineModule(List,List,List,ZZ) := Matrix => opts -> (V,F,E,r) -> (
    	AD := splineMatrix(V,F,E,r,opts);
	K := ker AD;
	b := #F;
    	image submatrix(gens K, toList(0..b-1),)
)

------------------------------------------
--Inputs: 
------------------------------------------
--V = list of vertices
--F = list of facets
--r = desired continuity of splines
--
--    OR!!!!
--
--V = list of pairs of adjacent faces
--F = list of linear forms defining codim 1 faces.
--r = desired continuity of splines
------------------------------------------
--Outputs:
--Spline module S^r(Delta)
------------------------------------------
splineModule(List,List,ZZ) := Matrix => opts -> (V,F,r) -> (
    	AD := splineMatrix(V,F,r,opts);
	K := ker AD;
	b := #F;
	if opts.InputType==="ByLinearForms" then (
		b = #(unique flatten V)
		);
    	image submatrix(gens K, toList(0..b-1),)
)
------------------------------------------
-------------------------------------------
-------------------------------------------
splineDimTable=method(Options => {
	symbol InputType => "ByFacets"
	}
    )
-------------------------------------------
-----Inputs:
-------------------------------------------
----- a= lower bound of dim table
----- b= upper bound of dim table
----- M= module
--------------------------------------------
------ Outputs:
--------------------------------------------
------ A net with the degrees between a and b on top row
------ and corresponding dimensions of graded pieces
------ of M in bottom row
-------------------------------------------

splineDimTable(ZZ,ZZ,Module):=Net=>opts->(a,b,M)->(
    r1:=prepend("Degree",toList(a..b));
    r2:=prepend("Dimension",apply(toList(a..b),i->hilbertFunction(i,M)));
    netList {r1,r2}
    )

-------------------------------------------
-----Inputs:
-------------------------------------------
----- a= lower bound of range
----- b= upper bound of range
----- L= list {V,F,E}, where V is a list of vertices, F a list of facets, E a list of codim 1 faces
----- r= degree of desired continuity
-------------------------------------------
-----Outputs:
-------------------------------------------
-------A table with the dimensions of the graded pieces
------ of the spline module in the range (a,b)
-------------------------------------------

splineDimTable(ZZ,ZZ,List,ZZ):= Net=>opts->(a,b,L,r)->(
    M := splineModule(L_0,L_1,L_2,r);
    splineDimTable(a,b,M)
    )

-------------------------------------------
-----Inputs:
-------------------------------------------
----- a= lower bound of range
----- b= upper bound of range
----- L= list {V,F}, where V is list of vertices, F a list of facets
-------------
-----OR!!
-------------------------------------------
----- a= lower bound of range
----- b= upper bound of range
----- L= list {V,F}, where V is a list of adjacent facets, F a list of forms
-----------defining codim 1 faces along which adjacent facets meet
-------------------------------------------
-------Outputs:
-------------------------------------------
-------A table with the dimensions of the graded pieces
------ of the spline module in the range (a,b)
-------------------------------------------

splineDimTable(ZZ,ZZ,List,ZZ):= Net => opts->(a,b,L,r)->(
    M := splineModule(L_0,L_1,r,opts);
    splineDimTable(a,b,M)
    )


-------------------------------------------
-------------------------------------------
posNum=method()
-------------------------------------------
-----Inputs:
-------------------------------------------
----- M, a graded module
--------------------------------------------
------ Outputs:
--------------------------------------------
------ The postulation number (largest integer 
------ for which Hilbert function and polynomial 
------ of M disagree).
--------------------------------------------
posNum(Module):= (N) ->(
    k := regularity N;
    while hilbertFunction(k,N)==hilbertPolyEval(k,N) do	(k=k-1);
    k
    )

------------------------------------------
-----------------------------------------

hilbertCTable=method()
-------------------------------------------
-----Inputs:
-------------------------------------------
----- a= an integer, lower bound
----- b= an integer, upper bound
----- M= graded module over a polynomial ring
--------------------------------------------
------ Outputs:
--------------------------------------------
------ A table whose top two rows are the same as
------ the output of splineDimTable and whose 
------ third row compares the first two to the
------ Hilbert Polynomial
--------------------------------------------

hilbertCTable(ZZ,ZZ,Module):= (a,b,M) ->(
    r1:=prepend("Degree",toList(a..b));
    r2:=prepend("Dimension",apply(toList(a..b),i->hilbertFunction(i,M)));
    r3:=prepend("HilbertPoly",apply(toList(a..b),i->hilbertPolyEval(i,M)));
    netList {r1,r2,r3}
    )
---------------------------------------------

hilbertPolyEval=method()
---------------------------------------------
-------------------------------------------
-----Inputs:
-------------------------------------------
----- i= integer at which you will evaluate the Hilbert polynomial
----- M= module
--------------------------------------------
------ Outputs:
--------------------------------------------
------ An Hilbert polynomial of the module M
------ evaluated at i.
--------------------------------------------

hilbertPolyEval(ZZ,Module):=(i,M)->(
    P:=hilbertPolynomial(M,Projective=>false);
    sub(P,(vars ring P)_(0,0)=>i)
    )

------------------------------------------

generalizedSplines = method(Options=>{
	symbol RingType => "Ambient"
	}
    )
------------------------------------------
------------------------------------------
-- This method computes the generalized spline module
-- associated to a graph whose edges are labeled by ideals.
------------------------------------------
--Inputs: 
------------------------------------------
--E = list of edges. Each edge is a list with two vertices.
----The set of vertices must be the integers 0..n-1.
--ideals = list of ideals that label the edges. Ideals must 
--be entered in same order as corresponding edges in E.
--If RingType is an integer,n, then this can be a list 
--of integers that label the edges.
-------------------------------------------------
--if RingType is an integer, n, then the ring is defined 
--internally as ZZ/n.
--If RingType is "Ambient", an ambient ring must already 
--be defined so that ideals can be entered.
------------------------------------------
--Outputs:
------------------------------------------
--Module of generalized splines on the graph given by the edgelist.
------------------------------------------

generalizedSplines(List,List) := Module => opts -> (E,ideals) ->(
    if opts.RingType === "Ambient" then(
    S := ring first ideals;
    --make sure ideals all lie in same ambient ring--
    ideals = apply(ideals, I->sub(I,S));
    ) else if toString(class(opts.RingType))==="ZZ" then(
       m := opts.RingType;
       if m==0 then(
	   S = ZZ;
	   ideals = apply(ideals,i-> ideal(i_S));
       ) else (
       	   R := ZZ[{}];
	   S = R/ideal(m_R);
	   ideals = apply(ideals,i-> ideal(i_S));
       )
    );
    --assume vertices are 0,...,n-1
    vertices := unique flatten E;
    n := #vertices;
    T := directSum(apply(ideals,I->coker gens I));
--Boundary Map from Edges to vertices (this encodes spline conditions)
    M := matrix apply(E,
	e->apply(n,
	    v->if(v===first e) then 1
	    else if(v===last e) then -1
	    else 0));
    ker(map(T,S^n,sub(M,S)))
)


------------------------------------------
simpBoundary = method()
------------------------------------------
--Input:
--F = list of codim i faces
--E = list of codim i+1 faces
------------------------------------------
--Output:
--B = boundary map matrix between codim i and codim i+1 faces
------------------------------------------
--Example:
--F = {{0,1,2},{0,1,3},{1,3,4},{1,2,4},{2,4,5},{0,2,5},{0,3,5}}
--E = {{0,1},{1,2},{0,2},{3,0},{1,3},{1,4},{2,4},{2,5},{0,5},{3,4},{4,5}}
--V = {{0},{1},{2},{4}}
------------------------------------------
simpBoundary(List,List) := Matrix => (F,E) -> (
    F = apply(F, f-> sort f);
    E = apply(E, e-> sort e);
    tempLF := {};
    rowList := {};
    apply(F, f-> (
	    tempLF = hashTable apply(#f, v-> position(E,e-> e == drop(f,{v,v})) => (-1)^v);
	    rowList = append(rowList,apply(#E, j->if member(j,keys tempLF) then tempLF#j else 0));
	    )
	);
    transpose matrix rowList
    )

------------------------------------------
orient = method()
------------------------------------------
--Input:
--I = the ideal generated by forms vanishing
--on codim one faces containing a face P
------------------------------------------
--Output:
--L = list of length dim P, smallest variables
--(w.r.t. standard lex order)
--which are dependent on the affine span of P
------------------------------------------

orient(Ideal):=List=> I->(
    LT:= leadTerm I;
    reverse sort flatten apply(flatten entries LT,f-> support f)
    )

------------------------------------------
polyBoundaryPair=method()
------------------------------------------
----This method is for use when there are precomputed
----lists of ideals and orientations
------------------------------------------
--Input:
--V = vertices of complex
--L1 = {G,IG,OG}
--L2 = {H,IH,OH}
--G = codim (i+1) face (list of indices of vertices) 
--IG = ideal of G
--OG = orientation of G
--H = codim i face (list of indices of vertices)
--IH = ideal of H
--OH = orientation of H
------------------------------------------
--Output:
--O =  +1,-1, or 0, encoding compatibility
--of orientations between G and H
------------------------------------------

polyBoundaryPair(List,List,List):=ZZ=>(V,L2,L1)->(
    G:=L1_0;
    H:=L2_0;
    --if G is not a face of H then return 0 --
    if not subsetL(G,H) then (
	ort:=0) else(	
	--get ambient ring--	
	S:=ring (L1_1);
	--make sure ideals and variables sit inside same ambient ring--
	IG:=sub(L1_1,S);
	OG:=apply(L1_2,var->sub(var,S));
	IH:=sub(L2_1,S);
	OH:=apply(L2_2,var->sub(var,S));
	--get a vertex of H that is not a vertex of G--
	testVindex :=select(1,H,v->not member(v,G));
	testV := transpose matrix({prepend(1,flatten V_testVindex)});
	--get a generator of IG that is not a generator of IH--
	outVect :=first select(1,flatten entries gens IG,f->(f%IH!=0));
	--get row vector whose entries are coefficients of outVect
	matV :=transpose jacobian ideal outVect;
	--modify outVect so vector corresponding to outVect points 
	--outward from H
	if sub(((matV*testV)_(0,0)),ZZ)<0 then(outVect=-outVect);
	--reduce outVect modulo IH--
	outVect =(outVect%IH);
	--Get independent variables on the affine span of H and G--
	indH :=select(flatten entries vars S,v->not member(v,OH));
	indG :=select(flatten entries vars S,v-> not member(v,OG));
	--There is one independent variable on H that is dependent on G.  
	--Reduce it mod IG-- 
	pos :=position(indH,v->not member(v,indG));
	redHG :=(indH_pos)%IG;
	--Get coefficients of contraction of outVect with indH--
	inCont :=matrix{apply(length indH,i->(
		    (-1)^i*coefficient(indH_i,outVect)
		    ))};
	--Get list of coefficients for pulling back inCont to G--
        conv :=matrix{apply(length indH,i->(
		omit := drop(indH,{i,i});
		if member(indH_pos,omit) then(
		    posN := position(omit,v->(v==(indH_pos)));
		    pmvar := position(indG,v->(v==(indH_i)));
		    missvar :=indG_pmvar;
		    coef :=(-1)^posN*coefficient(missvar,redHG)*(-1)^pmvar
		    ) else (coef=1);
		coef))};
	sgn :=(inCont*(transpose conv))_(0,0);
	if (sgn>0) then(ort=1) else (ort=-1);
	);
    ort
    )

------------------------------------------
polyBoundary=method()
------------------------------------------
----This method is for use when there are precomputed
----lists of ideals and orientations
------------------------------------------
--Input:
--V = vertices of complex
--L1 = {G,IG,OG}
--L2 = {H,IH,OH}
--G = codim (i+1) faces (list of indices of vertices) 
--IG = ideals of G
--OG = orientations of G
--H = codim i faces (list of indices of vertices)
--IH = ideals of H
--OH = orientations of H
------------------------------------------
--Output:
--O =  +1,-1, or 0, encoding compatibility
--of orientations between G and H
------------------------------------------

polyBoundary(List,List,List):=Matrix=>(V,L2,L1)->(
    G:=L1_0;
    H:=L2_0;
    IG:=L1_1;
    OG:=L1_2;
    IH:=L2_1;
    OH:=L2_2;
    matrix table(#G,#H,(i,j)->(
	    L1n :={G_i,IG_i,OG_i};
	    L2n :={H_j,IH_j,OH_j};
	    polyBoundaryPair(V,L2n,L1n)
	    ))
    )

------------------------------------------

boundaryComplex = method()
------------------------------------------
--Input:
--F= list of facets of a simplicial complex
----which is a pseudomanifold (Important!)
------------------------------------------
--Output:
--A list of codim one faces on the boundary
------------------------------------------
boundaryComplex(List) := List => F -> (
    n := #F;
    d := #(F_0);
    codim1faces := unique flatten apply(n,i-> subsets(F_i,d-1));
    select(codim1faces, f-> number(F, g-> all(f, v-> member(v,g))) === 1)
    )

------------------------------------------------

cellularComplex = method(
    	Options =>{
	    symbol InputType => "Polyhedral",
	    symbol Homogenize => true, 
	    symbol VariableName => getSymbol "t",
	    symbol CoefficientRing => QQ
	    }
    )
------------------------------------------------
---This method computes the cellular chain complex of a simplicial or
---complex with coefficients in a polynomial
---ring, modulo the boundary
------------------------------------------------
---Inputs (if simplicial): A list of facets
------------------------------------------------
---Outputs: The cellular chain complex whose homology
--- is the homology of the simplicial complex relative
--- to its boundary.
--------------------------------------------------
cellularComplex(List) := ChainComplex => opts -> (F) -> (
    if opts.InputType === "Polyhedral" then (
	print "Need a List of Vertices";
	chain := null;
	);
    if opts.InputType === "Simplicial" then (
	d := (# first F)-1;
	if opts.Homogenize then (
	    t := opts.VariableName;
	    S := (opts.CoefficientRing)[t_0..t_d];
	    ) else (
	    t = opts.VariableName;
	    S = (opts.CoefficientRing)[t_1..t_d];
	    );
	boundaryF := boundaryComplex(F);
	C := apply(d+1, i-> getCodimDFacesSimplicial(F,i));
	boundaryC := join({{}},apply(d, i-> getCodimDFacesSimplicial(boundaryF,i)));
    	intC := apply(#C, i -> select(C_i, f -> not member(f,boundaryC_i)));
    	chain = chainComplex(reverse apply(#intC-1, c-> simpBoundary(intC_c,intC_(c+1))))**S
	);
    chain
    )


------------------------------------------------
---Inputs: 
-- V= list of vertices
-- F= list of facets
------------------------------------------------
---Outputs: The cellular chain complex whose homology
--- is the homology of the simplicial or polyhedral complex relative
--- to its boundary.
--------------------------------------------------

cellularComplex(List,List) := ChainComplex => opts -> (V,F) -> (
    d := (# first V);
    if opts.Homogenize then (
	t := opts.VariableName;
	S := (opts.CoefficientRing)[t_0..t_d];
	) else (
	t = opts.VariableName;
	S = (opts.CoefficientRing)[t_1..t_d];
	);
    if issimplicial(V,F) then (
	boundaryF := boundaryComplex(F);
	C := apply(d+1, i-> getCodimDFacesSimplicial(F,i));
	boundaryC := join({{}},apply(d, i-> getCodimDFacesSimplicial(boundaryF,i)));
    	intC := apply(#C, i -> select(C_i, f -> not member(f,boundaryC_i)));
    	chain := chainComplex(reverse apply(#intC-1, c-> simpBoundary(intC_c,intC_(c+1))))
	) else (
	--Construct list whose ith element is intersections of codim i--
	current :=F;
	intC ={current};
	scan(d,i->(
		current=getCodim1Intersections(current);
		intC =append(intC,current)
	));
    	--get the forms defining codimension 1 faces--
	fList :=formsList(V,intC_1,0,opts);
    	--create a list whose ith element is ideals of codim i faces--
	idList :={apply(F,f->ideal(0_S)),apply(fList,f->ideal f)};
	scan(d-1,i->(
		idList=append(idList,apply(intC_(i+2),G->(
			    ind:=codim1Cont(intC_1,G);
			    sub(ideal fList_ind,S)
		)))
        ));
        --set orientations of all faces (lexicographically lowest dependent
	--variables on each face)
	orList := apply(idList, L->apply(L,I->orient I));
	--set up the chain complex
	chain = chainComplex(reverse apply(#intC-1, c-> (
		    L1 := {intC_(c+1),idList_(c+1),orList_(c+1)};
		    L2 := {intC_c,idList_c,orList_c};
		    polyBoundary(V,L2,L1)
		    )))
	);
    chain**S
    )



------------------------------------------
idealsComplex=method(Options=>{
	symbol Homogenize => true, 
	symbol VariableName => getSymbol "t",
	symbol CoefficientRing => QQ
    }
    )
------------------------------------------
--This function computes the Schenck-Stillman chain complex
--of ideals for a simplicial or polyhedral complex
------------------------------------------
--Inputs:
--V: vertex list
--F: facet list
--r: desired order of smoothness
------------------------------------------
--Outputs: The Schenck-Stillman complex of ideals
------------------------------------------

idealsComplex(List,List,ZZ):=ChainComplex => opts -> (V,F,r)->(
    d := #(first V);
    if opts.Homogenize then (
	t := opts.VariableName;
	S := (opts.CoefficientRing)[t_0..t_d];
	) else (
	t = opts.VariableName;
	S = (opts.CoefficientRing)[t_1..t_d];
	);
    if issimplicial(V,F) then (
	--list of interior faces in order of increasing codimension--
	boundaryF := boundaryComplex(F);
	C := apply(d+1, i-> getCodimDFacesSimplicial(F,i));
	boundaryC := join({{}},apply(d, i-> getCodimDFacesSimplicial(boundaryF,i)));
	intC := apply(#C, i -> select(C_i, f -> not member(f,boundaryC_i)));
	--list of forms defining codim 1 interior faces
	intformslist := formsList(V,intC_1,r,opts);
	--list of modules which will define chain complex--
	fullmodulelist:= apply(#intC,i->directSum apply(intC_i,e->(
		CE := positions(intC_1,f->subsetL(e,f));
		sub(module ideal (intformslist_CE),S)
		)));
	--defining the chain complex
	CCSS :=chainComplex(reverse apply(#intC-1, c-> (
		    inducedMap(fullmodulelist_(c+1),fullmodulelist_c,(simpBoundary(intC_c,intC_(c+1)))**S)
		    ))
	    )
    	) else (
	current :=F;
	intC ={current};
	scan(d,i->(
		current=getCodim1Intersections(current);
		intC =append(intC,current)
	));
    	--get the forms defining codimension 1 faces--
	fList :=formsList(V,intC_1,0,opts);
    	--create a list whose ith element is ideals of codim i faces--
	idList :={apply(F,f->ideal(0_S)),apply(fList,f->ideal f)};
	scan(d-1,i->(
		idList=append(idList,apply(intC_(i+2),G->(
			    ind:=codim1Cont(intC_1,G);
			    sub(ideal fList_ind,S)
		)))
        ));
        --set orientations of all faces (lexicographically lowest dependent
	--variables on each face)
	orList := apply(idList, L->apply(L,I->orient I));
	--set up list of forms to r+1 power--
	intformslist =apply(fList,f->f^(r+1));
	--list of modules which will define chain complex--
	fullmodulelist = apply(#intC,i->directSum apply(intC_i,e->(
		CE := codim1Cont(intC_1,e);
		sub(module ideal (intformslist_CE),S)
		)));
	--set up the chain complex
	CCSS = chainComplex(reverse apply(#intC-1, c-> (
		    L1 := {intC_(c+1),idList_(c+1),orList_(c+1)};
		    L2 := {intC_c,idList_c,orList_c};
		    M := polyBoundary(V,L2,L1)**S;
		    inducedMap(fullmodulelist_(c+1),fullmodulelist_c,M)
		    )))
	);
    CCSS
    )


------------------------------------------
splineComplex=method(Options=>{
	symbol Homogenize => true, 
	symbol VariableName => getSymbol "t",
	symbol CoefficientRing => QQ
    }
    )

------------------------------------------
--This function computes the Schenck-Stillman chain complex
--of quotient rings for a simplicial or polyhedral complex,
--called the spline complex
------------------------------------------
--Inputs:
--V: vertex list
--F: facet list
--r: desired order of smoothness
------------------------------------------
--Outputs: The Schenck-Stillman spline complex
------------------------------------------

splineComplex(List,List,ZZ):=ChainComplex => opts -> (V,F,r)->(
    d := #(first V);
    if opts.Homogenize then (
	t := opts.VariableName;
	S := (opts.CoefficientRing)[t_0..t_d];
	) else (
	t = opts.VariableName;
	S = (opts.CoefficientRing)[t_1..t_d];
	);
    if issimplicial(V,F) then (
	--list of interior faces in order of increasing codimension--
	boundaryF := boundaryComplex(F);
	C := apply(d+1, i-> getCodimDFacesSimplicial(F,i));
	boundaryC := join({{}},apply(d, i-> getCodimDFacesSimplicial(boundaryF,i)));
	intC := apply(#C, i -> select(C_i, f -> not member(f,boundaryC_i)));
	--list of forms defining codim 1 interior faces
	intformslist := formsList(V,intC_1,r,opts);
	--list of modules which will define chain complex--
	fullmodulelist:= {S^(#F)};
	scan(#intC-1,i->(
		newMod := directSum apply(intC_(i+1),e->(
		CE := codim1Cont(intC_1,e);
		coker sub(gens ideal (intformslist_CE),S)
		));
	    	fullmodulelist=append(fullmodulelist,newMod)
	));
	--defining the chain complex
	CCSS :=chainComplex(reverse apply(#intC-1, c-> (
		    map(fullmodulelist_(c+1),fullmodulelist_c,(simpBoundary(intC_c,intC_(c+1)))**S)
		    ))
	    );
    	) else (
	current :=F;
	intC ={current};
	scan(d,i->(
		current=getCodim1Intersections(current);
		intC =append(intC,current)
	));
    	--get the forms defining codimension 1 faces--
	fList :=formsList(V,intC_1,0,opts);
    	--create a list whose ith element is ideals of codim i faces--
	idList :={apply(F,f->ideal(0_S)),apply(fList,f->ideal f)};
	scan(d-1,i->(
		idList=append(idList,apply(intC_(i+2),G->(
			    ind:=codim1Cont(intC_1,G);
			    sub(ideal fList_ind,S)
		)))
        ));
        --set orientations of all faces (lexicographically lowest dependent
	--variables on each face)
	orList := apply(idList, L->apply(L,I->orient I));
	--set up list of forms to r+1 power--
	intformslist =apply(fList,f->f^(r+1));
	--list of modules which will define chain complex--
	fullmodulelist= {S^(#F)};
	scan(#intC-1,i->(
		newMod := directSum apply(intC_(i+1),e->(
		CE := codim1Cont(intC_1,e);
		coker sub(gens ideal (intformslist_CE),S)
		));
	    	fullmodulelist=append(fullmodulelist,newMod)
	));
	--set up the chain complex
	CCSS = chainComplex(reverse apply(#intC-1, c-> (
		    L1 := {intC_(c+1),idList_(c+1),orList_(c+1)};
		    L2 := {intC_c,idList_c,orList_c};
		    M := polyBoundary(V,L2,L1)**S;
		    map(fullmodulelist_(c+1),fullmodulelist_c,M)
		    )))
	);
    CCSS
    )


------------------------------------------
splineComplexMap=method(Options=>{
	symbol Homogenize => true, 
	symbol VariableName => getSymbol "t",
	symbol CoefficientRing => QQ
    }
    )
------------------------------------------
--Builds the three complexes idealsComplex,
--cellularComplex, splineComplex, along with
--chain maps between these complexes
------------------------------------------
--Inputs:
--V: vertex list
--F: facet list
--r: desired order of smoothness
------------------------------------------
--Outputs: The three chain complexes with
--maps between them
------------------------------------------

splineComplexMap(List,List,ZZ):=List=>opts->(V,F,r)->(
    )
	    

------------------------------------------
codim1Cont=method()
------------------------------------------
----Inputs:
----E = list of codim one intersections
----G = a face of the complex (V,F)
------------------------------------------
----Output: Indices of E corresponding to
----- codim one faces containing G
-----------------------------------------

codim1Cont(List,List):=List=> (E,G)->(
    --positions of codim one faces containing G
    positions(E,e->subsetL(G,e))
    )


------------------------------------------
issimplicial = method()
------------------------------------------
-- Assumes that the inputted complex is pure
------------------------------------------
--Inputs:
-- V = vertex coordinates of Delta
-- F = list of facets of Delta
------------------------------------------
--Outputs:
--Boolean, if Delta is simplicial,
--checking that each facet is a simplex
--of the appropriate dimension.
------------------------------------------
issimplicial(List,List) := Boolean => (V,F) ->(
    n := #first(V);
    all(F,f->#f==(n+1))
)

-----------------------------------------

------------------------------------------
------------------------------------------
-- Documentation
------------------------------------------
------------------------------------------

beginDocumentation()

-- Front Page
doc ///
    Key
        AlgebraicSplines
    Headline
        a package for building splines and computing bases
    Description
        Text
            This package provides methods for computations with piecewise polynomial functions (splines) over
	    polytopal complexes.
	Text
	    @SUBSECTION "Definitions"@
	Text
	    Let $\Delta$ be a partition (simplicial,polytopal,cellular,rectilinear, etc.) of a space $\RR^n$.
	    The spline module $S_d^{r}(\Delta)$ is the module of all functions $f\in C^r(\Delta)$ such that
	    $f$ is a polynomial of degree $d$ when restricted to each face $\sigma\in\Delta$.
	Text
	
	Text
	    This package computes the @TO splineModule@ and @TO splineMatrix@ of $\Delta$, as well
	    as defining new types @TO Splines@ and @TO Spline@ that contain geometric data 
	    for $\Delta$ (if entered) and details on the associated spline module $S_d^r(\Delta)$.
        Text
	    @SUBSECTION "Other acknowledgements"@
            --
            Methods in this package borrows heavily from code written by Hal Schenck
	    and Mike DiPasquale.
///

------------------------------------------
-- Data type & constructor
------------------------------------------

-- Spline Matrix Constructor
doc ///
    Key
        splineMatrix
	(splineMatrix,List,List,ZZ)
	(splineMatrix,List,List,List,ZZ)
	InputType
	CheckHereditary
	ByFacets
	ByLinearForms
    Headline
        compute matrix giving adjacent regions and continuity level
    Usage
    	S = splineMatrix(V,F,E,r)
	S = splineMatrix(B,L,r)
    Inputs
    	V:List
	    list of coordinates of vertices of Delta
	F:List
	    list of facets of Delta
	E:List
	    list of edges of Delta
	r:ZZ
	    degree of desired continuity
	InputType=>String
	    either "ByFacets", or "ByLinearForms"
	CheckHereditary=>Boolean
	    either "true" or "false", depending on if you want
	    to check if Delta is hereditary before attempting 
	    to compute splines.
    Outputs
    	S:Matrix
	  resulting spline module
    Description
        Text
	    This creates the basic spline matrix that has splines as
	    its kernel.
	Example
	    V = {{0,0},{1,0},{1,1},{-1,1},{-2,-1},{0,-1}};-- the coordinates of vertices
            F = {{0,2,1},{0,2,3},{0,3,4},{0,4,5},{0,1,5}};  -- a list of facets (pure complex)
            E = {{0,1},{0,2},{0,3},{0,4},{0,5}};   -- list of edges in graph
    	    splineMatrix(V,F,E,1)
        Text
            Alternately, spline matrices can be created directly from the
	    dual graph (with edges labeled by linear forms).  Note: This way of
	    entering data requires the ambient polynomial ring to be defined.
	Example
	    R = QQ[x,y]
	    B = {{0,1},{1,2},{2,3},{3,4},{4,0}}
	    L = {x-y,y,x,y-2*x,x+y}
	    splineMatrix(B,L,1,InputType=>"ByLinearForms")

///

doc ///
    Key
        isHereditary
	(isHereditary,List,List)
	(isHereditary,List)
    Headline
    	checks if a complex $\Delta$ is hereditary
    Usage
    	B = isHereditary(F,E)
	B = isHereditary(F)
    Inputs
    	F:List
	    list of facets of F
	E:List
	    list of codimension 1 faces of F
    Outputs
    	B:Boolean
	    returns true if F is hereditary
    Description
        Text
	    A complex $\Delta$ is hereditary if it is a pseudomanifold (all 
	    codimention 1 faces are contained in two facets), and the link of 
	    each vertex is connected.
	
	Text
	    The hereditary check can take both facets and codimension 1 faces:
	
	Example
	    F = {{1,2,3},{2,3,4},{3,4,5},{4,5,6}}
	    E = {{2,3},{3,4},{4,5},{5,6}}
	    isHereditary(F,E)
	    
	Example
	    F = {{1,2,3},{2,3,4},{3,4,5},{5,6,7}}
	    E = {{2,3},{3,4},{4,5}}
	    isHereditary(F,E)
	    
	Text
	    Alternately, if the complex is simplicial, codimension 1 faces can
	    be computed automatically.
	    
	Example
	    F = {{1,2,3},{2,3,4},{3,4,5},{4,5,6}}
	    isHereditary(F)
    SeeAlso
        splineMatrix
	
/// 

--<<<<<<< HEAD
--=======
doc ///
    Key
        splineModule
	(splineModule,List,List,List,ZZ)
	(splineModule,List,List,ZZ)
    Headline
        compute the module of all splines on partition of a space
    Usage
        M = splineModule(V,F,E,r)
	M = splineModule(V,F,r)
    Inputs
        V:List
	    V = list of coordinates of vertices
	F:List
	    F = list of facets
	E:List
	    E = list of codimension 1 faces (interior or not)
	r:ZZ
	    r = desired degree of smoothness
    Outputs
        M:Module
	    M = module of splines on $\Delta$
    Description
        Text
	    This is some text.
	Example
	    V = {{0,0},{1,0},{1,1},{0,1}}
	    F = {{0,1,2},{0,2,3}}
	    E = {{0,1},{0,2},{0,3},{1,2},{2,3}}
	    splineModule(V,F,E,1)
    Caveat
        I'm not sure if this is fully documented yet.
	
	
///

doc ///
    Key
        splineDimTable
	(splineDimTable,ZZ,ZZ,Module)
	(splineDimTable,ZZ,ZZ,List,ZZ)
    Headline
        a table with the dimensions of the graded pieces of a graded module
    Usage
        T=splineDimTable(a,b,M)
	T=splineDimTable(a,b,L,r)
    Inputs
        a:ZZ
	    a= lowest degree in the table
	b:ZZ
	    b= largest degree in the table
	N:Module
	    M= graded module
	L:List
	    L= a list {V,F,E} of the vertices, faces and edges of a polyhedral complex
	r:ZZ
	    r= degree of smoothnes 

    Outputs
        T:Table
	    T= table with the dimensions of the graded pieces of M in the range a,b
    Description
        Text
	    The output table gives you the dimensions of the graded pieces
	    of the module M where the degree is between a and b. 
	Example
	    V = {{0,0},{1,0},{1,1},{0,1}}
	    F = {{0,1,2},{0,2,3}}
	    E = {{0,1},{0,2},{0,3},{1,2},{2,3}}
	    M=splineModule(V,F,E,2)
	    splineDimTable(0,8,M)
	Text
	    You may instead input the list L={V,F,E} of the vertices, faces and edges of the spline.
	Example
	    L = {V,F,E};
	    splineDimTable(0,8,L,2)
	
      
///

doc ///
    Key
        hilbertPolyEval
	(hilbertPolyEval,ZZ,Module)
    Headline
        a function to evaluate the hilbertPolynomial of a graded module at an integer
    Usage
        v = hilbertPolyEval(a,M)
    Inputs
        a:ZZ
	    a= integer at which you will evaluate the hilbertPolynomial of the graded module M
	M:Module
	    M= graded module
    Outputs
        v:ZZ
	    v= hilbertPolynomial of the graded module M evaluated at a
    Description
        Text
            For any graded module M and any integer a, you may evaluate the hilberPolynomial of M
	    at a.
	Example
	    V = {{0,0},{1,0},{1,1},{0,1}};
	    F = {{0,1,2},{0,2,3}};
	    E = {{0,1},{0,2},{0,3},{1,2},{2,3}};
	    M = splineModule(V,F,E,2)
	    hilbertPolyEval(2,M)
	    
///

doc ///
    Key
        posNum
	(posNum,Module)
    Headline
        computes the largest degree at which the hilbert function of the graded module M is not equal to the hilbertPolynomial
    Usage
        v = posNum(M)
    Inputs
        M:Module
	    M= graded module
    Outputs
        v:ZZ
	    v= largest degree at which the hilbert function of the graded module M is not equal to the hilbertPolynomial
    Description
        Text
	    This function computes the postulation number of M which is defined as the
	    largest degree at which the hilbert function of the graded module M is not equal to the hilbertPolynomial
	Example
	    V = {{0,0},{1,0},{1,1},{0,1}};
	    F = {{0,1,2},{0,2,3}};
	    E = {{0,1},{0,2},{0,3},{1,2},{2,3}};
	    M = splineModule(V,F,E,2)
	    posNum(M)
	    
///
        

doc ///
    Key
        hilbertCTable
	(hilbertCTable,ZZ,ZZ,Module)
    Headline
        a table to compare the values of the hilbertFunction and hilbertPolynomial of a graded module
    Usage
        T = hilbertCTable(a,b,M)
    Inputs
        a:ZZ
	    a= lowest degree in the  table
	b:ZZ
	    b= largest degree in the table
	M:Module
	    M= graded module
    Outputs        
	T:Table
	    T= table with the degrees and values of the hilbertFunction and hilbertPolynomial
    Description
        Text
	    The first row of the output table contains the degrees, the second row contains the 
	    values of the hilbertFunction, the third row contains the values of the hilbertPolynomial
	Example
	    V = {{0,0},{1,0},{1,1},{0,1}}
	    F = {{0,1,2},{0,2,3}}
	    E = {{0,1},{0,2},{0,3},{1,2},{2,3}}
	    hilbertCTable(0,8,splineModule(V,F,E,1))

///

doc ///
    Key
        Splines
	VertexCoordinates
	Regions
	SplineModule
    Headline
    	a class for splines (piecewise polynomial functions on subdivisions)
    Description
    	Text
	    This class is a type of @TO "HashTable"@ that stores information on
	    a subdivision $\Delta$ of ${\mathbb R}^n$, given by a set of vertex
	    coordinates and a list of facets (and possibly edges), along with a
	    module of all splines on $\Delta$ of continuity $r$.
	Example
	    V = {{0,0},{1,0},{1,1},{-1,1},{-2,-1},{0,-1}};
	    F = {{0,2,1},{0,2,3},{0,3,4},{0,4,5},{0,1,5}};
	    E = {{0,1},{0,2},{0,3},{0,4},{0,5}};
	    S = splines(V,F,E,1) -- splines in R^2 with smoothness 1
    SeeAlso
        splines
	Spline
	spline
///
-->>>>>>> origin/master

TEST ///
V = {{0,0},{1,0},{1,1},{-1,1},{-2,-1},{0,-1}}
F = {{0,2,1},{0,2,3},{0,3,4},{0,4,5},{0,1,5}}
E = {{0,1},{0,2},{0,3},{0,4},{0,5}}
assert(splineMatrix(V,F,E,0) == matrix {{1, 0, 0, 0, -1, t_2, 0, 0, 0, 0}, {1, -1, 0, 0, 0, 0, t_1-t_2,
      0, 0, 0}, {0, 1, -1, 0, 0, 0, 0, t_1+t_2, 0, 0}, {0, 0, 1, -1, 0, 0, 0,
      0, t_1-2*t_2, 0}, {0, 0, 0, 1, -1, 0, 0, 0, 0, t_1}})
assert(splineMatrix(V,F,E,0,Homogenize=>false) == matrix {{1, 0, 0, 0, -1, t_2, 0, 0, 0, 0}, {1, -1, 0, 0, 0, 0, t_1-t_2,
      0, 0, 0}, {0, 1, -1, 0, 0, 0, 0, t_1+t_2, 0, 0}, {0, 0, 1, -1, 0, 0, 0,
      0, t_1-2*t_2, 0}, {0, 0, 0, 1, -1, 0, 0, 0, 0, t_1}})
assert(splineMatrix(V,F,E,1) == matrix {{1, 0, 0, 0, -1, t_2^2, 0, 0, 0, 0}, {1, -1, 0, 0, 0, 0,
      t_1^2-2*t_1*t_2+t_2^2, 0, 0, 0}, {0, 1, -1, 0, 0, 0, 0,
      t_1^2+2*t_1*t_2+t_2^2, 0, 0}, {0, 0, 1, -1, 0, 0, 0, 0,
      t_1^2-4*t_1*t_2+4*t_2^2, 0}, {0, 0, 0, 1, -1, 0, 0, 0, 0, t_1^2}})
assert(isHereditary(F,E) === true)
///

TEST ///
V={{-5,0},{-3,0},{-1,-4},{-1,4},{-1,-2},{-1,2},{0,-1},{0,1},{1,-2},{1,2},{1,-4},{1,4},{3,0},{5,0}}
F={{0, 1, 4, 2}, {0, 1, 5, 3}, {8, 10, 13, 12}, {9, 11, 13, 12}, {1, 4, 6, 7, 5}, {2, 4, 6, 8, 10}, {3, 5, 7, 9, 11}, {6, 7, 9, 12, 8}}
E={{0, 1}, {0, 2}, {0, 3}, {1, 4}, {1, 5}, {2, 4}, {2, 10}, {3, 5}, {3, 11}, {4, 6}, {5, 7}, {6, 7}, {6, 8}, {7, 9}, {8, 10}, {8, 12}, {9, 11}, {9, 12}, {10, 13}, {11, 13}, {12, 13}}
assert(splineMatrix(V,F,E,0) == matrix {{1, -1, 0, 0, 0, 0, 0, 0, t_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0}, {1, 0, 0, 0, -1, 0, 0, 0, 0, 3*t_0+t_1+t_2, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 3*t_0+t_1-t_2, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0,
      t_0+t_1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, -1, 0, 0,
      0, 0, 0, t_0+t_1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, -1, 0,
      0, 0, 0, 0, 0, 0, t_0-t_1+t_2, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0,
      1, 0, -1, 0, 0, 0, 0, 0, 0, 0, t_0-t_1-t_2, 0, 0, 0, 0, 0, 0, 0, 0}, {0,
      0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, t_1, 0, 0, 0, 0, 0, 0, 0}, {0,
      0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, t_0+t_1+t_2, 0, 0, 0, 0, 0,
      0}, {0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_0+t_1-t_2, 0,
      0, 0, 0, 0}, {0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      t_0-t_1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 3*t_0-t_1+t_2, 0, 0, 0}, {0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, t_0-t_1, 0, 0}, {0, 0, 0, 1, 0, 0, 0, -1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3*t_0-t_1-t_2, 0}, {0, 0, 1, -1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_2}})
assert(splineMatrix(V,F,E,0,Homogenize=>false) == matrix {{1, -1, 0, 0, 0, 0, 0, 0, t_2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0}, {1, 0, 0, 0, -1, 0, 0, 0, 0, t_1+t_2+3, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0}, {0, 1, 0, 0, -1, 0, 0, 0, 0, 0, t_1-t_2+3, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0}, {1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, t_1+1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, t_1+1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0,
      t_1-t_2-1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0,
      0, 0, 0, t_1+t_2-1, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, -1, 0,
      0, 0, 0, 0, 0, 0, t_1, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, -1, 0,
      0, 0, 0, 0, 0, 0, 0, t_1+t_2+1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 1,
      -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_1-t_2+1, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0,
      -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_1-1, 0, 0, 0, 0}, {0, 0, 1, 0,
      0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_1-t_2-3, 0, 0, 0}, {0, 0,
      0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_1-1, 0, 0}, {0,
      0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t_1+t_2-3,
      0}, {0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      t_2}})
assert(isHereditary(F,E) === true)
///


end

