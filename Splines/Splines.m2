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
    needsPackage "Polyhedra",
    needsPackage "Graphs"
    )

newPackage select((
    "Splines",
        Version => "0.1.0", 
        Date => "27. May 2015",
        Authors => {
            {Name => "Mike DiPasquale", Email => "midipasq@gmail.com", HomePage => "http://illinois.edu/~dipasqu1"},
            {Name => "Gwyn Whieldon", Email => "whieldon@hood.edu", HomePage => "http://cs.hood.edu/~whieldon"}
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
   "splineMatrix",
   "splineModule",
   "InputType",
   "ByFacets",
   "ByLinearForms",
   "isHereditary",
   "CheckHereditary",
   "Homogenize",
   "VariableName",
   "getCodimIFacesPolytope",
   "getCodimIFacesSimplicial",
   "interiorFaces"
    }

------------------------------------------
------------------------------------------
-- Methods
------------------------------------------
------------------------------------------


------------------------------------------
------------------------------------------
isHereditary= method()
------------------------------------------
------------------------------------------
-- This method checks if the polyhedral
-- complex with facets and edges (F,E)
-- is hereditary.
------------------------------------------
--Inputs: 
------------------------------------------
--F = ordered lists of facets
--E = list of edges
------------------------------------------
--Outputs:
--Boolean, if complex is hereditary
------------------------------------------

isHereditary(List,List) := Boolean => (F,E) -> (
    V := unique flatten join F;
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
--E' = list of interior edges
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
getCodim1FacesPolytope = method()
------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--F = list of facets of a polytope
------------------------------------------
--Outputs:
-----------------------------------------
--E = list of (interior) codim 1 faces
-----------------------------------------

getCodim1FacesPolytope(List) := List => F ->(
    --This function ASSUMES that the polytopal 
    --complex considered is hereditary.
    n := #F;
    --For each pair of facets, take their intersection:
    intersectFacets := unique flatten apply(#F-1, i-> 
	apply(toList(i+1..#F-1), 
	    j-> sort select(F_i, 
		v-> member(v,F_j))));
    --Remove any non-maximal faces in this intersections:
    select(intersectFacets, f -> (
    	(number(intersectFacets, g-> all(f, j-> member(j,g)))) === 1
    ))
)

------------------------------------------
------------------------------------------
getCodimIFacesPolytope = method()
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
getCodimIFacesPolytope(List,ZZ) := List => (F,d) ->(
    Fcodim := F;
    --Compute interior codime
    apply(d, i-> Fcodim = getCodim1FacesPolytope(Fcodim));
    Fcodim
    )

------------------------------------------
------------------------------------------
getCodimIFacesSimplicial = method()
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
getCodimIFacesSimplicial(List,ZZ) := List => (F,i) -> (
    d := getSize(F);
    unique flatten apply(F, f-> subsets(f,d-i))
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
splineMatrix(List,ZZ) := Matrix -> opts -> (L,r) -> (
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
	BM := matrix apply(facetEdgeH, i-> apply(#F, j-> if (j === first i) then 1 else if (j===last i) then -1 else 0));
	--To homogenize, we append a 1 as the final coordinate of each vertex coord in list V.
	--If not homogenizing, leave vertex coordinates V as is.
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
	minorList := apply(E, e-> gens gb minors(mM,M_e|varCol));
	if any(minorList, I-> ideal I === ideal 1) then (
	    error "Some vertices on entered face are not in codimension 1 face."
	    );
	T := diagonalMatrix(flatten apply(minorList, m -> (m_(0,0))^(r+1)));
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
	    E := getCodim1Intersections(F);
	    SM := splineMatrix(V,F,E,r,opts)  
	    )
	else(
	    E = getCodim1FacesPolytope(F);
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
    	submatrix(gens K, toList(0..b-1),)
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
    	submatrix(gens K, toList(0..b-1),)
)

------------------------------------------
------------------------------------------
getCodim1Intersections = method();
------------------------------------------
------------------------------------------
--Code to Compute Codim 1 Intersections:
------------------------------------------
--Inputs: 
------------------------------------------
--F = facets of a pure simplicial complex 
--(as lists of vertices)
------------------------------------------
--Outputs:
------------------------------------------
--the codimension-1 (interior) intersections
------------------------------------------
getCodim1Intersections(List) := List => F ->(
    n := #F;
    d := #(F_0);
    --For each non-final facet, construct all codimension 1 subsets.
    codim1faces := apply(n-1, i -> subsets(F_i,d-1));
    --Check if a codimension 1 subset is contained in another facet,
    --store it as a codim 1 intersection.
    sort flatten apply(#codim1faces, i -> 
	select(codim1faces_i, 
	    s -> any(F_{i+1..n-1}, 
		f-> all(s, v-> member(v,f)))))
)

------------------------------------------
getSize = method();
------------------------------------------
--Input: L = List of Lists
------------------------------------------
--Output: If all lists in L are same size,
-- the length of each individial list in L
------------------------------------------
getSize(List) := ZZ => L ->(
    if all(L, v-> #v == #(L_0)) then #L_0 else null
)


------------------------------------------
issimplicial = method();
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
    n := getSize(V);
    f := getSize(F);
    if not instance(n, Nothing) and not instance(f,Nothing) and n + 1 == f then true
    else(
	if instance(n, Nothing) then print "Vertices have inconsistent dimension."
	else false
    )
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
        Splines
    Headline
        a package for building splines and computing bases
    Description
        Text
            This package computing topological boundary maps and piecewise 
	    continuous splines on polyhedral complexes.
	Text
	    @SUBSECTION "Definitions of Splines"@
	    Given a complex $\Delta\subset{\mathbb R}^d$, a spline 
	    $f \in S^{r+1}_d(\Delta)$ is a piecewise polynomial function on $\Delta$, 
	    such that $f$ is of smoothness $r$ ($f \in C^r$) and $f$ restricts to a degree $d$ polynomial on each face 
	    $\sigma\in\Delta$.
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
	    list of oriented facets of Delta
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
	    its kernel. Note that the ambient ring of the appropriate
	    dimension needs to be defined.
	Example
	    R = QQ[x,y,z]
	    V = {{0,0},{1,0},{1,1},{-1,1},{-2,-1},{0,-1}};-- the coordinates of vertices
            F = {{0,2,1},{0,2,3},{0,3,4},{0,4,5},{0,1,5}};  -- a list of facets (pure complex)
            E = {{0,1},{0,2},{0,3},{0,4},{0,5}};   -- list of edges in graph
    	    splineMatrix(V,F,E,1)
        Text
            Alternately, spline matrices can be created directly from the
	    dual graph (with edges labeled by linear forms).
	Example
	    R = QQ[x,y]
	    B = {{0,1},{1,2},{2,3},{3,4},{4,0}}
	    L = {x-y,y,x,y-2*x,x+y}
	    splineMatrix(B,L,1,InputType=>"ByLinearForms")

///

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
