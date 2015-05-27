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
    needsPackage "Polyhedra"
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
        if version#"VERSION" > "1.4" then PackageExports => {}
        ), x -> x =!= null)

if version#"VERSION" <= "1.4" then (
    needsPackage "Polyhedra"
    )

export {
   "splineMatrix",
   "splineModule",
   "InputType",
   "ByFacets",
   "ByLinearForms",
   "CheckHereditary",
   "Homogenize",
   "VariableName"
    }

------------------------------------------
------------------------------------------
-- Methods
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
isHereditary(List,List) := Boolean => (F,E) -> (
    V := unique flatten join F;
    dualV := toList(0..#F-1);
    dualE := apply(#E, e-> positions(F, f-> all(E_e,v-> member(v,f))));
    if not all(dualE,e-> #e===2) then (
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

<<<<<<< Updated upstream
splineMatrix = method(Options => {symbol InputType => "ByFacets", symbol CheckHereditary => false})

=======
splineMatrix = method(Options => {
	symbol InputType => "ByFacets", 
	symbol CheckHereditary => false, 
	symbol Homogenize => true, 
	symbol VariableName => getSymbol "t",
	symbol CoefficientRing => QQ})
>>>>>>> Stashed changes
------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--verts = list of coordinates of vertices
--facets = ordered lists of facets
--edges = list of edges
--r = degree of desired continuity
------------------------------------------

splineMatrix(List,ZZ) := Matrix -> opts -> (L,r) -> (
    if opts.InputType === "ByFacets" then (
	splineMatrix(L_0,L_1,L_2,r)
	)
)

splineMatrix(List,List,List,ZZ) := Matrix => opts -> (V,F,E,r) -> (
    if opts.InputType === "ByFacets" then (
	if opts.CheckHereditary === true then (
	    --put hereditary check here.
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
	if opts.Homogenize === true then (
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
	    error return "Some vertices on entered face are not in codimension 1 face."
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
-- splineMatrix "ByLinearForms"
------------------------------------------
-- This method (essentially) inputs your
-- simplicial complex by the dual graph
-- (with edges labeled by linear forms
-- separating regions.)
--Inputs: 
------------------------------------------
--B = list of regions that are adjacent
--L = list of ordered linear forms that separate
--regions, given in order of input B
--r = degree of desired continuity
------------------------------------------


splineMatrix(List,List,ZZ) := Matrix => opts -> (B,L,r) ->(
    --Warn user if they are accidentally using ByFacets method with too few inputs.
    if opts.InputType === "ByFacets" then (
	--Function should compute E automatically, pretending it's simplicial or polytopal
	
	--Write function to compute E (given S or P complexes) here.
	--splineMatrix(B,L,E,r)
	print "'ByFacets' option not implemented yet for inputs (V,F,r)."
	);
    --If user DOES want to define complex by regions and dual graph.
    if opts.InputType === "ByLinearForms" then (
    m := max flatten B;
    A := matrix apply(B, i-> apply(toList(0..m), j-> if (j=== first i) then 1 else if (j===last i) then -1 else 0));
    D := matrix apply(#L, i-> apply(#L, j-> if i===j then L_i^(r+1) else 0));
    A|D
    )
)

splineModule = method(Options => {symbol InputType => "ByFacets", symbol CheckHereditary => false})

splineModule(List,List,List,ZZ) := Matrix => opts -> (verts,facets,edges,r) -> (
    )


--Interior Methods used in SplineMatrix--

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
            @SUBSECTION "Other acknowledgements"@
            --
            Methods in this package are put together from code written by Hal Schenck
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


end

    K := ker AD;
    b := max flatten B;
    submatrix(gens K, toList(0..b),)

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