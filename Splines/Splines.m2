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
   "CheckHereditary"
    }

------------------------------------------
------------------------------------------
-- Methods
------------------------------------------
------------------------------------------

splineMatrix = method(Options => {symbol InputType => "ByFacets", symbol CheckHereditary => false})
------------------------------------------
------------------------------------------
--Inputs: 
------------------------------------------
--verts = list of coordinates of vertices
--facets = ordered lists of facets
--edges = list of edges
--r = degree of desired continuity
------------------------------------------

splineMatrix(List,List,List,ZZ) := Matrix => opts -> (verts,facets,edges,r) -> (
    if opts.InputType === "ByFacets" then (
	if opts.CheckHereditary === true then (
	    --put hereditary check here.
	    )
	--put remainder of code for ByFacets Here
	)
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
	print "Need list of vertices, facets and edges, along with continuity r."
	);
    --If user DOES want to define complex by regions and dual graph.
    if opts.InputType === "ByLinearForms" then (
    m := max flatten B;
    A := matrix apply(B, i-> apply(toList(0..m), j-> if (j=== first i) then 1 else if (j===last i) then -1 else 0));
    D := matrix apply(#L, i-> apply(#L, j-> if i===j then L_i^(r+1) else 0));
    A|D
    )
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