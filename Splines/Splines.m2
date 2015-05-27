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
--Inputs: 
------------------------------------------
--B = list of regions that are adjacent
--L = list of ordered linear forms that separate
--regions, given in order of input B
--r = degree of desired continuity
------------------------------------------


splineMatrix(List,List,ZZ) := Matrix => opts -> (B,L,r) ->(
    if opts.InputType === "ByFacets" then (
	print "Need list of vertices, facets and edges, along with continuity r."
	);
    if opts.InputType === "ByLinearForms" then (
    m := max flatten B;
    A := matrix apply(B, i-> apply(toList(0..m), j-> if (j=== first i) then 1 else if (j===last i) then -1 else 0));
    D := matrix apply(#L, i-> apply(#L, j-> if i===j then L_i^(r+1) else 0));
    A|D
    )
)

end

    K := ker AD;
    b := max flatten B;
    submatrix(gens K, toList(0..b),)
