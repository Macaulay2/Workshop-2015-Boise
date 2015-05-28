needsPackage("NormalToricVarieties");

-----------
--NEW TYPES
-----------

-- defining the new type ToricMap
ToricMap = new Type of HashTable 
ToricMap.synonym = "toric map"
globalAssignment ToricMap

--------------
--CONSTRUCTORS
--------------

-- PURPOSE: construct a map between two (normal) toric varieties
-- INPUT: (Y,X,M) X and Y NormalToricVarities, M a ZZ-matrix
-- OUTPUT: a ToricMap from X to Y
-- COMMENTS: The matrix M describes a fan-compatible linear map from the one-parameter subgroup lattice 
-- of X to the one-parameter subgroup lattice of Y

toricMap = method(Options => {})


-- toricMap still needs a check for whether the map of lattices is
-- compatible with the fans defining the toric varieties.

toricMap (NormalToricVariety, NormalToricVariety, Matrix) := opts -> (Y,X,M) -> (

		m := dim Y;
		n := dim X;
		if not ((numRows M == m) and (numColumns M == n)) then (
			error ("expected a "| toString m | "x" | toString n | " matrix.");
			)
		else(
			new ToricMap from {
			symbol target => Y,
			symbol source => X,
			symbol matrix => M
			} 
			)
		)

source ToricMap := NormalToricVariety => f -> f.source
target ToricMap := NormalToricVariety => f -> f.target
matrix ToricMap := Matrix => o -> f -> f.matrix

compose = method(Options => {})

-- composing maps
compose (ToricMap, ToricMap) := opts -> (f,g) -> (
		
		if (not target g === source f) then error "unmatched domains"
		else return toricMap(source g, target f, (matrix f)*(matrix g))

		)
-- @@ operator (should it be * instead of @@? I don't think so)
ToricMap @@ ToricMap := ToricMap => (f,g) -> compose(f,g)

-- isIsomorphism: (just check if an inverse exists. is M^-1 a ZZ matrix
					and )