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
-- COMMENTS:

toricMap = method(Options => {})

toricMap (NormalToricVariety, NormalToricVariety, Matrix) := opts -> (Y,X,M) -> (
		new ToricMap from {
		symbol Target => Y,
		symbol Source => X,
		symbol Matrix => M
		}
		)