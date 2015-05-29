newPackage ("ToricMaps",
	Authors => {
		"Nathan Bliss",
		"Nathan Fieldsteel",
		"Robert Walker",
		"Jeff Poskin"
		},
	PackageExports => {
		"NormalToricVarieties"
		}
	)

--needsPackage("NormalToricVarieties");
load "coneContain.m2";

export{"ToricMap","checkCompatibility","toricMap","pullback","isIsomorphism","inverse"
		}

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

toricMap = method(Options => {checkCompatibility => true})


-- toricMap still needs a check for whether the map of lattices is
-- compatible with the fans defining the toric varieties.

toricMap (NormalToricVariety, NormalToricVariety, Matrix) := ToricMap => opts -> (Y,X,M) -> (

		n := dim X;
		m := dim Y;
		if not ((numRows M == m) and (numColumns M == n)) then (
			error ("expected a "| toString m | "x" | toString n | " matrix.");
			)
		else if (opts.checkCompatibility and not isCompatible(Y,X,M)) then (
			error "Lattice map not compatible with fans."
			)
		else(
			new ToricMap from {
			symbol target => Y,
			symbol source => X,
			symbol matrix => M
			} 
			)
		)

net ToricMap := f -> (
	m := net matrix f;
	w := width m+1;
	line := concatenate(apply((1..w), t -> "-"));
	arr := net target f | " <"|line|" "|net source f;
	w2 := width (net target f) + 2;
	sp := concatenate(apply((1..w2), t -> " "));
	return arr||(sp|m);
	)

source ToricMap := NormalToricVariety => f -> f.source
target ToricMap := NormalToricVariety => f -> f.target
matrix ToricMap := Matrix => o -> f -> f.matrix

compose := method()

-- composing maps
compose (ToricMap, ToricMap) := opts -> (f,g) -> (
		
		if (not target g === source f) then error "unmatched domains"
		else return toricMap(source g, target f, (matrix f)*(matrix g))

		)
-- @@ operator (should it be * instead of @@? I do not think so)

ToricMap @@ ToricMap := ToricMap => (f,g) -> compose(f,g)

cartierCoefficients := method()


-- Taken from NormalToricVarieties.m2 (not exported)
cartierCoefficients ToricDivisor := List => D -> (
	X := variety D;
	V := matrix rays X;
	a := matrix vector D;
	return apply(max X, s -> a^s // V^s)
	)

pullback = method()

pullback (ToricMap, ToricDivisor) := ToricDivisor => (f, D) -> (

		--check is D = divisor on target f?

		X2 := variety D;
        if X2 =!= target f then (
            error "variety of "|D|" is not the same as the target of"|f;
        );
		cdat := cartierCoefficients(D);
		maximalCones := max X2;
		numCones := length maximalCones;
        l := toList(apply(0..(numCones-1), i -> (maximalCones_i,cdat_i)));
		cartierDict := hashTable(l);
        pullbackDict := new MutableHashTable;
        for C in max source f do (
            coneC := posHull ((matrix f)*(transpose(matrix ((rays(source f))_C))));
            imC := null;
            for imCone in max target f do (
                polyCone := posHull(transpose(matrix ((rays(target f))_imCone)));
                if contains(polyCone,coneC) then (
                    imC = imCone;
                    break;
                );
            );
            --apply the transpose of matrix f to the cart data
            --associated to that maxl cone.
            pullbackDict#C = (transpose(matrix f))*(cartierDict#imC);
        );
        listContains := (l,x) -> (
            for i in l do (if x==i then return true);
            return false;
        );

		--then use that cart data to get the divisor coeffs
		--for the rays of source f.
        pullbackCoeffs := {};
        numRays := #(rays(source f)) -1;
        for i in 0..numRays do (
            j:=0;
            for k in keys(pullbackDict) do (
                if listContains(k,i) then (
                    pullbackCoeffs = append(pullbackCoeffs,((matrix({(rays(source f))_i})*(pullbackDict#k)))_(0,0));
                    break;
                );
            );
        );
        return toricDivisor(pullbackCoeffs,source f);
    );

ToricMap ^* := f -> D -> pullback(f,D)

isIsomorphism = method()

isIsomorphism (ToricMap) := Boolean => f -> (
	m := matrix f;
	d := det m;
	if not (d == 1 or d == -1) then (
		return false;
		)
	else (
		minv := inverse m;
		X := source f;
		Y := target f;
		if not (isCompatible(X,Y,minv)) then (
			return false;
			);
	return true;
		);
	)

inverse = method()

inverse (ToricMap) := (ToricMap) => f -> (
	if not isIsomorphism(f) then (
		error "map is not invertible."
		)
	else (
		return toricMap(source f, target f, inverse matrix f);
	)
	)

end

