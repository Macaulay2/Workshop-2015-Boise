loadPackage("ToricMaps")
--NormalToricVarieties already exists, but doesn't have a type for maps between toric varieties.
--our goal was to implement a ToricMaps type with some reasonable functionality.

--To show what we've done so far, let's look at some examples. Take the blowup of P2 at its three torus-fixed points
BlownUpP2 = normalToricVariety({{1,0},{1,1},{0,1},{-1,0},{-1,-1},{0,-1}},{{0,1},{1,2},{2,3},{3,4},{4,5},{0,5}})
	
--This maps to P1 x P1
P1crossP1 = normalToricVariety({{1,0},{0,1},{-1,0},{0,-1}},{{0,1},{1,2},{2,3},{0,3}})
Bl = toricMap(P1crossP1,BlownUpP2,matrix{{1,0},{0,1}})

--and P1 x P1 maps to P3 by the Segre embedding
P3 = projectiveSpace(3)
segre = toricMap(P3,P1crossP1,matrix{{1,1},{1,0},{0,1}})

--we can compose maps of toric varieties
compose(segre, Bl)

--(or, alternatively)
segre @@ Bl

--Let D be the divisor of a hyperplane section on P3
D = toricDivisor({1,0,0,0},P3)
--We can pull back along the segre embedding
pullback(segre,D)
--(or, alternatively)
segre ^* D

isAmple segre^*D

--you could also pullback all the way through the composition to get a divisor on the blown up P2
(segre @@ Bl)^* D
--(or, by functoriality)
Bl^* segre^* D

--We can also test whether a toricMap is an isomorphism, and if it is we can get its inverse
--Here's P2
P2 = projectiveSpace(2)

--Here's some other 2-dimensional toric variety
mysteryToricVariety = normalToricVariety({{0,1},{-1,0},{1,-1}},{{0,1},{1,2},{0,2}})
	
--I can find a map between them:
f = toricMap(mysteryToricVariety,P2,matrix{{-1,0},{0,1}})

--and observe:
isIsomorphismToricMap(f)
inverseToricMap(f)

--future / possibilities:

--we know how isProper should work, and have a running implementation, but we haven't thoroughly tested it

--There are many methods in the NTV package that could easily be modified to return maps,
--such as blowup, makeSmooth, makeSimplicial, etc. This shouldn't be hard

--and more.
