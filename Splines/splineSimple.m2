--To open Macaulay2, go to a terminal in SageCloud and type M2
--In the command line of M2, type load("splinetest.m2")
R = QQ[x,y]

----------------------------------------
-- Change the B and L for your examples.
----------------------------------------
--Sample Spline Example
--B is a list of adjacencies of regions
B = {{0,1},{1,2},{2,3},{3,4},{4,0}}
--L is the list of linear forms dividing regions, in order of B
L = {x-y,y,x,y-2*x,x+y}
---------------------------------------

splineMaker = method()

splineMaker(List,List) := (B,L) -> (
    m := max flatten B;
    A := matrix apply(B, i-> apply(toList(0..m), j-> if (j=== first i) then 1 else if (j===last i) then -1 else 0));
    D := matrix apply(#L, i-> apply(#L, j-> if i===j then L_i else 0));
    AD := A|D;
    K := ker AD;
    b := max flatten B;
    submatrix(gens K, toList(0..b),)
)

splineMaker(List,List,ZZ) := (B,L,r) -> (
    m := max flatten B;
    A := matrix apply(B, i-> apply(toList(0..m), j-> if (j=== first i) then 1 else if (j===last i) then -1 else 0));
    D := matrix apply(#L, i-> apply(#L, j-> if i===j then L_i^(r+1) else 0));
    AD := A|D;
    K := ker AD;
    b := max flatten B;
    submatrix(gens K, toList(0..b),)
)


S = splineMaker(B,L,1)
H = hilbertSeries(image S)

end


--Sample Spline Example
B = {{0,1},{1,2},{2,3},{3,0}}
L = {x,y,x,y}
S = splineMaker(B,L,1)
hilbertSeries(image S)

--Sample Spline Example
B = {{0,1},{0,2},{1,2},{1,3},{2,3}}
L = {y,y,x,x,y+z}
S = splineMaker(B,L,1)
H = hilbertSeries(image S)
