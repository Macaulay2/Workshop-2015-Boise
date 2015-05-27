-- Copyright 2014-2015: James Parson
-- You may redistribute this file under the terms of the GNU General Public
-- License as published by the Free Software Foundation, either version 2
-- of the License, or any later version.

------------------------------------------
------------------------------------------
-- Header
------------------------------------------
------------------------------------------

if version#"VERSION" <= "1.4" then (
    )

newPackage select((
    "Loci",
        Version => "0.0.1", 
        Date => "27. January 2015",
        Authors => {
            {Name => "James Parson", Email => "parson@hood.edu"},
            {Name => "Gwyn Whieldon", Email => "whieldon@hood.edu", HomePage => "http://www.hood.edu/Academics/Departments/Mathematics/Faculty/Gwyneth-Whieldon.html"}
        },
        Headline => "Package for computing open Loci for properties or rings and homomorphisms",
        Configuration => {},
        DebuggingMode => true,
        if version#"VERSION" > "1.4" then PackageExports => {}
        ), x -> x =!= null)

if version#"VERSION" <= "1.4" then (
    )

export {
    --
    -- methods in package
--    "presentGraph",
--    "extendScalars",
--    "genericFreeness",
--    "absFrob",
--    "inRad",
--    "isOpenSubLocus",
--    "equalOpenLoci",
--    "Omega1",
--    "naiveCotComplex",
--    "LSCotComplex",
--    "relativeFlatLocusINT",
--    "semiContinuity",
--    "isCM",
--    "kunzRegularLocus",
--    "equidimensionalINT",
--    "equidimensional",
--    "genericSaturate",
--    "primeFactors",
--    "regularLocusZZ",
--    "regularLocusFiber",
    "flatLocus",
    "relativeFlatLocus",
    "relativeDimension",
    "unramifiedLocus",
    "smoothLocus",
    "etaleLocus",
    "lciLocus",
    "CMLocus",
    "GorensteinLocus",
    "CILocus",
    "regularLocus",
    "kunzRegularLocus",
    "trueCodimension",
    "isNormalFixed",
    "codimPiece",
    "dimPiece",
    "isOpenSubLocus",
    "equalOpenLoci",
    "normalLocus",
    "isReduced",
    "reducedLocus",
    "isS",
    "isR"
    }

------------------------------------------
------------------------------------------
-- Methods
------------------------------------------
------------------------------------------

-- to do:
-- * catch errors in "presentGraph" : I know what I would like to do, but I am still having trouble
--     because of rings without coefficient rings.
-- * think about types for various methods. Perhaps "Ring" is not the right thing.
--     in the M2 source code, one often find Ring, PolynomialRing, QuotientRing, etc.
--     for a simple example, search source code for "isAffineRing" (in "M2/galois.m2" for some reason)
-- * add something to find reduced locus, normal locus.


-- * could we do something with constructible loci? image---or locus where dim has particular values
--     is constructible, as follows directly from generic freeness.
--     perhaps this should wait for another time, since the other functions are basically done.

------
-- General methods for manipulating rings
-----


presentGraph = method()

-- the error-catching business below does not work when "A" does not have a coefficientRing.
-- I don't know the proper way to handle this situation. Perhaps I should ask about it on the M2
-- mailing list. (But to do that, I will have to overcome my paralyzing fear of setting up a Google
-- account.)

presentGraph(RingMap) := f -> (
-- The function takes a homomorphism f:A-->B as input. The rings A and B must have the same coefficientRing R,
-- and the homomorphism f must be an R-algebra homomorphism.
-- The function returns an A-algebra D, which is B viewed as A-algebra via f, along with isomorphisms ii:B-->D and phi:D-->B.
-- To do: add code checking to make certain that the hypotheses hold.
-- See tensor product code for how to make missing checks on R and f.
--    if B =!= ring M
--        then error "expected target of RingMap to match ring of module";
--    R := coefficientRing B; -- assume that this is common coefficient ring of A and B and that f is an R-algebra homomorphism.
--    if R =!= coefficientRing A
--        then error "expected source and target of RingMap to have same coefficientRing";
    A := source f;
    B := target f;
    R := coefficientRing B;
--    if not ( R === coefficientRing A ) then error "expected source and target of ring map to have same coefficient ring";
--    gensR := generators(R, CoefficientRing => ZZ);
--    if not all(gensR, x -> promote(x,B) == f promote(x,A)) then error "expected ring map to be identity on coefficient ring";
    (B1,s) := flattenRing(B, CoefficientRing => R);
    C := A monoid ambient B1;
    i := map(C, ambient B1, gens C); -- important here to add the "gens C"; to see why try a homomorphism f from A to A such as frobenius
    j := map(C,A);
    I := i(ideal B1);
    II := ideal(apply(gens A, r->j(r) - i(lift(s(f(r)),ambient B1))));
    D := C/(I + II);
    ii := map(D,B1); -- this homomorphism is an isomorphism, but the coefficientRing of D is A; ** does this work? correct to be like "i" above?
    phi := map(B1,D,apply(gens B, r->s(r))|apply(gens A, r->s(f(r)))); -- this homomorphism is the inverse of ii
    (D,ii*s,(s^-1)*phi)
)

extendScalars = method()

extendScalars(Ring, RingMap) := (B,phi) -> (
-- input: an A-algebra B and a ring map A->C. Returns C-algebra obtained by extension of scalars and
--        extension of scalars homomorphism j:B-->this C-algebra
    A := source phi; -- we assume that A is also the coefficientRing of B; we should check this!
    C := target phi;
    (B1,r) := flattenRing(B, CoefficientRing => A);
    P := ambient B1;
    I := ideal B1;
    Q := C (monoid P);
    i := map(Q, P, gens Q | apply(gens A, t->(phi(t))));
    F := Q/i(I);
    j := map(F, B1, gens Q | apply (gens A, t->(phi(t))));
    (F,j*r)
)

---The next two methods are for flattening rings constructed with fraction fields
---fRing gives back a ring C over ZZ, QQ, or ZZ/p along with a map C->A that is a localization.
---To compute things about A, we may compute things about C and push forward along the localization map.
---This procedure is necessary for computing the regular locus when the coefficientRing of A
---is a field that is not perfect, since then regularity and smoothness over the coefficientRing do
---not coincide.

fRing = method()

fRing(Ring) := A -> (
    (B,f) := flattenRing A;
    if not instance(coefficientRing B, FractionField) then (B,f^-1)
    else (
	K := coefficientRing B;
	R := last K.baseRings;
	P := ambient B;
	I := ideal B;
	Q := R monoid P;
	bad := map(Q,P);
	J := ideal(apply(flatten entries gens I, x->bad(clearDenom x))); -- I read somewhere that I_* does the same as "fl ent gens I"
	g := map(B, Q/J);
	(C,i) := fRing(Q/J);
	(C, (f^-1)*g*i)
	)
)

clearDenom = method()

clearDenom(RingElement) := x -> (
    P := ring x;
    K := coefficientRing P;
    (M,D) := coefficients x;
    d := product apply(flatten entries D, a->denominator lift(a,K));
    d*x
)

genericFreeness = method()

genericFreeness(Ideal) := K -> ( -- We no longer use this method, but it reminds us how the module method should look
    G := flatten entries gens gb K;
    apply(G,leadCoefficient)
)


genericFreeness(Module) := M -> (
-- module/polynomial ring B with coefficient ring A: returns (possibly empty) list {f_1,...,f_r} of nonzero elements of B with M_f free over A_f
-- if the list is empty, then M is free over A.
-- we should just be able to use "leadCoefficient", but the method does not seem to work properly
    m := gens gb presentation M;
    n := numgens source m - 1;
    is := apply(toList(0..n),i->leadCoefficient m_((leadComponent m)_i,i));
    hs := toList set is; -- remove duplicates to cut down on recursive calls
    select(hs, f -> f!=1) -- remove 1's from the list to save silly recursive calls
)

genericFreenessFP = method()

genericFreenessFP(Module) := M -> (
-- module over B that is finitely presented over A: returns list as above
-- something bad is going on here, since I copied this code from the relativeFlatLocus code below
    B := ring(M);
    A := coefficientRing(B);
    C := ambient(B); -- C is a polynomial ring over A
    J := ideal B; -- We have C/J = B, presenting B as A-algebra
    m := presentation M; -- Next we lift M to a module over C
    N := (cokernel lift(m,C))/J; -- This C-module N is M restricted to C along C-->B
    genericFreeness N
)

absFrob = method()

-- returns the absolute Frobenius endomorphism of A, if A has characteristic p for some p)

absFrob(Ring) := A -> (
    p := char A;
    if not isPrime p then
        error "expected a ring with prime characteristic"
    else map(A,A,apply(gens(A,CoefficientRing=>ZZ),x->x^p))
)

---
-- partially ordered set of open loci
---

inRad = method()

inRad(RingElement, Ideal) := (a,I) -> (
--check whether the ring element a is in the radical of I without computing the radical.
--one should check to see that a and I have the same ring.
--t := symbol t; --This worked! Hooray! This is not the way I've used this in the past...
A := ring a;
B := A (monoid [Variables=>1]);
promote(I,B)+ideal(promote(a,B)*B_0-1) == B
)

isOpenSubLocus = method()

isOpenSubLocus(Ideal,Ideal) := (I,J) -> (
-- assume I and J are ideals in the same ring A---one should test for this!
-- check whether the open locus defined by I is a subset of the open locus defined by J
if (ring I) =!= (ring J) then error "expected ideals in the same ring";
all(flatten entries gens I, x->inRad(x,J)) -- perhaps "I_*" in place of "flatten entries ... "
)

equalOpenLoci = method()

equalOpenLoci(Ideal,Ideal) := (I,J) -> isOpenSubLocus(I,J) and isOpenSubLocus(J,I)




---
-- Flat locus methods
----

-- perhaps add a "verbose" flag to activate printing of evidence for flatness.

-- Almost everything breaks if we try to apply it to a module over a fraction field B, since fraction fields have no ambient rings.
-- One should catch this possibility and give an error message instead of crashing into the debugger.

flatLocus = method()

flatLocus(Module) := M -> trim ann Ext^1(trim M, trim image presentation M)
-- I put the "trim" in, since I was getting an error "wrong number of rows or columns"
-- when I called this function while trying to compute the regular locus for ZZ[x]/(x^3 - 2)

relativeFlatLocusINT = method()

relativeFlatLocusINT(Module, Ideal) := (M,I) ->(
-- this function can be horribly slow!
-- The recursive calls often seem to calculate the same things over and over.
-- Is there a way to cache values to speed up the recursion? Surely there must be a way to do this...
--    print I; -- remove this comment for an amusing look at the horrors of recursion
    B := ring(M); -- we are assuming B is a polynomial ring in what follows.
    A := coefficientRing(B); -- This line does not work, if B is ZZ,QQ,ZZ/p, GF(p,n),RR, CC, etc.
    if I == A then (ideal(1_B))
    else (
        C := (A/I) monoid B; -- is this the right way to do it?
        fs := apply(genericFreeness(C**M), f -> lift(f,A));
        xs := apply(fs, f -> trim relativeFlatLocusINT(M, trim(I+ideal(f)))); -- probably too much "trim"ming!
        L1 := if xs=={} then (ideal(1_B)) else trim (intersect(xs)); -- empty intersection needs to be handled separately
        L2 := trim ann(HH_1((chainComplex(map(A^1, module I, gens I))**B)**M)); -- The homology is Tor_1^A(A/I,M)
        trim intersect(L1,L2)
    )
)

relativeFlatLocusZZ = method()

--because Macaulay 2 only handles ZZ/n for n prime, we need a special
--method to handle flatness over ZZ.

relativeFlatLocusZZ(Module) := M ->(
--internal method: we assume that M is a module over a polynomial ring over ZZ
    B := ring(M);
    if 0_B == 1_B then ideal(1_B);
    f := product(genericFreeness(M)); --nonzero element of ZZ such that M is free over ZZ[1/f]
    ann(ker(map(M,M,{{f}}))) -- M is flat over ZZ if and only if it is torsion free.
                             -- The only torsion is f-power torsion. As long as the f-torsion
                             -- submodule vanishes, M is flat
)


relativeFlatLocus = method()

relativeFlatLocus(Module) := M ->(
    B := ring(M);
    A := coefficientRing(B);
    (D,f) := flattenRing(B,CoefficientRing => A); --we do this for safety
    C := ambient(D); -- C is a polynomial ring over A
    J := ideal D; -- We have C/J = D, presenting D as A-algebra
    m := presentation (f**M); -- Next we lift M to a module over C
    N := (cokernel lift(m,C))/J; -- This C-module N is M restricted to C along C-->B
    if A === ZZ then
        trim ((f^-1)(promote(relativeFlatLocusZZ(N),D)))
    else
        trim ((f^-1)(promote(relativeFlatLocusINT(N,ideal(0_A)),D)))
-- The internal flatLocusINT routine works with modules over polynomial rings. We promote back to B
)

relativeFlatLocus(Ring) := B -> relativeFlatLocus(B^1)

relativeFlatLocus(Ideal) := J -> lift(relativeFlatLocus(ring(J)^1/J), ring(J))

relativeFlatLocus(RingMap, Module) := (f,M) -> (
-- We expect here that the target of f is the ring of M. (We should check that and print and error, if it is not the case.)
-- The check on these hypotheses can occur in "presentGraph", which does the grunt work of setting up the proper rings.
-- The function returns the flat locus for M considered as a module over the source A of f.
-- See tensor product code for how to make missing checks on R and f.
    if target f =!= ring M then error "expected target of ring map to be ring of module";
    (D,i,j) := presentGraph(f);
    trim j(relativeFlatLocus(i**M))
)



relativeFlatLocus(RingMap) := f -> relativeFlatLocus(f,(target f)^1)
-- Returns flat locus for target of ring homomorphism f over its source.
-- We assume that source and target have the same coefficient ring and that f is a homomorphism over this ring.

-----
-- Semi-continuity of fiber dimension
-----


relativeDimensionINT = method()

relativeDimension = method()




relativeDimension(ZZ, RingMap) := (n,phi) -> (
-- work in progress: This function implements "semicontinuity of fiber dimension."
-- it yields the open locus (in Spec(B)) where the map f:Spec(B)-->Spec(A) has relative dimension < n.
-- For example, if we take n = 1, we get the quasi-finite locus of the morphism.
-- test example: A = QQ[x,y], B = QQ[u,v], phi = map(B,A,{u,u*v}). The when n = 0, we get ideal(0_B);
-- when n = 1, we get ideal(u); when n > 1, we get ideal(1_B).
    A := source phi;
    B := target phi;
    relativeDimensionINT(n, phi, ideal(0_A), ideal(0_B))
)

relativeDimensionINT(ZZ,RingMap,Ideal,Ideal) := (n,phi,I,J) -> (
--The internal function to be called recursively.
--The function returns the locus in Spec(B/J) over which the relative dimension of Spec(B/J)-->Spec(A/I) is less than n.
--I think if one looks carefully, then one can see Noetherian induction just on B.
    A := source phi; -- I should be an ideal of A
    B := target phi; -- J should be an ideal of B
    if J==B then J
    else(
    if not isPrime(J) then
        trim (intersect apply(minimalPrimes(J),P->relativeDimensionINT(n, phi, I, P)))
    else (
        if preimage(phi, J) != I then -- this check does not work if A does not have a coefficientRing, e.g. if A is QQ
            trim relativeDimensionINT(n, phi, preimage(phi, J), J)
        else (
            (D,i,j) := presentGraph(phi);
            r := map(frac (A/I),A);
            (F,k) := extendScalars(D/(i(J)), r); -- fiber of Spec(B/J) over generic point of Spec(A/I)
            e := dim F; -- dimension of generic fiber of Spec(B/J)->Spec(A/I); is this just the difference of the dimensions?
            if n <= e then J
            else (
                (R,s) := extendScalars(D,map(A/I,A));
                M := s**(i**(B^1/J));
                ls := apply(genericFreenessFP(M), r->lift(r,A));
                f := product ls;
                trim relativeDimensionINT(n, phi, I + ideal(f), J + ideal(phi(f)))
                )
            )
    )
)
)

----------------------------
-- Cotangent complex methods
-- The methods below only work when the ambient ring of the input ring is a polynomial ring.
-- The functions that call them always call (I hope) with input of this form.
----------------------------

Omega1 = method()

-- This module is the truncation of the full cotangent complex at degree 0.
-- The version with "Ring" input only works when the ambient ring of the input
-- is a polynomial ring.
Omega1(Ring) := B -> trim (B**(coker jacobian ideal B))

Omega1(RingMap) := f -> (
    (D,i,j) := presentGraph(f);
    trim (j**Omega1(D))
)

naiveCotComplex = method()

naiveCotComplex(Ring) := B -> (
-- naive cotangent complex for B/A, where A = coefficientRing B.
-- If B = P/J, where P is the ambient polynomial ring over A, then the complex is
-- J/J^2 --> Omega^1_{P/A}\otimes_P B,
-- where the differential comes from the derivative d.
-- This complex is the truncation of the full cotangent complex at degree -1.
-- This function assumes that the ambient ring P is a polynomial ring.
    J := ideal B;
    P := ambient B;
    n := numgens P;
    d := map(P^n, module J, jacobian J);
    C := chainComplex(d)**B
)

{*naiveCotComplex(RingMap) := f -> (
-- This does not work. Is it possible to extend scalars for a chain complex along a ring homomorphism?
    (D,i,j) := presentGraph(f);
    j(naiveCotComplex(D)) -- I don't know what to do here: j(chain complex) doesn't work and neither does any tensor product that I tried
    -- perhaps one could extract the components of the complex and extends scalars manually.
)
*}

LSCotComplex = method()

LSCotComplex(Ring) := B -> (
-- Lichtenbaum-Schlessinger cotangent complex, which is the truncation of the full cotangent complex
-- at degree -2. One writes it as
-- Rel/TrivRel --> F\otimes B --> Omega^1_{P/A}\otimes B
-- doesn't quite work yet ** but maybe it does now
-- in any case, the function assumes that the ambient ring P is a polynomial ring.
    J := ideal B;
    P := ambient B;
    n := numgens J;
    d1 := jacobian J;
    m := presentation module J;
    Rel := image m;
    TrivRel := image koszul(2, gens module J);
    d2 := map(P^n, Rel/TrivRel, m);
    C := chainComplex(P);
    C.dd_1 = d1; C.dd_2 = d2; -- I tried to do this with the "chainComplex" command, but it did not work
    C**B
)







----
-- Smooth, etale, unramified, relative complete intersection: Loci computed using cotangent complex; all relative properties
-- other things to add: syntomic (flat plus lci) locus.
---

unramifiedLocus = method()

unramifiedLocus(Ring) := B -> trim fittingIdeal(0, Omega1(B))
-- one could also use annihilator; unramified means vanishing of relative Omega^1

unramifiedLocus(RingMap) := f -> trim fittingIdeal(0, Omega1(f))
-- one could also use annihilator; unramified means vanishing of relative Omega^1


smoothLocus = method()

smoothLocus(Ideal) := J -> trim lift(smoothLocus ((ring J)/J), ring J)


smoothLocus(Ring) := B -> (
-- smooth locus for B over its coefficientRing
-- doesn't quite work yet ** in fact, I think it it working!
    C := naiveCotComplex(B);
    trim intersect(ann HH_1(C), flatLocus HH_0(C))
)

smoothLocus(RingMap) := f -> (
-- same hypotheses on f as in flatLocus(RingMap). Returns smooth locus of f.
    (D,i,j) := presentGraph(f);
    trim j(smoothLocus(D))
)



etaleLocus = method()

etaleLocus(Ring) := B -> trim (smoothLocus(B)*unramifiedLocus(B)) -- etale means smooth and unramified

etaleLocus(RingMap) := f -> trim (smoothLocus(f)*unramifiedLocus(f))

lciLocus = method()

lciLocus(Ring) := B-> (
-- relative complex intersection locus for B over its coefficientRing: we need T_2(B,A,B) = 0 (HH_2 of LS-cot. complex) and I/I^2 locally free,
-- where B = P/I with P a polynomial ring over A.
    I := module ideal B;
    C := LSCotComplex(B);
    L1 := HH_2(C);
    L2 := Ext^1(B**I, trim image presentation (B**I)); -- obstruction to I/I^2 = B**I locally free
    trim intersect(ann L1, ann L2)
)

lciLocus(RingMap) := f -> (
-- same hypotheses on f as in flatLocus(RingMap). Returns the relative CI locus of f.
    (D,i,j) := presentGraph(f);
    trim j(lciLocus(D))
)

----
-- Methods for absolute properties: CM, Gorenstein, regular, etc
-- possibly add: reduced locus, normal locus, S_n locus, R_n locus
-- expose and document R_n, S_n, reduced checks
-- extend methods using primary decomposition to affine rings by finding
-- models over QQ or ZZ/p and extending scalars as in regularLocus function
---

codimPiece = method()

--Given a closed subset cut out by an ideal, returns the open locus where the closed subset has codimension > n
--only works for affine rings because of dimensions and primary decomposition.
--This function realizes the lower semicontinuity of codim_x(Y,X) in the variable x for a closed subset Y of X.
--Note that for x in X-Y, the codimension is +infty. In general, the larger the requested codimension, the smaller
--the open set. (Only high codimension irreducible components of Y are included.)

codimPiece(ZZ,Ideal) := (n, I) -> (
    ps := select(minimalPrimes(I), P -> ( trueCodimension(P) <= n ));
    if ps == {} then
        ideal(1_(ring I))
    else
        trim intersect ps
)

dimPiece = method()

-- realizes upper semicontinuity of dim_x(X) in x: returns open locus of x where dim_x(X)<n.

dimPiece(ZZ, Ring) := (n, A) -> (
    ps := select(minimalPrimes(ideal(0_A)), P -> (dim(A/P) >= n));
    if ps == {} then
        ideal(1_A)
    else
        trim intersect ps
)

trueCodimension2 = method()

trueCodimension2(Ideal) := I -> (
-- This function gives the codimension of Z(I) in Spec(A), where A is the ring of I
-- the built-in function may give the wrong result, if Spec(A) is not biequidimensional
-- Since this function uses primary decomposition, it only works when
-- the base ring of A is QQ or ZZ/p
    A := ring I;
    if not isAffineRing(A) then
        error "expected affine ring"
    else (
        min(apply(minimalPrimes(ideal(0_A)), P->dim(A/P)-dim(A/(I+P))))
        )
)

trueCodimension = method()

trueCodimension(Ideal) := I -> (
-- This function does the same thing as above but without using the full
-- computation of minimal primes containing I. The same function would work
-- for finitely presented ZZ-algebras, if one could get "dim" to work for
-- such rings.
-- The key here is that equidimensional affine rings are biequidimensional, which lets us
-- compute codimension as below.
    A := ring I;
    if not isAffineRing(A) then
        error "expected affine ring"
    else (
        min(apply(equidimensional(ideal(0_A)), P->dim(A/P)-dim(A/(I+P))))
        )
)

equidimensional = method()

equidimensional(Ideal) := I -> (
-- Returns a list of equidimensional ideals whose product has the same radical as I,
-- i.e. that cuts out the same subset of Spec(A) as I.
-- The function is handy, since it does not use primary decomposition.
-- In place of primary decomposition, it uses homological methods from
-- Eisenbud-Huneke-Vasconcelos. It will only work for an ideal in a ring
-- finitely presented over a regular ring.
-- despite my hopes, this function does not yet work over ZZ
-- ** In fact, the E-H-V approach is already in "topComponents"
    A := ring I;
    (B,i) := flattenRing A;
    P := ambient B;
-- We will assume that the coefficient ring here is ZZ or a field -- but ZZ does not yet work
    J := lift (i^-1 I, P);
    apply(equidimensionalINT J, K -> i(promote(K,B)))
    )

equidimensionalINT = method()

equidimensionalINT(Ideal) := I -> (
-- This function assumes that the given ideal is in an ideal in a regular domain.
    P := ring I;
    if I == P then {}
    else
        (
        e := codim I; -- this does not work over ZZ!
        J := ann Ext^e(P^1/I,P);
-- alternative would be to use
--      J := topComponents(I);
        prepend(J,equidimensionalINT(saturate(I,J)))
        )
    )

regularLocus = method()

regularLocus(Ring) := A -> (
-- this function seems to work now for rings finitely presented over a field.
-- the tricky thing is to handle char = p, since to use the smooth locus to detect regularity
-- we need to work over a perfect field. To reduce to the case of a perfect base field,
-- we find a ring B that is finitely presented over ZZ/p such that A is a localization of B.
-- (The fRing function does the necessary calculation.)
-- a useful check to separate easier cases from harder cases is "isAffineRing", which checks whether
-- a ring is a quotient of a polynomial ring over a field.
    if isAffineRing(A) then (
        if char A == 0 then smoothLocus(A)
        else (
            (B,j) := fRing(A); j(smoothLocus(B))
            )
         )
    else if coefficientRing(A) === ZZ then regularLocusZZ(A)
    else error "only affine rings and rings finitely presented over ZZ implemented."
)

kunzRegularLocus = method()

-- An alternative way to calculate the regular locus for a ring of char. p

kunzRegularLocus(Ring) := A -> (
    (B,i) := fRing(A); i(relativeFlatLocus absFrob B)
-- we first extract a finitely presented ring over ZZ, QQ, or ZZ/p to work with;
-- this step is necessary, since the flatness function for a homomorphism requires that homomorphism
-- be linear over the coefficientRings of the source and target.
-- a theorem of Kunz shows that for a ring of char p, flatness of absolute Frobenius is equivalent to regularity
-- note that if B does not have characteristic p, then "absFrob" returns an error
)

primeFactors = method()

primeFactors(ZZ) := n -> (
-- returns a list of the prime factors of n
    apply(toList factor n, p -> p # 0)
)

regularLocusZZ = method()

regularLocusZZ(Ring) := A -> (
-- We assume without checking that the coefficientRing of A is ZZ for this internal function.
-- The goal is to define a function that returns the regular locus
-- of a finitely presented ZZ algebra
-- The function is just a toy, since M2 cannot handle ZZ/p when we have p > 32767 = 2^15 - 1
-- Such primes turn up incidentally in innocent-looking computations.
-- Consider, for example, regularLocusZZ (ZZ[x]/(x^7-5*x+5)), where the prime 590263 appears as a factor
-- of the discriminant and thus stops the computation.
-- Singular looks a bit better, since it works to 2^31 - 1.
    I := smoothLocus(A); -- calculate the smooth locus over ZZ: I cuts out the singular locus
    J := genericSaturate(I); -- calculate the closure of the singular locus of the fiber over QQ
    K := saturate(I,J); -- remove the components of J from I
    n := char (A/K); -- this characteristic cannot be 0, since we removed the part of the non-smooth locus over (0_ZZ)
    ps := primeFactors(n);
    if ps == {} then J -- the list "ps" is empty precisely when n = 1, i.e. when I and J have the same components.
    -- in this case, both cut out the singular locus; we can return either one.
    else (
        trim intersect(J, intersect apply(ps, p->regularLocusFiber(p, A, K)))
        -- the singular locus will be the union of J and the singular loci in the fibers that meet K
    )

)

regularLocusFiber = method()

regularLocusFiber(ZZ, Ring, Ideal) := (p,A,I) ->(
-- take a prime p, a finitely presented ZZ-algebra A, and an ideal I in A
-- Returns the locus in Spec(A) of primes containing I and p where A is regular
-- This function is quite awkward, since M2 has weak support for finitely presented ZZ-algebras
-- We compute some things in the fiber over p and then lift back to A
-- We use the local criterion for regularity: if A/I is regular then A is regular along I
-- precisely where A-->A/I is an lci morphism, i.e. where we have a regular immersion.
    Ap := A/ideal(p_A);
    Ip := promote(I,Ap);
    (F,phi,psi) := presentGraph(map(Ap,ZZ/p)); -- phi:Ap-->F and psi:F-->Ap isomorphism.
    if Ip == ideal(1_Ap) then ideal(1_A) 
    else
        (
        ps := minimalPrimes(phi(Ip));
        if phi(Ip) != ps_0 then trim intersect(apply(ps, P->regularLocusFiber(p,A,lift(psi(P),A))))
        -- if phi(Ip) (I + <p>) is not prime, may call the function recursively on the minimal primes.
        -- in the second part of the function, we may thus assume that phi(Ip) is prime.
        -- that recursion is probably not the cleverest thing to do, but it works.
        else
            (
            R := A/(I+ideal(p_A));
            J := lciLocus(map(R,A)); -- this tells us where Spec(A/(I+<p>))->Spec(A) is a regular immersion
            if J == ideal(0_R) then -- if it fails to be a regular immersion everywhere, along I, then
            -- A is not regular at the generic point of I+<p> and hence is not regular along all of I+<p>
            -- since the regular locus is open.
                (I + ideal(p_A))
            else
                regularLocusFiber(p,A,lift(psi(smoothLocus(phi(Ip))),A)*lift(J,A))
                -- A is regular at the places where A/(I + <p>) is regular and A->A/(I+<p>) is lci.
                -- we continue recursively, which is safe to do, since the regular locus of A/(I+<p>)
                -- (same as smooth locus over ZZ/p) is non-empty open, as A/(I+<p>) is a finitely
                -- presented domain over ZZ/p.
            )
        )

)

genericSaturate = method()

genericSaturate(Ideal) := I -> (
-- Hidden method: assume I is an ideal in a ring A finitely presented over ZZ
-- find a new ideal cutting out the closure of the generic fiber of V(I) in Spec(A)
    A := ring I; -- assume here that coefficientRing is ZZ
    d := product apply(flatten entries gens gb I, leadCoefficient);
--    P := ambient A;
--    J := lift(I,P);
--    d := product apply(flatten entries gens gb J, leadCoefficient);
-- It seems that Groebner bases over quotient rings
-- are set up so that we do not need to lift to polynomial rings
-- I am a bit nervous, but it seems to work.
    saturate(I, ideal(d_A))
)

CMLocus = method()

CMLocus(Ring) := A -> (
-- examples of non-CM rings: QQ[x,y]/(x^2,x*y) and QQ[x,y,z]/(x*y,x*z)
-- if A is an affine ring or a finitely presented ZZ algebra, the function returns its CM locus
    (R, f) := flattenRing A;
    k := coefficientRing R;
    if (isField k) or (k===ZZ) then (
        C := ambient R;
        I := ideal R;
        g := map(R,C);
        n := dim C;
        J := sum(0..n, i->product(select(0..n,j->j!=i),j->ann(Ext^j(C^1/I,C^1))));
        trim f^-1(g(J))
        )
    else error "expected affine ring or finitely presented ZZ-algebra"
)

GorensteinLocus = method()

GorensteinLocus(Ring) := A -> (
-- nice test case: take coimage of QQ[u,v,w]->Q[t] by u->t^3, v->t^4, w->t^5. CM but not Gorenstein.
-- Looks at Bruns-Herzog for more examples. Find something nice to show Gorenstein but not CI
    (R, f) := flattenRing A;
    k := coefficientRing R;
    if (isField k) or (k===ZZ) then (
        C := ambient R;
        I := ideal R;
        g := map(R,C);
        n := dim C;
        L := f^-1(g(product(0..n,j->fittingIdeal(1,Ext^j(C^1/I,C^1)))));
        trim intersect(L,CMLocus(A))
        )
   else error "expected affine ring or finitely presented ZZ-algebra"
)

CILocus = method()
-- calculates the absolute CI locus for affine rings and finitely presented ZZ algebras
-- using the fact that the CI locus is the same as the lci locus over the base field
-- or over ZZ.

CILocus(Ring) := A -> (
    (R, f) := flattenRing A;
    k := coefficientRing R;
    if (isField k) or (k===ZZ) then
        trim f^-1(lciLocus(map(R,k)))
    else
        error "expected affine ring or finitely presented ZZ-algebra"
)

isCM = method()

isCM(Ring) := A -> (
-- tests whether the ring A is CM.
    CMLocus(A) == A
)


isS = method()

isS(ZZ, Ring) := (n,A) -> (
-- tests whether the ring A satisfies the S_n condition.
    (trueCodimension(CMLocus(A))) > n -- since we use "true codim", this function works w/o equidimensional hypothesis
)

isR = method()

isR(ZZ, Ring) := (n,A) -> (
--tests whether the ring satisfies the R_n condition
--currently everything works for affine rings
    (trueCodimension(regularLocus(A))) > n -- since we use "true codim", this function works w/o equidimensional hypothesis
)

isNormalFixed = method()

isNormalFixed(Ring) := A -> isS(2,A) and isR(1,A)

normalLocus = method()

normalLocus(Ring) := A -> codimPiece(1,regularLocus A)*codimPiece(2,CMLocus A)

isReduced = method()

isReduced(Ring) := A -> isS(1,A) and isR(0,A)

reducedLocus = method()

reducedLocus(Ring) := A -> intersect minimalPrimes ideal(0_A)


beginDocumentation()

-- Front Page
doc ///
    Key
        Loci
    Headline
        A package for computing open loci corresponding to various properties of rings and ring maps
    Description
        Text
            Let $A$ be a commutative ring that is finitely presented over $\mathbf{Z}$ or over a field.
            For many properties $\Pi$ of local rings, the set $U$ of primes $P$ in ${\rm Spec}(A)$ such that $A_P$
            has the property $\Pi$ is open in ${\rm Spec}(A)$. This holds, for example, if $\Pi$ is the Gorenstein property,
            the Cohen-Macaulay property, or the regularity property.
            
            Let $\phi:A\to B$ be a ring homomorphism, and let $f:{\rm Spec}(B)\to {\rm Spec}(B)$ be the corresponding
            morphism of spectra. For some properties of scheme morphisms, the set $U$ of points $P$ in ${\rm Spec}(B)$ such that $f$ has the
            given property at $P$ is open in ${\rm Spec}(B)$. This openness holds, for example, for the properties of being flat,
            unramified, {\'e}tale, and smooth.
            
            This package provides functions that compute such open loci where various properties hold.
            Each function constructing such a locus returns an ideal. The elements of the desired locus
            are those primes that do not contain the ideal. (Two ideals thus describe the same locus when
            they are contained in the same prime ideals, i.e. when they have the same radical.)
///


doc ///
    Key
        isOpenSubLocus
        (isOpenSubLocus, Ideal, Ideal)
    Headline
        calculates order relation on open subsets of spectra
    Usage
        i=isOpenSubLocus(I,J)
    Inputs
        I:Ideal
        J:Ideal
            assumed to be in the same ring as I
    Outputs
        i:Boolean
            indicating whether the open locus defined by I is contained in the open locus defined by J
    Description
        Text
            Let $I$ and $J$ be finitely generated ideals in a ring $A$. These ideals correspond
            to open subspaces $U$ and $V$ of ${\rm Spec}(A)$. We have $U\subset V$ precisely when
            for each generator $f$ of $I$, the monoid generated by $f$ meets $J$. In other words
            we have $U\subset V$ precisely when some power of $f$ is in $J$ or, equivalently,
            when $f$ is in the radical of $J$. We check this condition by checking whether
            $1-ft$ and $J$ generate the unit ideal in $A[t]$.
    SeeAlso
        equalOpenLoci
///

doc ///
    Key
        equalOpenLoci
        (equalOpenLoci, Ideal, Ideal)
    Headline
        determines if the open loci defined by ideals are equal
    Usage
        I = equalOpenLoci(I,J)
    Inputs
        I:Ideal
        J:Ideal
            assumed to be in the same ring as I
    Outputs
        i:Boolean
            indicating whether the open loci defined by I and J are equal
    Description
        Text
            Let $I$ and $J$ be finitely generated ideals in a ring $A$. These ideals correspond
            to open subspaces $U$ and $V$ of ${\rm Spec}(A)$. Using the function isOpenSubLocus,
            we determine if these subspaces are equal.
    SeeAlso
        isOpenSubLocus
///

-- flatLocus

doc ///
    Key
        flatLocus
        (flatLocus, Module)
    Headline
        calculates the flat locus of a finitely presented module over a ring
    Usage
        I = flatLocus(M)
    Inputs
        M:Module
    Outputs
        I:Ideal
            defining the (open) flat locus of $M$
    Description
        Text
            Let $A$ be the ring of $M$. For a prime ideal $P$ of $A$, the localization
            $M_P$ is flat over $A$ if and only if $P$ does not contain $I$, i.e. the flat
            locus of $M$ in ${\rm Spec}(A)$ is complement of the closed set defined by $I$.
        Example
            A = QQ[x,y];
            M = module ideal(x,y);
            flatLocus M
        Text
           The method works as follows: let $0\to K\to F\to M\to 0$ be an exact sequence with
           $F$ projective. Then $M$ is locally free precisely when ${\rm Ext}^1(M,K)$ vanishes.
           Since the formation of this ${\rm Ext}$ module is compatible with localization, the
           open locally free locus is the complement of the support of the ${\rm Ext}$ module.
           The ideal returned is an ideal cutting out this support.
    SeeAlso
        relativeFlatLocus
///

-- relativeFlatLocus

doc ///
    Key
        relativeFlatLocus
        (relativeFlatLocus, Module)
        (relativeFlatLocus, Ring)
        (relativeFlatLocus, Ideal)
        (relativeFlatLocus, RingMap)
        (relativeFlatLocus, RingMap, Module)
    Headline
        calculates the relative flat locus of a finitely presented module over an algebra
    Usage
        I = relativeFlatLocus M
        I = relativeFlatLocus B
        I = relativeFlatLocus J
        I = relativeFlatLocus(phi, M)
    Inputs
        M:Module
        B:Ring
        J:Ideal
        phi:RingMap
    Outputs
        I:Ideal
            defining the (open) relative flat locus of M
            over the coefficient ring of B or over the source of phi
    Description
        Text
            Let $M$ be a finitely presented module over a commutative ring $B$. Suppose
            that $B$ is a finitely presented $A$-algebra, either given as a quotient of
            a polynomial ring over $A$ or provided with its $A$-module structure via a
            ring homomorphism $\phi:A\to B$. Under these hypotheses, there is an open set
            $U$ in ${\rm Spec}(B)$ such that for all primes $P$ of $B$, the localization
            $M_P$ is flat over $A$ if and only if $P$ is in $U$.
            
            This method returns an ideal $I$ of $B$ with the property that for all primes
            $P$ of $B$, the localization $M_P$ is flat over $A$ if and only if $P$ does not
            contain $I$, i.e. the relative flat locus for $M$ over $A$ is the complement in
            ${\rm Spec}(B)$ of the closed set defined by $I$.
            
            When there is no module specified, the method returns the relative flat locus
            for $B$ viewed as an $A$-module over $A$.
        Example
            A = QQ[x,y]/(x^3 - y^2);
            B = QQ[t];
            phi = map(B, A, {x_A=>t^2, y_A=>t^3});
            relativeFlatLocus phi
        Text
            The method uses a variant of the construction implicit in Grothendieck's proof of
            the openness of the flat locus for a finitely presented algebra from EGA IV part 3,
            Theorem 11.1.1.
    Caveat
        If the ring $B$ has no coefficient Ring (e.g. if $B$ is $\mathbf{Q}$ or $\mathbf{Z}/p$),
        then the function will fail. The proper way to calculate the absolute flat locus is to
        use the command flatLocus
    SeeAlso
        flatLocus
///

--relativeDimension

doc ///
    Key
        relativeDimension
        (relativeDimension, ZZ, RingMap)
    Headline
        calculates the locus where the relative dimension is less than an integer n
    Usage
        I = relativeDimension(n, phi)
    Inputs
        n:ZZ
        phi:RingMap
    Outputs
        I:Ideal
            belonging to the target of phi and defining the (open) locus where the relative dimension
            of the map on spectra corresponding to phi has relative dimension less than n
    Description
        Text
            Let $\phi:A\to B$ be a homomorphism of rings that makes $B$ a finitely
            presented $A$-algebra. Let $f:{\rm Spec(B)\to {\rm Spec(A)}$ be the corresponding map on spectra,
            and let $n$ be an integer. The locus of points $x$ in ${\rm Spec}(B)$ such that
            ${\rm dim}_x(f^{-1}(f(x)))<n$ is open in ${\rm Spec(B)}$.
            (One calls this fact ``semicontinuity of fiber dimension.'') Note, for example, that
            if we take $n=1$, this locus is quasi-finite locus of $f$.
            
            This method returns an ideal $I$ in $B$ such that for all $P\in{\rm Spec}(B)$, the fiber
            dimension at $P$ is less than $n$ if and only if $P$ does not contain $I$, i.e. the
            open locus where the local fiber dimension is less than $n$ is the complement of the
            closed subset of ${\rm Spec}(B)$ defined by $I$.
        Example
            A = QQ[u,v];
            B = QQ[x,y];
            phi = map(B, A, {u_A=>x, v_A=>x*y});
            relativeDimension(1, phi)
        Text
            The method uses the construction implicit in Grothendieck's proof of semicontinuity of
            fiber dimension ("Chevalley's Theorem") in EGA IV part 3, Theorem 13.1.3.
    Caveat
        The rings $A$ and $B$ must be a rings over QQ or over ZZ/p because the algorithm
        uses minimal primes, which have only been implemented in Macaulay 2 for such
        rings, and dimensions, which have only been implemented in Macaulay 2 over
        affine rings.
///

--unramifiedLocus
doc ///
    Key
        unramifiedLocus
        (unramifiedLocus, Ring)
        (unramifiedLocus, RingMap)
    Headline
        calculates the locus where an algebra is unramified
    Usage
        I = unramifiedLocus B
        I = unramifiedLocus phi
    Inputs
        B:Ring
            finitely presented over its coefficient ring A
        phi:RingMap
    Outputs
        I:Ideal
            in the A-algebra B (or in the target of phi) defining the (open) locus where
            B is unramified over A
            (or where the target of phi is unramified over the source of phi)
    Description
        Text
            Let $B$ be a finitely presented $A$-algebra, either given as a quotient of
            a polynomial ring over $A$ or provided with its $A$-module structure via a
            ring homomorphism $\phi:A\to B$. Let $f:{\rm Spec(B)}\to {\rm Spec(A)}$ be
            the corresponding map on spectra. The locus of points $P$ in ${\rm Spec}(B)$
            such that $f$ is unramified at $P$ is open in ${\rm Spec}(B)$.
            
            This method returns an ideal $I$ in $B$ such that for all $P\in{\rm Spec}(B)$, the
            morphism $f$ is unramified at $P$ if and only if $I$ does not contain $P$, i.e. the
            open locus where $f$ is unramified is the complement of the closed subset of
            ${\rm Spec}(B)$ defined by $I$.
        Example
            A = QQ[x,y];
            B = A/(x^2, x*y, y^2);
            unramifiedLocus(map(B,A))
        Text
            The unramified locus for $B$ over $A$ is precisely the complement of the support
            of $\Omega^1_{B/A}$. The method returns an ideal cutting out this support.
    SeeAlso
        etaleLocus
        smoothLocus
        lciLocus
///

--etaleLocus
doc ///
    Key
        etaleLocus
        (etaleLocus, Ring)
        (etaleLocus, RingMap)
    Headline
        calculates the locus where an algebra is etale
    Usage
        I = etaleLocus B
        I = etaleLocus phi
    Inputs
        B:Ring
            finitely presented over its coefficient ring A
        phi:RingMap
    Outputs
        I:Ideal
            in the A-algebra B (or in the target of phi) defining the (open) locus where
            B is etale over A (or where the target of phi is etale over the source of phi)
    Description
        Text
            Let $B$ be a finitely presented $A$-algebra, either given as a quotient of
            a polynomial ring over $A$ or provided with its $A$-module structure via a
            ring homomorphism $\phi:A\to B$. Let $f:{\rm Spec(B)}\to {\rm Spec(A)}$ be
            the corresponding map on spectra. The locus of points $P$ in ${\rm Spec}(B)$
            such that $f$ is etale at $P$ is open in ${\rm Spec}(B)$.
            
            This method returns an ideal $I$ in $B$ such that for all $P\in{\rm Spec}(B)$, the
            morphism $f$ is etale at $P$ if and only if $I$ does not contain $P$, i.e. the
            open locus where $f$ is etale is the complement of the closed subset of
            ${\rm Spec}(B)$ defined by $I$.
        Example
            B = ZZ[x]/(x^2 + 1);
            etaleLocus B
        Text
            The etale locus for $B$ over $A$ is the intersection of the unramified locus and
            the smooth locus.
    SeeAlso
        unramifiedLocus
        smoothLocus
        lciLocus
///

--smoothLocus
doc ///
    Key
        smoothLocus
        (smoothLocus, Ring)
        (smoothLocus, RingMap)
    Headline
        calculates the locus where an algebra is smooth
    Usage
        I = smoothLocus B
        I = smoothLocus phi
    Inputs
        B:Ring
            finitely presented over its coefficient ring A
        phi:RingMap
    Outputs
        I:Ideal
            in the A-algebra B (or in the target of phi) defining the (open) locus where
            B is smooth over A (or where the target of phi is smooth over the source of phi)
    Description
        Text
            Let $B$ be a finitely presented $A$-algebra, either given as a quotient of
            a polynomial ring over $A$ or provided with its $A$-module structure via a
            ring homomorphism $\phi:A\to B$. Let $f:{\rm Spec(B)}\to {\rm Spec(A)}$ be
            the corresponding map on spectra. The locus of points $P$ in ${\rm Spec}(B)$
            such that $f$ is smooth at $P$ is open in ${\rm Spec}(B)$.
            
            This method returns an ideal $I$ in $B$ such that for all $P\in{\rm Spec}(B)$, the
            morphism $f$ is smooth at $P$ if and only if $I$ does not contain $P$, i.e. the
            open locus where $f$ is smooth is the complement of the closed subset of
            ${\rm Spec}(B)$ defined by $I$.
        Example
            B = ZZ[x,y]/(x*y + 5);
            smoothLocus B
        Text
            One can calculate the smooth locus using a presentation of $B$ as $A$-algebra. It
            is the locus where the associated 2-term truncated cotangent complex has vanishing ${\rm H}_1$
            and locally free ${\rm H}_0$. (This condition is called ``the Jacobian criterion for smoothness.'')
            The method computes the locus where these homology modules have
            the desired properties.
    SeeAlso
        unramifiedLocus
        etaleLocus
        lciLocus
///

--lciLocus
doc ///
    Key
        lciLocus
        (lciLocus, Ring)
        (lciLocus, RingMap)
    Headline
        calculates the locus where an algebra is a local complete intersection
    Usage
        I = lciLocus B
        I = lciLocus phi
    Inputs
        B:Ring
            finitely presented over its coefficient ring A
        phi:RingMap
    Outputs
        I:Ideal
            in the A-algebra B (or in the target of phi) defining the (open) locus where
            B is a local complete intersection over A
            (or where the target of phi is a local complete
            intersection over the source of phi)
    Description
        Text
            Let $B$ be a finitely presented $A$-algebra, either given as a quotient of
            a polynomial ring over $A$ or provided with its $A$-module structure via a
            ring homomorphism $\phi:A\to B$. Let $f:{\rm Spec(B)}\to {\rm Spec(A)}$ be
            the corresponding map on spectra. The locus of points $P$ in ${\rm Spec}(B)$
            such that $f$ is a local complete intersection at $P$ is open in ${\rm Spec}(B)$.
            
            This method returns an ideal $I$ in $B$ such that for all $P\in{\rm Spec}(B)$, the
            morphism $f$ is a local complete intersection at $P$ if and only if $I$ does not
            contain $P$, i.e. the open locus where $f$ is a local complete intersection
            is the complement of the closed subset of ${\rm Spec}(B)$ defined by $I$.
        Example
            A = QQ[x,y];
            B = A[t,z]/(x*y-z*t,x*z-t^2,y*t-z^2);
            lciLocus B
        Text
            One can calculate the local complete intersection locus using a presentation of
            $B$ as an $A$-algebra. It is the locus where the associated 3-term truncated
            cotangent complex (Lichtenbaum-Schlessinger cotangent complex)
            has vanishing ${\rm H}_2$ and locally free ${\rm H}_1$. The method computes the
            locus where these homology modules have the desired properties.
    SeeAlso
        unramifiedLocus
        etaleLocus
        smoothLocus
        CILocus
///

--regularLocus
doc ///
    Key
        regularLocus
        (regularLocus, Ring)
    Headline
        calculates the regular locus of a ring
    Usage
        I = regularLocus A
    Inputs
        A:Ring
            assumed to be either an affine ring (finitely presented over a field) or
            a finitely presented Z-algebra
    Outputs
        I:Ideal
            defining the locus where A is regular
    Description
        Text
            Let $A$ be a commutative ring and consider the locus of $P\in{\rm Spec}(A)$ such
            that $A_P$ is a regular local ring. For many sorts of rings, this locus is open.
            The openness holds, in particular, if $A$ is an affine ring (finitely presented over
            a field) or if $A$ is finitely presented over $\mathbf{Z}$. In either of these two cases
            this function returns an ideal $I$ of $A$ such that for all $P\in{\rm Spec}(A)$, the
            localization $A_P$ is regular if and only if $P$ does not contain $I$. In other words
            the (open) regular locus is the complement in ${\rm Spec}(A)$ of the closed set
            defined by $I$.
            
            As a first example, we consider a finitely presented algebra over $\mathbf{Z}$ that is
            not smooth over $\mathbf{Z}$ but which is regular.
        Example
            A = ZZ[x,y]/(x*y + 5);
            regularLocus A
        Text
            A second example is an inseparable field extension $B/F$. The field $B$ is a regular
            ring, but it is not smooth as an algebra over $F$.
        Example
            F = frac (ZZ/5[a]);
            B = F[x]/(x^5 - a);
            regularLocus B
        Text
            The following example comes from the introduction to Zariski's paper
            "The concept of a simple point on an abstract algebraic variety". It is regular
            but not smooth over $K$.
        Example
            K = frac (ZZ/3[b]);
            C = K[x,y]/(x^3 + y^3 - b);
            regularLocus C
        Text
            In the two previous examples, the difference between regularity and smoothness
            comes from the fact that the coefficient field is ``too small.'' Both algebras are
            smooth over larger fields than their given coefficient fields. The next example
            (due to Chevalley) also comes from the introduction to
            "The concept of a simple point on an abstract algebraic variety". It is notable
            since the associated scheme is regular and geometrically reduced and irreducible
            over $K$ but not smooth over $K$ at $y=0$.
        Example
            K = frac (ZZ/3[b]);
            D = K[x,y]/(y^2 + x^3 - b);
            regularLocus D
        Text
            If $A$ is an affine algebra whose base field $k$ is perfect, the regular locus for $A$
            is the same as the smooth locus of $A$ over $k$, which one may compute using the
            homology of the cotangent complex.
            
            If the base field $k$ is not perfect, then the regular locus may be larger than the
            smooth locus. We assume that the base field $k$ is finitely generated over its
            prime field $\mathbf{Z}/p\mathbf{Z}$. To compute the regular locus, we find a model
            $A'$ for $A$ over $\mathbf{Z}/p\mathbf{Z}$, whose regular locus we compute using the
            homology of the cotangent complex. Extending scalars along $A'\to A$, we find the regular
            locus of $A$. One generally refers to the method used in these affine-ring cases as
            the Jacobian criterion for regularity.
            
            If $A$ is a finitely presented $\mathbf{Z}$-algebra, then we follow Nagata's proof
            of the openness of the regular locus to construct the regular locus. The construction
            works recursively (corresponding to the Noetherian induction in the proof), using two key tools:
            
            * the Jacobian criterion for regularity and
            
            * the fact that if $A\to B$ is a surjection of rings that makes $B$ a local complete intersection
              over $A$ (regular immmersion) and if $B$ is regular, then $A$ is regular.
    Caveat
        Fairly large primes may occur when computing the regular locus of simple-looking finitely
        presented $\mathbf{Z}$-algebras. Since Macaulay 2 currently only handles primes up to $2^{15}-1$,
        one may encounter mysterious errors while using this function. For example, computing the regular
        locus for ZZ[x]/(x^7 – 5*x +5) generates an error, since it requires working modulo the prime 590263.
    SeeAlso
        smoothLocus
        CILocus
        CMLocus
        GorensteinLocus
///

doc ///
    Key
        kunzRegularLocus
        (kunzRegularLocus, Ring)
    Headline
        calculates the regular locus for a ring of prime characteristic p>0
    Usage
        I = kunzRegularLocus A
    Inputs
        A:Ring
            assumed to have prime characteristic p > 0
    Outputs
        I:Ideal
            defining the locus where A is regular
    Description
        Text
            Let $A$ be a commutative ring and consider the locus of $P\in{\rm Spec}(A)$ such
            that $A_P$ is a regular local ring. For many sorts of rings, this locus is open. This
            function computes the locus when the ring $A$ has characteristic $p$, assuming that
            Macaulay 2 can recognize the ring as a localization of a finitely presented
            $\mathbf{Z}$-algebra of characteristic $p$ or of a finintely presented
            $\mathbf{Z}/p\mathbf{Z}$-algebra.
            
            By a theorem of Kunz, the regular locus is the flat locus for the absolute Frobenius
            ($p$-power) endomorphism of $A$. This function returns an ideal $I$ of $A$ such that
            for all $P\in{\rm Spec}(A)$, the localization $A_P$ is regular if and only if $P$
            does not contain $I$. In other words the (open) regular locus is the complement in
            ${\rm Spec}(A)$ of the closed set defined by $I$.
        Text
            Our first example $B$ is an inseparable field extension of a field $F$ of characteristic $5$.
            The field $B$ is a regular ring, but it is not smooth as an algebra over $F$.
        Example
            F = frac (ZZ/5[a]);
            B = F[x]/(x^5 - a);
            kunzRegularLocus B
        Text
            The second example comes from the introduction to Zariski's paper
            "The concept of a simple point on an abstract algebraic variety". It is regular
            but not smooth over $K$.
        Example
            K = frac (ZZ/3[b]);
            C = K[x,y]/(x^3 + y^3 - b);
            kunzRegularLocus C
        Text
            In the two previous examples, the difference between regularity and smoothness
            comes from the fact that the coefficient field is ``too small.'' Both algebras are
            smooth over larger fields than their given coefficient fields. The next example
            (due to Chevalley) also comes from the introduction to
            "The concept of a simple point on an abstract algebraic variety". It is notable
            since the associated scheme is regular and geometrically reduced and irreducible
            over $K$ but not smooth over $K$ at $y=0$.
        Example
            K = frac (ZZ/3[b]);
            D = K[x,y]/(y^2 + x^3 - b);
            kunzRegularLocus D
        Text
            The function regularLocus calculates the same locus (and works for rings not
            of characteristic $p$), but it uses an entirely different method based on the
            connection between smoothness over $\mathbf{Z}/p\mathbf{Z}$ and regularity.
            The two functions may thus yield different descriptions of the same locus, as
            in the next example.
        Example
            A = (ZZ/3)[x,y]/(x^3 - y^2);
            regularLocus A
            kunzRegularLocus A
    SeeAlso
        relativeFlatLocus
        regularLocus
///


doc ///
    Key
        CMLocus
        (CMLocus, Ring)
    Headline
        calculates the Cohen-Macualay locus of a ring
    Usage
        I = CMLocus A
    Inputs
        A:Ring
            assumed to be either an affine ring (finitely presented over a field) or
            a finitely presented Z-algebra
    Outputs
        I:Ideal
            defining the locus where A is Cohen-Macaulay
    Description
        Text
            Let $A$ be a commutative ring and consider the locus of $P\in{\rm Spec}(A)$ such
            that $A_P$ is a Cohen-Macaulay local ring. For many sorts of rings, this locus is open.
            The openness holds, in particular, if $A$ is an affine ring (finitely presented over
            a field) or if $A$ is finitely presented over $\mathbf{Z}$. In either of these two cases
            this function returns an ideal $I$ of $A$ such that for all $P\in{\rm Spec}(A)$, the
            localization $A_P$ is Cohen-Macaulay if and only if $P$ does not contain $I$. In other words
            the (open) Cohen-Macaulay locus is the complement in ${\rm Spec}(A)$ of the closed set
            defined by $I$.
            
            The following example corresponds geometrically to a plane pierced by a line. This structure
            fails to be Cohen-Macaulay at the intersection point of the plane and the line.
        Example
            A = QQ[x,y,z]/(x*y,x*z);
            CMLocus A
        Text
            The ring in the following example (a semigroup ring) is Cohen-Macaulay but not Gorenstein. For much
            more on properties of semigroup rings, see Chapter 6 of "Cohen-Macaulay Rings" by Bruns and Herzog.
        Example
            A = coimage(map(QQ[t],QQ[u,v,w],{t^3,t^4,t^5}));
            CMLocus A
        Text
            Note that one can use the output of this function to check Serre's property $S_n$: a ring has
            the $S_n$ property precisely when the codimension of its non-Cohen-Macaulay locus is greater than
            $n$. For an affine ring $A$, one can compute this codimension in Macaulay 2. (The necessary dimension
            function for finitely presented $\mathbf{Z}$ algebras is not available at the time of this writing.)
///

doc ///
    Key
        GorensteinLocus
        (GorensteinLocus, Ring)
    Headline
        calculates the Gorenstein locus of a ring
    Usage
        I = GorensteinLocus A
    Inputs
        A:Ring
            assumed to be either an affine ring (finitely presented over a field) or
            a finitely presented Z-algebra
    Outputs
        I:Ideal
            defining the locus where A is Gorenstein
    Description
        Text
            Let $A$ be a commutative ring and consider the locus of $P\in{\rm Spec}(A)$ such
            that $A_P$ is a Gorenstein local ring. For many sorts of rings, this locus is open.
            The openness holds, in particular, if $A$ is an affine ring (finitely presented over
            a field) or if $A$ is finitely presented over $\mathbf{Z}$. In either of these two cases
            this function returns an ideal $I$ of $A$ such that for all $P\in{\rm Spec}(A)$, the
            localization $A_P$ is Gorenstein if and only if $P$ does not contain $I$. In other words
            the (open) Gorenstein locus is the complement in ${\rm Spec}(A)$ of the closed set
            defined by $I$.

            The ring in the following example (a semigroup ring) is Cohen-Macaulay but not Gorenstein. For much
            more on properties of semigroup rings, see Chapter 6 of "Cohen-Macaulay Rings" by Bruns and Herzog.
        Example
            A = coimage(map(QQ[t],QQ[u,v,w],{t^3,t^4,t^5}));
            GorensteinLocus A
    SeeAlso
        CILocus
        CMLocus
        regularLocus
///

doc ///
    Key
        CILocus
        (CILocus, Ring)
    Headline
        calculates the complete intersection locus of a ring
    Usage
        I = CILocus A
    Inputs
        A:Ring
            assumed to be either an affine ring (finitely presented over a field) or
            a finitely presented Z-algebra
    Outputs
        I:Ideal
            defining the locus where A is a complete intersection
    Description
        Text
            Let $A$ be a commutative ring and consider the locus of $P\in{\rm Spec}(A)$ such
            that $A_P$ is a complete intersection local ring. For many sorts of rings, this locus is open.
            The openness holds, in particular, if $A$ is an affine ring (finitely presented over
            a field) or if $A$ is finitely presented over $\mathbf{Z}$. In either of these two cases
            this function returns an ideal $I$ of $A$ such that for all $P\in{\rm Spec}(A)$, the
            localization $A_P$ is a complete intersection if and only if $P$ does not contain $I$. In other words
            the (open) complete intersection locus is the complement in ${\rm Spec}(A)$ of the closed set
            defined by $I$.

            To compute the complete intersection locus, we use the fact that for the types of rings $A$ in question
            (finitely presented over a field or over $\mathbf{Z}$), the complete interesection locus for $A$ is
            the same as the local complete intersection locus for $A$ viewed as an algebra over its base ring.
    SeeAlso
        lciLocus
        CMLocus
        GorensteinLocus
        regularLocus
///

doc ///
    Key
        trueCodimension
        (trueCodimension, Ideal)
    Headline
        calculates the codimension (height) of an ideal in an affine ring
    Usage
        n = trueCodimension I
    Inputs
        I:Ideal
            assumed to belong to an affine ring (finitely presented over a field)
    Outputs
        n:ZZ
            the codimension (height) of the ideal I
    Description
        Text
            The Macaulay 2 function "codim" applied to an ideal I in an affine ring A
            returns dim(A)-dim(A/I). This difference of dimensions is only guaranteed
            to be the codimension of I, if the ring A is biequidimensional, a property that
            follows for affine rings from equidimensionality. This function returns the
            actual codimension in all cases, even if A is not equidimensional.
        Example
            A = QQ[x,y];
            I = ideal(x,y)*ideal(x-1);
            B = A/I;
            codim ideal(x,y)
            trueCodimension ideal(x,y)
        Text
            Calculating the actual codimension is useful, if one wishes to check, for example,
            normality of an affine ring A using Serre's criterion. This criterion characterizes
            normality using the codimensions of the non-regular (codimension > 1) and non-CM
            (codimension > 2) loci.
    Caveat
        Since Macaulay 2 only handles dimension computations universally
        for affine rings, this function only works for ideals in affine rings.
///

doc ///
    Key
        isNormalFixed
        (isNormalFixed, Ring)
    Headline
        decides whether an affine ring is normal
    Usage
        i = isNormalFixed A
    Inputs
        A:Ring
            assumed to be an affine ring (finitely presented over a field)
    Outputs
        i:Boolean
            indicating whether A is normal
    Description
        Text
            This function uses Serre's criterion to determine whether an affine ring
            A is normal. Since it uses the trueCodimension function and the regularLocus
            function implemented in this package, it works even if A is not equidimensional
            and even if the coefficient ring of A is not a perfect field.
            
            The following example comes from a Macaulay 2 bug report filed by Karl Schwede
            on 6/30/2014. The ring B below is regular (and hence normal) but not equidimensional.
        Example
            A = QQ[x,y,z];
            I = ideal(x,y)*ideal(x-1);
            B = A/I;
            isNormal B
            isNormalFixed B
        Text
            The second example is one of the examples from the regularLocus function. The ring C
            is a field and hence normal. Since F is an inseparable extension of its coefficient field,
            one cannot use the Jacobian criterion for smoothness to detect regularity directly. The
            Macaulay 2 function "isNormal" thus does not detect normality.
        Example
            F = frac(ZZ/5[a]);
            C = F[x]/(x^5 - a);
            isNormal C
            isNormalFixed C
    Caveat
        Since Macaulay 2 only handles dimension computations uniformly for affine rings, this function
        only works for affine rings.
    SeeAlso
        CMLocus
        regularLocus
        normalLocus
///

doc ///
    Key
        dimPiece
        (dimPiece, ZZ, Ring)
    Headline
        Calculates the open locus where the spectrum of an affine ring has dimension less than n
    Usage
        I = dimPiece(n, A)
    Inputs
        n:ZZ
        A:Ring
            over QQ or ZZ/p
    Outputs
        I:Ideal
            defining the open locus in Spec(A) where the local dimension is less than n
    Description
        Text
            Let A be a commutative ring, and let X be Spec(A). The function $x\mapsto {\rm dim}_x(X)$ is upper
            semicontinuous on $X$, meaning that for each $n$, the locus of $x$ such that ${\rm dim}_x(X)<n$ is
            open in $X$: it is the complement of the irreducible components of dimension greater than or equal
            to $n$. This function returns an ideal whose corresponding closed set in $X$ is this union of
            irreducible components. The open subset of $X$ defined by the ideal is thus the
            open set where $X$ has dimension less than $n$.
            
            The function relativeDimension provides a relative version of this construction,
            usually known as "semicontinuity of fiber dimension."
    Caveat
        Since this function uses primary decomposition, it only works for rings over QQ and rings over ZZ/p.
    SeeAlso
        relativeDimension
        trueCodimension
///

doc ///
    Key
        isS
        (isS, ZZ, Ring)
    Headline
        decides whether an affine ring satisfies Serre's S_n criterion
    Usage
        i = isS(n,A)
    Inputs
        n:ZZ
        A:Ring
            assumed to be an affine ring (finitely presented over a field)
    Outputs
        i:Boolean
            indicating whether A satisfies the S_n condition
    Description
        Text
            This function determines whether a ring $A$ satisfies Serre's $S_n$
            condition, i.e. it determines whether the non-CM locus for $A$ has
            codimension greater than $n$.
    Caveat
        Since Macaulay 2 only handles dimension computations uniformly for affine rings, this function
        only works for affine rings.
    SeeAlso
        isReduced
        isNormal
        isR
        CMLocus
///

doc ///
    Key
        isR
        (isR, ZZ, Ring)
    Headline
        decides whether an affine ring satisfies Serre's R_n condition
    Usage
        i = isR(n, A)
    Inputs
        n:ZZ
        A:Ring
            assumed to be an affine ring (finitely presented over a field)
    Outputs
        i:Boolean
            indicating whether A satisfies the R_n condition
    Description
        Text
            This function determines whether a ring $A$ satisfies Serre's $R_n$
            condition, i.e. it determines whether the non-regular (singular) locus
            for $A$ has codimension greater than $n$.
    Caveat
        Since Macaulay 2 only handles dimension computations uniformly for affine rings, this function
        only works for affine rings.
    SeeAlso
        isReduced
        isNormal
        isS
        regularLocus
        kunzRegularLocus
///

doc ///
    Key
        isReduced
        (isReduced, Ring)
    Headline
        decides whether an affine ring is reduced
    Usage
        i = isReduced A
    Inputs
        A:Ring
            assumed to be an affine ring (finitely presented over a field)
    Outputs
        i:Boolean
            indicating whether A is reduced
    Description
        Text
            This function uses the characterization "$R_0$ and $S_1$" of reduced local rings
            to determine whether an affine ring A is reduced. Since it uses the trueCodimension
            function and the regularLocus function implemented in this package, it works even if
            A is not equidimensional and even if the coefficient ring of A is not a perfect field.
            
            One could also check that A is reduced by computing the radical of A, but currently
            the radical function (and other functions related to primary decomposition) in Macaulay
            2 require the base field to be $\mathbf{Q}$ or $\mathbf{Z}/p\mathbf{Z}$ for a prime $p$.
    Caveat
        Since Macaulay 2 only handles dimension computations uniformly for affine rings, this function
        only works for affine rings.
    SeeAlso
        CMLocus
        regularLocus
        isNormalFixed
        reducedLocus
///

doc ///
    Key
        reducedLocus
        (reducedLocus, Ring)
    Headline
        calculates the reduced locus of a ring
    Usage
        I = reducedLocus(A)
    Inputs
        A:Ring
            assumed to be finitely presented over QQ or ZZ/p
    Outputs
        I:Ideal
            defining the reduced locus in Spec(A)
    Description
        Text
            Let $A$ be a commutative ring, and consider the locus of $P\in{\rm Spec}(A)$ such that
            $A_P$ is reduced. This locus is the open subset of ${\rm Spec}(A)$ of primes $P$ that
            do not contain the nilradical of $A$. The function returns the nilradical, computed
            as the intersection of the minimal primes of $A$.
    Caveat
        Since the function uses primary decomposition directly to compute the radical, it
        only works for rings finitely presented over $\mathbf{Q}$ or $\mathbf{Z}/p\mathbf{Z}$
        for a prime $p$.
    SeeAlso
        isReduced
///

doc ///
    Key
        normalLocus
        (normalLocus, Ring)
    Headline
        calculuates the normal locus of a ring
    Usage
        I = normalLocus(A)
    Inputs
        A:Ring
            assumed to be finitely presented over QQ or ZZ/p
    Outputs
        I:Ideal
            defining the normal locus in Spec(A)
    Description
        Text
            Let $A$ be a commutative ring, and consider the locus of $P\in{\rm Spec}(A)$ such
            that $A_P$ is a normal local ring. For many sorts of rings, this locus is open.
            If $A$ is finitely presented over $\mathbf{Q}$ or over $\mathbf{Z}/p\mathbf{Z}$ for
            a prime $p$, this openness holds (as it does in many more general cases). In this case
            the function returns an ideal $I$ such that for all $P\in{\rm Spec}(A)$, the
            localization $A_P$ is normal if and only if $P$ does not contain $I$. In other words
            the (open) normal locus is the complement in ${\rm Spec}(A)$ of the closed set
            defined by $I$.
            
            To calculate the normal locus, the function uses the CM locus (CMLocus) and the
            regular locus (regularLocus): by Serre's criterion, the normal locus is the intersection
            of the (open) subsets of ${\rm Spec}(A)$ where the CM locus has codimension greater
            than 1 and where the regular locus has codimension. The function calculates these
            open subsets using primary decomposition and the trueCodimension function in this
            package.
    Caveat
            Since the function uses primary decomposition, it only works for rings that are
            finitely presented over $\mathbf{Q}$ or $\mathbf{Z}/p\mathbf{Z}$ for a prime $p$.
            Note that the function isNormalFixed works for affine rings, without restriction on the
            base field.
    SeeAlso
            regularLocus
            CMLocus
            codimPiece
            isNormalFixed
            isS
            isR
            trueCodimension
///

doc ///
    Key
        codimPiece
        (codimPiece, ZZ, Ideal)
    Headline
        calculates the open locus where the codimension of a closed set is greater than n
    Usage
        I = codimPiece(n, J)
    Inputs
        n:ZZ
        J:Ideal
            in a ring over QQ or ZZ/p
    Outputs
        I:Ideal
            defining the open locus where Z(J) has codimension greater than n
    Description
        Text
            Let A be a commutative Noetherian ring, and let X be Spec(A). Let $Y$ be a closed subset of $X$.
            The function $x\mapsto {\rm codim}_x(Y,X)$ is lower semicontinuous on $X$, meaning that for each $n$,
            the locus of $x$ such that ${\rm codim}_x(Y,X)>n$ is open in $X$. This function returns an ideal whose
            corresponding open set in $X$ is the locus where the codimension of $Y$ is greater than $n$.
            
            For a discussion of the relevant dimension theory, one may see Section 14.2 of EGA Chapter 0
            (found at the beginning of the first piece of EGA IV). It is useful to keep in mind, in particular,
            that if we have $m<n$, then the open set where $Y$ has codimension greater than $m$ contains
            the open set where $Y$ has codimension greater than $n$. In the extreme cases, the set where
            the codimension is greater than $-1$ is all of $X$, and the set where the codimension is greater
            than the dimension of $X$ is $X\setminus Y$. (By definition if $x$ is not in $Y$, then $Y$ has
            infinite codimension in $X$ at $x$.)
    Caveat
        Since this function uses primary decomposition, it only works for rings over QQ and rings over ZZ/p.
    SeeAlso
        trueCodimension
///



--Package Tests

TEST ///
A = QQ[x];
assert(inRad(x,ideal"x2"));
///

TEST ///
A = QQ[x,y];
I = ideal"x2y";
J = ideal"xy2";
K = ideal"x2,y2"
assert(equalOpenLoci(I,J));
assert(not equalOpenLoci(I,K));
assert(inRad(x*y,I)); -- why three "=" signs above and yet only two here?
///

TEST ///
A = ZZ[x,y]/(x^3-y^2);
B = ZZ[t];
f = map(B,A,{t^2,t^3});
I = ideal(t);
assert(equalOpenLoci(relativeFlatLocus(f),I));
assert(equalOpenLoci(unramifiedLocus(f),I));
///

TEST ///
A = ZZ[x,y]/(x^3-y^2);
B = ZZ[t];
f = map(B,A,{t^2,t^3});
I = ideal(t);
assert(equalOpenLoci(relativeFlatLocus(f),I));
assert(equalOpenLoci(unramifiedLocus(f),I));
///

TEST ///
F = frac (ZZ/5[a]);
A = F[x]/(x^5-a);
assert(regularLocus(A) == A);
///

TEST ///
A = QQ[x,y];
B = QQ[x,y,z,t]/(x*y-z*t,x*z-t^2,y*t-z^2);
f = map(B,A);
isCM B == true;
equalOpenLoci(relativeFlatLocus f, B) == true;
equalOpenLoci(etaleLocus f, ideal(x*y)) == true;
equalOpenLoci(smoothLocus f, ideal(x*y)) == true;
equalOpenLoci(lciLocus f, ideal(x,y,z,t)) == true;
///


end;


More tests to add: look at examples mentioned in comments on CMLocus and GorensteinLocus

The following test is very slow! (It takes nearly a minute to run.)

TEST ///
B = QQ[p,q,r];
C = QQ[a,b,c,d]/(a^3 - b*(b-a));
phi = map(B,C,{p^2-p,p^3-p^2,p*q,r});
assert(equalOpenLocus(relativeFlatLocus phi, ideal(p*(p-1))));
///
