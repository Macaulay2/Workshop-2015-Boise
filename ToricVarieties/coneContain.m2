--input: cone inputCone, fan F such that inputCone is contained in a face of F
--output: a cone, the smallest cone of F containing C
needsPackage("Polyhedra");
smallestContainingCone = (F,inputCone) -> (
    recurser = (bigCone,checkCone) -> (
        for C in faces(1,bigCone) do (
            if contains(C,checkCone) then (
                return recurser(C,checkCone);
            );
        );
        --checkCone is not in any of the smaller cones
        return bigCone; 
    );
    for C in maxCones(F) do (
        if contains(C,inputCone) then (
            return recurser(C,inputCone);
        );
    );
    {*  if we've gotten this far, then inputCone isn't contained
        in any of the maximal cones of F, i.e. bad input  *}
    error "-- input cone is not contained in any cone of fan";
);

-- EXAMPLE
F = normalFan hypercube 2;
C = posHull matrix {{1,0},{1,1}};
C2 = posHull matrix {{1,-1},{1,1}};
--smallestContainingCone(F,C)
--smallestContainingCone(F,C2)

--input: fan1, source fan; fan2, target fan; M, matrix map of lattices
--output: list of tuples {c,imC} where imC is the smallest cone in fan2
--        containing the image of c (a cone in fan1) under the map M
minImageCones = (fan2,fan1,M) -> (
    return apply(maxCones(fan1),c->{c,smallestContainingCone(fan2,posHull(M*rays(c)))})
);


maxContainingCone = (F,inputCone) -> (
    for C in maxCones(F) do (
        if contains(C,inputCone) then (
            return C;
        );
    );
    {*  if we've gotten this far, then inputCone isn't contained
        in any of the maximal cones of F, i.e. bad input  *}
    error "-- input cone is not contained in any cone of fan";
);


--For each maximal cone C of fan1, returns a maximal cone from fan2
--that contains the image of C under M
maxImageCones =  (fan2,fan1,M) -> (
    return apply(maxCones(fan1),c->maxImageCones(fan2,posHull(M*rays(c))))
);


--input: M, a matrix; X and Y, source and target normal toric varieties
--output: b, a boolean value, true iff M respects the fans of X and Y
isCompatible = (Y,X,M) -> (
    for Cx in maxCones(fan(X)) do (
        xConeContained = false;
        imCx = posHull(M*rays(Cx));
        for Cy in maxCones(fan(Y)) do (
            if contains(Cy,imCx) then (
                xConeContained = true;
                break;
            );
        );
        if not xConeContained then return false;
    );
    return true;
);
