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




{* an attempt at non-recursion
smallestContainingCone = (F,inputCone) -> (
    firstTime = true;
    currentFan = F;
    while true do (
        for C in maxCones(currentFan) do (
            if contains(C,inputCone) then (
                currentFan = fan(C);
                break;
            )
        )
        firstTime = false
    )
    if firstTime then error "-- input cone is not contained in any cone of fan";
*}



--input: M, a matrix; X and Y, source and target normal toric varieties
--output: b, a boolean value, true iff M respects the fans of X and Y
isCompatible = (M,X,Y) -> (
    for Cx in maxCones(fan(X)) do (
        xConeContained = false;
        for Cy in maxCones(fan(Y)) do (
            imCx = cone(M*rays(Cx));
            if contains(Cy,imCx) then (
                xConeContained = true;
                break;
            );
        );
        if not xConeContained then return false;
    );
    return true;
);
