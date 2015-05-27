loadPackage "SchurRings"
R = schurRing(QQ,s,4)
f = s_{3,2,1}^6

getCoefficient = (s, f) ->  (
    tms := select(listForm f, x -> x#0 == s);
    if #tms == 0 then 0_(coefficientRing ring f)
    else tms#0#1)

getCoefficient({9,9,9,9},f)
getCoefficient({10,10,10,10},f)

-- The following doesn't work:
coefficient(s_{9,9,9,9},f)
-- it outputs:
stdio:12:1:(3): error: expected polynomial ring
