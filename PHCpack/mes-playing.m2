restart
needsPackage "PHCpack"

R = CC[x,y,z]
system = {y-x^2,z-x^3,x+y+z-1}
solns = solveSystem(system)

restart
needsPackage "PHCpack"
R = CC[s,t]
rand = (R) -> (kk := coefficientRing R; map(kk,R,for x in gens R list random kk))
phi1 = rand R
phi2 = rand R
P = {s^3-t-1, s^2*t^3-s-t^2-1, s*t-1}
Q = drop(P,-1)
P1 = P - (P/phi1)
P2 = P - (P/phi2)
Q1 = Q - (Q/phi1)
Q2 = Q - (Q/phi2)
sP1 = solveSystem P1
sP2 = solveSystem P2
use R
sQ1 = solveSystem Q1
sQ2 = solveSystem Q2
use R
-- trackPaths(P1, P2, sP1) -- fails, since overdetermined
use R
trackPaths(Q1, Q2, sQ1)
trackPaths(Q1, Q2, sP1)

polys = {s^3-t-1, s^2*t^3-s-t^2-1, s*t-1}
polys = drop(polys,-1)
polys0 = polys - (polys/phi1)
polys1 = polys - (polys/phi2)
(solveSystem polys0)
netList oo
errorDepth=0
S1 = drop(polys0,-1)
S2 = drop(polys1,-1)
solveSystem S1
trackPaths(drop(polys0,-1), drop(polys1,-1), solveSystem polys0)
trackPaths(polys0, polys1, solveSystem polys0)


R = QQ[s,t]
polys = {s^3-t-1, s^2*t^3-s-t^2-1, s*t-1}
drop(oo,-1)
jacobian ideal oo
det oo
