-- This script applies PHCpack.m2 to generate several instances 
-- for random values of the parameters for a system defined by
-- a polynomial map.  The blackbox solver is used to solve a first
-- instance and the solutions of this first instance are then the
-- start solutions for the path tracker to the second instance.

loadPackage("PHCpack");
makeIdentifiabilitySystem = method()
makeIdentifiabilitySystem (List) := (somemap) -> (
  R := ring ideal somemap;
  ValueList = for v in R.gens list (v => random(CC));
  phi := map(CC, R, ValueList);
  return (somemap/(f -> f - phi f), ValueList)
)

doOneLoop = method()
doOneLoop (List, List, List) := (FirstSystem, SecondSystem, sols) -> (
  sols1 = trackPaths(SecondSystem, FirstSystem, sols);
  sols2 = trackPaths(FirstSystem, SecondSystem, sols1);
  return sols2
)

doMonodromy = method(Options => {NumLoops => 5, Tolerance => 4})
doMonodromy (List) := o -> (System) -> (
  (FirstSystem, FirstSolution) = makeIdentifiabilitySystem(System);
  print FirstSolution;
  Sol = {point({for v in FirstSolution list v#1})};
  for i from 1 to o.NumLoops do (
    NewSol = doOneLoop(FirstSystem, (makeIdentifiabilitySystem(System))#0, Sol);
    for i from 0 to (length (NewSol#0#Coordinates) - 1) do (
      PointDiff = (((Sol#0)#Coordinates)#i) - (((NewSol#0)#Coordinates)#i);
      if ((round(o.Tolerance, realPart PointDiff) != 0)
      or (round(o.Tolerance, imaginaryPart PointDiff) != 0))
      then (
        print "Different solutions found for given variable ordering:";
        print gens ring (System#0);
        print "Solution one:";
        print NewSol;
        print "";
        print "Solution two:";
        print Sol;
        return;
      );
    );
    Sol = NewSol;
  );
  print "One solution found for given variable ordering:";
  print gens ring (System#0);
  print Sol;
)

doBlackbox = method()
doBlackbox (List) := (System) -> (
  print "Blackbox solution(s) for given variable ordering:";
  print gens ring (System#0);
  return solveSystem((makeIdentifiabilitySystem(System))#0);
)

R = CC[a11,a12,a13,a21,a32];
m = {a11 - a12 - a13 - a32,
     a11*a12 + a11*a13 - a12*a13 + a12*a21 + a11*a32 - a13*a32,
     a11*a12*a13 + a12*a13*a21 + a11*a13*a32 + a13*a21*a32,
     a12 + a13 + a32, a12*a13 + a13*a32};
{*
R=CC[a11, a12, a13, a21, a22, a32];
m = {a11 - a13 + a22, a11*a13 + a12*a21 - a11*a22 + a13*a22, 
  a12*a13*a21 - a11*a13*a22 + a13*a21*a32, -a13*a22, a21, 
  a13*a21}
*}
doMonodromy(m);
print "";
print doBlackbox(m);
