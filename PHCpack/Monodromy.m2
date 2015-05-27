-- This script applies PHCpack.m2 to generate several instances 
-- for random values of the parameters for a system defined by
-- a polynomial map.  The blackbox solver is used to solve a first
-- instance and the solutions of this first instance are then the
-- start solutions for the path tracker to the second instance.

loadPackage("PHCpack");
makeIdentifiabilitySystem = method()
makeIdentifiabilitySystem (List) := (somemap) -> (
  -- IN: a list of polynomials that represent some map
  -- OUT: the map evaluated at random values for the parameters, and
  --      the random values generated for the parameters.
  R := ring ideal somemap;
  -- Create an option list of random complex values
  ValueList = for v in R.gens list (v => random(CC));
  -- Create a map from the ring to the complex numbers
  -- with respect to our random complex values
  phi := map(CC, R, ValueList);
  -- Apply the map to find values
  return (somemap/(f-> f - phi f), ValueList)
)

doOneLoop = method()
doOneLoop (List, List, Point) := (FirstSystem, SecondSystem, sols) -> (
  -- IN: 
  -- OUT:
  stdio << "tracking paths from first to second ..." << endl;
  sols1 = trackPaths(SecondSystem, FirstSystem, sols);
  stdio << "tracking paths from second to first ..." << endl;
  sols2 = trackPaths(FirstSystem, SecondSystem, sols1);
  return sols2
)

doOneLoop (List, List, List) := (FirstSystem, SecondSystem, sols) -> (
  -- IN: 
  -- OUT:
  stdio << "tracking paths from first to second ..." << endl;
  sols1 = trackPaths(SecondSystem, FirstSystem, sols);
  stdio << "tracking paths from second to first ..." << endl;
  sols2 = trackPaths(FirstSystem, SecondSystem, sols1);
  return sols2
)

doMonodromy = method(Options => {NumLoops => 5, UseBlackbox => False})
doMonodromy (List) := o -> (System) -> (

  (FirstSystem, FirstSolution) = makeIdentifiabilitySystem(System);
  (SecondSystem, SecondSolution) = makeIdentifiabilitySystem(System);
  FirstSolution = for v in FirstSolution list v#1;
  SecondSolution = for v in SecondSolution list v#1;

  if o.UseBlackbox == False then 
    sols = solveSystem(FirstSystem)
  else
    sols = {point({FirstSolution})};
  doOneLoop(FirstSystem, SecondSystem, sols);
  for i from 1 to o.NumLoops do sols = doOneLoop(FirstSystem, SecondSystem, sols);
  return sols;
)

R = CC[a11,a12,a13,a21,a32];
m = {a11 - a12 - a13 - a32,
     a11*a12 + a11*a13 - a12*a13 + a12*a21 + a11*a32 - a13*a32,
     a11*a12*a13 + a12*a13*a21 + a11*a13*a32 + a13*a21*a32,
     a12 + a13 + a32, a12*a13 + a13*a32}

print doMonodromy(m)
