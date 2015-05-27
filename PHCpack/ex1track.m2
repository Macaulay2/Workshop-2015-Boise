-- This script applies PHCpack.m2 to generate several instances 
-- for random values of the parameters for a system defined by
-- a polynomial map.  The blackbox solver is used to solve a first
-- instance and the solutions of this first instance are then the
-- start solutions for the path tracker to the second instance.

makeIdentifiabilitySystem = method()
makeIdentifiabilitySystem (List) := (somemap) -> (
  -- IN: a list of polynomials that represent some map
  -- OUT: the map evaluated at random values for the parameters, and
  --      the random values generated for the parameters.
  R := ring ideal somemap; 
  parameters = R.gens;
  subargs = {};
  valpars = {};
  local randomvalue;
  for v in parameters do (
      randomvalue = random(CC);
      valpars = append(valpars, randomvalue);
      subargs = append(subargs, v=>randomvalue);
  );
  result = {};
  for pol in m do (
      result = append(result, pol - sub(pol, subargs))
  );
  return (result, valpars)
)

R = CC[a11,a12,a13,a21,a32];
m = {a11 - a12 - a13 - a32,
     a11*a12 + a11*a13 - a12*a13 + a12*a21 + a11*a32 - a13*a32,
     a11*a12*a13 + a12*a13*a21 + a11*a13*a32 + a13*a21*a32,
     a12 + a13 + a32, a12*a13 + a13*a32}

(tosolvefirst, asolutionfirst) = makeIdentifiabilitySystem(m)
(tosolvesecond, asolutionsecond) = makeIdentifiabilitySystem(m)

stdio << "a solution : " << asolutionfirst << endl;

loadPackage("PHCpack")
stdio << "solving the first system ..." << endl;
sols = solveSystem(tosolvefirst)
sols/print

stdio << "tracking paths from first to second ..." << endl;

sols2 = trackPaths(tosolvesecond, tosolvefirst, sols)
sols2/print

stdio << "tracking paths from second to first ..." << endl;

sols3 = trackPaths(tosolvefirst, tosolvesecond, sols2)
sols3/print
