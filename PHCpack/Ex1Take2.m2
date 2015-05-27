R = CC[a11,a12,a13,a21,a32];
m = {a11 - a12 - a13 - a32,
  a11*a12 + a11*a13 - a12*a13 + a12*a21 + a11*a32 - a13*a32,
  a11*a12*a13 + a12*a13*a21 + a11*a13*a32 + a13*a21*a32,
  a12 + a13 + a32, a12*a13 + a13*a32}

-- Make a list of options
VarList = R.gens;
ValueLists = new MutableList;
for i from 0 to length VarList - 1 do 
	ValueLists#i = VarList#i => random(CC);
ValueLists = new List from ValueLists

-- Substitute values and make list of equations for PHC
EquationList = new MutableList;
for i from 0 to length m - 1 do EquationList#i = m#i - sub(m#i, ValueLists)
EquationList = new List from EquationList

print EquationList
loadPackage("PHCpack")
sols = solveSystem(EquationList)
#sols
sols/print
