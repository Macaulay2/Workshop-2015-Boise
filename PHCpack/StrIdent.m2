newPackage(
  "StrIdent",
  Version => "0.1.0", 
  Date => "29 May 2015",
  Authors => {
    {Name => "Taylor Brysiewicz"},
    {Name => "Nikki Meshkat"},
    {Name => "Jeff Sommars"},
    {Name => "Jan Verschelde"}
  },
  Headline => "Structural Identifiability",
  DebuggingMode => false
)

export{
"characteristicPoly",
"doBlackbox",
"doMonodromy",
"doMultiMonodromy",
"getCoefficients",
"makeMultiIdentifiabilitySystem",
"maxAlgInd",
"NumLoops",
"pullCoefficients",
"restrictRing",
"strIdentPolys",
"Tolerance",
"tryEvalSol"
}

--needsPackage("PHCpack");

loadPackage "PHCpack"
characteristicPoly = method()
characteristicPoly (Matrix) := (M) -> (	
  --IN: a square matrix (not a mutable matrix)
  --OUT: the characteristic polynomial of a matrix
  x := symbol x;
  S := CC[gens ring M][x];
  phi := map(S, ring M, gens coefficientRing S);
  newM := phi M;
  --create the identity matrix
  IdM := matrix(mutableIdentity(S, numgens source M));
  --find characteristic polynomial
  return(det(x*IdM - newM));
)

pullCoefficients = method()
pullCoefficients (RingElement) := (f) -> (
  --IN: a polynomial f in x
  --OUT: a list of the coefficients of x^i, excluding occurences  of 0,-1, and 1
  CoeffList := new MutableList from {};
  --Create a list to accumulate coefficients
  k := 0;
  t := 0;
  x := (gens ring f)#0;
  --loop through powers of x
  while k <= degree(x,f) do (
    --find coefficient of x^k
    C := coefficient(x^k,f);
    --if its a bad value (0,1,-1), don't include it
    if not (C == 0 or C==1 or C==-1) then (
      CoeffList#t = C;
      t=t+1;
    );
    k = k+1;
  );
  return(CoeffList);
)

restrictRing = method()
restrictRing (Matrix) := (A) -> (
  SbRing := CC[support A];    
  L := indices A;
  i := 0;
  t := 0;
  mapList := new MutableList from {};
  while i<numgens ring A do (
    if member(i,L) then 
      (mapList#i = (gens SbRing)#t; t=t+1;)
    else 
      mapList#i=0;
      i = i+1; 
  );
  MapList := new List from mapList;
  phi := map(SbRing,ring A,MapList);
  return(phi(A));
)

restrictRing (RingElement) := (A) -> (
  SbRing := CC[support A];    
  L := indices A;
  t := 0;
  mapList := new MutableList from {};
  for i from 0 to (numgens ring A - 1) do (
    if member(i,L) then 
      (mapList#i = (gens SbRing)#t; t=t+1;)
    else 
      mapList#i=0;
  );
  MapList := new List from mapList;
  phi := map(SbRing,ring A,MapList);
  return(phi(A));
)

getCoefficients = method()
getCoefficients (Matrix, List, List) := (B,I,J) -> (
  --IN: a weighted adjacency matrix corresponding to a compartment model, A
  --    a list of input indices, I
  --    a list of output indices, J
  --OUT: a list of polynomials describing a map on the variables a_(i,j) 
  R := ring B;
  AA := restrictRing(B);
  y := symbol y;
  SS := CC[join(gens ring AA,{y})];
  Phi := map(SS,ring AA, delete(y,gens SS));
  A := Phi AA;
  chA := characteristicPoly(A);
  --find characteristic polynomial
  allCoeff := pullCoefficients(chA);
  --begin a list to collect the relevant polynomials and include the coeff of chA
  for jindex from 0 to (#J - 1) do
   --scroll through inputs
  (
    j := J#jindex;
    for iindex from 0 to (#I - 1) do
    --scroll through outputs
    (
      i := I#iindex;
      --if i=j then collect relevant polynomials quickly via a submatrix
      if j==i then (
        chAjj := characteristicPoly(submatrix'(A,{j},{j})); 
        allCoeff = join(allCoeff,pullCoefficients(chAjj));
      )
      else if j!=i then ( --otherwise use a derivative trick
        AnewMM := mutableMatrix(A);
        use SS;
        AnewMM_(i,j) = y;
        Anew := matrix(AnewMM);
        chAnew := characteristicPoly(Anew);
        allCoeff = join(allCoeff,pullCoefficients(diff(y,chAnew)));
      );
    );
  );
  AllCoeff := new List from allCoeff;
  AllCoeff2 := new MutableList from {};
  oldRing := ring AllCoeff#0;    
  LL := delete((gens(oldRing))#(numgens oldRing-1),gens(ring AllCoeff#0));  
  Ring3 := CC[LL];
  for i from 0 to (length AllCoeff - 1) do (
    PHI := map(Ring3,ring AllCoeff#i,join(gens Ring3,{0}));
    AllCoeff2#i = PHI(AllCoeff#i);
  );
  finalCoeff := new List from AllCoeff2;
  return(finalCoeff);
)

makeIdentifiabilitySystem = method()
makeIdentifiabilitySystem (List) := (somemap) -> (
  R := ring ideal somemap;
  ValueList := for v in R.gens list (v => random(CC));
  phi := map(CC, R, ValueList);
  return (somemap/(f -> f - phi f), ValueList)
)

doOneLoop = method()
doOneLoop (List, List, List) := (FirstSystem, SecondSystem, sols) -> (
  sols1 := trackPaths(SecondSystem, FirstSystem, sols);
  sols2 := trackPaths(FirstSystem, SecondSystem, sols1);
  return sols2
)

doMonodromy = method(Options => {NumLoops => 5, Tolerance => 4})
doMonodromy (List) := o -> (System) -> (
  (FirstSystem, FirstSolution) := makeIdentifiabilitySystem(System);
  print FirstSolution;
  Sol := {point({for v in FirstSolution list v#1})};
  for i from 1 to o.NumLoops do (
    NewSol := doOneLoop(FirstSystem, (makeIdentifiabilitySystem(System))#0, Sol);
    for i from 0 to (length (NewSol#0#Coordinates) - 1) do (
      PointDiff := (((Sol#0)#Coordinates)#i) - (((NewSol#0)#Coordinates)#i);
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
    Sol := NewSol;
  );
  print "One solution found for given variable ordering:";
  print gens ring (System#0);
  print Sol;
)

tryEvalSol = method()
tryEvalSol (List, List, ZZ) := (NewSol, ExtraPolys, Tolerance) -> (
  ValueMap := for i from 0 to (length NewSol - 1) list ((gens(ring (ExtraPolys#0)))#i => NewSol#i);
  for Poly in ExtraPolys do (
    TestValue := sub(Poly, ValueMap);
    if ((round(Tolerance, realPart TestValue) != 0)
    or (round(Tolerance, imaginaryPart TestValue) != 0)) then
      return false;
  );
  return true;
)

makeMultiIdentifiabilitySystem = method()
makeMultiIdentifiabilitySystem (List, List) := (IndSystem, ExtraPolys) -> (
  R := ring ideal IndSystem;
  ValueList := for v in R.gens list (v => random(CC));
  phi := map(CC, R, ValueList);
  return (IndSystem/(f -> f - phi f), ValueList, ExtraPolys/(f -> f - phi f))
)

removeDups = method()
removeDups (List, ZZ) := (Sols, Tolerance) -> (
  CleanSols := new MutableList;
  ShortSols := new MutableList;
  NextIndex := 0;
  for TestSol in Sols do (
    ShortenedTestSol := for Index in (TestSol#0)#Coordinates list round(Tolerance, realPart Index) + ii*round(Tolerance, imaginaryPart Index);
    for i from 0 to NextIndex + 1 do (
      if i == NextIndex then (
        ShortSols#NextIndex = ShortenedTestSol;
        CleanSols#NextIndex = TestSol;
        NextIndex = NextIndex + 1;
        break;
      );
      if ShortenedTestSol == ShortSols#i then
        break;
    );
  );
  return(delete(null,new List from CleanSols));
)

doMultiMonodromy = method(Options => {NumLoops => 20, Tolerance => 4})
doMultiMonodromy (List, List) := o -> (IndSystem, ExtraPolys) -> (
  (FirstSystem, FirstSolution, EvalExtraPolys) := makeMultiIdentifiabilitySystem(IndSystem, ExtraPolys);
  Sol := {point({for v in FirstSolution list v#1})};
  Sols := new MutableList;
  for i from 1 to o.NumLoops do (
    NewSol := doOneLoop(FirstSystem, (makeIdentifiabilitySystem(IndSystem))#0, Sol);
    Sols#(i-1) = NewSol;
    Sol := NewSol;
  );
  TestSols := removeDups(new List from Sols, o.Tolerance);
  
  --NEED TO CLEAN UP THE CODE SO THAT MutableLists AREN'T USED
  GoodSols := new MutableList;
  for i from 0 to length TestSols - 1 do (
    if tryEvalSol(((TestSols#i)#0)#Coordinates, EvalExtraPolys, o.Tolerance) then
      GoodSols#i = (TestSols#i)#0;
  );
  return(delete(null,new List from GoodSols));
)

doBlackbox = method()
doBlackbox (List) := (System) -> (
  print "Blackbox solution(s) for given variable ordering:";
  print gens ring (System#0);
  return solveSystem((makeIdentifiabilitySystem(System))#0);
)

loadPackage "EliminationMatrices"
maxAlgInd = method()
maxAlgInd (List) := (F) -> (
  M := matrix {F};    
  JM := jacobian(M);
  mCol := maxCol(JM);
  Indices := mCol#1;
  linIndlist := for i in Indices list M_(0,i); 
  return(linIndlist);   
)

strIdentPolys = method()
strIdentPolys (Matrix, List, List) := (A,I,J) ->(
  F:=getCoefficients(A,I,J);
  sup := for f in F list(support f);
  flatsup :=unique flatten sup;
  if #(flatsup) == #F then 
    return({F,{}})
  else if #(flatsup) > #F then 
    return("This System is Underdetermined")
  else (
    oldPolys:=toList(set F - set maxAlgInd(F));
    return({maxAlgInd(F),oldPolys});
  );
)
--##########################################################################--
-- TESTS
--##########################################################################

TEST/// 
  R = CC[a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44]
  Ex2 = matrix{{a11,a12,a13},{a21,-(a12+a32),0},{0,a32,-a13}}
  System = getCoefficients(Ex2,{0},{0})
  sol=doMonodromy(System)
///;

TEST/// 
  R = CC[a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44]
  Ex1=matrix{{a11,a12,a13},{a21,a22,0},{0,a32,-a13}}
  System = getCoefficients(Ex1,{0},{0,1})
  -- Need to kill equation because PHCpack doesn't deal well with overdetermined
  -- systems yet.
  System = delete(System#(4),System)
  sol=doMonodromy(System)
///;

TEST///
  R = CC[a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44]
  Ex2=matrix{{a11,a12,a13},{a21,-(a12+a32),0},{0,a32,-a13}}
  stri := strIdentPolys(Ex2,{0},{0})
  print doBlackbox(stri#0)
///;

TEST///
  R := CC[a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44];
  Ex1:=matrix{{a11,a12,a13},{a21,a22,0},{0,a32,-a13}};
  STRI := strIdentPolys(Ex1,{0},{0,1});
  print doMultiMonodromy (STRI#0,STRI#1)
///;
end
