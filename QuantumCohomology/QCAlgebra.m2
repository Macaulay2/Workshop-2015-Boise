newPackage("QCAlgebra",
     Headline => "Quantum Cohomology for the modern world",
     Version => "0.01",
     Date => "May 27, 2015",
     Authors => {
	  {Name => "Dave Swinarski",
	   HomePage => "",
	   Email => ""},
	  {Name => "Anna Kazanova",
	   HomePage => "",
	   Email => ""},
	  {Name => "Robert Williams",
	   HomePage => "",
	   Email => ""},
	  {Name => "Corey Harris",
	   HomePage => "",
	   Email => ""}},
     --AuxiliaryFiles => true,
     AuxiliaryFiles => false,
     DebuggingMode => true,
     CacheExampleOutput =>true
     )


export { QCRing, QCQuotientRing, QCPolynomialRing,
         QCRingMap, QCRingElement,
         QCGroebnerBasis, ncGroebnerBasis,
	 Basis, qcMap, qcIdeal, qcGroebnerBasis, CheckPrefixOnly, weights, monList, qcMatrix, generatorSymbols
}

MAXDEG = 40
MAXSIZE = 1000

QCRing = new Type of Ring
QCPolynomialRing = new Type of QCRing
QCRingElement = new Type of HashTable
QCMonomial = new Type of HashTable
QCRingMap = new Type of HashTable

globalAssignment QCRing

---------------------------------------------------------------
--- Helpful general-purpose functions
---------------------------------------------------------------

removeNulls = xs -> select(xs, x -> x =!= null)

removeZeroes = myHash -> select(myHash, c -> c != 0)

minUsing = (xs,f) -> (
   n := min (xs / f);
   first select(1,xs, x -> f x == n)
)

sortUsing = (xs,f) -> (sort apply(xs, x -> (f x, x))) / last
reverseSortUsing = (xs, f) -> (reverse sort apply(xs, x -> (f x, x))) / last

-- get base m representation of an integer n
rebase = (m,n) -> (
   if (n < 0) then error "Must provide an integer greater than or equal to 1";
   if n == 0 then return {};
   if n == 1 then return {1};
   k := floor (log_m n);
   loopn := n;
   reverse for i from 0 to k list (
      digit := loopn // m^(k-i);
      loopn = loopn % m^(k-i);
      digit
   )   
)

-- list-based mingle
functionMingle = method()
functionMingle (List, List, FunctionClosure) := (xs,ys,f) -> (
   -- This function merges the lists xs and ys together, but merges them according to the function f.
   -- If f(i) is true, then the function takes the next unused element from xs.  Otherwise, the next
   -- unused element of ys is taken.  If either list is 'asked' for an element and has run out,
   -- then the list built so far is returned.
   xlen := #xs;
   ylen := #ys;
   ix := -1;
   iy := -1;
   for ians from 0 to xlen+ylen-1 list (
      if f ians then (
         ix = ix + 1;
         if ix >= xlen then break;
         xs#ix
      )
      else (
         iy = iy + 1;
         if iy >= ylen then break;
         ys#iy
      )
   )
)

-------------------------------------------
--- QCRing methods ------------------------
-------------------------------------------
coefficientRing QCRing := A -> A.CoefficientRing

generators QCRing := opts -> A -> (
    if A.?generators then A.generators else {}
)

numgens QCRing := A -> #(A.generators)

use QCRing := A -> (scan(A.generatorSymbols, A.generators, (sym,val) -> sym <- val); A)

setWeights = method()
setWeights (QCRing,List) := (A,weightList) -> (
   gensA := A.generatorSymbols;
   A#(symbol weights) = new HashTable from apply(#gensA, i -> (gensA#i,weightList#i));
   A
)

promoteHash = (termHash,A) -> (
   hashTable apply(pairs termHash, p -> (qcMonomial(p#0#monList,A),p#1))
)

-------------------------------------------
--- QCPolynomialRing methods --------------
-------------------------------------------

new QCPolynomialRing from List := (QCPolynomialRing, inits) -> new QCPolynomialRing of QCRingElement from new HashTable from inits

Ring List := (R, varList) -> (
   -- get the symbols associated to the list that is passed in, in case the variables have been used earlier.
   if #varList == 0 then error "Expected at least one variable.";
   if #varList == 1 and class varList#0 === Sequence then varList = toList first varList;
   varList = varList / baseName;
   A := new QCPolynomialRing from {(symbol generators) => {},
                                   (symbol generatorSymbols) => varList,
                                   (symbol degreesRing) => degreesRing 1,
				   (symbol CoefficientRing) => R,
                                   (symbol cache) => new CacheTable from {},
				   (symbol baseRings) => {ZZ}
                                   };
   newGens := apply(varList, v -> v <- putInRing({v},A,1));

   A#(symbol generators) = newGens;
   
   setWeights(A, toList (#(gens A):1));
      
   --- all these promotes will need to be written between this ring and all base rings.
   promote (ZZ,A) := (n,A) -> putInRing({},A,promote(n,R));


   promote (QQ,A) := (n,A) ->  putInRing({},A,promote(n,R));
        
   promote (R,A) := (n,A) -> putInRing({},A,n);
      
   promote (A,A) := (f,A) -> f;
   
   addVals := (c,d) -> (
      e := c+d;
      if e == 0 then continue else e
   );

   multVals := (c,d) -> c*d;
      
   --multKeys := (m,n) -> m | n;
   multKeys := (m,n) -> (
       mtally := new HashTable from tally (m#monList | n#monList);
       monList := flatten (for k in keys mtally list (
           toList (mtally#k : k)
       ));
       qcMonomial(monList,m.ring)   
   );

   A + A := (f,g) -> (
      newHash := removeZeroes merge(f.terms,g.terms,addVals);
      if newHash === hashTable {} then newHash = (promote(0,f.ring)).terms;
      new A from hashTable {(symbol ring, f.ring),
                            (symbol cache, new CacheTable from {("isReduced",false)}),
                            (symbol terms, newHash)}   
   );

   A ? A := (f,g) -> (
      m := first pairs (leadMonomial f).terms;
      n := first pairs (leadMonomial g).terms;
      m ? n
   );

   A * A := (f,g) -> (
      -- new way
      -- time not very predictable...
      newHash := removeZeroes combine(f.terms,g.terms,multKeys,multVals,addVals);
      -- old way
      --newHash := new MutableHashTable;
      --for t in pairs f.terms do (
      --   for s in pairs g.terms do (
      --      newMon := t#0 | s#0;
      --      newCoeff := (t#1)*(s#1);
      --      if newHash#?newMon then newHash#newMon = newHash#newMon + newCoeff else newHash#newMon = newCoeff;
      --   );
      --);
      --newHash = removeZeroes hashTable pairs newHash;
      if newHash === hashTable {} then newHash = (promote(0,f.ring)).terms;
      new A from hashTable {(symbol ring, f.ring),
                            (symbol cache, new CacheTable from {("isReduced",false)}),
                            (symbol terms, newHash)}
   );

   A ^ ZZ := (f,n) -> product toList (n:f);

   R * A := (r,f) -> promote(r,A)*f;
   A * R := (f,r) -> r*f;
   QQ * A := (r,f) -> promote(r,A)*f;
   A * QQ := (f,r) -> r*f;
   ZZ * A := (r,f) -> promote(r,A)*f;

   A * ZZ := (f,r) -> r*f;
   A - A := (f,g) -> f + (-1)*g;
   - A := f -> (-1)*f;
   A + ZZ := (f,r) -> f + promote(r,A);
   ZZ + A := (r,f) -> f + r;
   A + QQ := (f,r) -> f + promote(r,A);
   QQ + A := (r,f) -> f + r;
   A + R  := (f,r) -> f + promote(r,A);
   R + A  := (r,f) -> f + r;
   
   A ? A := (f,g) -> (leadQCMonomial f) ? (leadQCMonomial g);

   A == A := (f,g) -> #(f.terms) == #(g.terms) and (sort pairs f.terms) == (sort pairs g.terms);
   A == ZZ := (f,n) -> (#(f.terms) == 0 and n == 0) or (#(f.terms) == 1 and ((first pairs f.terms)#0#monList === {}) and ((first pairs f.terms)#1 == n));
   ZZ == A := (n,f) -> f == n;
   A == QQ := (f,n) -> (#(f.terms) == 0 and n == 0) or (#(f.terms) == 1 and ((first pairs f.terms)#0#monList === {}) and ((first pairs f.terms)#1 == n));
   QQ == A := (n,f) -> f == n;
   A == R := (f,n) -> (#(f.terms) == 0 and n == 0) or (#(f.terms) == 1 and ((first pairs f.terms)#0#monList === {}) and ((first pairs f.terms)#1 == n));
   R == A := (n,f) -> f == n;

   A
)

net QCRing := A -> (
    hasAttribute := value Core#"private dictionary"#"hasAttribute";
    getAttribute := value Core#"private dictionary"#"getAttribute";
    ReverseDictionary := value Core#"private dictionary"#"ReverseDictionary";
    if hasAttribute(A,ReverseDictionary) then toString getAttribute(A,ReverseDictionary)
    else net A.CoefficientRing | net A.generators
)

-------------------------------------------
--- QCMonomial functions ------------------
-------------------------------------------

qcMonomial = method()
qcMonomial (List,QCRing) := (monL,B) -> (
   newMon := new QCMonomial from {(symbol monList) => monL,
                                  (symbol ring) => B};
   newMon
)

QCMonomial | List := (mon,symList) -> qcMonomial(mon#monList | symList, mon.ring)
List | QCMonomial := (symList,mon) -> qcMonomial(symList | mon#monList, mon.ring)

QCMonomial | QCMonomial := (mon1,mon2) -> qcMonomial(mon1#monList | mon2#monList, mon1.ring)

net QCMonomial := mon -> (
   mon = mon#monList;
   if mon === {} then return net "";
   myNet := net "";
   tempVar := first mon;
   curDegree := 0;
   for v in mon do (
      if v === tempVar then curDegree = curDegree + 1
      else (
          myNet = myNet | (net tempVar) | if curDegree == 1 then (net "") else ((net curDegree)^1);
          tempVar = v;
          curDegree = 1;
      );
   );
   myNet | (net tempVar) | if curDegree == 1 then (net "") else ((net curDegree)^1)
)

toString QCMonomial := mon -> (
   mon = mon#monList;
   if mon === {} then return "";
   myNet := "";
   tempVar := first mon;
   curDegree := 0;
   for v in mon do (
      if v === tempVar then curDegree = curDegree + 1
      else (
          myNet = myNet | (toString tempVar) | if curDegree == 1 then "*" else "^" | curDegree | "*";
          tempVar = v;
          curDegree = 1;
      );
   );
   myNet | (toString tempVar) | if curDegree == 1 then "" else "^" | curDegree
)

degree QCMonomial := mon -> sum apply(#(mon#monList), i -> ((mon.ring).weights)#(mon#monList#i))

tensorLength = method()
tensorLength QCMonomial := mon -> #(mon#monList)

putInRing = method()
putInRing (QCMonomial, ZZ) := 
putInRing (QCMonomial, QQ) :=
putInRing (QCMonomial, RingElement) := (mon,coeff) ->
      new (mon.ring) from {(symbol ring) => mon.ring,
			   (symbol cache) => new CacheTable from {("isReduced",false)},
                           (symbol terms) => new HashTable from {(mon,promote(coeff,coefficientRing (mon.ring)))}}
putInRing (List, QCRing, ZZ) := 
putInRing (List, QCRing, QQ) :=
putInRing (List, QCRing, RingElement) := (monList,A,coeff) -> (
    mon := qcMonomial(monList,A);
    new A from {(symbol ring) => A,
                (symbol cache) => new CacheTable from {("isReduced",false)},
                (symbol terms) => new HashTable from {(mon,promote(coeff,coefficientRing A))}}
)

QCMonomial _ List := (mon,substr) -> qcMonomial((mon#monList)_substr,mon.ring)

findSubstring = method(Options => {CheckPrefixOnly => false})
findSubstring (QCMonomial,QCMonomial) := opts -> (lt, mon) -> (
   mon = mon#monList;
   lt = lt#monList;
   deg := length lt;
   if opts#CheckPrefixOnly and take(mon, deg) === lt then
      return true
   else if opts#CheckPrefixOnly then
      return false;
   if not isSubset(lt,mon) then return null;
   substrIndex := null;
   for i from 0 to #mon-1 do (
      if #mon - i < deg then break;
      if lt === mon_{i..i+deg-1} then (
         substrIndex = i;
         break;
      );
   );
   if substrIndex =!= null then (take(mon,substrIndex),take(mon,-#mon+deg+substrIndex)) else null
)

QCMonomial ? QCMonomial := (m,n) -> if (m#monList) === (n#monList) then symbol == else
                                    if #m < #n then symbol < else
                                    if #m > #n then symbol > else
                                    if (#m == #n and (toIndexList m) < (toIndexList n)) then symbol > else symbol <

QCMonomial == QCMonomial := (m,n) -> (m#monList) === (n#monList)

toIndexList = method()
toIndexList QCMonomial := mon -> (
   apply(mon#monList, x -> position(mon.ring.generatorSymbols, y -> x === y))
)

-----------------------------------------

-----------------------------------------
--- QCRingElement methods ---------------
-----------------------------------------

net QCRingElement := f -> (
   if #(f.terms) == 0 then return net "0";
   
   firstTerm := true;
   myNet := net "";
   for t in sort pairs f.terms do (
      tempNet := net t#1;
      printParens := ring t#1 =!= QQ and
                     ring t#1 =!= ZZ and
		     (size t#1 > 1 or (isField ring t#1 and 
			               numgens coefficientRing ring t#1 > 0 and
				       size sub(t#1, coefficientRing ring t#1) > 1));
      myNet = myNet |
              (if not firstTerm and t#1 > 0 then
                 net "+"
              else 
                 net "") |
              (if printParens then net "(" else net "") | 
              (if t#1 != 1 and t#1 != -1 then
                 tempNet
               else if t#1 == -1 then net "-"
               else net "") |
              (if printParens then net ")" else net "") |
              (if t#0#monList === {} and (t#1 == 1 or t#1 == -1) then net "1" else net t#0);
      firstTerm = false;
   );
   myNet
)

toStringMaybeSort := method(Options => {Sort => false})
toStringMaybeSort QCRingElement := opts -> f -> (
   sortFcn := if opts#Sort then sort else identity;
   if #(f.terms) == 0 then "0" else (
      firstTerm := true;
      myNet := "";
      for t in sortFcn pairs f.terms do (
         tempNet := toString t#1;
         printParens := ring t#1 =!= QQ and ring t#1 =!= ZZ and size t#1 > 1;
         myNet = myNet |
                 (if not firstTerm and t#1 > 0 then
                    "+"
                 else 
                    "") |
                 (if printParens then "(" else "") |
                 (if t#1 != 1 and t#1 != -1 and t#0#monList =!= {} then
                    tempNet | "*"
                  else if t#1 == -1 and t#0#monList =!= {} then "-"
                  else if t#0#monList === {} then tempNet
		  else "") |
                 (if printParens then ")" else "") |
                 (if t#0#monList =!= {} then toString t#0 else "");
         firstTerm = false;
      );
      myNet
   )
)

clearDenominators = method()
clearDenominators QCRingElement := f -> (
   if coefficientRing ring f =!= QQ then (f,f) else (
      coeffDens := apply(values (f.terms), p -> if class p === QQ then denominator p else 1);
      myLCM := lcm coeffDens;
      (f*myLCM,myLCM)
   )
)

baseName QCRingElement := x -> (
   A := class x;
   pos := position(gens A, y -> y == x);
   if pos === null then error "Expected a generator";
   A.generatorSymbols#pos
)

ring QCRingElement := QCRing => f -> f.ring

flagReduced := f -> new (ring f) from {(symbol ring) => ring f,
                                       (symbol cache) => new CacheTable from {("isReduced",true)},
                                       (symbol terms) => f.terms}


coordinates = method(Options =>{Basis=>null})

coordinates QCRingElement := opts -> f -> (
  if opts#Basis === null then
     sparseCoeffs({f})
  else
    coordinates({f},Basis=>opts#Basis)
)

coordinates List := opts -> L -> (
   if opts#Basis === null then
      sparseCoeffs L
   else
      bas := opts#Basis;
      R := ring bas#0;
      d := degree bas#0;
      if not all(L, m-> (isHomogeneous(m) and ((degree m)== d))) then 
	error "Expected homogeneous elements of the same degree.";
      mons := flatten entries basis(d,R);
      M := sparseCoeffs(bas, Monomials=>mons);
      N := sparseCoeffs(L, Monomials=>mons);
      if rank (M|N) != rank M then error "Expected elements in the span of given basis.";
      I := id_(target M);
      T := I // M;
      T*N
)


sparseCoeffs = method(Options => {Monomials => null})
sparseCoeffs QCRingElement := opts -> f -> (
  if opts#Monomials === null then
     sparseCoeffs({f})
  else
    sparseCoeffs({f},Monomials=>opts#Monomials)
)

sparseCoeffs List := opts -> L -> (
  d := if all(L, m -> m == 0) then 0 else L#(position(L,m->m!=0));
  if not all(L, m-> (isHomogeneous(m) and (m == 0 or (degree m)==(degree d)))) then 
	error "Expected homogeneous elements of the same degree.";
  B := (L#0).ring;
  R := coefficientRing B;
  mons := if opts#Monomials === null then (
              unique flatten apply(L, e-> flatten entries monomials e)) 
          else opts#Monomials;
  
  m := #mons;
  
  mons =  (mons / (m -> first keys m.terms));
  mons =  hashTable apply(m, i -> (mons#i,i));
   
  termsF := pairs (L#0).terms;
  
  coeffs := if (L#0)!=0 then (apply(termsF, (k,v) -> if v!=0 then (mons#k,0) => promote(v,R))) else {};

  l:=length L;
  if l>1 then
     for i from 1 to l-1 do (
        if not isHomogeneous L#i then error "Expected a homogeneous element.";
        if (L#i) !=0 then (
        termsF = pairs (L#i).terms;
        newCoeffs := (apply(termsF, (k,v) -> if v!=0 then (mons#k,i) => promote(v,R)));
       
	coeffs = coeffs | newCoeffs;);
     ); 
   map(R^m , R^l, coeffs)
)

monomials QCRingElement := opts -> f -> (
    qcMatrix {apply(sort keys f.terms, mon -> putInRing(mon,1))}
)
toString QCRingElement := f -> toStringMaybeSort(f,Sort=>true)
degree QCRingElement := f -> (
    fkeys := (keys f.terms);
    max (fkeys / degree)
)
tensorLength QCRingElement := f -> (
    ltf := leadTerm f;
    ltfkeys := (keys (ltf.terms));
    first (ltfkeys / tensorLength)
)
size QCRingElement := f -> #(f.terms)
leadTerm QCRingElement := f -> (
    if f.cache.?leadTerm then return f.cache.leadTerm;
    degf := degree f;
    lt := first sort select((pairs f.terms), p -> degree p#0 == degf);
    retVal := putInRing(lt#0,lt#1);
    f.cache.leadTerm = retVal;
    retVal
)
leadMonomial QCRingElement := f -> putInRing(leadQCMonomial f,1);
leadQCMonomial = f -> first keys (leadTerm f).terms;
leadMonomialForGB = f -> (
   R := coefficientRing ring f;
   if isField R then
      (leadQCMonomial f, 1_R)
   else
      (leadQCMonomial f, leadMonomial leadCoefficient f)
)
leadCoefficient QCRingElement := f -> if size f == 0 then 0 else (pairs (leadTerm f).terms)#0#1;
isConstant QCRingElement := f -> f.terms === hashTable {} or (#(f.terms) == 1 and f.terms#?(qcMonomial({},ring f)))
isHomogeneous QCRingElement := f -> (
    B := ring f;
    if f == promote(0,B) then true
    else (
    fTerms := keys f.terms;
    degf := degree first fTerms;
    all(fTerms, g -> degree g == degf))
)
terms QCRingElement := f -> (
    for p in pairs (f.terms) list (
       putInRing(p#0,p#1)
    )
)
support QCRingElement := f -> (
   varSymbols := unique flatten apply(pairs (f.terms), (m,c) -> unique toList m#monList);
   apply(varSymbols, v -> putInRing({v},f.ring,1))
)

QCRingElement * List := List => (f,xs) -> apply(xs, x -> f*x);
List * QCRingElement := List => (xs,f) -> apply(xs, x -> x*f);

toM2Ring = method(Options => {SkewCommutative => false})
toM2Ring QCRing := opts -> B -> (
   gensB := gens B;
   R := coefficientRing B;
   gensI := gens ideal B;
   abB := R [ gens B , SkewCommutative => opts#SkewCommutative];
   if gensI == {} then
      abB
   else (
      phi := ambient qcMap(abB,B,gens abB);
      abI := ideal flatten entries mingens ideal ((gensI) / phi);
      if abI == 0 then
         abB
      else
         abB/abI
   )
)

{*
toQCRing = method()
toQCRing Ring := R -> (
   --isComm := isCommutative R;
   --isExter := isExterior R;
   --if not isComm and not isExter then error "Input ring must be either strictly (-1)-skew commutative or commutative.";
   --- generate the (skew)commutivity relations
   Q := coefficientRing R;
   A := Q (gens R);
   phi := qcMap(A,ambient R,gens A);
   --commRelations := apply(subsets(gens A,2), s-> s_0*s_1+(-1)^(if isComm then -1 else 0)*s_1*s_0);
   extRelations := if isExter then apply(gens A, s -> s^2) else {};
   --- here is the defining ideal of the commutative algebra, inside the tensor algebra
   I := qcIdeal (commRelations | extRelations | ((flatten entries gens ideal R) / phi));
   commIgb := gb ideal R;
   maxDeg := (((flatten entries gens commIgb) / degree) | {{0}}) / sum // max;
   Igb := qcGroebnerBasis(I, DegreeLimit => 2*maxDeg);
   A/I
)
*}

isNormal QCRingElement := f -> (
   if not isHomogeneous f then error "Expected a homogeneous element.";
   all(gens ring f, x -> findNormalComplement(f,x) =!= null)
)

normalAutomorphism = method()
normalAutomorphism QCRingElement := f -> (
   B := ring f;
   normalComplements := apply(gens B, x -> findNormalComplement(f,x));
   if any(normalComplements, f -> f === null) then error "Expected a normal element.";
   qcMap(B, B, normalComplements)
)

findNormalComplement = method()
findNormalComplement (QCRingElement,QCRingElement) := (f,x) -> (
   B := ring f;
   if B =!= ring x then error "Expected elements from the same ring.";
   if not isHomogeneous f or not isHomogeneous x then error "Expected homogeneous elements";
   n := degree f;
   m := degree x;
   leftFCoeff := sparseCoeffs(f*x,Monomials=>flatten entries basis(n+m,B));
   rightMultF := rightMultiplicationMap(f,m);
   factorMap := (leftFCoeff // rightMultF);
   if rightMultF * factorMap == leftFCoeff then
      first flatten entries (basis(m,B) * factorMap)
   else
      null
)

getMinMaxDegrees = gensList -> (
   minQCGBDeg := minCoeffDeg := infinity;
   maxQCGBDeg := maxCoeffDeg := -infinity;
   scan(gensList, f -> (degf := tensorLength f;  --- had to change this to tensor length for GB reduction...
                        degLeadCoeff := if isField coefficientRing ring f then 0 else first degree leadCoefficient f;
                        if degf > maxQCGBDeg then maxQCGBDeg = degf;
                        if degf < minQCGBDeg then minQCGBDeg = degf;
                        if degLeadCoeff > maxCoeffDeg then maxCoeffDeg = degLeadCoeff;
                        if degLeadCoeff < minCoeffDeg then minCoeffDeg = degLeadCoeff;));
   (minQCGBDeg,maxQCGBDeg,minCoeffDeg,maxCoeffDeg)
)



leftMultiplicationMap = method()
leftMultiplicationMap(QCRingElement,ZZ) := (f,n) -> (
   B := f.ring;
   m := degree f;
   if m === -infinity then m = 0;
   nBasis := flatten entries basis(n,B);
   nmBasis := flatten entries basis(n+m,B);
   leftMultiplicationMap(f,nBasis,nmBasis)
)

leftMultiplicationMap(QCRingElement,ZZ,ZZ) := (f,n,m) -> (
   B := f.ring;
   nBasis := flatten entries basis(n,B);
   mBasis := flatten entries basis(m,B);
   leftMultiplicationMap(f,nBasis,mBasis)
)

leftMultiplicationMap(QCRingElement,List,List) := (f,fromBasis,toBasis) -> (
   local retVal;
   A := ring f;
   R := coefficientRing A;
   if not isHomogeneous f then error "Expected a homogeneous element.";
   if fromBasis == {} and toBasis == {} then (
      retVal = map(R^0,R^0,0);
      retVal
   )
   else if fromBasis == {} then (
      retVal = map(R^(#toBasis), R^0,0);
      retVal
   )
   else if toBasis == {} then (
      retVal = map(R^0,R^(#fromBasis),0);
      retVal
   )
   else (
      sparseCoeffs(f*fromBasis, Monomials=>toBasis)
   )
)

{*
TEST ///
restart
needsPackage "QCAlgebra"
A = QQ{a,b}
I = qcIdeal {a*a*a,a*a*b,a*b*a,a*b*b,b*a*a,b*a*b,b*b*a,b*b*b}
B = A/I
basis(0,B)
basis(1,B)
basis(2,B)
basis(3,B)
leftMultiplicationMap(a,-1,0)
leftMultiplicationMap(a,-1,0)
leftMultiplicationMap(a,-1,0)
///
*}

rightMultiplicationMap = method()
rightMultiplicationMap(QCRingElement,ZZ) := (f,n) -> (
   B := f.ring;
   m := degree f;
   if m === -infinity then m = 0;
   nBasis := flatten entries basis(n,B);
   nmBasis := flatten entries basis(n+m,B);
   rightMultiplicationMap(f,nBasis,nmBasis)
)

rightMultiplicationMap(QCRingElement,ZZ,ZZ) := (f,n,m) -> (   
   if f != 0 and degree f != m-n then error "Expected third argument to be the degree of f, if nonzero.";
   B := f.ring;
   nBasis := flatten entries basis(n,B);
   mBasis := flatten entries basis(m,B);
   rightMultiplicationMap(f,nBasis,mBasis)
)

rightMultiplicationMap(QCRingElement,List,List) := (f,fromBasis,toBasis) -> (
   local retVal;
   A := ring f;
   R := coefficientRing A;
   if not isHomogeneous f then error "Expected a homogeneous element.";
   if fromBasis == {} and toBasis == {} then (
      retVal = map(R^0,R^0,0);
      retVal
   )
   else if fromBasis == {} then (
      retVal = map(R^(#toBasis), R^0,0);
      retVal
   )
   else if toBasis == {} then (
      retVal = map(R^0,R^(#fromBasis),0);
      retVal
   )
   else (
      sparseCoeffs(fromBasis*f, Monomials=>toBasis)
   )
)

{*
---------------------------------------
----QCRingMap Commands -----------------
---------------------------------------

qcMap = method(Options => {Derivation=>false})
--- qcMap from Ring to QCRing not implemented.
qcMap (Ring,QCRing,List) := 
qcMap (QCRing,Ring,List) := 
qcMap (QCRing,QCRing,List) := opts -> (B,C,imageList) -> (
   genCSymbols := (gens C) / baseName;
   if opts#Derivation and B=!=C then error "Source and target of a derivation must be the same.";
   if not all(imageList / class, r -> r === B) then error "Expected a list of entries in the target ring.";
   new QCRingMap from hashTable {(symbol fuqctionHash) => hashTable apply(#genCSymbols, i -> (genCSymbols#i,imageList#i)),
                                 (symbol source) => C,
                                 (symbol target) => B,
				 (symbol Derivation) => opts#Derivation,
				 (symbol cache) => new CacheTable from {}}
)

source QCRingMap := f -> f.source
target QCRingMap := f -> f.target
matrix QCRingMap := opts -> f -> (
     if member(QCRing, aqcestors class f.target) then
        qcMatrix {(gens source f) / f}
     else
        matrix {(gens source f) / f}
)
--id _ QCRing := B -> qcMap(B,B,gens B)

QCRingMap QCRingElement := (f,x) -> (
   if x == 0 then return promote(0, target f);
   if ring x =!= source f then error "Ring element not in source of ring map.";
   C := ring x;
   if f.Derivation then
      sum for t in pairs x.terms list(   
         mon := t#0#monList;
         monImage := sum apply(# mon, j-> 
	  product apply(replace(j,f.fuqctionHash#(mon_j),mon), s-> 
		 if class s===Symbol or class s===IndexedVariable then putInRing({s},C,1) else s
		 )
	    );
         promote(sub(t#1,coefficientRing target f)*monImage,C) -- an empty sum is 0, which we need to promote 
      )
   else
      sum for t in pairs x.terms list (
         monImage := promote(product apply(t#0#monList, v -> f.fuqctionHash#v),target f);
         sub(t#1,coefficientRing target f)*monImage
      )
)


QCRingMap RingElement := (f,x) -> (
   if x == 0 then return promote(0, target f);
   if ring x =!= source f then error "Ring element not in source of ring map.";
   C := ring x;

   if f.Derivation then
      sum for t in terms x list (
         coeff := leadCoefficient t;
         mon := leadMonomial t;
         supp := support leadMonomial t;
         monImage := sum apply(# supp, j -> (
		    d:=degree(supp_j,mon);
		    product apply(# supp, k-> 
			 if k==j then  
			 -- chain rule
                            d*(f.fuqctionHash#(baseName supp_k))^(d-1)
		      	 else (supp_k)^(degree(supp_k,mon))))
	            );
         promote(coeff*monImage,C) -- an empty sum is 0, which we need to promote
      )    
   else
      sum for t in terms x list (
         coeff := leadCoefficient t;
         mon := leadMonomial t;
         monImage := promote(product apply(support mon, y -> (f.fuqctionHash#(baseName y))^(degree(y,mon))),target f);
         coeff*monImage
      )
)

QCRingMap QCMatrix := (f,M) -> (
   newMatr := applyTable(M.matrix, x -> f x);
   newM := qcMatrix newMatr;
   if isHomogeneous M and isHomogeneous f then
      assignDegrees(newM,M.target,M.source)
   else
      assignDegrees newM;
   newM
)

List / QCRingMap := (xs, f) -> apply(xs, x -> f x)

net QCRingMap := f -> (
   net "QCRingMap " | (net target f) | net " <--- " | (net source f)
)

ambient QCRingMap := f -> (
--- check whether the source or target is an QCRing
   C := source f;
   ambC := ambient C;
   genCSymbols := (gens C) / baseName;
   qcMap(target f, ambC, apply(genCSymbols, c -> f.fuqctionHash#c))
)

isWellDefined QCRingMap := f -> (
--- check whether the source or target is an QCRing
   defIdeal := ideal source f;
   liftf := ambient f;
   all(gens defIdeal, x -> liftf x == 0)
)

isHomogeneous QCRingMap := f -> (
   gensB := gens source f;
   gensC := gens target f;
   all(gensB, x -> degree (f x) == degree x or f x == 0)
)

QCRingMap _ ZZ := (f,n) -> (
   if not isHomogeneous f then error "Expected a homogeneous QCRingMap.";
   if f.cache#?"DegreeMatrices" and f.cache#"DegreeMatrices"#?n then
      return f.cache#"DegreeMatrices"#n;
   B := source f;
   C := target f;
   srcBasis := flatten entries basis(n,B);
   tarBasis := flatten entries basis(n,C);
   imageList := srcBasis / f;
   if #(unique (select(imageList, g -> g != 0) / degree)) != 1 then
      error "Expected the image of degree " << n << " part of source to lie in single degree." << endl;
   retVal := sparseCoeffs(imageList,Monomials=> tarBasis);
   if not f.cache#?"DegreeMatrices" then f.cache#"DegreeMatrices" = new CacheTable from {};
   f.cache#"DegreeMatrices"#n = retVal;
   retVal
)

QCRingMap ? QCRingMap := (f,g) -> (matrix f) ? (matrix g)

QCRingMap @@ QCRingMap := (f,g) -> (
   if target g =!= source f then error "Expected composable maps.";
   qcMap(target f, source g, apply(gens source g, x -> f g x))
)

QCRingMap QCIdeal := (phi, I) -> qcIdeal ((gens I) / phi)

QCRingMap QCGroebnerBasis := (phi,Igb) -> (
   I := qcIdeal gens Igb;
   Iphi := phi I;
   qcGroebnerBasis(Iphi, InstallGB => true)
)

QCRingMap + QCRingMap := (f,g) -> (
    if source f =!= source g and target f =!= target g then
       error "Expected maps between the same source and target.";
    qcMap(target f, source f, apply(gens source f, x -> f x + g x))
)

ZZ * QCRingMap := 
QQ * QCRingMap := 
RingElement * QCRingMap := (a,f) -> (
    qcMap(target f, source f, apply(gens source f, x -> a*(f x)))
)

QCRingMap ^ ZZ := (f,n) -> (
   if source f =!= target f then
      error "Expected a ring endomorphism.";
   A := source f;
   if n < 0 then (
      M := f_1;
      if (rank M != numgens A) then 
         error "Expected an invertible ring map.";
      gensA := qcMatrix {gens A};
      g := qcMap(A, A, flatten entries (gensA*(M^(-1))));
      fold(abs(n):g, (a,b) -> a @@ b)
   )
   else if n == 0 then qcMap(A, A, gens A)
   else fold(n:f, (a,b) -> a @@ b)
)

kernelComponent = method()
kernelComponent (ZZ,QCRingMap) := (d,f) -> (
   -- computes the kernel of a homogeneous ring map in a specified degree
   if not isHomogeneous f then error "Expected degree 0 map.";
   R := source f;
   bas := basis(d,R);
   K := mingens ker f_d;
   if K == 0 then return qcMatrix{{promote(0,R)}} else bas*K
)

gddKernel = method()
gddKernel (ZZ,QCRingMap) := (d,f) -> (
  -- computes a generating set for the kernel of a homogeneous ring map up to a specified degree
  K := {};
  for i from 1 to d do (
     << "Computing kernel in degree " << i << endl;
     K = K | flatten entries kernelComponent(i,f);
  );
  minimizeRelations(select(K,r-> r!=0))
)
*}
-------------------------------------------------------------
------- end package code ------------------------------------
{*
-------------------- timing code ---------------------------
wallTime = Command (() -> value get "!date +%s.%N")  --- %N doesn't work on Mac
wallTiming = f -> (
    a := wallTime(); 
    r := f(); 
    b := wallTime();  
    << "wall time : " << b-a << " seconds" << endl;
    r);
------------------------------------------------------------
*}
--- include the documentation
--load (currentFileDirectory | "QCAlgebra/QCAlgebraDoc.m2")

end

restart
uninstallPackage "QCAlgebra"
installPackage "QCAlgebra"
needsPackage "QCAlgebra"
QQ{x,y}
x.terms
x*x
x*y*x*y*y
f = x
g = y
removeZeroes combine(f.terms,g.terms,multKeys,multVals,addVals);

--viewHelp "QCAlgebra"

{*
--- arithmetic benchmark
restart
needsPackage "QCAlgebra"
A = QQ{x,y,z}
f = x+y+z
time(f^12);
time(f^6*f^6);
*}