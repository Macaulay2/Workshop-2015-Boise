newPackage("QuantumCohomology",
     Headline => "Quantum Cohomology for the modern world",
     Version => "0.01",
     Date => "May 27, 2015",
     Authors => {
	  {Name => "Corey Harris",
	   HomePage => "http://coreyharris.name",
	   Email => "charris@math.fsu.edu"},
   	  {Name => "Anna Kazanova",
	   HomePage => "https://sites.google.com/site/annakazanova/",
	   Email => "kazanova@uga.edu"},    
	  {Name => "Dave Swinarski",
	   HomePage => "http://faculty.fordham.edu/dswinarski/",
	   Email => "dswinarski@fordham.edu"},
	  {Name => "Robert Williams",
	   HomePage => "",
	   Email => ""}
     },
     --AuxiliaryFiles => true,
     AuxiliaryFiles => false,
     DebuggingMode => true,
     CacheExampleOutput =>true
     )

export {
    generatorSymbols,
    qcRing,
    sName
}

QCRing = new Type of Ring
QCPolynomialRing = new Type of QCRing
QCRingElement = new Type of HashTable

globalAssignment QCRing

removeZeroes = myHash -> select(myHash, c -> c != 0)

possibleTableaux := (r,l) ->(
    MyList:= toList((r+1):0)..toList((r+1):l);
    desc := L -> (
	if #L < 2 then return true;
    	if L#0 < L#1 then (
	    return false
    	    )
    	else (
	    return desc(drop(L,1))
    	    )
	);
    NewList:= for i in MyList list(
    	if desc(i) then (
	    i
    	    )
    	else (
	    "drop me"
   	    )
	);
    NewList=delete("drop me", NewList);
    NewList=  for i in NewList list (
    	delete(0,i)
	);
    NewList
)


pieriProduct = (p, r, l, yt) -> (
    if #yt == 0 then (
	return {{p}}
    );
    if r == 0 then (
        return {{yt#0+p}}	
    );
    b:=0;
    if #yt==r+1 then b=yt#r;
    ilist := reverse(toList(max(0,p-yt#0+b)..min(p,l-yt#0)));
    sublist := apply(ilist, i -> (
        M :=  pieriProduct(p-i,r-1,yt#0,drop(yt,1));
	apply(M,a -> prepend(yt#0+i,a))
    ));
    return flatten sublist
)

----------------------------------------------------------------
----------------------------------------------------------------
----------------------------------------------------------------
-- Quantum tableaux multiplication
----------------------------------------------------------------
----------------------------------------------------------------
----------------------------------------------------------------

-- These functions are copied and pasted from the files
-- PieriSums.m2 and
-- RingfreeGiambelli.m2
-- in the Boise workshop GitHub repository

{* 
Implement the quantum Pieri rule

We do the first and second sums in the quantum Pieri rule separately

*}

{*
Ideas for future work: 

If memory is a problem, try iterating over
the solutions instead of recursing
*}

{* 
Note: we tried using the Polyhedra package to enumerate the possible
mu's in the first sum.  It was extremely slow.  We did not even try 
to do this for the second sum.  Actually we think enumerating the points 
was OK, what took so long was creating the polytopes.  
*}


{*
For the first sum in quantum Pieri: 
Our strategy is recursive
For each allowed of mu_1, compute the possible 
diagrams with the new bounding rectangle of width
lambda_1 and adding p-(mu_1-lambda_1) boxes
*}

-- This function was originally called pieritest
pieriFirstSum = (p, r, l, yt) -> (
    if #yt == 0 then (
	return {{p}}
    );
    if r == 0 then (
        return {{yt#0+p}}	
    );
    b:=0;
    if #yt==r+1 then b=yt#r;
    ilist := reverse(toList(max(0,p-yt#0+b)..min(p,l-yt#0)));
    sublist := apply(ilist, i -> (
        M :=  pieriFirstSum(p-i,r-1,yt#0,drop(yt,1));
	apply(M,a -> prepend(yt#0+i,a))
    ));
    return flatten sublist
)


{* For the second sum in quantum Pieri, our strategy is:

0. If there is no full column, then the second sum is zero
1. Remove a full column
2. Construct the complementary diagram cyt
3. Apply pieritest to cyt
4. Take complement of all these

*}

complementaryDiagram = (r,l,yt) -> (
    while #yt <= r do yt = append(yt,0);    
    reverse apply(#yt, i -> l-yt_i)
);

-- corrected bug May 29; when forming C, use yt#0 as l
-- (not cyt#0 like we had previously)
pieriSecondSum = (p,r,l,yt) -> (
    if (#yt < r+1) or (#yt == r+1 and yt_r == 0) then return {};    
    yt=apply(#yt, i -> yt_i -1); 
    cyt:=complementaryDiagram(r,yt#0,yt);
    C:=pieriFirstSum(l-p,r,yt#0,cyt);
    apply(#C, i-> complementaryDiagram(r,yt#0,C_i))  
);

{*
Now put them together.

Let Y be the name in our package for the quantum coefficient ring, probably either ZZ[q] or QQ[q]

Then our format to represent an element in the quantum cohomology ring 
in a ringfree way is as a list of pairs
{partition as a list of integers,coefficient as an element in Y}

Example: (1+q)*s_{2,1,1} - 3*s_{3,2} is represented as 
L= { {{2,1,1},1+q}, {{3,2},-3} }

*}


quantumPieriOneTerm = (p,r,l,T,Y) -> (
    P1:=pieriFirstSum(p,r,l,T_0);    
    P2:=pieriSecondSum(p,r,l,T_0);
    P1=apply(#P1, i -> {delete(0,P1_i),T_1});
    P2=apply(#P2, i -> {delete(0,P2_i),(Y_0)*(T_1)});
    return flatten {P1,P2}
);

simplify = (L,Y) -> (
    H:=partition(first,L);
    L=apply(keys H, k -> { k, (1_Y)*(last sum(H#k))});
    select(L, k -> last(k) != 0_Y)
);

--Do each term, then combine them
quantumPieri = (p,r,l,L,Y) -> (
    if p==0 then return simplify(L,Y);
    simplify(flatten apply(#L, i -> quantumPieriOneTerm(p,r,l,L_i,Y)),Y)
);


{* 
Ringfree Giambelli

Implement Giambelli without putting things in a ring first

*}

--load "PieriSums.m2";


-- This function takes an exponent vector of a monomial to 
-- a product of generators written out that many times
-- e.g. the exponent vector of x_2^2*x_3*x_4^3 is {0,2,1,3,0} 
-- and we map this to {2,2,3,4,4,4} 
expToList = (v) -> (
    answer:={};
    j:=0;
    for i from 0 to #v-1 do (     
        answer=append(answer,apply(v_i, j -> i+1)) 
    );
    flatten answer
);

-- We compute the determinant of the Giambelli matrix 
-- in an abstract ring that the user should never see
-- and then use it to write out the Pieri multiplications that we need to do
giambelliDet = (l,yt) -> (
    s := getSymbol "s";
    S := QQ(monoid[s_1..s_l]);
    M:=matrix apply(#yt, i -> apply(#yt, j -> if yt_i+j-i > l then 0_S else if yt_i+j-i< 0 then 0_S else if yt_i+j-i == 0 then 1_S else S_(s_(yt_i+j-i))));
    L:=listForm(det(M));
    apply(#L, i -> { expToList(L_i_0),L_i_1})
);

-- This function multiplies an element g represented by the list L
-- by several horizontal strips, given in plist
iteratedPieri = (plist,r,l,L,Y) -> (
    plist = reverse plist;
    for i from 0 to #plist-1 do (
        L=quantumPieri(plist_i,r,l,L,Y)
    );    
    L
);

-- The main function.  Install as multiplication in the quantum cohomology ring
-- The inputs are: r and l indicate we're in the Grassmannian Gr(r+1,r+1+l)
-- yt1 and yt2 are the two diagrams that we want to multiply
-- Y should be either ZZ[q] or QQ[q]
quantumMultiplication = (r,l,yt1,yt2,Y) -> (
    W:=giambelliDet(l,yt2);
   -- print concatenate("giambelliDet(l,yt2) = ",toString(W)) << endl;
    L:=apply(#W, i -> iteratedPieri(W_i_0,r,l,{ {yt1,1}},Y));
   -- for i from 0 to #L-1 do (print concatenate(toString(W_i_1)," ",toString(L_i)) << endl);
    L=flatten apply(#L, i-> apply(#(L_i), j -> {L_i_j_0,(W_i_1)*(L_i_j_1)}));
    simplify(L,Y)
)



new QCPolynomialRing from List := (QCPolynomialRing, inits) -> new QCPolynomialRing of QCRingElement from new HashTable from inits

qcRing = method()
qcRing (ZZ,ZZ,String,String) := (r,l,s,q) -> (
   varList := possibleTableaux(r,l);
   R:=QQ(monoid[getSymbol q]);
   -- get the symbols associated to the list that is passed in, in case the variables have been used earlier.
   if #varList == 0 then error "Expected at least one variable.";
   if #varList == 1 and class varList#0 === Sequence then varList = toList first varList;
   --varList = varList / baseName;
   if r < 0 then error "Expected a non-negative subdimension";
   if l < 1 then error "Expected a possitive dimension for the vector space";
   if s == q then error "s and q should be different";
   A := new QCPolynomialRing from {(symbol generators) => {},
                                   (symbol generatorSymbols) => varList,
				   (symbol CoefficientRing) => R,
                                   (symbol cache) => new CacheTable from {},
				   (symbol baseRings) => {ZZ},
				   (symbol sName) => getSymbol s
                                   };
   newVars := for v in varList list (
       new A from {(symbol terms) => new HashTable from {v => 1},
	       	   (symbol ring) => A,
		   (symbol cache) => new CacheTable from {}
	           }
   );
   
   --newGens := apply(varList, v -> v <- putInRing({v},A,1));

   A#(symbol generators) = newVars;
   
   promote(ZZ,A) := (z,A) -> (
       new A from hashTable {(symbol ring, A),
                            (symbol terms, new HashTable from {{} => sub(z,R)})}
   );
   
   promote(QQ,A) := (z,A) -> (
       new A from hashTable {(symbol ring, A),
                            (symbol terms, new HashTable from {{} => sub(z,R)})}
   );

  addVals := (c,d) -> (
      e := c+d;
      if e == 0 then continue else e
  ); 

   multVals := (c,d) -> c*d;
   

   A + A := (f,g) -> (
      newHash := removeZeroes merge(f.terms,g.terms,addVals);     
      if newHash === hashTable {} then newHash = (promote(0,f.ring)).terms;
      new A from hashTable {(symbol ring, f.ring),
                            (symbol terms, newHash)}   
   );
   
   ZZ * A := (z,a) -> (
       new A from hashTable {(symbol ring, a.ring),
	                     (symbol terms, applyValues(a.terms,x -> z*x))}
   );
   A * ZZ := (a,z) -> z*a;
   A - A := (c,d) -> c + (-1)*d;
   - A := (a) -> (-1)*a;
   A + ZZ := (a,z) -> a + promote(z,A);
   ZZ + A := (z,a) -> a + z;
   A / ZZ := (a,z) -> (1/z) * a;

   QQ * A := (z,a) -> (
       new A from hashTable {(symbol ring, a.ring),
	                     (symbol terms, applyValues(a.terms,x -> z*x))}
   );
   A * QQ := (a,z) -> z*a;
   A + QQ := (a,z) -> a + promote(z,A);
   QQ + A := (z,a) -> a + z;
   A / QQ := (a,z) -> (1/z) * a;
   
   R * A := (z,a) -> (
       new A from hashTable {(symbol ring, a.ring),
	                     (symbol terms, applyValues(a.terms,x -> z*x))}
   );
   A * R := (a,z) -> z*a;
   A + R := (a,z) -> a + promote(z,A);
   R + A := (z,a) -> a + z;

   A == A := (a,b) -> a.terms === b.terms;
   
   A * A := (a,b) -> (
       --<< a << " *  " << b << endl;
       {*
       if #(a.terms) == 1 and #(b.terms) == 1 and #(first keys a.terms) == 1 then (
       	   coeff := (first values a.terms) * (first values b.terms);
	   termkeys := pieriProduct(first first keys a.terms, r, l, first keys b.terms);
	   putInRing(termkeys,coeff,A)
       ) else if #(a.terms) == 1 and #(b.terms) == 1 and #(first keys b.terms) == 1 then (
	   a*b
       ) else if #(a.terms) > 1 then (
       	   sum ( for t in keys a.terms list (putInRing({t},(a.terms)#t,A))*b )
       ) else if #(b.terms) > 1 then (
       	   sum ( for t in keys b.terms list a*(putInRing({t},(b.terms)#t,A)) )
       ) else promote(1,R)
       *}
       if #(a.terms) > 1 then (
       	   sum ( for t in pairs a.terms list (putInRing({t},(a.terms)#(t#0),A))*b )
       ) else if #(b.terms) > 1 then (
       	   sum ( for t in pairs b.terms list a*(putInRing({t},(b.terms)#(t#0),A)) )
       ) else (
       	   (at,ac) := first pairs a.terms;
       	   (bt,bc) := first pairs b.terms;
           putInRing(quantumMultiplication(r,l,at,bt,R),ac*bc,A)
       )
   );
   
   A ^ ZZ := (a,n) -> product toList (n:a);

   A
)
   
use QCRing := A -> (
    scan(A.generatorSymbols, A.generators, (l,val) -> (A.sName)_l <- val);
    getSymbol(toString (A.CoefficientRing)_0) <- (A.CoefficientRing)_0;
    A
)

putInRing = method()
putInRing (List,QCRing) := (lst, A) -> (
    termlist := new HashTable from lst;
    new A from {(symbol ring) => A,
	        (symbol cache) => new CacheTable from {},
		(symbol terms) => termlist}
)

putInRing (List,ZZ,QCRing) :=
putInRing (List,RingElement,QCRing) := (lst, z, A) -> (
    termlist := new HashTable from apply(lst, l -> (l#0, sub(z,A.CoefficientRing)*sub(l#1,A.CoefficientRing)));
    new A from {(symbol ring) => A,
	        (symbol cache) => new CacheTable from {},
		(symbol terms) => termlist}
)


QCPolynomialRing _ List := (A,l) -> A.generators#(position(A.generatorSymbols,a -> a === l));

net QCRing := A -> (
    hasAttribute := value Core#"private dictionary"#"hasAttribute";
    getAttribute := value Core#"private dictionary"#"getAttribute";
    ReverseDictionary := value Core#"private dictionary"#"ReverseDictionary";
    if hasAttribute(A,ReverseDictionary) then toString getAttribute(A,ReverseDictionary)
    else net A.CoefficientRing | net A.generators
)


---------------------------------------------
-- QCRingElement methods
---------------------------------------------

QCRingElement _ List := (q,l) -> (
    q.terms#l
)

toString QCRingElement := q -> (
    A := q.ring;
    C := A.CoefficientRing;
    concatenate between("+",apply(keys q.terms, k -> toString(q.terms#k)|"*"|toString(A.sName_k)))
)

net QCRingElement := q -> (
    s := q.ring.sName;
    sum ( for key in sort(keys q.terms) list (
    	if #key == 1 then (q.terms)#key * (hold s)_(first key)
	else (
    	    coeff := (q.terms)#key;
    	    (coeff)*((hold s)_(toSequence key))
	)
    ))
)

-- leadTerm QCPolynomialRing := 

end

restart
uninstallPackage "QuantumCohomology"
--installPackage "QuantumCohomology"
debug needsPackage "QuantumCohomology"
QH = qcRing(3,4,"s","q")
e = 3*q*s_{2,1}+(43/3)*(q^4*s_{4,2,1})
s_{1} * s_{2}
(5*s_{1}) * s_{2}
(s_{1} + s_{2,1}) * s_{2}
(s_{2}) * (s_{1} + s_{2,1}) 
s_{4,4,4,4} * s_{4,4,4,4}

toString e
--QH.generators
3 * QH_{3,2}
QH_{3,2} * 3
QH_{3,2,2} + QH_{3,2,1} + QH_{3,2,1} - QH_{3,2,2}
QH_{4,2} + 86
7 + QH_{4,2}
3/2 * QH_{3,2}
QH_{3,2} * (3/7)
QH_{3,2,2} + QH_{3,2,1} + QH_{3,2,1} - QH_{3,2,2}
QH_{4,2} + 86/3
7/3 + QH_{4,2}
QH_{3,1} / 2
3*QH_{3,1} / (2/7)
try QH_{4} / 0
2*QH_{4} == QH_{4}+QH_{4}
3*QH_{4} == QH_{4}+QH_{4}

e = 3*q*s_{2,1}+(43/3)*(q^4*s_{4,2,1})
(4/3*s_{2,1})*(3*q*s_{1})
(4/3*s_{1})*(3*q*s_{2,1})

s_{1} * s_{2}
(s_{1} + s_{2}) * (s_{1})
(s_{1}) * (s_{2} + s_{1})
(s_{1} + s_{2}) * (s_{2} + s_{1})


((s_{1} + s_{2}) * (s_{2} + s_{1}))^3
