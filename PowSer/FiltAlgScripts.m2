FamilyOfIdeals := new Type of GradedModule  

familyOfIdeals = method()

-- in the following want to assume input as list of ideals
-- {I_1,I_2,dots, I_l} where I_i > I_(i+1) and we don't assume I_1 = R}


familyOfIdeals(List) := FamilyOfIdeals => (L) -> (
    --- assume all ideals are defined over the same ring ---
    H := new FamilyOfIdeals;
    H.ring = ring first L;
    apply(#L, i-> H#(i+1) = L#i);
    return H
)



--the following script inductively creates a family of ideals 
FamilyOfIdeals^ZZ := Ideal =>(H,n) -> (
    N:=n;
    L := select(keys H, i-> (class i) === ZZ);
    if (N-1) > max L then error"N-1 not yet made";
    --I := H#(max L);
    genlist := {};
    if N==0 then return ideal(1_(H.ring));
    if N>0 and L#?(N-1) then return H#N
     else  ( for i from 1 to N-1 do (
	  for j from 1 to N-1 do (
	     (if i+j==N then genlist = append(genlist, (H#i)*(H#j)););
	     );
	 );
     return H#N= trim ideal genlist	
    )
)


--the following creates a finite type algebra given a family of ideals
-- with certain defining properties
-- want the the family of ideals to have the proerties that 
-- (a) I_i > I_(i+1)
-- (b) I_i I_j < I_(i+j) for all i,j

finiteTypeFilteredAlgebra = method()

finiteTypeFilteredAlgebra(List) := (L) -> (
    x:=getSymbol"x"; 
    t:=getSymbol"t";
    R:=ring first L;
    S:=R[t];
    genlist := { };
    degreelist := { };
    for i to #L-1 do ( I = L#i; 
	genlist = 
	append(genlist,
	    apply(I_*,f->promote(f,S)*t_S^(i+1)));
	)   ;
    myList := flatten genlist;
    degreeList := apply(myList, i-> degree i);
    n:= #myList;
    T:= R[x_1..x_n, Degrees => degreeList, Join=> false];
    A := QQ[gens T, Degrees => degreeList/first];
    f := map(S,T,myList);
    (ker f, A, f)
    )

-- want to write script to check if a family of ideals {I_0,..,I_l} satisfies the conditions:
-- (a) I_i > I_(i+1)
-- (b) I_i I_j < I_(i+j) for all i,j
-- to do later

--
--

--
-- construct the "associated graded object"
grIdeal = method()
grIdeal (Ring,RingMap,ZZ,FamilyOfIdeals) := (A,F,n,I) ->(
    In := I^n;
    In1 := I^(n+1);
    S:= source F;
    B := sub (basis (n,A),S);
    B' := F B;
    RT:=ring B';
    R := coefficientRing RT;
    B'' := sub (B',RT_0 => 1_R);
    (B,B'');
    C := (B'')//(gens In);
    D := map(In/In1, source C, C);
    E := sub (gens ker D,S);
    trim ideal (B*E)
    )

-- to do later 

associatedGraded = method()

associatedGraded(List,ZZ) :=(L,n) ->(
    I := familyOfIdeals(L);
    (J,A,F):=finiteTypeFilteredAlgebra(L);
    S:= ring J;
    K:=sum(for i from 0 to n list grIdeal(A,F,i,I));
    trim K
    )

end

restart
load"FiltAlgScripts.m2"


R=QQ[X,Y]

I1=ideal(X^2,Y)
I2=ideal(X^2,Y^2)
L:={I1,I2}
associatedGraded(L,2)
associatedGraded(L,3)
associatedGraded(L,4)
associatedGraded(L,5)
associatedGraded(L,10)
minimalPresentation oo
(gens ring oo)/degree
H=familyOfIdeals(L)
L
H^0
H^1
H^2
H^3
H^4
H^3

(J,A,F)=finiteTypeFilteredAlgebra(L)
F
finiteTypeFilteredAlgebra(L)
grIdeal(A,F,1,H)
grIdeal(A,F,2,H)
grIdeal(A,F,3,H)
compactMatrixForm = false
grIdeal(A,ring J,2,H)
gens o21

A
describe A
isHomogeneous(J)
netList J_*

S= ring J / J -- the "rees Hilbert series"
hilbertSeries(S)

--
R = QQ[x,y]
S=R[a,b,Degrees=>{{1,1},{2,1}}, Join=>false]
degree a
degree x_S

