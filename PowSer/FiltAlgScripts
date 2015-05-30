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
    if N==0 then return H.ring;
    if N>0 and L#?(N-1) then return H#N
     else  ( for i from 1 to N-1 do (
	  for j from 1 to N-1 do (
	     (if i+j==N then genlist = append(genlist, (H#i)*(H#j)););
	     );
	 );
     return H#N= ideal genlist	
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
       T:= R[x_1..x_n, Degrees => degreeList];
       f := map(S,T,myList);
       return Q := T/ker f
    )

-- want to write script to check if a family of ideals {I_0,..,I_l} satisfies the conditions:
-- (a) I_i > I_(i+1)
-- (b) I_i I_j < I_(i+j) for all i,j
-- to do later

--
--

--
-- construct the "associated graded object"
-- to do later 

