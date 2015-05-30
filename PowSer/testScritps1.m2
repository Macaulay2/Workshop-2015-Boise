
--  some code here is modified from filtered complex code


FilteredVectorSpace = new Type of HashTable

filteredVectorSpace = method()

filteredVectorSpace(List) := HashTable => (L) -> (
    maps = L;
       V = target maps#0;-- By default the ambient vspace is target of first map.           
 P := {0 => V} | apply (#maps,  p -> p+1 => image maps#p);
  new FilteredVectorSpace from reverse(P)
   )


FilteredVectorSpace^ ZZ := Module => (V,j) -> (
    Keys := keys V;
   Max := max Keys;
   Min := min Keys;
   if j >= Min and j <= Max  then return V#j else (
       if j < Min then (return V#Min) else( return image(0*id_(V#Min)))
       );
    )


gr = method()
gr(FilteredVectorSpace,ZZ) := (V,n) -> (
    return V^n/V^(n+1)
    )

hfgr = method()

hfgr(FilteredVectorSpace,ZZ,ZZ):= (V,i,j) -> ( 
    hilbertFunction(j,V^i/V^(i+1))
    )
 
---
---
---

-- given ideals I,J,K with IJ < K want to construct the bilinear map I x J --> K
multiplyIdeals = method()

multiplyIdeals(Ideal,Ideal,Ideal) :=(I,J,K) ->(
    if isSubset(I*J,K)==false then error"IJ not in K" 
    else(
    i :=  module I;
    j := module J;
    k := module K;
    r := I.ring; -- assum I,J,K are ideals in the same ring
    f := inducedMap(r^1, k, id_k);
    g := inducedMap(r^1,i,id_i) ** inducedMap(r^1,j,id_j);  
    return g//f);    
)



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

spots := (H) -> ( select(keys H, i-> (class i) === ZZ))


-- the following will extend the family of ideals inductively by a suitable future rule
powerFamilyOfIdeals = method()

powerFamilyOfIdeals(FamilyOfIdeals,ZZ) := Ideal => (H,N) -> (
    L := select(keys H, i-> (class i) === ZZ);
    if (N-1) > max L then error"N-1 not yet made";
    --I := H#(max L);
    genlist := {};
    if N==0 then return H.ring;
    if N>0 and L#?N then return H#N
     else  ( for i from 1 to N-1 do (
	  for j from 1 to N-1 do (
	     (if i+j==N then genlist = append(genlist, (H#i)*(H#j)););
	     );
	 );
     return H#N= ideal genlist	
    )
)

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


--
--

