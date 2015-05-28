
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
    H.ring = (L#0).ring;
    H#0 = H.ring; -- want H#0 to be our ambient ring
    apply(#L, i-> H#(i+1) = L#i);
    return H
)


spots = method()
spots(HashTable) := (H) ->( 
    select(keys H, i-> (class i) === ZZ)
    )


-- the following will extend the family of ideals by a suitable power
-- of the last ideal of the list
powerFamilyOfIdeals = method()

powerFamilyOfIdeals(FamilyOfIdeals,ZZ) := Ideal => (H,i) -> (
    L := spots(H);
    I= H#(max L);
    if i>=0 and L#?i then return H#i
     else  return I^i
    )






