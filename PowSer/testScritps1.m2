

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


----
----

multiplyIdeals(I,I,I)
multiplyIdeals(I,I,I^2)

----


R=QQ[x,y,z]


I = ideal(x,y,z)

gens I

I^2

--
-- Want to construct the R-bilinear map I x I --> I^2, defined by (a,b) \mapsto ab 

I_{0}
I_{1}

I_{2}

I_(0)
I_(1)
I_(2)


image I_{0}**I_{1}

id_(module I)

I = module ideal vars R

J = module ideal vars R

K = module (ideal vars R)^2


f = inducedMap(R^1,K,id_K)

a = id_I

b = id_J

g = inducedMap(R^1,I,a)**inducedMap(R^1,I,b)

target g

g//f








I_(1)*I_(2)

matrix({{I_(1)*I_(2)}})

map(K,R^1,matrix({{I_(1)*I_(2)}}))


K=I^2

map(K,R^1,matrix({{I_(1)*I_(2)}}))



K_{4}
K_(4)




K_{0}
K_{1}
K_(0)
K_(3)
K_(4)
K_(5)

gens K

---
---

--AssociatedGradedObject = new Type of HashTable

--associatedGradedObject = method()

--associatedGradedObject(FilteredVectorSpace) := HashTable =>(V) -> (
-V--spots = keys V;
--gr := new HashTable from apply(spots, i-> i=> V^(i)/V^(i+1));
--return gr
--    )


--AssociatedGradedObject^ ZZ := Module =>(Gr, j) -> (
--    spots = keys Gr;
--    if spots#?j then return Gr#j else return $V
--    )

--
--
-- scratch test

keys V

k=QQ
V=k^4
f0=id_V
V0=image f0
f1=map(V,V,matrix(k,{{0,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,0}}))
V1=image f1

L={f0,f1}


apply(#L,i->i)

V = filteredVectorSpace({f0,f1})
l = keys V

V^2

V^3

V^4

V^1

V^(-1)

V^(-2)

V^0

gr(V,0)
gr(V,-1)
gr(V,3)
hilbertFunction(3,gr(V,2))
hilbertFunction(3,V)

hilbertFunction(V)

R=k[x,y]

W=R^{-1}++R^{-2}

hilbertFunction(3,W)

g0=id_W
W0=image g0
g1=map(W,W,matrix(R,{{0,0},{0,1}}))
W1=image g1

L={g0,g1}

WW=filteredVectorSpace({g0,g1})

gr(WW,1)

hilbertFunction(4,gr(WW,1))

hfgr(WW,2,2)

hilbertSeries(gr(WW,1))

hilbertSeries(gr(WW,3))

R
vars R

S=QQ[a,b]

hfgr(WW,2,2)*a^2*b^2


-----


