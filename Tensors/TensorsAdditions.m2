newPackage(
	"TensorsAdditions",
    	Version => "0.1", 
    	Date => "May 27, 2015",
    	Authors => {
	     {Name => "Hirotachi Abo"},
	     {Name => "Roberto Barrera"},
	     {Name => "Robert Krone"},
	     {Name => "Benjamin Reames"},
	     {Name => "Zach Teitler"}
	     },
    	Headline => "tensors",
	AuxiliaryFiles => true,
	PackageExports => {
	  "Tensors"
	},
    	PackageImports => {
	    "PHCpack"
	    },
        DebuggingMode => true 
)

export {
  symmetrize,
  tensorEigenvectors,
  tensorToPolynomial,
  tensorToMultilinearForm,
  multiplicationTensor,
  eigenDiscriminant,
  tensorEigenvectorsCoordinates,
  isSymmetric,
  polynomialToTensor,
  tensorFromSlices,
  flattenTensor,
  factorMap
}

symmetrize = method()
symmetrize (Tensor) := (T) -> symmetrize(T,toList (0..#(tensorDims T)-1))
symmetrize (Tensor,List) := (T,L) -> (
    d := #(tensorDims T);
    R := ring T;
    S := apply(permutations L, p->(
	    ind := new MutableList from (0..d-1);
	    scan(#L, j->(ind#(L#j) = p#j));
	    T@(toList ind)
	    ));
    (1_R/((#L)!))*(sum S)
    )

isSymmetric = method()
isSymmetric Tensor := T -> (
    D := tensorDims T;
    if not all(#D, i->D#i == D#0) then return false;
    S := permutations(#D);
    all(S, s->T@s == T)
    )

contract (Tensor,List,List) := (T,K,L) -> (
    D := tensorDims T;
    KD := apply(K, k->D#k);
    LD := apply(L, k->D#k);
    if KD != LD then error "dimension mismatch";
    Tslices := apply((#K:0)..<(toSequence KD), i-> (
	    sliceList := new MutableList from (#D:null);
	    scan(#K, j->(sliceList#(K#j) = i#j; sliceList#(L#j) = i#j));
	    T_(toList sliceList)
	    ));
     sum toList Tslices
     );
contract (Tensor,Number,Number) := (T,k,l) -> contract(T,{k},{l})
contract (Tensor,Tensor,List,List) := (T,U,K,L) -> (
    Td := tensorDims T;
    Ud := tensorDims U;
    KD := apply(K, k->Td#k);
    LD := apply(L, k->Ud#k);
    if KD != LD then error "dimension mismatch";
    slices := apply((#K:0)..<(toSequence KD), i-> (
	    TsliceList := new MutableList from (#Td:null);
	    UsliceList := new MutableList from (#Ud:null);
	    scan(#K, j->(TsliceList#(K#j) = i#j; UsliceList#(L#j) = i#j));
	    Tslice := T_(toList TsliceList);
	    Uslice := U_(toList UsliceList);
	    if #K == #Td or #K == #Ud then Tslice*Uslice else Tslice**Uslice
	    ));
     sum toList slices
     );
contract (Tensor,Tensor,Number,Number) := (T,U,k,l) -> contract(T,U,{k},{l})

tensorEigenvectors = method()
tensorEigenvectors (Tensor,Number,Ring,RingElement) := (T,k,S,x) -> (
    R := ring T;
    d := tensorDims T;
    n := d#0;
    xpos := position(gens S, y->y==x);
    v := new MutableList from (n:0_S);
    for ind in (#d:0)..(#d:n-1) do (
	monList := toList apply(#ind, j->(if j != k then S_(xpos + ind#j) else 1_S));
	mon := product monList;
	v#(ind#k) = v#(ind#k) + sub(T_ind,S)*mon;
	);
    minors(2, matrix{toList v, take(gens S,{xpos,xpos + numgens S - 1})})
    );
tensorEigenvectors (Tensor,Number,Ring) := (T,k,S) -> tensorEigenvectors (T,k,S,S_0)
tensorEigenvectors (Tensor,Number,Symbol) := (T,k,x) -> (
    R := ring T;
    n := (tensorDims T)#0;
    S := R[apply(n,i->x_i)];
    tensorEigenvectors(T,k,S,S_0)
    );

tensorToPolynomial = method()
tensorToPolynomial (Tensor,Symbol) := (T,x) -> (
    R := ring T;
    D := tensorDims T;
    if not all(#D, i->(D#i == D#0)) then error "tensor is not square";
    n := D#0;
    S := R[apply(n,i->x_i)];
    tensorToPolynomial(T,S,S_0)
    );
tensorToPolynomial (Tensor,Ring) := (T,S) -> tensorToPolynomial(T,S,S_0)
tensorToPolynomial (Tensor,Ring,RingElement) := (T,S,x) -> (
    R := ring T;
    D := tensorDims T;
    n := D#0;
    xpos := position(gens S, y->y==x);
    xTen := makeTensor(take(gens S,{xpos,xpos+n-1}));
    U := xTen ^** #D;
    L := toList (0..#D-1);
    contract(sub(T,S),U,L,L)
    );

tensorToMultilinearForm = method()
tensorToMultilinearForm (Tensor,Ring) := (T,S) -> tensorToMultilinearForm(T,S,S_0)
tensorToMultilinearForm (Tensor,Ring,RingElement) := (T,S,x) -> (
    D := tensorDims T;
    xpos := position(gens S, y->y==x);
    U := null;
    for n in D do (
	u := makeTensor(take(gens S, {xpos, xpos + n -1}));
	if U === null then U = u else U = U**u;
	xpos = xpos + n;
	);
    L := toList (0..#D-1);
    contract(sub(T,S),U,L,L)
    );
tensorToMultilinearForm (Tensor,Symbol) := (T,x) -> (
    R := ring T;
    D := tensorDims T;
    varList := flatten apply(#D, i->apply(D#i, j->x_(i,j)));
    S := R[varList];
    tensorToMultilinearForm(T,S)
    )
    

tensorModule Tensor := T -> class T

tensor Matrix := o -> M -> makeTensor entries M

tensor Vector := o -> V -> makeTensor entries V

multiplicationTensor = method()
multiplicationTensor Ring := R -> (
    Bmatrix := basis R;
    B := flatten entries Bmatrix;
    K := coefficientRing R;
    V := tensorModule(K, {#B});
    L := for i from 0 to #B-1 list (
	for j from 0 to #B-1 list (
	    pVect := sub(last coefficients(B#i * B#j, Monomials=>Bmatrix), K);
	    pTens := makeTensor flatten entries pVect;
	    V_(1:i) ** V_(1:j) ** pTens
	    )
	);
    sum flatten L
    )


eigenDiscriminant = method()
eigenDiscriminant (Number,Number,Ring) := (n,d,Sa) -> (
    K := coefficientRing Sa;
    x := symbol x;
    vx := toList apply(n,i->x_i);
    vs := gens Sa;
    S := K[vs,vx];
    vx = take(gens S, -n);
    vs = take(gens S, n^d);
    T := genericTensor(S,toList (d:n));
    I := tensorEigenvectors(T,0,S,first vx);
    jj := diff(transpose matrix{vx},gens I);
    singI := minors(n-1,jj)+I;
    J := saturate(singI,ideal vx);
    sub(eliminate(vx,J),Sa)
    )

tensorEigenvectorsCoordinates = method()
tensorEigenvectorsCoordinates (Tensor,Number,Symbol) := (T,k,x) -> (
    n := (tensorDims T)#0;
    I := tensorEigenvectors(T,k,x);
    R := ring I;
    S := CC[toSequence entries vars R];
    J := sub(I,S);
    rr := (vars S | matrix{{1_S}})*transpose random(CC^1,CC^(n+1));
    L := J + ideal rr;
    F := first entries gens L;
    solveSystem(F)
    )


polynomialToTensor = method()
polynomialToTensor (RingElement ) := (f) -> (
    d := degree f;
    n := numgens ring f;
    M := tensorModule(QQ, toList(d_0:n));
    T := apply(terms f, m -> (
	    c := leadCoefficient m;
	    j := toSequence(flatten apply(#(exponents m)_0, 
		i -> toList(((exponents m)_0)_i:i)
		));
	    c*M_j
	    ));
    symmetrize(sum T)
    )

substitute (Tensor,Thing) := (T,S) -> (
    E := entries T;
    E = apply(E, e->sub(e,S));
    M := tensorModule(ring first E,tensorDims T);
    tensor(M,E)
    )

Tensor ^** Number := (T,n) -> (
    if n == 0 then return 1_(ring T);
    U := T;
    for i from 1 to n-1 do U = U**T;
    U
    )

tensorFromSlices = method()
tensorFromSlices List := S -> (
    Sf := ultimate(flatten,S);
    D1 := initialDimensions S;
    fSf := first Sf;
    D2 := if instance(fSf,Tensor) then tensorDims fSf else {};
    R := ring fSf;
    M := tensorModule(R,D1|D2);
    tensor(M, if instance(fSf,Tensor) then flatten apply(Sf,entries) else Sf)
    )

flattenTensor = method()
flattenTensor (Tensor,List) := (T,L) -> (
   D := tensorDims T;
   LD := apply(L, k->D#k);
   Tslices := apply((#L:0)..<(toSequence LD), i-> (
	   sliceInd := new MutableList from (#D:null);
	   scan(#L, j->(sliceInd#(L#j) = i#j));
	   T_(toList sliceInd)
	   ));
   U := tensorFromSlices toList Tslices;
   d := #(tensorDims U);
   k := min L;
   w := toList apply(d,i->if i < k then i+1 else if i == k then 0 else i);
   U@w
   )

factorMap = method()
factorMap (Tensor,Matrix,Number) := (T,M,k) -> (
    D := tensorDims T;
    if D#k != numcols M then error "dimension mismatch";
    slices := toList apply(numrows M, i-> contract(T,makeTensor first entries M^{i},k,0));
    U := tensorFromSlices slices;
    w := toList apply(#D,i->if i < k then i+1 else if i == k then 0 else i);
    print(U,w);
    U@w
    )

initialDimensions=method()
initialDimensions List := L -> (d:={};
     while instance(L,List) do (d=d|{#L},L=L_0);
     return d)

end
