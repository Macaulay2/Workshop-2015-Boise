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
        DebuggingMode => true 
)

export {
  symmetrize,
  tensorEigenvectors,
  tensorToPolynomial,
  multiplicationTensor
}

symmetrize = method()
symmetrize (Tensor) := (T) -> (
    d := #(tensorDims T);
    R := ring T;
    L := apply(permutations d, p->T@p);
    (1_R/(d!))*(sum L)
);

contract (Tensor,Number,Number) := (T,k,l) -> (
    d := #(tensorDims T);
    n := (tensorDims T)#k;
    m := (tensorDims T)#l;
    assert(n == m);
    Tslices := apply((0,0)..(n-1,n-1), ij-> (
	sliceList := toList apply(d, p->(if p==k then ij#0 else if p==l then ij#1 else null));
	T_sliceList
	));
    sum toList Tslices
    );

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
    minors(2, matrix{toList v, gens S})
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
    n := D#0;
    S := R[apply(n,i->x_i)];
    f := 0_S;
    for ind in (#D:0)..(#D:n-1) do (
	mon := product toList apply(#ind, j->S_(ind#j));
	f = f + sub(T_ind,S)*mon;
	);
    f
    );
tensorToPolynomial (Tensor,List) := (T,X) -> (
    R := ring T;
    D := tensorDims T;
    n := D#0;
    S := R[flatten apply(#D,i->(apply(D#i,j -> (X#i)_j)))];
    f := 0_S;
    for ind in (#D:0)..(#D:n-1) do (
	mon := product toList apply(#ind, j->S_(ind#j));
	f = f + sub(T_ind,S)*mon;
	);
    f
    );
    

tensorModule Tensor := T -> tensorModule(ring T, tensorDims T)

tensor Matrix := o -> M -> makeTensor entries M

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

end
