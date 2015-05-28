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

contract (Tensor,List,List) := (T,K,L) -> (
    D := tensorDims T;
    KD := apply(K, k->D#k);
    LD := apply(L, k->D#k);
    if KD != LD then error "dimension mismatch";
    Tslices := apply((#K:0)..<(toSequence KD), i-> (
	    sliceList := new MutableList from (#D:null);
	    scan(#K, j->(sliceList#(K#j) = i#j; sliceList#(L#j) = i#j));
	    print T_(toList sliceList);
	    T_(toList sliceList)
	    ));
    sum toList Tslices
    );
contract (Tensor,Number,Number) := (T,k,l) -> contract(T,{k},{l})

tensorEigenvectors = method()
tensorEigenvectors (Tensor,Number,Symbol) := (T,k,x) -> (
    R := ring T;
    d := tensorDims T;
    n := d#0;
    S := R[apply(n,i->x_i)];
    v := new MutableList from (n:0_S);
    for ind in (#d:0)..(#d:n-1) do (
	monList := toList apply(#ind, j->(if j != k then S_(ind#j) else 1_S));
	mon := product monList;
	v#(ind#k) = v#(ind#k) + sub(T_ind,S)*mon;
	);
    minors(2, matrix{toList v, gens S})
    );

tensorToPolynomial = method()
tensorToPolynomial (Tensor,Symbol) := (T,x) -> (
    R := ring T;
    d := tensorDims T;
    n := d#0;
    S := R[apply(n,i->x_i)];
    f := 0_S;
    for ind in (#d:0)..(#d:n-1) do (
	mon := product toList apply(#ind, j->S_(ind#j));
	f = f + sub(T_ind,S)*mon;
	);
    f
    );

tensorModule Tensor := T -> class T

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
