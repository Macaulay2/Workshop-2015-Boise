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
	}
)

export {
  symmetrize
}

symmetrize = method()
symmetrize (Tensor) := (T) -> (
  error "not implemented yet";
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
    sum Tslices
    );

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
    )

end
