restart
loadPackage "Tensors"

R = QQ[a..h]
T = genericTensor(R,{2,2,2})
R = QQ
T = randomTensor(R,{3,3,3})
T = makeTensor{{{1_R,2},{3,4}},{{5,6},{7,8}}}
tensorDims T

polynomialVector = (T,k) -> (
    R := ring T;
    d := tensorDims T;
    n := d#0;
    S := R[x_1..x_n];
    v := new MutableList from (n:0_S);
    for ind in (#d:0)..(#d:n-1) do (
	monList := toList apply(#ind, j->(if j != k then S_(ind#j) else 1_S));
	mon := product monList;
	v#(ind#k) = v#(ind#k) + sub(T_ind,S)*mon;
	);
    toList v
    )



V = polynomialVector(T,2)
S = ring first V
M = matrix {V,gens S}
I = minors(2,M)
primaryDecomposition I
hilbertPolynomial I

polynomialToTensor = p -> (
    d := degree p;
    n := numgens ring p;
    makeTensor(toList(d:n),
    
    )
R = QQ
p = permutations 3
T = makeTensor({3,3,3}, ind->if member(toList ind,p) then 1_R else 0_R) 
