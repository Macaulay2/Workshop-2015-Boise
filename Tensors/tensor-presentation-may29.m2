loadPackage "TensorsAdditions"
--Constructing tensors
R = QQ
T = randomTensor(R,{3,3,3})
tensorDims T
M = tensorModule(R,{3,3})
N = tensorModule(R,{2})
P = M**N
tensorDims P
P_(1,0,0) + P_(0,1,0)

R = QQ[a..h]
T = genericTensor(R,{2,2,2})

--Entries and slices
T_{0,1,1}
T_{0,0,}
T_{,0,}
T_{,1,}

--Eigenvectors
A = makeTensor{{0,1},{-2,3}}
R = QQ[x,y]    
right = tensorEigenvectors(A,0,R)
left = tensorEigenvectors(A,1,R)
primaryDecomposition left, primaryDecomposition right

T = makeTensor{{{1,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,1,0},{0,0,0}},{{0,0,0},{0,0,0},{0,0,1}}}
R = QQ[x,y,z]
I = tensorEigenvectors(T,0,R)
hilbertPolynomial I
netList primaryDecomposition I

T = makeTensor{{{1,1},{0,0}},{{0,0},{1,1}}}
R = QQ[x,y]
I = tensorEigenvectors(T,0,R)


T = makeTensor{{{1,0,0},{0,0,-1/2},{0,-1/2,0}},{{0,0,-1/2},{0,1,0},{-1/2,0,0}},{{0,-1/2,0},{-1/2,0,0},{0,0,1}}}
isSymmetric T

-- isSymmetric? 
R = QQ[x,y,z]
p = tensorToPolynomial (T,R)
polynomialToTensor p == T
S = QQ[a_0..a_2,b_0..b_2,c_0..c_2]
tensorToMultilinearForm(T,S)

I = tensorEigenvectors(T,0,R)
netList primaryDecomposition I

d = 3
n = 3
T = randomTensor(QQ,toList (d:n))
I = tensorEigenvectors(T,0,symbol x)
primaryDecomposition I
--netList (s = tensorEigenvectorsCoordinates(T,0,symbol x))
R = ring I
S = CC[toSequence entries vars R]
J = sub(I,S)
rr = (vars S | matrix{{1_S}})*transpose random(CC^1,CC^(n+1))
L = J + ideal rr
F = first entries gens L
--s = solveSystem(F)
s = bertiniZeroDimSolve(F)

ne = (n,d)->sum for i from 0 to n-1 list (d-1)^i
ne(n,d) == #s

T = QQ[a_(d:0)..a_(d:n-1)]
eigenDiscriminant(2,3,T)
