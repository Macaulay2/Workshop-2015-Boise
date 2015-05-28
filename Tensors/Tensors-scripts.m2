restart
loadPackage "TensorsAdditions"

R = QQ[a..h]
T = genericTensor(R,{2,2,2})
R = QQ
T = randomTensor(R,{3,3,3})

I = tensorEigenvectors(T,2)
hilbertPolynomial I

contract(T,0,1)

polynomialToTensor = p -> (
    d := degree p;
    n := numgens ring p;
    makeTensor(toList(d:n),
    
    )
R = QQ
p = permutations 3
T = makeTensor({3,3,3}, ind->if member(toList ind,p) then 1_R else 0_R) 
