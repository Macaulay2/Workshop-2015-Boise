restart
loadPackage "TensorsAdditions"

--Creating Tensors--
R = QQ[a..h]
T = genericTensor(R,{2,2,2})
R = QQ
T = randomTensor(R,{3,3,3})
--makeTensor
--TensorModule

--Operations on Tensors--
symmetrize T
tensorToPolynomial(T,symbol x)
--polynomialToTensor
--etc


--Tensor Eigenvectors--
R = QQ[x]
S = R/ideal{x^3-1}
T = multiplicationTensor S
I = tensorEigenvectors(T,2,symbol y)
primaryDecomposition I

