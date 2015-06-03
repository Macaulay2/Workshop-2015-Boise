beginDocumentation()

------------
-- Tensor --
------------

doc ///
Key
  Tensor
  (symbol _,Tensor,Sequence)
Headline
  The class of all tensors
Description
  Text
   In Macaulay2, a tensor is stored as a vector which is 
   a member of a @TO2{"TensorModule","tensor module"}@.

  Example
    R=QQ[a..h]
    T=makeTensor({2,2,2},a..h)
    class T
    vector T

  Text
    Tensor products of tensors have the 
    appropriate dimensions.

  Example
    X=makeTensor{{a,b},{c,d}}
    Y=makeTensor{{1_R,2},{3,4}}
    X**Y

  Text
   Tensors can be manipulated similarly to vectors.

  Example
    U=makeTensor({2,2,2},{h,g,f,e,d,c,b,a})
    T+2*U
     
///

------------------
-- TensorModule --
------------------
doc ///
Key
  TensorModule
  (symbol ^,TensorModule,ZZ)
  (symbol **,TensorModule,TensorModule)
Headline
  The class of all tensor modules
Description
  Text
   A tensor module is a module that is a tensor product of smaller modules, which "remembers"
   that it is a tensor product.  Mathematically, one could define a tensor module as a module M
   augmented with a list of other modules M_1...M_n and a choice of isomorphism to M'=M_1**...**M_n.
   In Macaulay2, this isomoprhism is implicit in that M and M' are == as modules, and the isomorphism
   is accessible as inducedMap(M,M').
   
  Example
    R=QQ[x]
    M=R^3 ** R^3 ** R^4 -- doesn't remember it's a tensor product, but
    N=tensorModule(R,{3,3,4}) -- does.
    (class M,class N)
    M==N -- they are equal as modules, 
    M'===M -- but not as typed objects, and
    O=tensorModule(R,{4,3,4})
    not M==O -- tensor modules with different factors are considered different
    N.factors
    O.factors
    
///

----------------
-- makeTensor --
----------------
doc ///
Key
  makeTensor
  (makeTensor,List)
  (makeTensor,VisibleList,VisibleList)
  (makeTensor,VisibleList,Function)
Headline
  Constructor for making tensors from their entries.
Usage   
  makeTensor(rect)
  makeTensor(dims,ents)
  makeTensor(dims,func)
Inputs
  rect:List
    a rectangular nested list of entries of the tensor.
  dims:VisibleList
    dimensions of the tensor
  ents:VisibleList
    entries of the tensor
  func:Function
    specifying each entry of the tensor as 
    a function of its position key (a sequence)
Description
  Text
   A tensor can be built like a matrix, 
   by entering a rectangular nested list
   of its entries,

  Example
    makeTensor{{{1,2},{3,4}},{{5,6},{7,8}}}

  Text
   or by entering its dimensions followed by
   a non-nested visible list of entries

  Example
    makeTensor({2,2,2},{1,2,3,4,5,6,7,8})
    makeTensor({3,3,3},1..27)

  Text
   or by entering its dimenions followed by
   a function which gives the entries of 
   tensor as a function of their position keys,
   
  Example
    R=QQ[x,y,z]
    T=makeTensor({3,3,3},(i,j,k)->x^i*y^j*z^k)
    T_(1,2,2)==x^1*y^2*z^2

Caveat
  The constructor @TO"tensor"@ is used to 
  return tensors from tensor-related types 
  defined in the Tensors package.

///

------------------
-- tensorModule --
------------------
--a.c. incomplete
doc ///
Key
  tensorModule
  (tensorModule,Ring,List)
  (tensorModule,Module)
  (tensorModule,Module,List)
  (tensorModule,TensorModule)
Headline
  Constructor for making modules whose elements are tensors
Usage   
  tensorModule(R,dims)
  tensorModule(M)
  tensorModule(M,dims)
  tensorModule(N)
Inputs
  R:Ring
  dims:List
    specifying dimensions for the tensors in the module
  M:Module
    a free module or quotient module.
  N:TensorModule
Consequences
  Item
    Returns a module whose elements are interpreted 
    and displayed as tensors, which "remembers" if
    it is a tensor product of smaller modules.
  Item
    {\tt tensorModule(R,dims)} 
Description
  Example
    R=QQ[x,y]
    M=tensorModule(R,{4,3,2})
    M_(2,0,1) -- same as M_13
    M.factors
    M.dimensions
  Example
    I = ideal(x*y)
    M=R^4/I
    N=tensorModule(M,{2,2})
    N_(0,0)
    oo*x*y==0_N
///

-------------------------
-- tensorModuleProduct --
-------------------------
--doc ///
--Key

--Headline
 
--Usage
 
-- Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

----------------
-- tensorDims --
----------------
--doc ///
--Key
--  (tensorDims, Module)
--  (tensorDims, TensorModule)
--  tensorDims
--Headline
 
--Usage
 
-- Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

---------------------------------
-- genericTensor --
---------------------------------
--doc ///
--Key
--  (genericTensor, Ring, ZZ, List)
--  (genericTensor, Ring, List)
--  (genericTensor, Ring, RingElement, List)
--  genericTensor
--Headline
 
--Usage
 
--Inputs
 
--Outputs
 
--Description
 
--  Example

--///  

------------------
-- randomTensor --
------------------
--doc ///
--Key
--  (randomTensor, Ring, List)
--  (randomTensor, Ring, Thing, List)
--  randomTensor
--Headline
 
--Usage
 
--Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

--------------------------
-- indexedTensorProduct --
--------------------------
--doc ///
--Key

--Headline
 
--Usage
 
-- Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

-----------------
-- einsteinSum --
-----------------
--doc ///
--Key

--Headline
 
--Usage
 
-- Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

----------------
-- symmetrize --
----------------
--doc ///
--Key
--  (symmetrize, Tensor)
--  (symmetrize, Tensor, List)
--  symmetrize
--Headline
 
--Usage

--Inputs
 
--Outputs
 
--Description
  
--  Example

--///

-----------------
-- isSymmetric --
-----------------
--doc ///
--Key
--  (isSymmetric, Tensor)
--  isSymmetric
--Headline
 
--Usage
 
--Inputs
 
--Outputs
 
-- Description
 
--  Exmaple
  
--///

------------------------
-- tensorToPolynomial --
------------------------
--doc ///
--Key
--  (tensorToPolynomial, Tensor, Symbol)
--  (tensorToPolynomial, Tensor, Ring)
--  (tensorToPolynomial, Tensor, Ring, RingElement)
--  tensorToPolynomial
--Headline
 
--Usage
 
--Inputs
 
--Outputs
 
-- Description
 
--  Example
  
--///

-----------------------------
-- tensorToMultilinearForm --
-----------------------------
--doc ///
--Key
--  (tensorToMultilinearForm, Tensor, Ring)
--  (tensorToMultilinearForm, Tensor, Ring, RingElement)
--  (tensorToMultilinearForm, Tensor, Symbol)
--  tensorToMultilinearForm
--Headline
 
--Usage
 
--Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

------------------------
-- polynomialToTensor --
------------------------
--doc ///
--Key
--  (polynomialToTensor, RingElement)
--  polynomialToTensor
--Headline
 
--Usage
 
--Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

--------------------------
-- multiplicationTensor --
--------------------------
--doc ///
--Key

--Headline
 
--Usage
 
-- Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

--------------------------------
-- matrixMultiplicationTensor --
--------------------------------
--doc ///
--Key

--Headline
 
--Usage
 
-- Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

------------------------
-- tensorEigenvectors --
------------------------
--doc ///
--Key
--  (tensorEigenvectors, Tensor, ZZ, Ring, RingElement)
--  (tensorEigenvectors, Tensor, ZZ, Ring)
--  (tensorEigenvectors, Tensor, ZZ, Symbol)
--  tensorEigenvectors
--Headline
 
--Usage
 
--Inputs
 
--Outputs
 
--Description
  
--  Example
  
--///

-----------------------
-- eigenDiscriminant --
-----------------------
--doc///
--Key
--  (eigenDiscriminant, ZZ, ZZ, Ring)
--  eigenDiscriminant
--Headline
 
--Usage
 
-- Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

-----------------------------------
-- tensorEigenvectorsCoordinates --
-----------------------------------
--doc///
--Key
--  (tensorEigenvectorsCoordinates, Tensor, ZZ, Symbol)
--  tensorEigenvectorsCoordinates
--Headline
 
--Usage
 
-- Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

----------------------
-- tensorFromSlices --
----------------------
--doc///
--Key
--  (tensorFromSlices, List)
--  tensorFromSlices
--Headline
 
--Usage
 
-- Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

------------------
--flattenTensor --
------------------
--doc///
--Key
--  (flattenTensor, Tensor, List)
--  flattenTensor
--Headline
 
--Usage
 
-- Inputs
 
--Outputs
 
--Description
 
--  Example

--///

---------------
-- factorMap --
---------------
--doc///
--Key
--  (factorMap, Tensor, Matrix, ZZ)
--  factorMap
--Headline
 
--Usage
 
--Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

---------------------
-- permutationSign --
---------------------
--doc///
--Key
--  (permutationSign, List)
--  permutationSign
--Headline
 
--Usage
 
--Inputs
 
--Outputs
 
--Description
 
--  Example
  
--///

---------------------------------
-- associativeCartesianProduct --
---------------------------------
--doc///
--Key
--  (associativeCartesianProduct, VisibleList)
--  associativeCartesianProduct
--Headline

--Usage

--Inputs

--Outputs

--Description

--  Example

--///

end

1/0

restart
debug loadPackage"Tensors"

restart
debug loadPackage("Tensors",DebuggingMode=>true)

restart
path = unique path|{"/Users/andrew/Dropbox/Macaulay2","/Users/andrew/SVN/WFU-2012/Tensors-Package/"};
uninstallPackage"Tensors"
installPackage"Tensors"
viewHelp"Tensors"
