beginDocumentation()
{*     doc ///
        Key
	  TensorArray
        Headline
	  The class of all tensor arrays
        Usage
        Inputs
        Outputs
        Consequences
         Item
        Description
         Text
         Code
         Pre
         Example
        Subnodes
        Caveat
        SeeAlso
     ///
*}

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
   augmented with a list of other modules M_1...M_n and a choice of isomorphism to M'=M_1**...**M_n
   and M.  In Macaulay2, this isomoprhism is implicit in that M and M' are == as modules, and the isomorphism
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


end

1/0

restart
debug loadPackage"Tensors"

restart
debug loadPackage("Tensors",DebuggingMode=>true)

restart
uninstallPackage"Tensors"
installPackage"Tensors"
viewHelp"Tensors"
