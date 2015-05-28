----------------------------
--Part 3: Indexed Tensors
--should only depend on part 2
----------------------------
--needsPackage"Tensors"
export{IndexedTensor,indexedTensor}
IndexedTensor = new Type of HashTable
IndexedTensor.synonym="indexed tensor"
subscriptNet :=method(Dispatch=>Thing)
subscriptNet VisibleList := inds -> toString(inds_0)|concatenate(
     (take(inds,{1,#inds}))/(i->","|toString(i)))
net IndexedTensor := A -> net (hold A.cache#(gs"name"))_(subscriptNet A#(cs"indices"))
noname="[unnamed IndexedTensor]"
it=
indexedTensor=method()
it (Tensor,VisibleList) := (t,inds) -> (
     c:=new CacheTable from {(gs"name") => noname};
     new IndexedTensor from hashTable{
     	  cache => c,
       	  symbol indices => toSequence inds,
     	  symbol tensor => t}
     )

tensor IndexedTensor := opts -> t -> t.tensor
indices IndexedTensor := t -> t.indices
entries IndexedTensor := entries@@tensor

IndexedTensor.GlobalAssignHook = (sym,val) -> (
     if val.cache#(gs"name") === noname then val.cache#(gs"name") = sym;
     )
IndexedTensor#{Standard,AfterPrint} = T -> (
     << endl;				  -- double space
     t:= T.tensor;
     << concatenate(interpreterDepth:"o") << lineNumber << " : "
     << net class t
     << "-IndexedTensor with indices "
     << subscriptNet indices T
     << endl
     )

rit=
renameIndexedTensor=method()
rit (IndexedTensor,Symbol) := (t,s) -> t.cache#(gs"name") = s

-----------------------
--Indexed tensors from
--subscripting tensors
----------------------
--a.c. needs error checking
parseIndexedTensor=method()
parseIndexedTensor (Tensor,Sequence) := (T,s) -> (
     dims:=tensorDimensions T;
     if not #dims==#s then error "tensor subscripted with the wrong number of indices";
     syms:=toList select(s,isSymbolic);
     syms':=unique syms;
     --if a list of unique symbols is given...
     if #syms'==#s then return (
	  indexedTensor(T,syms'));
     p:=hashTable apply(syms,i->i=>position(syms',j->j===i));
     f := l -> apply(s,i->if instance(i,ZZ) then i 
	  else l#(p#i));
     firstsyms:=apply(syms',i->position(s,j->i===j));
     --a.c. need error check that repeated indices
     --have the same ranges...
     odims:=dims_firstsyms;
     ents:=toList apply(tensorKeys odims,i->fastTensorAccess(T,f i));
     M:=class T;
     factors:=M#(gs"factors")_firstsyms;
     M':=tensorModuleProduct factors;
     T':=tensor(M',ents);
     indexedTensor(T',syms')
     )

Tensor_Sequence := (T,s) -> (
     if allInstances(s,ZZ) then return (
	  tensorAccess(T,s));
     if not allInstances(s,
	  {ZZ,Symbol,IndexedVariable}) then 
          error "Tensor_Sequence expected a sequence 
          of integers, symbols, or indexed variables";
     parseIndexedTensor(T,s)
     )
Tensor_Symbol := (T,s) -> T_(1:s)
Tensor_IndexedVariable := (T,v) -> T_(1:v)

----------------------------------
-- Basic Indexed Tensor operations
----------------------------------
tensorDimensions IndexedTensor := t -> tensorDimensions tensor t

IndexedTensor_Sequence := (t,s) -> tensorAccess(tensor t,s)

indexPosition:=method()
indexPosition (Symbol,IndexedTensor):=(i,t) -> position(indices t,j->j===i)
dim (Symbol,IndexedTensor) := (i,t) -> (
     p:=indexPosition(i,t);
     if p===null then 0 else dim(p,tensor t)
     )

-------------------------------------------
-- Permuting the axes of indexed tensors
-------------------------------------------
IndexedTensor @ List := (t,l) -> (
     indexedTensor((tensor t)@l,(indices t)@l))

IndexedTensor @ Sequence := (t,s) -> (
    if not allInstances(s,{Symbol,IndexedVariable}) then error "
    IndexedTensor @ Sequence expected a sequence of Symbols 
    or IndexedVariables";
    inds:=toList indices t;
    inds':=sort inds;
    n:=#inds;
    perm:=toList apply(0..<n,i->find(s_i,inds'));
    assert(inds'_perm===toList s);
    indexedTensor((tensor t)@perm,s)    
    )

sort IndexedTensor := opts -> t -> t@(toSequence sort toList indices t)


---------------------------------------
--"Hadamard" products of indexed tensors
---------------------------------------
itprod=
indexedTensorProduct = method()

itprod List := its -> (
     assertInstances(its,IndexedTensor,"indexedTensorProduct(List)");
     inds := unique splice apply(its,indices);
     n := #inds;
     ranges := hashTable apply(inds,i->i=>apply(its,t->dim(i,t)));
     --check index ranges match
     for i in inds do 
       if not #unique delete(0,ranges#i) == 1 then 
       error("dimension mismatch for index "|toString(i));
     dims := apply(inds,i->max ranges#i);
     keys':=apply(tensorKeys dims,toList);
     subs := (keyList,indSeq) -> (
	  h:=hashTable toList apply(0..<n,i->inds#i=>keyList#i);
	  apply(indSeq,i->h#i)
	  );
     entry := keyList -> product apply(its,t->
	  fastTensorAccess(tensor t,subs(keyList,indices t))
	  );	  
     ents := apply(keys',entry);
     T:=makeTensor(dims,ents);
     indexedTensor(T,inds)
     )

IndexedTensor*IndexedTensor := (t,u) -> itprod{t,u}
--note that itprod is faster than * iterated by folding

--Symbolic marginalization

sum(List,IndexedTensor):=(tosum,t)->(
     inds:=toList indices t;
     n:=#inds;
     tosum':=toList select(0..<n,i->member(inds_i,tosum));
     T:=marginalize(tensor t,tosum');
     indexedTensor(T,inds-set(tosum))
     )
sum(Symbol,IndexedTensor):=(s,t)->sum({s},t)
sum(IndexedVariable,IndexedTensor):=(s,t)->sum({s},t)


-------------------------------------------------
--Einstein summation of lists of indexed tensors
-------------------------------------------------
einsum=
einsteinSum=method(Dispatch=>Thing)
einsteinSum VisibleList := l -> (
     tosum:=repeatedEntries flatten(apply(l,toList@@indices));
     sum(tosum,indexedTensorProduct l)
     )


end

restart
debug loadPackage"Tensors"

restart
debug loadPackage("Tensors",DebuggingMode=>true)

restart
uninstallPackage"Tensors"
installPackage"Tensors"
viewHelp"Tensors"