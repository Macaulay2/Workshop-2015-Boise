newPackage(
	"Tensors",
    	Version => "0.2.0", 
    	Date => "May 30, 2015",
    	Authors => {
	     {Name => "Andrew Critch", Email => "critch@math.berkeley.edu", HomePage => "http://www.acritch.com/"},
	     {Name => "Claudiu Raicu", Email => "claudiu@math.berkeley.edu", HomePage => "http://math.berkeley.edu/~claudiu/"},
	     {Name => "Hirotachi Abo"},
	     {Name => "Roberto Barrera"},
	     {Name => "Robert Krone"},
	     {Name => "Benjamin Reames"},
	     {Name => "Zach Teitler"}
	     },
    	Headline => "tensors",
	AuxiliaryFiles => true
    	)
 --Macaulay2-1.4/share/Macaulay2/Core/matrix1.m2 
 --needs to replaced for this package to work 

----------------------------------------
--Searchable comment legend:
--a.c. : andrew critch
----------------------------------------
export{Tensor,TensorModule,
     makeTensor,tensorModule,tensorModuleProduct,
     tensorDims,genericTensor,
     randomTensor,
     indexedTensorProduct,einsteinSum}
export{associativeCartesianProduct}

-------------------------
--Symbol methods
-------------------------

gs = getSymbol

cs=coreSymbol = method()     
cs String := nam -> (
     getGlobalSymbol(Core#"private dictionary",nam))

isSymbolic = x -> instance(x,Symbol) or instance(x,IndexedVariable)

--these are currently used for einstein summation,
--which needs to be rewritten

-----------------------
--Error methods
-----------------------
assertInstances=method()
assertInstances (List,Type) := (L,T) -> if not all(L,i->instance(i,T)) then (
     error ("expected list entries to be instances of "|toString(T)|"s"))
assertInstances (List,Type,String) := (L,T,context) -> if not all(L,i->instance(i,T)) then (
     error (context|" expected list entries to be instances of "|toString(T)|"s"))

allInstances = method()
allInstances (VisibleList,List) := (things,types) -> (
     all(things,i->not all(types,j->not instance(i,j)))
     )
allInstances (VisibleList,HashTable) := (things,type) -> (
     allInstances(things,{type})
     )

--------------------------------------------
--Load part 1 (minimize dependence on this)
--------------------------------------------
load "./tensors/cartesian-list-methods.m2"

inserts=method()
inserts(VisibleList,VisibleList,VisibleList):=(locs,things,host)->(
     if not #locs===#things then error "#locations =!= #things to insert";
     for i in 0..<#locs do host=insert(locs#i,things#i,host);
     host
     )

find = method()
find (Thing,VisibleList) := (x,l) -> (
     position(l,i->i===x))

----------------------------------
--Part 2 of 3:
--Tensors and Tensor Modules
----------------------------------

----------------
--Tensor Modules
----------------
Tensor=new Type of Vector
Tensor.synonym="tensor"
Tensor.cache = new CacheTable
vector Tensor := t -> new Vector from t

TensorModule = new Type of Module
TensorModule.cache = new CacheTable
module TensorModule := M -> M.module
module Module := identity
------
--Using dimensions method previously defined for
--RNLs now for...
--if not class tensorDims === MethodFunction then (

tensorDims =
tensorDimensions = method()
tensorDims Module := M -> {rank ambient M}
tensorDims TensorModule := M -> M#(gs"dimensions")

tensorKeys = method(Dispatch=>Thing)
tensorKeys VisibleList := l ->  toList acp(apply(l,i->0..<i)) 
tensorKeys Tensor := t -> tensorKeys tensorDims t

--Printing TensorModules:
moduleSummary=M->(
     n := rank ambient M;
     if M.?generators then
     if M.?relations then << ", subquotient of " << ambient M
     else << ", submodule of " << ambient M
     else if M.?relations then << ", quotient of " << ambient M
     else if n > 0 then (
	  if not all(degrees M, d -> all(d, zero)) 
	  then << ", degrees " << if degreeLength M === 1 then flatten degrees M else degrees M;
	  );
     )

iftm=
isFreeTensorModule = method()
iftm TensorModule := M -> (
     try M.cache#(gs"isFree") else
     M.cache#(gs"isFree")=all(M#(gs"factors"),isFreeModule)
     )

TensorModule.synonym="tensor module"

net TensorModule := M -> (
     if isFreeTensorModule M then (
	  (net module M)|
	  "{"|(fold(M#(gs"dimensions")/toString,(i,j)->i|"x"|j))|"}"
	  ) else (
	  fold(apply(M#(gs"factors"),net@@module),(i,j)->i|" ** "|j)
	  )
     )

TensorModule#{Standard,AfterPrint} = M -> (
     << endl;				  -- double space
     n := rank ambient M;
     << concatenate(interpreterDepth:"o") << lineNumber << " : "
     << (if isFreeTensorModule M then "Free " else "")
     << ring M
     << "-TensorModule of order "|toString(#M#(gs"dimensions"))|
     ", dimensions "|toString(M#(gs"dimensions"));
     moduleSummary M;
     << endl;
     )

-------------------------------------
--Building tensor modules:
-------------------------------------
tm=
tensorModule = method()

--make a free module into a tensor module:
tensorModule (Ring,List) := (R,dims) -> (
     d:=product dims;
     new TensorModule of Tensor from (
	  new HashTable from (pairs R^d)|{
      	       gs"factors" =>  apply(dims,i->R^i),
     	       gs"dimensions" =>  dims,
	       symbol module => R^d}
     	  )
     )

--make a possibly non-free module into an order 1 tensor module, 
--for tensoring with other such modules to build higher-order
--non-free tensor modules:
tensorModule Module := M -> (
     if not isQuotientModule M then error "tensorModule(Module) expected a free module or quotient module";
      new TensorModule of Tensor from (
       	   new HashTable from (pairs M)|{
		gs"factors" =>  {M},
       	   	gs"dimensions" =>  {rank ambient M},
	        symbol module => M}
	   )
     )
tensorModule TensorModule := identity

--this is conceptually weird if M is not free and #L>1
tensorModule (Module,List) := (M,dims) -> (
     d:=product dims;
     if not rank ambient M == d then error "dimensions do not multiply to the number of entries";
     if not isQuotientModule M then error "tensorModule (Module,List) expected a quotient module";
     new TensorModule of Tensor from (
	   new HashTable from (pairs M)|{
	   	gs"factors" =>  {M},
       	   	gs"dimensions" =>  dims,
	        symbol module => M})
     )

--perhaps this should instead be
--t-> (classes := ancestors class t;
--     return classes#(position(classes,i->class i===TensorModule))
--     )

fm=--[INTERNAL]
factorModules=method()
factorModules TensorModule := T -> T#(gs"factors")
factorModules Module := M -> {M}

tensorDims Tensor := t -> tensorDims class t
dim (ZZ,Tensor) := (n,T) -> (tensorDims T)#n

--Tensor module from a list of modules to tensor product,
--which themselves may be tensor modules
tmp=
tensorModuleProduct=method(Dispatch=>Thing)
tensorModuleProduct Sequence := fctrs -> tensorModuleProduct toList fctrs
tensorModuleProduct List := fctrs -> (
     assertInstances(fctrs,Module,"tensorModuleProduct(List)");
     dims:=flatten(fctrs/tensorDims);
     f:=flatten(fctrs/factorModules);
     M:=fold(fctrs/module,(i,j)->i**j);
     T:=if all(fctrs,isFreeModule) then TensorModule else TensorModule;
      new T of Tensor from (
	   new HashTable from (pairs M)|{
	   	gs"factors" => f,
       	   	gs"dimensions" => dims,
	        symbol module => M})
     )

----------------------------
--Comparing tensor modules
----------------------------
TensorModule == TensorModule := (M,N) -> (M#(gs"factors") / module)==(N#(gs"factors") / module)

----------------------------
--TensorModule operations
----------------------------
TensorModule^ZZ := (M,n) -> tensorModuleProduct (n:M)
TensorModule**TensorModule := (M,N) -> tensorModuleProduct(M,N)

--permute the factors of a tensor module:
TensorModule @ List := (M,l) -> tensorModuleProduct M#(gs"factors")_l

-----------------------------
--Basic tensor methods
-----------------------------
--Get the ambient module of a tensor
module Tensor := t -> module class t;

--Convert a tensor back into a vector
vector Tensor := t -> new (module t) from t

--Extract an entry of a tensor
--by a multi-index

--fast access without error checking
tensorAccess = method()
tensorAccess (Tensor,Sequence) := (t,s) -> (
     dims := tensorDims t;
     if not #s == #dims then error "dimension mismatch";
     if not all(0..<#s,i->s#i<dims#i) then error "index out of range";
     ind := s#0;
     for i in 1..<#s do ind = ind*(dims#i) + s#i;
     t_ind
     )

Tensor _ Sequence := tensorAccess

fta=
fastTensorAccess = method()
fta (Tensor,Sequence) := (t,s) -> (
     dims := tensorDims t;
     ind := s#0;
     for i in 1..<#s do ind = ind*(dims#i) + s#i;
     t_ind
     )

------------------------------------------
--Making tensors without RNLs (previously TensorArrays)
------------------------------------------
tensor (TensorModule,Vector) := opts -> (M,v) -> new M from v
tensor (TensorModule,VisibleList) := opts -> (M,l) -> (
     new M from map(M,(ring M)^1,for i in l list {i}))
--
makeTensor=method()
--a.c. fix "new M from" here...
makeTensor (VisibleList,VisibleList):=(dims,ents)->(
     R:=commonRing toList ents;
     M:=tensorModule(R,dims);
     tensor(M,ents)
     )
makeTensor (VisibleList,Function):=(dims,f)->(
     ents:=apply(tensorKeys dims,f);
     makeTensor(dims,ents))

Ring**Tensor := (r,t) -> error "not implemented yet"

Tensor/Function := (t,f) -> tensor(class t,apply(entries t,f))

----------------------------
--Access to basis elements
--by multi-index
----------------------------
TensorModule _ Sequence := (M,s) -> (
     dims := M#(gs"dimensions");
     if not #s == #dims then error "dimension mismatch";
     if not all(0..<#s,i->s#i<dims#i) then error "index out of range";
     ind := s#0;
     for i in 1..<#s do ind = ind*(dims#i) + s#i;
     M_ind
     )

------------------------------
--Conversions between Tensors
--and RectangularRestedLists
------------------------------

makeTensor List := L -> (
     if not isrect(L) then error "makeTensor List expected a rectangular nested list";
     dims:=initialDimensions L;
     ents:=ultimate(flatten,L);
     makeTensor(dims,ents)
     )


rnl=
rectangularNestedList=method()
rnl(List,List):=(dims,L) -> (
     if not product dims == #L then error "dimensions mismatch";
     while #dims>1 do (
	  d:=last dims;
	  L = for i in 0..<round(#L/d) list take(L,{i*d,(i+1)*d-1});
	  dims = take(dims,{0,-2+#dims}));
     L)

tensorNet = method()
tensorNet Tensor := T -> (
     dims := tensorDims T;
     if #dims < 3 then return netList rnl(dims,entries T);
     colKeys := tensorKeys(remove(dims,0));
     rowKeys := 0..<dims#0;
     colWidth := j -> j => max apply(rowKeys,i->width net T_((1:i)|j));
     colWidths := hashTable apply(colKeys,colWidth);
     padding := I -> concatenate(colWidths#(remove(I,0)) - (width net T_I):" ");
     padEntry := I -> (net T_I)|(padding I);
     netList rnl(dims,apply(tensorKeys dims,padEntry))
     )

net Tensor := memoize tensorNet;

---------------------------
--Tensor operations
---------------------------
Tensor + Tensor := (v,w) -> (
     if not class v === class w then error "Tensor+Tensor not from the same TensorModule";
     tensor(class v,(vector v)+(vector w))
     )
Tensor - Tensor := (v,w) -> (
     if not class v === class w then error "Tensor-Tensor not from the same TensorModule";
     tensor(class v,(vector v)-(vector w))
     )
RingElement * Tensor := (r,w) -> (
     if not ring r === ring w then error "RingElement*Tensor not over the same ring";
     tensor(class w,r*(vector w))
     )
Tensor * RingElement := (w,r) -> r*w
- Tensor := w -> (-1)*w
Tensor ** Tensor := (v,w) -> (
     M:=(class v)**(class w);
     tensor(M,(vector v)**(vector w))
     )
Tensor ^ ZZ := (t,n) -> fold(for i in 0..<n list t,(i,j)->i**j)

--------------------------------
--Permuting the axes of a tensor
--------------------------------
Tensor @ List := (T,l) -> (
     assertInstances(l,ZZ,"tensor(Tensor,List)");
     dims:=tensorDims T;
     if not set l === set(0..<#dims) then error "
     Tensor @ List expected a permutation of 0..<#d, where
     d is the number of dimensions of the Tensor";
     l':=inversePermutation l;
     M:=(class T)@l;
     inds:=tensorKeys dims_l;
     ents:=apply(inds,i->T_(toSequence i_l'));
     tensor(M,ents)
     )

--------------------------------------
--Turn a free tensor into a function that
--accesses its entries
--------------------------------------
assertFreeTensor=method()
assertFreeTensor Tensor := t -> (
     if not isFreeModule class t then 
     error "expected a tensor in a free tensor module"
     )

tensorFunction = method()
tensorFunction Tensor := t -> (
     f:=method(Dispatch=>Thing);
     f Sequence := s -> t_s;
     f
     )

---------------------
--Tensor slices
---------------------
--use inserts function here!
Tensor_List := (t,l) -> (
     l':=toSequence select(l,i->not i===null);
     if #l'==#l then return t_(toSequence l);
     assertFreeTensor t;
     dims:=tensorDims t;
     blanks:=positions(l,i->i===null);
     odims:=dims_blanks;
     M:=class t;
     M':=tensorModuleProduct((M#(gs"factors"))_blanks);
     keylists:=toList \ tensorKeys odims;
     ents:=toList apply(keylists,i->t_(inserts(blanks,i,l')));
     tensor(M',ents)
     )
--------
TEST///


///

-------------------
--Tensor marginals
--------------------
marg=
marginalize=method()
marg(Tensor,List) := (T,tosum) -> (
     assertFreeTensor T;
     dims:=tensorDims T;
     n:=#dims;
     if not all(tosum,i->instance(i,ZZ) and i<n) then 
      error "marginalize(Tensor,List) expected a list of integers less than the dimensions of the tensor";
     if #tosum===n then return sum entries T;
     tokeep := toList(0..<n)-set(tosum);
     keepkeys:=tensorKeys dims_tokeep;
     sumkeys:=tensorKeys dims_tosum;
     f := l -> sum apply(sumkeys,i->T_(inserts(tosum,i,l)));
     ents:=toList apply(keepkeys,f);
     M:=tensorModuleProduct((class T)#(gs"factors")_tokeep);
     tensor(M,ents)
     )


TEST ///
R=QQ[x 1]
M=tm(R,{2,2})
N=tm(R,{4})
assert(M==R^4)--they are equal as modules
assert(not M===R^4)
assert(not M==N)
assert(not M===N)
h=new MutableHashTable
h#M=1
h#N==1--unfortunately
///


--a.c. SORT THIS UPWARD:
diff(Tensor,RingElement) := (t,r) -> t/(i->diff(i,r))

---------------------
--Load part 3
---------------------
load "./tensors/gentensors.m2"
load "./tensors/indexedtensors.m2"

--

TEST  ///


///

load "./tensors/tensors-documentation.m2"

end

restart
debug loadPackage"Tensors"

restart
debug loadPackage("Tensors",DebuggingMode=>true)

