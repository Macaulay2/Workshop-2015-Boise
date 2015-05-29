-----------------------------------------------
--PART 1 of 3:
--Methods for nested lists...
--MINIMIZE dependence on this 
--section in future versions
------------------------------------------------

{*The following cartesian product lists have sequences
as their entries, rather than lists.  this is intentional, 
both for consistency with Set**Set, and for planned later 
use of nested lists of sequences of indices.
*}
cp=
cartesianProduct=method()

cartesianProduct(VisibleList,VisibleList):= (L,M) -> flatten apply(L,l->apply(M,m->(l,m)))
cartesianProduct(Sequence,Sequence) := (L,M) -> join apply(L,l->apply(M,m->(l,m)))

acp=--INTERNAL ABBREVIATION
associativeCartesianProduct = method(Dispatch=>Thing)
associativeCartesianProduct VisibleList := L -> (
     p:=apply(fold(L,cp),deepSplice);
     if not class p_0 === Sequence then p=apply(p,i->1:i);
     p
     )

--Compute the initial dimensions of a list;
--if the list is nested and rectangular
--this equals its array dimension:
--INTERNAL METHOD:
initialDimensions=method()
initialDimensions List := L -> (d:={};
     while instance(L,List) do (d=d|{#L},L=L_0);
     return d)

---Recursive function to test if a nested list is rectangular
---
isrect=--INTERNAL ABBREVIATION
isRectangular = method()
isrect(Thing) := (x) -> true
isrect(List) := (L) -> (
     if not instance(L_0,List) then return all(L,i->not instance(i,List));
     all(L,i->instance(i,List) and isrect(i) and #i==#L_0)
     )

repeatedEntries = method(Dispatch=>Thing)
repeatedEntries VisibleList := l -> (
     select(unique l,i->#positions(l,j->j===i)>1)
     )
repeatedEntries{1,2,3,1,3,4,5}


end

restart
debug loadPackage"Tensors"

restart
debug loadPackage("Tensors",DebuggingMode=>true)

restart
uninstallPackage"Tensors"
installPackage"Tensors"
viewHelp"Tensors"