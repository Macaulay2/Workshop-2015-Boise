-- Given an r and an l, PossibleTableaux will output all Young Tableaux that
--    will fit in an r+1 by l box
PossibleTableaux := (r,l) ->(
    MyList= toList((r+1):0)..toList((r+1):l);
    desc := L -> (
	if #L < 2 then return true;
    	if L#0 < L#1 then (
	    return false
    	    )
    	else (
	    return desc(drop(L,1))
    	    )
	)
    NewList:= for i in MyList list(
    	if desc(i) then (
	    i
    	    )
    	else (
	    "drop me"
   	    )
	);
    NewList:=delete("drop me", NewList);
    NewList:=  for i in NewList list (
    	delete(0,i)
	);
    delete({},NewList)
)