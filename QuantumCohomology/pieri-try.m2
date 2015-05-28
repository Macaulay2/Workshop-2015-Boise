pieritest = (p, l, r, yt) -> (
    --<< "p,l,r,yt " << toString (p,l,r,yt) << endl;
    if #yt == 0 then (
	return {{p}}
    );
    ilist := reverse(toList(max(0,p-yt#0)..min(p,l-yt#0)));
    --<< "our ilist: " << toString ilist << endl;
    sublist := apply(ilist, i -> (
        M :=  pieritest(p-i,yt#0,r-1,drop(yt,1));
	--<< "M : " << toString M << endl;
	apply(M,a -> prepend(yt#0+i,a))
	--<< "x : " << toString x << endl;
	--x
    ));
    return flatten sublist
)

end

break
load "pieri-try.m2"
time pieritest(4,6,5,{4,2,2,1})
time pieritest(10,20,10,{9,8,8,6,4,4,3,2})
