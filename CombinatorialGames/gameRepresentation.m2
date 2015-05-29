loadPackage("SimplicialComplexes")
cgtJoin = method()
cgtJoin(List,HashTable) := String => (L,dummyH) -> (
    if #L === 0 then (
	"") else (
    fold(concatenate,apply(#L, i-> if i=!=#L-1 then concatenate(dummyH#(L_i),",") else dummyH#(L_i)))
    )
)

gameRepresentation = method()
gameRepresentation(SimplicialComplex,List,List) := String => (Delta,L,R) -> (
    S := ring Delta;
    L = apply(L,v-> sub(v,S));
    R = apply(R,v-> sub(v,S));
    V := S_*;
    if ((#L+#R) === #V) and (join(L,R) == V) then (
	S = QQ[L,R,Degrees=>join(apply(L,i->{1,0}),apply(R,i->{0,1}))];
	d := dim(Delta);
	Fvec := apply(toList(0..d+1), i->flatten entries faces(d-i,Delta));
	H := hashTable apply(Fvec_0, f-> f=>"{|}");
	E := hashTable {{}=>""};
	--print H;
	dummyH := merge(H,E,join);
	tempStringL ={};
	tempStringR ={};
	Hnew ={};
	for F in drop(Fvec,1) do (
	    tempStringL = apply(apply(F, m-> select(keys H, k-> (degree(k)-degree(m) == {1,0}) and (k % m == 0))), i-> if i=={} then {} else i);
    	    tempStringR = apply(apply(F, m-> select(keys H, k-> (degree(k)-degree(m) == {0,1}) and (k % m == 0))), i-> if i=={} then {} else i);
    	    Hnew = hashTable apply(#F, i-> F_i => concatenate("{",cgtJoin(tempStringL_i,dummyH),"|",cgtJoin(tempStringR_i,dummyH),"}"));
    	    H = merge(H,Hnew,join);
	    --print H;
    	    dummyH = merge(dummyH,Hnew,join);
	    );
	print H;
	H#(sub(1,ring Delta))
	)
    else "Variables of Delta not bi-partitioned."
    )
end

restart

loadPackage("SimplicialComplexes")
S = QQ[x_1,x_2,x_3,y_1]
L = {x_1,x_2,x_3}
R = {y_1}
S = QQ[L,R,Degrees=>join(apply(L,i->{1,0}),apply(R,i->{0,1}))]
Delta = simplicialComplex(apply({{x_1,x_2,x_3},{x_2,x_3,y_1}},product))

gameRepresentation(Delta,L,R)