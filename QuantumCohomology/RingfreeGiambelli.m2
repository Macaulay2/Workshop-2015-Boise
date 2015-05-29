{* 
Ringfree Giambelli

Implement Giambelli without putting things in a ring first,
if this can be done??

*}

--Old stuff
complementaryDiagram = (r,l,yt) -> (
    while #yt <= r do yt = append(yt,0);    
    reverse apply(#yt, i -> l-yt_i)
);

--this function was originally called pieritest
pieriFirstSum = (p, r, l, yt) -> (
    if #yt == 0 then (
	return {{p}}
    );
    if r == 0 then (
        return {{yt#0+p}}	
    );
    b:=0;
    if #yt==r+1 then b=yt#r;
    ilist := reverse(toList(max(0,p-yt#0+b)..min(p,l-yt#0)));
    sublist := apply(ilist, i -> (
        M :=  pieriFirstSum(p-i,r-1,yt#0,drop(yt,1));
	apply(M,a -> prepend(yt#0+i,a))
    ));
    return flatten sublist
)

pieriSecondSum = (p,r,l,yt) -> (
    if (#yt < r+1) or (#yt == r+1 and yt_r == 0) then return {};    
    yt=apply(#yt, i -> yt_i -1); 
    cyt:=complementaryDiagram(r,yt#0,yt);
    C:=pieriFirstSum(l-p,r,cyt#0,cyt);
    delete(0,apply(#C, i-> complementaryDiagram(r,yt#0,C_i))) 
);


{*

giambelliDet = (l,yt) -> (
    M:=apply(#yt, i -> apply(#yt, j -> if yt_i+j-i > l then null else if yt_i+j-i< 0 then null else yt_i+j-i));
    if #M==1 then return {M_0};
    
);
*} 
end

--Want: {0,2,1,3,0} |-> {2,2,3,4,4,4}
expToList = (v) -> (
    answer:={};
    j:=0;
    for i from 0 to #v-1 do (     
        answer=append(answer,apply(v_i, j -> i+1)) 
    );
    flatten answer
);

giambelliDet = (l,yt) -> (
    S:=QQ[s_1..s_l];  
    M:=matrix apply(#yt, i -> apply(#yt, j -> if yt_i+j-i > l then 0 else if yt_i+j-i< 0 then 0 else if yt_i+j-i == 0 then 1_S else s_(yt_i+j-i)));
    L:=listForm(det(M));
    apply(#L, i -> { expToList(L_i_0),L_i_1})
);


restart
break
load "RingfreeGiambelli.m2"
giambelliMatrix(5,{3,2})
M=giambelliMatrix(5,{3,2,1,1})
M_0_3 == null











