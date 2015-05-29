{* 
Ringfree Giambelli

Implement Giambelli without putting things in a ring first,
if this can be done??

*}


load "PieriSums.m2";


{*

giambelliDet = (l,yt) -> (
    M:=apply(#yt, i -> apply(#yt, j -> if yt_i+j-i > l then null else if yt_i+j-i< 0 then null else yt_i+j-i));
    if #M==1 then return {M_0};
    
);
*} 


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
    S := QQ(monoid[s_1..s_l]);
    M:=matrix apply(#yt, i -> apply(#yt, j -> if yt_i+j-i > l then 0_S else if yt_i+j-i< 0 then 0_S else if yt_i+j-i == 0 then 1_S else S_(s_(yt_i+j-i))));
    L:=listForm(det(M));
    apply(#L, i -> { expToList(L_i_0),L_i_1})
);

--quantumPieri input order: (p,r,l,L,Y)


iteratedPieri = (plist,r,l,L,Y) -> (
    plist = reverse plist;
    for i from 0 to #plist-1 do (
        L=quantumPieri(plist_i,r,l,L,Y)
    );    
    L
);


quantumMonomialMultiplication = (r,l,yt1,yt2,Y) -> (
    W:=giambelliDet(l,yt2);
    L:=apply(#W, i -> iteratedPieri(W_i_0,r,l,{ {1,yt1}},Y));
    for i from 0 to #L-1 do (print concatenate(toString(W_i_1)," ",toString(L_i)) << endl);
    L=flatten apply(#L, i-> apply(#(L_i), j -> {(W_i_1)*(L_i_j_0),L_i_j_1}));
    simplify(L,Y)
)
end

restart
break
load "RingfreeGiambelli.m2"
giambelliDet(5,{3,2})
M=giambelliMatrix(5,{3,2,1,1})
M_0_3 == null



--Test iterated Pieri
--First just do it twice to get the right answer
Y=QQ[q]
f=quantumPieri(2,2,5,{ {1,{2,2}} },Y)
quantumPieri(3,2,5,f,Y)
iteratedPieri({2,3},2,5,{ {1,{2,2}} },Y)
oo_0


--Test quantumMonomialMultiplication
Y=QQ[q]
quantumPieri(3,2,5,{ {1,{2,2}} },Y)
quantumMonomialMultiplication(2,5,{2,2},{3,2},Y)


quantumMonomialMultiplication(2,5,{2,2},{3,2,1},Y)
--Bug: our function gives {{1, {4, 4, 2}}, {-q, {1, 1}}, {1, {5, 4, 1}}, {1, {5, 3, 2}}, {1, {4, 3, 3}}}
--but Anders Buch's Maple program gives S[5, 4, 1] + S[5, 3, 2] + S[4, 4, 2] + S[4, 3, 3] 
-- :-(





