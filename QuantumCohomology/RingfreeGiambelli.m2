{* 
Ringfree Giambelli

Implement Giambelli without putting things in a ring first

*}

load "PieriSums.m2";


-- This function takes an exponent vector of a monomial to 
-- a product of generators written out that many times
-- e.g. the exponent vector of x_2^2*x_3*x_4^3 is {0,2,1,3,0} 
-- and we map this to {2,2,3,4,4,4} 
expToList = (v) -> (
    answer:={};
    j:=0;
    for i from 0 to #v-1 do (     
        answer=append(answer,apply(v_i, j -> i+1)) 
    );
    flatten answer
);

-- We compute the determinant of the Giambelli matrix 
-- in an abstract ring that the user should never see
-- and then use it to write out the Pieri multiplications that we need to do
giambelliDet = (l,yt) -> (
    S := QQ(monoid[s_1..s_l]);
    M:=matrix apply(#yt, i -> apply(#yt, j -> if yt_i+j-i > l then 0_S else if yt_i+j-i< 0 then 0_S else if yt_i+j-i == 0 then 1_S else S_(s_(yt_i+j-i))));
    L:=listForm(det(M));
    apply(#L, i -> { expToList(L_i_0),L_i_1})
);

-- This function multiplies an element g represented by the list L
-- by several horizontal strips, given in plist
iteratedPieri = (plist,r,l,L,Y) -> (
    plist = reverse plist;
    for i from 0 to #plist-1 do (
        L=quantumPieri(plist_i,r,l,L,Y)
    );    
    L
);

-- The main function.  Will export this in the package
quantumMonomialMultiplication = (r,l,yt1,yt2,Y) -> (
    W:=giambelliDet(l,yt2);
   -- print concatenate("giambelliDet(l,yt2) = ",toString(W)) << endl;
    L:=apply(#W, i -> iteratedPieri(W_i_0,r,l,{ {yt1,1}},Y));
   -- for i from 0 to #L-1 do (print concatenate(toString(W_i_1)," ",toString(L_i)) << endl);
    L=flatten apply(#L, i-> apply(#(L_i), j -> {L_i_j_0,(W_i_1)*(L_i_j_1)}));
    simplify(L,Y)
)
end

restart
break
load "RingfreeGiambelli.m2"


--Test quantumMonomialMultiplication
Y=QQ[q]

quantumMonomialMultiplication(2,2,{1},{2,1,1},Y)

quantumPieri(3,2,5,{ {1,{2,2}} },Y)
quantumMonomialMultiplication(2,5,{2,2},{3,2},Y)
-- Anders Buch's Maple program gives S[5, 4] + S[4, 3, 2] + S[4, 4, 1] + S[5, 2, 2] + S[5, 3, 1]

TEST ///
   Y=QQ[q];
   L={{{5, 4}, 1_Y}, {{5, 3, 1}, 1_Y}, {{5, 2, 2}, 1_Y}, {{4, 4, 1}, 1_Y}, {{4, 3, 2}, 1_Y}};
   assert(set(quantumMonomialMultiplication(2,5,{2,2},{3,2},Y)) === set(L) )
///

quantumMonomialMultiplication(2,5,{2,2},{3,2,1},Y)
--Anders Buch's Maple program gives S[5, 4, 1] + S[5, 3, 2] + S[4, 4, 2] + S[4, 3, 3] 

TEST ///
   Y=QQ[q];
   L={{{4, 4, 2}, 1_Y}, {{5, 4, 1}, 1_Y}, {{5, 3, 2}, 1_Y}, {{4, 3, 3}, 1_Y}};
   assert(set(quantumMonomialMultiplication(2,5,{2,2},{3,2,1},Y)) === set(L) )
///

