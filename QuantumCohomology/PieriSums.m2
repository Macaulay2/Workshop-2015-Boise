{* 
PieriSecondSum

Implement the second sum in the quantum Pieri rule 

*}

{*
Ideas for future work: 

If memory is a problem, try iterating over
the solutions instead of recursing
*}

{* 
We're not even going to try using Polyhedra package to enumerate
the possible mu's this time.  It was too slow for the first sum.
*}




{*First try:
Strategy:
0. If there is no full column, then the second sum is zero
1. Remove a full column
2. Construct the complementary diagram cyt
3. Apply pieritest to cyt
4. Take complement of all these


*}

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
    apply(#C, i-> complementaryDiagram(r,yt#0,C_i))  
);

{*
Format for an element: list of pairs
{coefficient as an element in QQ[q], partition as a list of integers}

L= { {1+q,{2,1,1}}, {-3,{3,2}} }

*}

--let Y be the name in our package for the quantum coefficient ring, probably either ZZ[q] or QQ[q]
quantumPieriOneTerm = (p,r,l,T,Y) -> (
    P1:=pieriFirstSum(p,r,l,T_1);    
    P2:=pieriSecondSum(p,r,l,T_1);
    P1=apply(#P1, i -> {T_0,delete(0,P1_i)});
    P2=apply(#P2, i -> {(Y_0)*(T_0),delete(0,P2_i)});
    return flatten {P1,P2}
);

simplify = (L,Y) -> (
    H:=partition(last,L);
    L=apply(keys H, k -> {first sum(H#k), k});
    select(L, k -> k_0 != 0_Y)
    --delete(null, apply(#L, i -> if L_i_0 != 0 then L_i))
);


    
--Do each term, then combine them
quantumPieri = (p,r,l,L,Y) -> (
    if p==0 then return simplify(L,Y);
    simplify(flatten apply(#L, i -> quantumPieriOneTerm(p,r,l,L_i,Y)),Y)
);

end

restart
break
load "PieriSums.m2";
Y = QQ[q];
L = { {1+q,{5,1,1}}, {-3,{3,2}} };
quantumPieriOneTerm(2,2,5,L_0,Y)
quantumPieri(2,2,5,L,Y)





Z = { {1+q,{5,1,1}}, {-3,{3,2}}, {-3+q^2,{5,1,1} }};
simplify(Z)
Z = { {1+q,{5,1,1}}, {-3,{3,2}}, {-1-q,{5,1,1} }};





complementaryDiagram(2,5,{4,3,2})
complementaryDiagram(3,2,{2,2,1})

pieriSecondSum(4,3,4,{3,3,2,1})
pieriSecondSum(3,3,4,{3,3,2,1})
pieriSecondSum(2,3,4,{3,3,2,1})

--Test pieri against Anders's program
callqcalc(2,4,"S[3]*S[3,2,1]")
o2 = HashTable{{1, 1} => {q}  }
               {2} => {q}
               {4, 3, 2} => {1}
	       
pieriFirstSum(3,2,4,{3,2,1})
pieriSecondSum(3,2,4,{3,2,1})

Works! :-D


callqcalc(5,8,"S[4]*S[6,4,2,2]")
S[8, 6, 2, 2] + S[8, 5, 3, 2] + S[8, 4, 4, 2] + S[7, 6, 3, 2] + S[7, 5, 4, 2]

     + S[6, 6, 4, 2] + S[8, 5, 2, 2, 1] + S[8, 4, 3, 2, 1] + S[7, 6, 2, 2, 1]

     + S[7, 5, 3, 2, 1] + S[7, 4, 4, 2, 1] + S[6, 6, 3, 2, 1]

     + S[6, 5, 4, 2, 1] + S[8, 4, 2, 2, 2] + S[7, 5, 2, 2, 2]

     + S[7, 4, 3, 2, 2] + S[6, 6, 2, 2, 2] + S[6, 5, 3, 2, 2]

     + S[6, 4, 4, 2, 2]

pieriFirstSum(4,5,8,{6,4,2,2})
pieriSecondSum(4,5,8,{6,4,2,2})


callqcalc(5,8,"S[4]*S[6,4,2,2,1,1]")
q S[3, 2, 1] + q S[4, 1, 1] + q S[3, 1, 1, 1] + S[6, 5, 4, 2, 2, 1]

     + S[6, 6, 3, 2, 2, 1] + S[6, 6, 4, 2, 1, 1] + S[7, 4, 4, 2, 2, 1]

     + S[7, 5, 3, 2, 2, 1] + S[7, 5, 4, 2, 1, 1] + S[7, 6, 2, 2, 2, 1]

     + S[7, 6, 3, 2, 1, 1] + S[8, 4, 3, 2, 2, 1] + S[8, 4, 4, 2, 1, 1]

     + S[8, 5, 2, 2, 2, 1] + S[8, 5, 3, 2, 1, 1] + S[8, 6, 2, 2, 1, 1]
     
pieriFirstSum(4,5,8,{6,4,2,2,1,1})
pieriSecondSum(4,5,8,{6,4,2,2,1,1})     

Works!
