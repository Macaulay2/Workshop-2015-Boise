{* 
Implement the quantum Pieri rule

We do the first and second sums in the quantum Pieri rule separately

*}

{*
Ideas for future work: 

If memory is a problem, try iterating over
the solutions instead of recursing
*}

{* 
Note: we tried using the Polyhedra package to enumerate the possible
mu's in the first sum.  It was extremely slow.  We did not even try 
to do this for the second sum.  Actually we think enumerating the points 
was OK, what took so long was creating the polytopes.  
*}


{*
For the first sum in quantum Pieri: 
Our strategy is recursive
For each allowed of mu_1, compute the possible 
diagrams with the new bounding rectangle of width
lambda_1 and adding p-(mu_1-lambda_1) boxes
*}

-- This function was originally called pieritest
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


{* For the second sum in quantum Pieri, our strategy is:

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

-- corrected bug May 29; when forming C, use yt#0 as l
-- (not cyt#0 like we had previously)
pieriSecondSum = (p,r,l,yt) -> (
    if (#yt < r+1) or (#yt == r+1 and yt_r == 0) then return {};    
    yt=apply(#yt, i -> yt_i -1); 
    cyt:=complementaryDiagram(r,yt#0,yt);
    C:=pieriFirstSum(l-p,r,yt#0,cyt);
    apply(#C, i-> complementaryDiagram(r,yt#0,C_i))  
);

{*
Now put them together.

Let Y be the name in our package for the quantum coefficient ring, probably either ZZ[q] or QQ[q]

Then our format to represent an element in the quantum cohomology ring 
in a ringfree way is as a list of pairs
{partition as a list of integers,coefficient as an element in Y}

Example: (1+q)*s_{2,1,1} - 3*s_{3,2} is represented as 
L= { {{2,1,1},1+q}, {{3,2},-3} }

*}


quantumPieriOneTerm = (p,r,l,T,Y) -> (
    P1:=pieriFirstSum(p,r,l,T_0);    
    P2:=pieriSecondSum(p,r,l,T_0);
    P1=apply(#P1, i -> {delete(0,P1_i),T_1});
    P2=apply(#P2, i -> {delete(0,P2_i),(Y_0)*(T_1)});
    return flatten {P1,P2}
);

simplify = (L,Y) -> (
    H:=partition(first,L);
    L=apply(keys H, k -> { k, (1_Y)*(last sum(H#k))});
    select(L, k -> last(k) != 0_Y)
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
pieriFirstSum(0,2,2,{2,1,1})
pieriSecondSum(1,2,2,{2,1,1})

quantumPieriOneTerm(1,2,2,{{2,1,1},1},Y)
quantumPieri(1,2,2,{{{2,1,1},1}},Y)


Y = QQ[q];
L = { {{5,1,1},1+q}, {{3,2},-3} };
quantumPieriOneTerm(2,2,5,L_0,Y)
quantumPieriOneTerm(2,2,5,L_1,Y)
quantumPieri(2,2,5,L,Y)

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
