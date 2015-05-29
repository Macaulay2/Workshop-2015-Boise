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


end

restart
break
load "PieriSecondSum.m2"
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
