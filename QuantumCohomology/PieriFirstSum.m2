{* 
PieriFirstSum

Implement the first sum in the quantum Pieri rule 
(basically the classical formula except partitions must
stay within bounding rectangle)

Format for monomials:
{coefficient, integer for power of q, partition as a list of integers}

Inputs:
1. a monomial whose partition has only one part
2. a monomial whose partition is arbitrary 
*}

{*
Ideas for future work: 

If memory is a problem, try iterating over
the solutions instead of recursing
*}

{* 
First try: use the Polyhedra package to enumerate
the possible mu's
*}


loadPackage("Polyhedra");

pieriFirstSum = (p,l,r,L) -> (
    M1:=matrix apply(r+1, i -> apply(r+1,j -> if i==j then 1 else 0));
    v1:=matrix apply(r+1, i -> if i==0 then {l} else if i<=#L then {L_(i-1)} else {0});
    M2:=matrix apply(r+1, i -> apply(r+1,j -> if i==j then -1 else 0));
    v2:=matrix apply(r+1, i -> if i<#L then {-1*(L_i)} else {0});
    M:=M1||M2;
    v:=v1||v2;
    N:=matrix {apply(r+1, i -> 1)};
    w:=matrix {{sum(L)+p}};
    P:=intersection(M,v,N,w);
    latticePoints(P)
);



{*
Examples: 
time pieriFirstSum(4,6,4,{4,2,2,1})
-- used 0.171507 seconds
time pieriFirstSum(10,20,10,{9,8,8,6,4,4,3,2})
-- used 6.08612 seconds
*}


{*Second try: 

A recursive function for the first term in quantum Pieri
For each allowed of mu_1, compute the possible 
diagrams with the new bounding rectangle of width
lambda_1 and adding p-(mu_1-lambda_1) boxes
*}