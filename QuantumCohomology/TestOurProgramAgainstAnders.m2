load "OurMapleInterface-Dave.m2";
load "RingfreeGiambelli.m2";

sl3weights=(l) -> (
    v:={};
    for i from 0 to l do (
    for j from i to l do (
    for k from j to l do(
    lambda:= {k,j,i};
    v=append(v, delete(0,lambda))
    ); ););
    return drop(v,1)
);

end


restart
load "TestOurProgramAgainstAnders.m2";



--Test 1 example
Y=QQ[q];
HA=qcalcMonomialMultiplication(2,5,{2,2},{3,2},Y)
HM2=new HashTable from quantumMonomialMultiplication(2,5,{2,2},{3,2},Y)
HA === HM2

HA=qcalcMonomialMultiplication(2,2,{2,2},{2,2,2},Y)
HM2=new HashTable from quantumMonomialMultiplication(2,2,{2,2},{2,2,2},Y)
HA === HM2


HM2=new HashTable from quantumMonomialMultiplication(2,2,{1},{2,1,1},Y)

--Test all the monomials for sl3 and level 2
Y=QQ[q];
l=4;
B=sl3weights(l);
B

answer = {};
for i from 0 to #B-1 do (
    for j from 0 to #B-1 do (
        HA=qcalcMonomialMultiplication(2,l,B_i,B_j,Y);
	HM2=new HashTable from quantumMonomialMultiplication(2,l,B_i,B_j,Y);
	if HA =!= HM2 then answer = append(answer, {i,j,B_i,B_j, pairs(HA),pairs(HM2),apply(values(HA), i-> class(i)), apply(values(HM2), i -> class(i))})
    )
)

--apply(#B, i -> apply(#B, j -> {B_i,B_j, qcalcMonomialMultiplication(2,l,B_i,B_j,Y) === new HashTable from quantumMonomialMultiplication(2,l,B_i,B_j,Y)}))

