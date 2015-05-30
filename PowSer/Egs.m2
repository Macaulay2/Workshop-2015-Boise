R=QQ[X,Y]

I1=ideal(X^2,Y)
I2=ideal(X^2,Y^2)
L:={I1,I2}
H=familyOfIdeals(L)

H^0
H^1
H^2
H^3
H^10
H^4
H^3

mingens H^3

finiteTypeFilteredAlgebra(L)
S=oo
hilbertSeries(S)

