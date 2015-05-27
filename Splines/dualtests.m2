V = {{0,0},{1,0},{1,1},{-1,1},{-2,-1},{0,-1}}
F = {{0,2,1},{0,2,3},{0,3,4},{0,4,5},{0,1,5}}
E = {{0,1},{0,2},{0,3},{0,4},{0,5}}

dualV = toList(0..#F-1)
H = hashTable apply(#E, e-> e=>positions(F, f-> all(E_e,v-> member(v,f))))
dualE = values H
