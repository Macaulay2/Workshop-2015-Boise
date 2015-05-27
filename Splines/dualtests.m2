restart
V = {{0,0},{1,0},{1,1},{-1,1},{-2,-1},{0,-1}}
F = {{0,2,1},{0,2,3},{0,3,4},{0,4,5},{0,1,5}}
E = {{0,1},{0,2},{0,3},{0,4},{0,5}}

installPackage("Graphs")
vi
dualV = toList(0..#F-1)
dualE = values edgeH
dualG = 

edgeH = hashTable apply(#E, e-> e=>positions(F, f-> all(E_e,v-> member(v,f))))
linkH = hashTable apply(#V, v-> v=>select(#F, f -> member(v,F_f)))

