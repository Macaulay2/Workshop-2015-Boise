V={{-5,0},{-3,0},{-1,-4},{-1,4},{-1,-2},{-1,2},{0,-1},{0,1},{1,-2},{1,2},{1,-4},{1,4},{3,0},{5,0}}
F={{0, 1, 4, 2}, {0, 1, 5, 3}, {8, 10, 13, 12}, {9, 11, 13, 12}, {1, 4, 6, 7, 5}, {2, 4, 6, 8, 10}, {3, 5, 7, 9, 11}, {6, 7, 9, 12, 8}}

E={{0, 1}, {0, 2}, {0, 3}, {1, 4}, {1, 5}, {2, 4}, {2, 10}, {3, 5}, {3, 11}, {4, 6}, {5, 7}, {6, 7}, {6, 8}, {7, 9}, {8, 10}, {8, 12}, {9, 11}, {9, 12}, {10, 13}, {11, 13}, {12, 13}}



d = INTgetSize(V);
M = transpose matrix V || matrix {apply(#V, i->1)}
G = select(select(unique(flatten(apply(#F-1,
    i->apply(toList(i+1..#F-1),
	j-> select(F_i,
	    l->member(l,F_j)))))),
    	    	k->#k > d-1), C -> rank M_C == d)



INTgetSize = method();
INTgetSize(List) := ZZ => vectors ->(
    n := #(vectors_0);
    if instance(position(vectors,v->#v != n),Nothing) then return n
    else(
	return null
    )
)
