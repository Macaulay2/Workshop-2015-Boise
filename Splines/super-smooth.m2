R = Q-Q[x,y]

getDimension = method();
--Interior method, get dimension of collection of vertices
--f = list of vertices forming a face (e.g. {0,2,4,5})
--V = list of coordinates of all vertices
getDimension(List,List) := ZZ => (f,V)->(
   M := (transpose matrix V) || (matrix {V/(i -> 1)}); -- add row of 1s to matrix of vertices
   rank M_f-1
)

getSize = method();
getSize(List) := ZZ => L ->(
    if all(L, v-> #v == #(L_0)) then #L_0 else null
)

--SmallExample:Star of Vertex 2D--
--V = {{0,0},{1,0},{1,1},{-1,1},{-2,-1},{0,-1}};
--F = {{0,2,1},{0,2,3},{0,3,4},{0,4,5},{0,1,5}};
--E = {{0,1},{0,2},{0,3},{0,4},{0,5}};
--forms = {y,x-y,x+y,x-2*y,x};

--Schlegel Diagram of Cube--
--V={{-1,-1},{-1,1},{1,1},{1,-1},{-2,-2},{-2,2},{2,2},{2,-2}};
--F={{0,1,2,3},{0,1,4,5},{1,2,5,6},{2,3,6,7},{0,3,4,7}};
--E={{0,1},{1,2},{2,3},{0,3},{0,4},{1,5},{2,6},{3,7}};
--forms = {1+x,1-y,1-x,1+y,x-y,x+y,x-y,x+y};

--Schlegel Diagram of Octahedron
V = {{0,0},{3/2,1},{3,10},{2,2},{3/2,3},{1,2}};
F = {{0,1,2},{0,1,5},{0,4,5},{1,2,3},{1,3,5},{2,3,4},{3,4,5}};
E = {{0,1},{0,5},{1,2},{1,3},{1,5},{2,4},{3,4},{3,5},{4,5}};
forms = {2*x-3*y,2*x-y,8-6*x+y,2-2*x+y,4-2*x-y,6-2*x-y,2-y,2-y,2*x-y};

r = 0;
s = 1;
d = getSize(V);


--The list E must ONLY contain the interior edges
formIdeals := forms/ideal;
edgeIdeals := for i from 0 to #E-1 list(
    I := formIdeals_i^(r+1);
    for j from 0 to #E-1 do(
	if i == j then continue;
	intersectionIdeal := formIdeals_i;
    	int := select(E_i, elt -> member(elt,E_j));
	if getDimension(int,V) == d-2 then(
	    intersectionIdeal = (intersectionIdeal + formIdeals_j)^(s+1);
	    I = intersect(I,intersectionIdeal);
	); -- end if
    ); -- end inner for
    I
); -- end outer for


