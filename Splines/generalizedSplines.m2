generalizedSplines = method();
--assume vertices are 0,...,n-1
generalizedSplines(List,List) := Module => (E,ideals) ->(
    S := ring first ideals;
    vertices := unique flatten E;
    n := #vertices;
    T := directSum(apply(ideals,I->coker gens I));
    M := matrix apply(E,
	e->apply(n,
	    v->if(v===first e) then 1
	    else if(v===last e) then -1
	    else 0));
   ker(map(T,S^n,sub(M,S)))
);

S = QQ[x,y,z,w];
E = {{1,2},{5,6},{4,7},{0,3},{0,4},{1,5},{2,6},{3,7},{0,1},{4,5},{6,7},{2,3}};
ideals = {ideal((x-w)^2),ideal((x-w)^2),ideal((x-w)^2),ideal((x-w)^2),
    ideal((z-w)^2),ideal((z-w)^2),ideal((z-w)^2),ideal((z-w)^2),
    ideal((y-w)^2),ideal((y-w)^2),ideal((y-w)^2),ideal((y-w)^2)}

V={{0,0,0},{0,0,1},{0,0,2},{0,1,0},{0,1,1},{0,1,2},{0,2,0},{0,2,1},{0,2,2},
    {1,0,0},{1,0,1},{1,0,2},{1,1,0},{1,1,1},{1,1,2},{1,2,0},{1,2,1},{1,2,2},
    {2,0,0},{2,0,1},{2,0,2},{2,1,0},{2,1,1},{2,1,2},{2,2,0},{2,2,1},{2,2,2}};

F={{1,2,4,5,10,11,13,14},{10,11,13,14,19,20,22,23},
    {2,3,5,6,11,12,14,15},{11,12,14,15,20,21,23,24},
    {4,5,7,8,13,14,16,17},{13,14,16,17,22,23,25,26},
    {5,6,8,9,14,15,17,18},{14,15,17,18,23,24,26,27}};
--must fix the fact that I enumerated from 1 instead of 0 because I'm stupid and careless
F = apply(F,f->apply(f,i->i-1));


adjfaces = E;
forms = {x-w,x-w,x-w,x-w,z-w,z-w,z-w,z-w,y-w,y-w,y-w,y-w};

M = generalizedSplines(E,ideals);
N = splineModule(V,F,1);
O = splineModule(adjfaces,forms,1,InputType=>"ByLinearForms");
h1 = hilbertSeries M
h2 = hilbertSeries N
h3 = hilbertSeries O
