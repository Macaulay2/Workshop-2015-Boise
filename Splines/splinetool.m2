-- Given a matrix listing the vertex coordinates) and the toList of
-- the complex of top dimension (specified as tuples of vertices), the
-- program builds the chain complexes corresponding to J, R, and R/J
-- An important note: Splines are distinguished by the fact
-- that they are syzygies on a certain matrix of the form Boundary|diag
-- where the entries in Boundary have a single +1 and a single -1  in
-- each row. This corresponds to a specific orientation on the max.
-- simplices, NOT the one where the vertices are lex ordered. However,
-- it turns out that if we modify the top boundary map by diag orient(s,v)
-- this gives us the SPLINES, but we don't really need this, because
-- what we get using the lex ordering is isomorphic to the module of
-- splines. Last modified to work with M2 version 0.9.2 HKS

diag = (f) -> map(R^(#f), #f,(i,j) -> if i == j then promote(f#i,R) else 0)

orient = (simp, verts) -> toList apply(#simp, i ->  1/det verts_(simp#i))

zeromat = (m,n) -> (matrix(R^m, n, (i,j) -> 0))

signedbd = (C, r) -> (
           m := # C.ifaces#(r-1);
           n := # C.ifaces#r;
           signs(C.ifaces#r);
           map(R^m,R^n, (i,j) -> if t#?(C.ifaces#r#j, C.ifaces#(r-1)#i) then
                                      t#(C.ifaces#r#j, C.ifaces#(r-1)#i)
else 0))
-- Produce the rth boundary map for REDUCED (mod boundary) homology.
-- Entry i, j corresponds to the boundary component of the jth r-face which
-- lies in the ith r-1 face, and to see the orientation of it, we look at t,
-- which is indexed on these faces HS 10/15/95


bd = (f) -> apply(#f, i -> drop(f, {i,i}))

faces = (f) -> unique flatten apply(f, bd)

topbfaces = (faces) -> (
    facets = flatten apply(faces, bd);
    u := new MutableHashTable;
    scan(facets, g -> if u#?g then u#g = false else u#g = true);
    select(keys u, f -> u#f == true))

signs = (f) -> (
        t = new MutableHashTable;
        u := new MutableHashTable;
        scan(0..# f - 1, i -> (
        u#i = bd(f#i);
        scan(0..# u#i - 1, j -> t#(f#i,u#i#j) = (-1)^j);)))
-- create an indexed object to reflect boundary orientation. i.e.
-- if f = {{0 1 2}{2 3 4}} then t#(f#0,u#0#1) = t#({0 1 2},{0 2})
-- should be -1 HS 10/15/95

------------------------------------------------------------------------------
--Spline portion uses parts from outside - signs, bd, faces, topbfaces.
--This just builds the matrix consisting of the top bd map and the diag
--matrix of hyperplanes, then computes the kernel and projects. For this,
--we don't need any of the fancy simplicial complex stuff, so ignore it.


spline = (sim, ver, r) ->(
    n := (#sim#0) - 1;
--    R = ZZ/31991[x_0..x_n];
R = QQ[x_0..x_n];
    submatrix(gens kernel(topbd(sim)*diag(orient1(sim,ver)) |
e1(sim,ver,r)), {0..#sim-1},))
--returns a matrix whose columns are the splines

topbd = (topfaces) -> (
          nfaces = toList apply(topfaces, i -> sort i);
          nminus1faces = faces nfaces;
          nminus1bfaces = topbfaces nfaces;
          nminus1ifaces = toList ((set  nminus1faces) - (set nminus1bfaces));
          m := # nminus1ifaces;
          n := # nfaces;
          signs(nfaces);
          map(R^m,R^n, (i,j) -> if t#?(nfaces#j, nminus1ifaces#i) then
                                   t#(nfaces#j, nminus1ifaces#i) else 0))
-- Produce the top boundary map for REDUCED (mod boundary) homology.
-- Remark - if we really want the SPLINES, then we have to modify the
-- topbd map by a matrix reflecting the real orientation of the top
-- simplices. In other words, we will multiply topbd by diag orient1
-- Modified HS 12/30/96

orient1 = (simp, verts) -> toList apply(#simp, i ->
                                    if (det verts_(simp#i) < 0) then -1
                                    else 1)

faceideal1 = (f, ver) ->
             minors(1+#f, (transpose vars R) | ver_f**R)
-- find the linear form vanishing on a top - 1 dimensional face

e1 = (sim, ver, r) -> (
    m := #nminus1ifaces;
    map(R^m, R^{m:(-r-1)}, (i,j) ->
       if i =!= j then 0
       else (generators(faceideal1(nminus1ifaces#i,ver)))_(0,0)^(r+1)))
-- build the diagonal matrix of r+1 powers of the linear forms
-- modified to work with version of M2 of 11/7/96 - changes in last else above.


-----------------------------------------------------------------------------
--Homology portion, wherein we build all the chain complexes.

CCR = (C) -> (
      D := new ChainComplex;
      D.ring = R;
      scan(-1..C.dim+1, i -> D#i = R^(# C.ifaces#i));
      scan(0..C.dim+1, i -> D.dd#i = signedbd(C,i));
      D)
-- Takes a simplicial complex C (given by it's top dimensional part), produces
-- the chain complex corresponding to the sheaf R.  HS 10/15/95


Splitup = (m) -> (
    n = transpose m;
    I = id_(target m);
    p = rank target m;
    degs = apply(p, i -> i*(p+1));
    transpose (I ** n)_degs)
-- Takes a matrix and makes each row a row with all non-row entries 0 in the
-- new object.


CCJ = (C, vertices, r) -> (
      J = new ChainComplex;
      J.ring = R;
      K = new ChainComplex;
      a = new ChainComplexMap;
      K = CCR(C);
      temp = eqns(C,vertices,r);
      --J#(C.dim-1) = coker presentation image temp;
      J#(C.dim-1) = image temp;
      a#(C.dim-1) = map(K_(C.dim-1), J_(C.dim-1), temp);
      scan(reverse(0..C.dim-2), i -> (
                     temp = mingens image Splitup(K.dd_(i+1)*temp); --changed
                     J#i = coker presentation image temp;
                     --J_i is obtained by doing the boundary map, and then
                     --"expanding" so that we get, on each i face, the sum
                     --of all the hyperplanes thru it
                     a#i = map(K_i,J_i, temp);));
      J#(C.dim+1) = J#(C.dim) = J#-1 = R^0;
      a#(-1) = map(K_-1, J_-1, 0);
      a#(C.dim) = map(K_(C.dim), J_(C.dim), 0);
      a#(C.dim+1) = map(K_(C.dim+1), J_(C.dim+1), 0);
      scan(reverse(0..C.dim), i -> (
        J.dd#i = map(J_(i-1), J_i, (K.dd_i*a_i) // a_(i-1));));  --changed
      J)
-- Takes a simplicial complex C, a list of the vertex locations (J depends on
-- embedding, R does not), and the order of differentiability desired, and
-- returns the chain complex J. Also creates a ChainComplexMap a: J -> R


qComplex = (U,tu) -> (
      K := new ChainComplex;
      b = new ChainComplexMap;
      scan(drop(sort keys U, 2), i -> (K#i = coker matrix tu_i));
      scan(0..max keys U, i -> (K.dd#i = map(K_(i-1), K_i, U.dd_i)));
      K.ring = R;
      K)
--Given a complex U, and a map of complexes tu, forms the cokernel complex

Init = (simp, vertices, r) -> (
       D = new SimplicialComplex from simp;
       B = CCR(D);
       A = CCJ(D, vertices, r);
       C = qComplex(B,a);)

SimplicialComplex <- new Type of MutableHashTable
SimplicialComplex.name = "SimplicialComplex"
new SimplicialComplex from List := (s,x) -> (
    C := new MutableHashTable;
    n := C.dim = (# x#0) - 1;
--    R = ZZ/31991[vars(0..n)];
R = QQ[vars(0..n)];
    C.faces = new MutableHashTable;
    C.bfaces = new MutableHashTable;
    C.ifaces = new MutableHashTable;
    C.faces#n = toList apply(x, i -> sort i);
    C.ifaces#n = C.faces#n;
    D := C.faces#(n-1) = faces C.faces#n;
    E := C.bfaces#(n-1) = topbfaces C.faces#n;
    C.ifaces#(n-1) = toList ((set D) - (set E));
    scan(reverse (0..n-2), i -> (
        D = C.faces#i = faces C.faces#(i+1);
        E = C.bfaces#i = faces C.bfaces#(i+1);
        C.ifaces#i = toList ((set D) - (set E));
        ));
    C.ifaces#-1 = C.faces#-1 = C.bfaces#-1 = {};
    C.ifaces#(n+1) = C.faces#(n+1) = C.bfaces#(n+1) = {};
    C)
-- 0 faces are vertices, and the -1, n+1 face is the 0 module (added in to
-- compute homology) Also, by sorting the top boundary faces, we avoid
-- any ambiguity  HS 10/15/95

-- faceideal = (C, f, vertices) ->
--    minors(1 + #f, (transpose vars R) | promote(vertices_f,R))
-- error 3/26/97 promote has been replaced by tensor

faceideal = (C, f, vertices) ->
    minors(1 + #f, (transpose vars R) | vertices_f**R)


--eqns = (C, vertices, r) -> (
--    n := C.dim;
--    m := #C.ifaces#(n-1);
--    map(R^m, R^{m:(-r-1)}, (i,j) ->
--       if i =!= j then 0
--       else ((faceideal(C,C.ifaces#(n-1)#i,vertices))_(0,0))^(r+1)))

eqns = (C, vertices, r) -> (
    n := C.dim;
    m := #C.ifaces#(n-1);
    map(R^m, R^{m:(-r-1)}, (i,j) ->
       if i =!= j then 0
       else
((generators(faceideal(C,C.ifaces#(n-1)#i,vertices)))_(0,0))^(r+1)))
--replaced old eqns (above) 11/20/96, due to change in M2. HKS

------------------------------------------------------------------------------
--Examples

line = {{0,1},{1,2},{2,3}}

lverts1 = transpose matrix(
{{0,1},{2,1},{3,1},{4,1}})

lverts2 = transpose matrix(
{{1,1},{3,1},{5,1},{7,1}})

lverts3 = transpose matrix(
{{1,1},{2,1},{4,1},{8,1}})

lverts4 = transpose matrix(
{{1,1},{2,1},{4,1},{7,1}})

line2= {{0,1},{1,2},{2,3},{3,4}}

l2verts = transpose matrix(
{{0,1},{2,1},{3,1},{4,1},{5,1}})

triv1= {{0,1,2},{0,2,3},{0,1,3}}

trivverts1 = transpose matrix(
{{0,0,1},
 {0,1,1},
 {1,0,1},
 {-1,-1,1}})

triv= {{0,1,2},{0,2,3},{0,3,4},{0,1,4}}

trivverts = transpose matrix(
{{0,0,1},
 {0,1,1},
 {1,0,1},
 {0,-1,1},
  {-1,0,1}})

tri1 = {{0,1,2},{0,1,4},{0,2,5},{1,2,3},{2,3,5},{1,3,4},{3,4,5}}

t1 =  transpose matrix(
{{0,4,1},
 {2,-2,1},
 {0,2,1},
 {-2,-2,1},
 {6,-4,1},
 {-6,-4,1}})
--t1 has a hidden symmetry!!!!!!

t1a =  transpose matrix(
{{0,4,1},
 {2,-2,1},
 {0,2,1},
 {-2,-2,1},
 {6,-5,1},
 {-6,-4,1}})

t2 =  transpose matrix(
{{0,5,1},
 {1,0,1},
 {-1,0,1},
 {0,-2,1},
 {6,-4,1},
 {-6,-4,1}})

t3 =  transpose matrix(
{{0,5,1},
 {1,0,1},
 {-1,0,1},
 {0,-3,1},
 {6,-4,1},
 {-6,-3,1}})

s1 = {{0,1,3},{1,3,4},{1,2,4},{2,4,7},{4,6,7},{3,4,6},{3,5,6},{0,3,5}}

sverts=transpose matrix(
{{0,2,1},
 {2,2,1},
 {4,2,1},
 {1,1,1},
 {3,1,1},
 {0,0,1},
 {2,0,1},
 {4,0,1}})

tri2 =
{{0,1,2},{0,1,4},{0,2,5},{2,3,5},{1,3,4},{3,5,6},{3,6,4},{4,6,5},{1,2,3}}

tri2g1 = {{0,1,2},{0,1,4},{0,2,5},{2,3,5},{1,3,4},{3,5,6},{3,6,4},{4,6,5}}

tri2vertices =  transpose matrix(
{{0,8,1},
 {2,-2,1},
 {0,2,1},
 {-2,-2,1},
 {8,-6,1},
 {-8,-6,1},
 {-1,-4,1}})

symmg1={{0,3,4},{0,1,4},{1,4,5},{1,5,6},{1,2,6},{2,6,7},{2,7,8},{0,2,8},{0,3,8}}

s2={{0,1,2},{0,3,4},{0,1,4},{1,4,5},{1,5,6},{1,2,6},{2,6,7},{2,7,8},{0,2,8},{0,3,8}}

s2verts = transpose matrix(
{{0,1,1},
 {1,-1,1},
 {-1,-1,1},
 {0,5,1},
 {2,1,1},
 {4,-3,1},
 {0,-3,1},
 {-4,-3,1},
 {-2,1,1}})


tetra =
{{0,2,5,6},{0,3,6,7},{0,7,5,1},{4,5,6,7},{0,1,2,3},{2,3,4,6},{1,2,4,5},{1,4,3,7}
,{2,4,5,6},{1,5,4,7},{3,4,6,7},{0,2,1,5},{0,3,2,6},{0,1,3,7},{1,2,3,4}}

tetrag1 =
{{0,2,5,6},{0,3,6,7},{0,7,5,1},{4,5,6,7},{2,3,4,6},{1,2,4,5},{1,4,3,7},{2,4,5,6}
,{1,5,4,7},{3,4,6,7},{0,2,1,5},{0,3,2,6},{0,1,3,7},{0,1,2,3}}

tetravertices = transpose matrix(
{{0,0,3,1},
 {1,1,1,1},
 {0,-2,1,1},
 {-1,1,1,1},
 {0,0,0,1},
 {6,-4,-1,1},
 {-6,-4,-1,1},
 {0,5,-1,1}})

tet1 = transpose matrix(
{{0,1,10,1},
 {1,1,1,1},
 {0,-2,1,1},
 {-1,1,1,1},
 {0,0,0,1},
 {7,-4,-1,1},
 {-6,-4,-1,1},
 {0,5,-1,1}})

tet2 = transpose matrix(
{{0,1,10,1},
 {1,1,1,1},
 {0,-2,1,1},
 {0,1,1,1},
 {0,0,0,1},
 {7,-4,-1,1},
 {-6,-4,-1,1},
 {0,5,-1,1}})



tetravertices1= transpose matrix(
{{0,0,5,1},
 {1,-1,1,1},
 {-1,-1,1,1},
 {0,1,1,1},
 {0,0,0,1},
 {6,-4,-1,1},
 {-6,-4,-1,1},
 {0,5,-1,1}})

tetravertices2= transpose matrix(
{{0,0,6,1},
 {1,1,1,1},
 {0,-2,1,1},
 {-1,1,1,1},
 {0,0,0,1},
 {6,-6,-1,1},
 {-6,-6,-1,1},
 {0,6,-1,1}})

tetravertices3= transpose matrix(
{{0,0,6,1},
 {1,1,1,1},
 {0,-1,1,1},
 {-1,1,1,1},
 {0,0,-1,1},
 {6,-6,-6,1},
 {-6,-6,-6,1},
 {0,6,-6,1}})

tetravertices4= transpose matrix(
{{1,0,5,1},
 {1,-1,1,1},
 {-1,-1,1,1},
 {0,1,1,1},
 {0,0,0,1},
 {16,-16,-5,1},
 {-12,-10,-4,1},
 {0,5,-1,1}})

oct =
{{0,1,2,5},{0,1,2,6},{0,1,4,5},{0,1,4,6},{0,2,3,5},{0,2,3,6},{0,3,4,5},{0,3,4,6}
}

oct1vertices = transpose matrix(
{{0,0,0,1},
 {-4,3,-1,1},
 {3,3,2,1},
 {3,-3,0,1},
 {-2,-3,1,1},
 {-2,1,10,1},
 {0,0,-10,1}})

oct2vertices = transpose matrix(
{{0,0,0,1},
 {1,0,0,1},
 {0,1,0,1},
 {-1,0,0,1},
 {0,-1,0,1},
 {0,0,1,1},
 {0,0,-1,1}})

oct3vertices = transpose matrix(
{{1,1,1,1},
 {4,0,0,1},
 {0,4,0,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,4,1},
 {0,0,-4,1}})

oct4vertices = transpose matrix(
{{1,1,1,1},
 {4,0,0,1},
 {0,4,0,1},
 {-4,0,0,1},
 {0,-4,1,1},
 {0,0,4,1},
 {0,0,-4,1}})

ogon =
{{0,1,2,3},{0,1,2,6},{0,1,3,5},{0,1,5,6},{0,2,3,4},{0,2,4,6},{0,3,4,5},{0,4,5,6}
}

o3vertices = transpose matrix(
{{0,0,0,1},
 {4,0,0,1},
 {0,4,0,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o4vertices = transpose matrix(
{{0,0,0,1},
 {3,0,1,1},
 {0,4,0,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o5vertices = transpose matrix(
{{0,0,0,1},
 {2,1,1,1},
 {0,4,0,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o6vertices = transpose matrix(
{{1,0,0,1},
 {4,0,0,1},
 {0,4,0,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o6avertices = transpose matrix(
{{0,0,0,1},
 {3,0,1,1},
 {0,3,-1,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o7vertices = transpose matrix(
{{0,0,0,1},
 {2,1,1,1},
 {0,3,1,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o8vertices = transpose matrix(
{{0,0,0,1},
 {2,1,1,1},
 {1,2,1,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o9vertices = transpose matrix(
{{1,0,1,1},
 {4,0,0,1},
 {0,4,0,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o10vertices = transpose matrix(
{{1,0,1,1},
 {3,0,1,1},
 {0,4,0,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o11vertices = transpose matrix(
{{1,0,1,1},
 {2,1,1,1},
 {0,4,0,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o12vertices = transpose matrix(
{{1,1,1,1},
 {4,0,0,1},
 {0,4,0,1},
 {0,0,4,1},
 {-4,0,0,1},
 {0,-4,0,1},
 {0,0,-4,1}})

o13vertices = transpose matrix(
{{0,0,0,1},
 {4,0,1,1},
 {0,4,1,1},
 {0,1,4,1},
 {-4,0,-1,1},
 {0,-4,-1,1},
 {0,-1,-4,1}})

o14vertices = transpose matrix(
{{0,1,0,1},
 {4,0,1,1},
 {0,4,1,1},
 {0,1,4,1},
 {-4,0,-1,1},
 {0,-4,-1,1},
 {0,-1,-4,1}})

isbundle = (cpx, verts,r)->(
                                         Init(cpx,verts, r);
                                         t := HH_3 C;
                                         i:=3;
                                         while i>0 do (
                                         print i;
                                         print dim Ext^i(t,R);
                                         i=i-1)
)

isbundle1 = (cpx, verts,r)->(
                                         Init(cpx,verts, r);
                                         print "dim H_2";
                                         print dim HH_2 C;
                                         print "dim H_1";
                                         print dim HH_1 C)

bpyr = {{0,1,2,3},{0,1,3,4},{0,2,3,4},{0,1,2,5},{0,1,4,5},{0,2,4,5}}

bverts = transpose matrix(
{{0,0,0,1},
 {0,3,0,1},
 {3,0,0,1},
 {0,0,3,1},
 {-2,-1,1,1},
  {1,2,-3,1}})  --free for r=1, nonfree for r=2 (and get all 9 hypers)

bverts0 = transpose matrix(
{{0,0,0,1},
 {0,3,0,1},
 {3,0,0,1},
 {0,0,3,1},
 {-2,-2,1,1},
  {1,0,-3,1}})  --free for  r <5, Nonfree r=5!, free r =6, nonfree r=7 

bverts1 = transpose matrix(
{{0,0,0,1},
 {0,3,0,1},
 {3,0,0,1},
 {0,0,3,1},
 {-2,-2,1,1},
  {0,0,-3,1}})  --free for all r

bverts2 = transpose matrix(
{{0,0,0,1},
 {0,3,0,1},
 {3,0,0,1},
 {0,0,3,1},
 {-2,-1,1,1},
  {1,0,-3,1}})  --free for  r <5, Nonfree r=5!, free r =6, nonfree r=7 

bverts3 = transpose matrix(
{{0,0,0,1},
 {0,3,0,1},
 {3,0,0,1},
 {0,0,3,1},
 {-2,-1,1,1},
  {0,1,-3,1}})  --free for  r <5, Nonfree r=5!, free r =6, nonfree r=7 

simps2 = {{0,1,2,4},{0,2,3,4},{2,3,4,5}}

verts2 = transpose matrix(
{{0,0,0,1},
 {3,0,0,1},
 {0,3,0,1},
 {0,0,3,1},
 {3,0,3,1},
 {0,3,3,1}}) --free for all r (dual graph a segment)

mtet= {{0,1,2,3},{0,1,2,4},{0,2,3,4},{0,1,3,4}}

vertsm = transpose matrix(
{{0,0,1,1},
 {-1,-1,0,1},
 {1,0,0,1},
 {0,1,0,1},
 {0,0,5,1}}) --free for all r (dual graph a segment)


