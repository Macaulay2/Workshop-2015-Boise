--Package for computing topological boundary maps and piecewise continuous splines on polyhedral complexes  --
--Making a small change

newPackage("Splines",DebuggingMode => true)

export {analyze,splines,bMap,bComp,relbMap,ref,splineComplex,splineComplexS,pComp,Graphics,IntGraphics,bGraphics,saveComplex,h0pres,h0pres0,splineSyz,bsplines,edgeAnn,combSplines,vertSplines,cvSplines,splineQuotient}

needsPackage "Polyhedra"

--Input: A list of vertices and a list of maximal faces.--
--Output: A polyhedral complex with those vertices and consisting of the polyhedra prescribed by the maximal face list--

pComp=method();
pComp(List, List):=(V,P)->(
	PList=apply(P,x->(convexHull(transpose matrix V_x)));
	polyhedralComplex PList)


--Input: A 2 or 3 dimensional Polyhedral Complex PC
--Output: A pair of lists (V,F) where V is a list of the vertices of PC and F is a list of codim 1 faces of P.  Can use these in mathematica to draw the (interior of the) polyhedral complex (first function) and the boundary of the polyhedral complex (second function)--

Graphics=method();
Graphics(PolyhedralComplex,String):=(PC,S)->(
	d:=dim PC;
	V:=vertices PC;
	V0:=apply(V,v->(flatten entries v));
	Facets:=polyhedra(d,PC);
	Faces:=polyhedra(d-1,PC);
	FV0:=apply(Facets,F->(vertices F) );
	FV1:=apply(FV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) )+ 1 ))));
	fV0:=apply(Faces,f->(vertices f) );
	fV1:=apply(fV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) )+ 1 ))));
	S<<"vert="<<toString(V0)<<endl<<"facets="<<toString(FV1)<<endl<<"faces="<<toString(fV1)<<endl<<close )

IntGraphics=method();
IntGraphics(PolyhedralComplex,String):=(PC,S)->(
	d:=dim PC;	
	V:=vertices PC;
	V0:=apply(V,v->(flatten entries v));
	Facets:=polyhedra(d,PC);
	Faces:=polyhedra(d-1,PC);
	IntFaces:=select(Faces,f->(
		vf:=apply(numcols vertices f,x->(vertices f)_x);
		L:=apply(Facets,F->(
			vF:=apply(numcols vertices F,y->(vertices F)_y);
			isSubset(set(vf),set(vF)) ));
		#select(L,j->(j==true))>1));
	IFPV0:=apply(IntFaces,p->vertices p);
	IFPV1:=apply(IFPV0,M->(apply(numcols M,i->( position(V,v->( (vector flatten entries transpose v)==M_i) ) +1 ))));
	S<<"vert="<<toString(V0)<<endl<<"faces="<<toString(IFPV1)<<endl<<close 
)

bGraphics=method();
bGraphics(PolyhedralComplex,String):=(PC,S)->(
	BC:=bComp(PC);
	d:=dim BC;
	V:=vertices BC;
	V0:=apply(V,v->(flatten entries v));
	Facets:=polyhedra(d,BC);
	FV0:=apply(Facets,F->(vertices F) );
	FV1:=apply(FV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) )+ 1 ))));
	S<<"vert="<<toString(V0)<<endl<<"faces="<<toString(FV1)<<endl<<close )

--shuffles the rows of halfspaces matrix so that they appear in same order as codim 1 faces of polyhedron --

halfsort=method(); 
halfsort(Sequence,List):=(x,y)->(
	r=numrows(x#0); 
	L=new MutableList from toList(1..r); 
	for i from 0 to r-1 do(
		j:=0; 
		V:={(x#1)_(i,0)}; 
		U:=toList( set(flatten entries ((x#0)^{i} * vertices(y_j))) ); 
		while U !=V do(
			j=j+1; 
			U=toList(set(flatten entries ((x#0)^{i} * vertices(y_j))) )  ); 
			L#j=flatten entries (x#0)^{i} ); 
		matrix toList L)

--Inputs: A pure, connected, hereditary polyhedral complex PC of dimension d in R^d (pseudomanifold?). (if a pair of facets X,Y of PC intersect nontrivially then there must be a chain of facets (F0,F1,..,Fn) of PC such that F0=X,Fn=Y, and each consecutive pair of facets Fi and F(i+1) intersect in a codim 1 face of both).
--Outputs: The boundary complex of PC.

bComp=method();
bComp(PolyhedralComplex):=PC->(
	d=dim PC;
	Facets=polyhedra(d,PC);
	Faces=polyhedra(d-1,PC);
	polyhedralComplex(
		select(Faces,f->(
			vf=apply(numcols vertices f,x->(vertices f)_x);
			L=apply(Facets,F->(
				vF=apply(numcols vertices F,y->(vertices F)_y);
				isSubset(set(vf),set(vF)) ));
			#select(L,j->(j==true))==1)) ))


--Inputs: A list of vertices, a list of top dim facets, a list of codim 1 faces, or a Polyhedral Complex
--Outputs: The matrix whose kernel is the module of continuous splines on the polyhedral complex.

splines=method();
splines(ZZ,List):=(r,L)->(
	V:=L_0;
	F:=L_1;
	f:=L_2;	
	d:=length(V_0);
	originPos=position(V,v->(v=={0,0,0}));
	if originPos===null then S:=QQ[t_0..t_d] else test:= all(F,f->member(originPos,f));
	if test=true then S:=QQ[t_1..t_d] else S:=QQ[t_0..t_d];
	BM0=for i from 0 to length(f)-1 list(
		Z:={};
		j:=0;
		h:=1;
		while j<length(F) do(
			if isSubset(set(f_i),set(F_j)) 
			then (Z=append(Z,h); 
				h=-1) 
			else (Z=append(Z,0));
			j=j+1);
		Z );
	Ipos=select(toList(0..#BM0-1),i->(
		#(select(BM0_i,j->(j!=0)) )==2));
	fint=f_(Ipos);
	BM=BM0_(Ipos);
	if #(flatten vars S)==d then M=transpose(matrix(S,V)) else M=transpose(matrix(S,V))||matrix(S,{apply(V,t->1)});
	DGList=for i from 0 to length(fint)-1 list(
		gens gb minors(numrows M,(transpose(vars S))|(M_(fint_i)) ) );
	if any(DGList,x->(x== matrix{{1}} **S )) then return "some vertices on entered face not in codim 1 face";
	T=diagonalMatrix apply(DGList,x->((x_(0,0))^(r+1)));
	(matrix BM)**S|T
)

splines(ZZ,PolyhedralComplex):=(r,PC)->(
	d:=dim PC;
	V:=vertices PC;
	V0:=entries(transpose V);
	Facets:=polyhedra(d,PC);
	Faces:=polyhedra(d-1,PC);
	FV0:=apply(Facets,F->(vertices F) );
	FV1:=apply(FV0,M->(apply(numcols M,i->(position(V0,v->( (vector v)==M_i ) ) ))));
	fV0:=apply(Faces,f->(vertices f) );
	fV1:=apply(fV0,M->(apply(numcols M,i->(position(V0,v->( (vector v)==M_i ) ) ))));
	splines(r,{V0,FV1,fV1}) )

analyze=method();
analyze(ZZ,List,String):=(r,L,S)->(
	sp=apply(L,l->splines(r,toSequence l));
	S<<toString(sp)<<endl<<close )

--Inputs: A list of vertices, a list of top dim facets, a list of codim 1 faces, or a Polyhedral Complex.

--Outputs: The matrix whose kernel is the module of continuous splines on the polyhedral complex, with boundary condition 0.

bsplines=method();
bsplines(ZZ,List):=(r,L)->(
	V:=L_0;
	F:=L_1;
	if length(L)<4 then ve={{}} else ve=L_3;
	f:=L_2-set(ve);
	d:=length(V_0);
	S:=QQ[t_0..t_d];
	BM:=for i from 0 to length(f)-1 list(
		Z={};
		j=0;
		h=1;  
		while j<length(F) do(
			if isSubset(set(f_i),set(F_j)) 
			then (Z=append(Z,h); 
				h=-1) 
			else (Z=append(Z,0));
			j=j+1);
		Z );
	M:=transpose(matrix(S,V))||matrix(S,{apply(V,t->1)});
	DGList:=apply(f,e->gens gb minors(d+1,(transpose(vars S))|(M_e)) );
	if any(DGList,x->(x== matrix{{1}} **S )) then return "some vertices on entered face not in codim 1 face";
	T:=diagonalMatrix apply(DGList,x->((x_(0,0))^(r+1)));
	(matrix BM)**S|T
)

bsplines(ZZ,PolyhedralComplex):=(r,PC)->(
	d:=dim PC;
	V:=vertices PC;
	V0:=entries transpose V;
	Facets:=polyhedra(d,PC);
	Faces:=polyhedra(d-1,PC);
	FV0:=apply(Facets,F->(vertices F) );
	FV1:=apply(FV0,M->(apply(numcols M,i->(position(V0,v->( (vector flatten entries transpose v)==M_i ) ) ))));
	fV0:=apply(Faces,f->(vertices f) );
	fV1:=apply(fV0,M->(apply(numcols M,i->(position(V0,v->( (vector flatten entries transpose v)==M_i ) ) ))));
	bsplines(r,{V0,FV1,fV1}) )

--Inputs: A list of vertices, a list of top dim facets, a list of codim 1 faces, or a Polyhedral Complex in two dimensions (with no cut edge, i.e. an edge both of whose vertices are on the boundary)
--Output: A matrix whose kernel is the nontrivial splines (shifted by 1 degree) in the form of simultaneous syzygies about vertices

splineSyz=method();
splineSyz(ZZ,List):=(r,L)->(
	V=L_0;
	F=L_1;
	f=L_2;	
	d=length(V_0);
	S=QQ[t_0..t_d];
	fint=select(f,i->(
		#(select(F,j->(isSubset(i,j))))==2));
	fb=f-set(fint);
	vb=set(flatten fb);
	vint=toList(0..length(V)-1)-vb;
	M=transpose(matrix(S,V))||matrix(S,{apply(V,t->1)});
	DGList=apply(fint,e->(gens gb minors(d+1,(transpose(vars S))|(M_e) ) )_(0,0));
	if any(DGList,x->(x==sub(1,S))) then return "Error: some vertices on entered face not in codim 1 face";	
	matrix apply(vint,i->(
		apply(length(fint),e->(
			if member(i,fint_e) 
			then (
				p:=position(fint_e,x->(x==i));
				return (-1)^p*(DGList_e)^(r+1))
			else return 0))))
	)

splineSyz(ZZ,PolyhedralComplex):=(r,PC)->(
	d:=dim PC;
	V:=vertices PC;
	V0:=apply(V,v->(flatten entries v));
	Facets:=polyhedra(d,PC);
	Faces:=polyhedra(d-1,PC);
	FV0:=apply(Facets,F->(vertices F) );
	FV1:=apply(FV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) ) ))));
	fV0:=apply(Faces,f->(vertices f) );
	fV1:=apply(fV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) ) ))));
	splineSyz(r,{V0,FV1,fV1}) 
	)

--Inputs:  A list of vertices, a list of top dim facets, a list of codim 1 faces in two dimensions
--Outputs:  Two matrices.  The image of the first is the symbolic module of combinatorial splines, the image of the second is the module of combinatorial splines.

combSplines=method()
combSplines(List):=L->(
	V=L_0;
	F=L_1;
	f=L_2;
	IntE=select(f,e->(
		(length select(F,a->(isSubset(e,a))))==2
		));
	n=(length IntE)-1;
	R=QQ[e_0..e_n,x,y,z];
	vs=flatten entries vars R;	
	fs=apply(F,f->(positions(IntE,e->(isSubset(e,f)))));
	dtab=apply(fs,j->(product apply(j,i->e_i)));
	etab=table(fs,length IntE,(f,i)->(
		if isSubset({i},f) then(
			T=select(fs,F->isSubset({i},F));
			ind=toList((set(T_0)+set(T_1))-((set T_0)*(set T_1)));
			product apply(ind,i->e_i)) else 0
		));
	Mat=transpose(matrix(R,V))||matrix(R,{apply(V,t->1)});
	DGList=apply(IntE,e->(gens gb minors(3,(transpose(matrix{{x,y,z}}))|(Mat_e) ) )_(0,0));
	if any(DGList,x->(x==sub(1,R))) then return "Error: some vertices on entered face not in codim 1 face";
	MG=(diagonalMatrix dtab)|(matrix etab)|(transpose matrix(R,{apply(length fs,i->1)}));
	rules=apply(n+1,i->(e_i=>DGList_i));
	MS0=sub(MG,rules);
	S=QQ[x,y,z];
	MS=sub(MS0,S);
	{MG,MS}
)

combSplines(PolyhedralComplex):=PC->(
	V:=vertices PC;
	V0:=apply(V,v->(flatten entries v));
	Facets:=polyhedra(2,PC);
	Faces:=polyhedra(1,PC);
	FV0:=apply(Facets,F->(vertices F) );
	FV1:=apply(FV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) ) ))));
	fV0:=apply(Faces,f->(vertices f) );
	fV1:=apply(fV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) ) ))));
	combSplines({V0,FV1,fV1})
)


--building vertex splines--

vertSplines=method();
vertSplines(List):=L->(
	vert=L_0;
	fs=L_1;
	es=L_2;
	IntE=select(es,e->(
		(length select(fs,a->(isSubset(e,a))))==2
		));
	R=QQ[x,y,z];
	Matr=transpose(matrix(R,vert))||matrix(R,{apply(vert,t->1)});
	FList=apply(IntE,e->(gens gb minors(3,(transpose(matrix{{x,y,z}}))|(Matr_e) ) )_(0,0));
	BndE=toList(set(es)-IntE);
	Bv=unique flatten BndE;
	Vn=toList(0..length(vert)-1);
	mlist=apply(Vn,i->matrix(
--list of edges containing vertex i--		
		pev=positions(es,e->isSubset({i},e));
		ev=es_pev;
--list of vertices adjacent to i--
		lvi=sort unique flatten ev;
		lv=vert_lvi;
--renumber edges as just between vertices adjacent to i--
		nev=apply(ev,e->(
			apply(e,n->position(lvi,i->(i==n)))));
--faces containing i--
		pfv=positions(fs,t->isSubset({i},t));
		fv=fs_pfv;
--list of all vertices in faces containing i--
		lvitot=sort unique flatten fv;
		lvtot=vert_lvitot;
--list of edges connecting edges adjacent to i (making link(i)) and renumbering--
		auxev=apply(fv,t->(
			Y=select(ev,e->(isSubset(e,t)));
			toList((set(Y_0)+set(Y_1))-set({i}))
			));
		nauxev=apply(auxev,e->(
			apply(e,n->position(lvi,i->(i==n)))));
--list of faces in triangulated star of vertex i and renumbering--
		auxfv=apply(length fv,n->(prepend(i,auxev_n)));
		nauxfv=apply(auxfv,e->(
			apply(e,n->position(lvi,i->(i==n)))));
--edges in star(i)--
		tev=join(ev,auxev);
		ntev=join(nev,nauxev);
--discard edges on the boundary of complex--
		bauxev=toList(set(tev)*set(BndE));
		nbauxev=apply(bauxev,e->(
			apply(e,n->position(lvi,i->(i==n)))));
--compute splines on star(i), just interested in the linear spline--
		M0=gens ker bsplines(0,{lv,nauxfv,ntev,nbauxev});
		phi=map(R,M0.ring,{x,y,z});
		M1=phi(M0);
		p=(id_(R^(numrows M1)))^(toList(0..length pfv -1));
		lin=(p*M1)_0;
--edges in faces containing i--
		esvtot=select(IntE,e->(
			test=false; 
			cnt=0; 
			while (test==false and cnt<(length fv)) do(test=isSubset(e,fv_cnt);cnt=cnt+1); 
			test));
--take the product of all forms corresponding to edges of faces containing i but not in star(i)--
		esforms=select(esvtot,e->(not isSubset(e,lvi)));
		pesforms=apply(esforms,e->position(IntE,j->(j==e)));
		P=product(FList_pesforms);
--build the continuous spline corresponding to vertex i--
		apply(length fs,n->(
			if isSubset({n},pfv) then (pos=position(pfv,x->(x==n)); return {P*(lin_pos)}) else return {0}
			))
		));
	mat=mlist_0;
	for i from 1 to length(mlist)-1 do(mat=mat|mlist_i);
	mat
	)


vertSplines(PolyhedralComplex):=PC->(
	V:=vertices PC;
	V0:=apply(V,v->(flatten entries v));
	Facets:=polyhedra(2,PC);
	Faces:=polyhedra(1,PC);
	FV0:=apply(Facets,F->(vertices F) );
	FV1:=apply(FV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) ) ))));
	fV0:=apply(Faces,f->(vertices f) );
	fV1:=apply(fV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) ) ))));
	vertSplines({V0,FV1,fV1})
)

cvSplines=method()
cvSplines(List):=L->(
	CS0=(combSplines(L))_1;
	VS=vertSplines(L);
	phi=map(VS.ring,CS0.ring,{x,y,z});
	CS=phi CS0;
	(CS|VS)
)

cvSplines(PolyhedralComplex):=PC->(
	CS0=(combSplines(PC))_1;
	VS=vertSplines(PC);
	phi=map(VS.ring,CS0.ring,{x,y,z});
	CS=phi CS0;
	(CS|VS)
)

--computes quotient module of C^0(\hat{\PC}) modulo CS(\hat{\PC})--

splineQuotient=method()
splineQuotient(List):=L->(
	K=gens ker splines(0,L);	
	numF=length L_1;
	numE=length L_2;
	pr=(id_((K.ring)^(numrows K)))^(toList(0..numF-1));
	C01=pr*K;
	CV=cvSplines(L);
	phi=map(CV.ring,K.ring,{x,y,z});
	C0=phi(C01);
	(image C0)/(image CV)
)

splineQuotient(PolyhedralComplex):=PC->(
	K=gens ker splines(0,PC);	
	numF=length L_1;
	numE=length L_2;
	pr=(id_((K.ring)^(numrows K)))^(toList(0..numF-1));
	C01=pr*K;
	CV=cvSplines(L);
	phi=map(CV.ring,K.ring,{x,y,z});
	C0=phi(C01);
	(image C0)/(image CV)
)

--Inputs:  A list of vertices, a list of top dim facets, a list of codim 1 faces in two dimensions
--Outputs:  A matrix presentation for the module H^0(A) (isomorphic to H^1(C) in the case of a disk)

h0pres=method();
h0pres(ZZ,List):=(r,L)->(
	V:=L_0;
	F:=L_1;
	f:=L_2;
	S=QQ[x,y,z];
	--boundary edges--
	eb=select(f,e->(
		cl=apply(length(F),j->isSubset(set(e),set(F_j)));
		#(select(cl,x->(x==true)))==1 ));
	--boundary vertices--
	vb=toList set flatten eb;
	--interior vertices		
	vint=toList(0..length(V)-1)-set(vb);
	--interior edges--
	ein=f-set(eb);
	--not totally interior edges--	
	enti=select(toList(0..length(ein)-1),i->(#(ein_i-set(vb))<=1));
	--equations for interior edges--
	M=transpose(matrix(S,V))||matrix(S,{apply(V,t->1)});
	DGList=for i from 0 to length(ein)-1 list(
		det( (transpose(vars S))|(M_(ein_i)) ) );
	--matrix recording positions of not totally interior edges in interior edge list--
	H0=S**(id_(ZZ^(length(ein))))_enti;
	--the matrix whose columns are syzygies	of the forms around each interior vertex, plus a generator for each not totally interior edge--
	for i from 0 to length(vint)-1 do(
		topH=matrix{apply(length(ein),j->(
			if isSubset(set({vint_i}),set(ein_j)) then (DGList_j)^(r+1) else 0*x))};
		evind=positions(ein,e->isSubset(set({vint_i}),set(e)));
		cevind=toList(0..length(ein)-1)-set(evind);
		bottomH=S**(transpose (id_(ZZ^(length(ein))))_cevind);
		H0=syz(topH||bottomH)|H0 );
		M=map((S^{-r-1})^(numrows H0),(S^{-2*r-2})^(numcols H0),H0);
	M)

h0pres(PolyhedralComplex):=PC->(
	d:=dim PC;
	V:=vertices PC;
	V0:=apply(V,v->(flatten entries v));
	Facets:=polyhedra(d,PC);
	Faces:=polyhedra(d-1,PC);
	FV0:=apply(Facets,F->(vertices F) );
	FV1:=apply(FV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) ) ))));
	fV0:=apply(Faces,f->(vertices f) );
	fV1:=apply(fV0,M->(apply(numcols M,i->(position(V,v->( (vector flatten entries transpose v)==M_i ) ) ))));
	h0pres({V0,FV1,fV1}) )

--Inputs: A list of vertices, a list of top dim facets, a list of codim 1 faces, and the vertices of a new edge
--Outputs:  The annihilator of element corresponding to the new edge in H_0(A)

edgeAnn=method();
edgeAnn(List,List):=(L,e)->(
	V:=L_0;
	F:=L_1;
	f:=L_2;
	S=QQ[x,y,z];
	sf0:= select(F,x->isSubset(e,x));
	if (length sf0)>=2 then return "Edge is already in complex";
	if (length sf0)==0 then return "Edge does not subdivide a face";
	sf:=flatten sf0;
	l:=det matrix {append(V_(e_0),1),append(V_(e_1),1),{x,y,z}};
	sf10:=select(sf,i->(
		ev:=sub(l,{x=>(V_i)_0,y=>(V_i)_1,z=>1});
		ev>0)
	    );
	sf20:=select(sf,i->(
		ev:=sub(l,{x=>(V_i)_0,y=>(V_i)_1,z=>1});
		ev<0)
	    );
	sf1:=join(sf10,e);
	sf2:=join(sf20,e);
	F1:=join(F-set({sf}),{sf1,sf2});
	f1:=append(f,e);
	L1:={V,F1,f1};
	H:=coker h0pres L;
	M1:=h0pres L1;
	r:=numrows M1-1;
	pos:=positions(flatten entries M1^{r},i->(i!=0));
	N:=M1_pos;
	c0:=positions(toList(0..r-1),i->(N_(i,0)!=0));
	if (length c0)==0 then c0={0,0};
	if (length c0)==1 then c0=append(c0,(c0_0+1)%r);
	c1:=positions(toList(0..r-1),i->(N_(i,1)!=0));
	if (length c1)==0 then c1={0,0};
	if (length c1)==1 then c1=append(c1,(c1_0+1)%r);
	a0:=sub(N_(r,0),QQ);
	a1:=sub(N_(c0_0,0),QQ);
	a2:=sub(N_(c0_1,0),QQ);
	b0:=sub(N_(r,1),QQ);
	b1:=sub(N_(c1_0,1),QQ);
	b2:=sub(N_(c1_1,1),QQ);	
	m:=image(a1/a0*H_{c0_0}+a2/a0*H_{c0_1}-b1/b0*H_{c1_0}-b2/b0*H_{c1_1});
	{annihilator m,L1}
)

--Returns the portion of the presentation matrix of H^0(A) that has degree 0--

h0pres0=method();
h0pres0(List):=L->(
	V=L_0;
	F=L_1;
	f=L_2;
	S=QQ[x,y,z];
	eb=select(f,e->(
		cl=apply(length(F),j->isSubset(set(e),set(F_j)));
		#(select(cl,x->(x==true)))==1 ));
	vb=toList set flatten eb;
	vint=toList(0..length(V))-set(vb);
	ein=f-set(eb);
	enti=select(toList(0..length(ein)-1),i->(#(ein_i-set(vb))<=1));
	M=transpose(matrix(S,V))||matrix(S,{apply(V,t->1)});
	DGList=for i from 0 to length(ein)-1 list(
		det( (transpose(vars S))|(M_(ein_i)) ) );
	H00=S**(id_(ZZ^(length(ein))))_enti;
	for i from 0 to length(vint)-1 do(
		topH=matrix{apply(length(ein),j->(
			if isSubset(set({vint_i}),set(ein_j)) then DGList_j else 0*x))};
		evind=positions(ein,e->isSubset(set({vint_i}),set(e)));
		cevind=toList(0..length(ein)-1)-set(evind);
		bottomH=S**(transpose (id_(ZZ^(length(ein))))_cevind);
		SYZ=syz(topH||bottomH);
		cls=toList(0..numcols SYZ -2);
		H00=SYZ_cls|H00);
	H00)

--Computes an echelon form of a matrix, returning a list of columns with leading 1s (code mostly copied from Polyhedra).  Note: The rows of the matrix can be viewed as coefficients of linear forms whose common vanishing locus is a hyperplane.  Interpreted thus, the list of columns with leading 1s is a list of variables which can be considered as dependent variables on the hyperplane.--

ref=method();
ref(Matrix):=M->(
	n=numcols M;
	s=numrows M;
	m=min(n,s);
	N=mutableMatrix(M**QQ);
	i:=0;
	stopper:=0;
	leading:={};
	while i<m and stopper < n do(
		j:=select(1,toList(i..s-1),k->N_(k,stopper) != 0);
		if j != {} then (
			j=j#0;
			leading=append(leading,stopper);
			scan((j+1)..(s-1),k->(
				if N_(k,i) != 0 then(
					a:=N_(j,i);
					b:=N_(k,i);
					N=rowAdd(N,k,-b/a,j) ) ) );
			if i != j then(N=rowSwap(N,i,j) );
			i=i+1);
		stopper=stopper+1);
	leading)

--Inputs: an integer i and a polyhedron P. 
--Outputs: Cellular boundary map delta_i, from codim i faces to codim i+1 faces. 
--For a pair (l1,l2) consisting of a (d-i) dimensional face and a (d-i-1) dimensional face the function determines whether +1,-1, or 0 should be assigned as follows.  The orientation on each face (cell) is assigned to be the wedge of all local coordinates (as determined by ref) in order of increasing index.  If l2 is not a face of l1, 0 is assigned.  If l2 is a face of l1, +1 is assigned if the orientation of l1 induces the chosen orientation on l2, -1 is assigned otherwise.

bMap=method(); 
bMap(ZZ,Polyhedron):=(i,P)->(
	L1:=faces(i,P); 
	L2:=faces(i+1,P); 
	n:=P#"ambient dimension"; 
	d:=dim P;

--sets indices of variables considered dependent on each face.  The orientation on each face is chosen to be the wedge of independent variables (technically the wedge of differentials dx_i for independent variables x_i) in order of increasing index--
	O1:=apply(L1,f1->(
		H:=(hyperplanes(f1))#0;
		ref(H)) );
	O2:=apply(L2,f2->(
		H:=(hyperplanes(f2))#0;
		ref(H)) );

	--the top differential--

	if (n==d and i==0) then return transpose matrix apply(L1,l1->(
		HS:=halfsort(halfspaces(l1),L2) ** QQ; 
		apply(L2,l2->(
			p2=position(L2,u->(u==l2)); 
			nm=flatten entries(HS^{p2}); 
			p3=position(nm,g->g!=0); 
			H2=(hyperplanes(l2))#0 ** QQ; 
			M2=inverse(H2_(O2_p2))*H2; 
			S1=sort(toList(set(0..n-1)-set({p3}) )); 
			S2=sort(toList(set(0..n-1)-set(O2_p2))); 
			D=det(matrix table(length S1,length S2,(i,j)->(
				if S1_i==S2_j then return 1/1; 
				if member(S1_i,set(S2)) and S1_i != S2_j then return 0/1; 
				if member(S1_i,O2_p2) then return -(M2_(S2))_(i,j) ) ) ); 
			if (-1)^(p3)*nm_p3*D<0 then return -1; 
			if (-1)^(p3)*nm_p3*D>0 then return 1) ) ) ); 

	--the middle differential--

	if (0<i and i<d-1) or (n != d and i==0) then return transpose matrix apply(L1,l1->(
		f1=faces(1,l1); 
		H1=(hyperplanes(l1))#0 ** QQ; 
		HS=halfsort(halfspaces(l1),f1) ** QQ; 
		p1=position(L1,v->(v==l1)); 
		--dependent variables of l1 are written in terms of independent variables determined in ref.  Relations encoded in M1.--
		M1=inverse(H1_(O1_p1))*H1; 
		apply(L2,l2->( 
			if not isFace(l2,l1) then return 0; 
			p2=position(L2,u->(u==l2)); 
			pf=position(f1,q->(q==l2)); 
			N=mutableMatrix (M1||HS^{pf}); 
			for c from 0 to (n-d+i-1) do rowAdd(N,n-d+i,-HS_(pf,O1_p1_c),c); 
		--nm is the outward normal to l2 in the local coordinates of l1--			
			nm=flatten entries((matrix N)^{n-d+i}); 
		--induced local coordinates on l2 are all local coordinates on l1 except the one indexed by p3--
			p3=position(nm,g->g!=0); 
			H2=(hyperplanes(l2))#0 ** QQ;
		--dependent variables of l2 are written in terms of independent variables determined in ref.  Relations encoded in M2.--
			M2=inverse(H2_(O2_p2))*H2; 
			S=sort(toList(set(0..n-1)-set(O1_p1))); 
			S1=sort(toList(set(0..n-1)-(set(O1_p1)+set({p3})) )); 
			S2=sort(toList(set(0..n-1)-set(O2_p2)));
		--D is determinant of matrix changing from induced local coordinates to chosen local coordinates--
			D=det(matrix table(length S1,length S2,(i,j)->(
				if S1_i==S2_j then return 1/1; 
				if member(S1_i,set(S2)) and S1_i != S2_j then return 0/1; 
				if member(S1_i,O2_p2) then return -(M2_S2)_(i,j) ) ) );
		 --The consistency of the induced orientation with the chosen one is determined by the sign of the first nonvanishing entry of nm, the position in which this entry appears in the local coordinates of l1, and the sign of D.--
			if (-1)^(position(S,t->(t==p3)))*nm_p3*D<0 then return -1; 
			if (-1)^(position(S,t->(t==p3)))*nm_p3*D>0 then return 1) ) ) );

	--the bottom differential--

	if i==d-1 then return transpose matrix apply(L1,l1->(
		f1=faces(1,l1); 
		H1=(hyperplanes(l1))#0 ** QQ; 
		HS=halfsort(halfspaces(l1),f1) ** QQ; 
		p1=position(L1,v->(v==l1)); 
		M1=inverse(H1_(O1_p1))*H1; 
		apply(L2,l2->( 
			if not isFace(l2,l1) then return 0; 
			p2=position(L2,u->(u==l2)); 
			pf=position(f1,q->(q==l2)); 
			N=mutableMatrix (M1||HS^{pf}); 
			for c from 0 to (n-d+i-1) do rowAdd(N,n-d+i,-HS_(pf,O1_p1_c),c); 
			nm=flatten entries((matrix N)^{n-d+i}); 
			p3=position(nm,g->g!=0); 
			if nm_p3<0 then return -1; 
			if nm_p3>0 then return 1) ) ) ) )

bMap(ZZ,PolyhedralComplex):=(i,PC)->(
	n:=PC#"ambient dimension"; 
	d:=dim PC;	
	L1:=polyhedra(d-i,PC); 
	L2:=polyhedra(d-i-1,PC);  
	O1:=apply(L1,f1->(
		H:=(hyperplanes(f1))#0;
		ref(H) ) );
	O2:=apply(L2,f2->(
		H:=(hyperplanes(f2))#0;
		ref(H) ) );
	if (n==d and i==0) then return transpose matrix apply(L1,l1->(
		f1=faces(1,l1);
		HS:=halfsort(halfspaces(l1),f1) ** QQ; 
		apply(L2,l2->(
			if not isFace(l2,l1) then return 0;
			p2=position(L2,u->(u==l2));
			pf=position(f1,q->(q==l2));
			nm=flatten entries(HS^{pf}); 
			p3=position(nm,g->g!=0); 
			H2=(hyperplanes(l2))#0 ** QQ; 
			M2=inverse(H2_(O2_p2))*H2; 
			S1=sort(toList(set(0..n-1)-set({p3}) )); 
			S2=sort(toList(set(0..n-1)-set(O2_p2))); 
			D=det(matrix table(length S1,length S2,(i,j)->(
				if S1_i==S2_j then return 1/1; 
				if member(S1_i,set(S2)) and S1_i != S2_j then return 0/1; 
				if member(S1_i,O2_p2) then return -(M2_(S2))_(i,j) ) ) ); 
			if (-1)^(p3)*nm_p3*D<0 then return -1; 
			if (-1)^(p3)*nm_p3*D>0 then return 1) ) ) ); 
	if (0<i and i<d-1) or (n != d and i==0) then return transpose matrix apply(L1,l1->(
		f1=faces(1,l1); 
		H1=(hyperplanes(l1))#0 ** QQ; 
		HS=halfsort(halfspaces(l1),f1) ** QQ; 
		p1=position(L1,v->(v==l1)); 
		M1=inverse(H1_(O1_p1))*H1; 
		apply(L2,l2->( 
			if not isFace(l2,l1) then return 0; 
			p2=position(L2,u->(u==l2)); 
			pf=position(f1,q->(q==l2)); 
			N=mutableMatrix (M1||HS^{pf}); 
			for c from 0 to (n-d+i-1) do rowAdd(N,n-d+i,-HS_(pf,O1_p1_c),c); 
			nm=flatten entries((matrix N)^{n-d+i}); 
			p3=position(nm,g->g!=0); 
			H2=(hyperplanes(l2))#0 ** QQ; 
			M2=inverse(H2_(O2_p2))*H2; 
			S=sort(toList(set(0..n-1)-set(O1_p1))); 
			S1=sort(toList(set(0..n-1)-(set(O1_p1)+set({p3})) )); 
			S2=sort(toList(set(0..n-1)-set(O2_p2))); 
			D=det(matrix table(length S1,length S2,(i,j)->(
				if S1_i==S2_j then return 1/1; 
				if member(S1_i,set(S2)) and S1_i != S2_j then return 0/1; 
				if member(S1_i,O2_p2) then return -(M2_S2)_(i,j) ) ) ); 
			if (-1)^(position(S,t->(t==p3)))*nm_p3*D<0 then return -1; 
			if (-1)^(position(S,t->(t==p3)))*nm_p3*D>0 then return 1) ) ) );
	if i==d-1 then return transpose matrix apply(L1,l1->(
		f1=faces(1,l1); 
		H1=(hyperplanes(l1))#0 ** QQ; 
		HS=halfsort(halfspaces(l1),f1) ** QQ; 
		p1=position(L1,v->(v==l1)); 
		M1=inverse(H1_(O1_p1))*H1; 
		apply(L2,l2->( 
			if not isFace(l2,l1) then return 0; 
			p2=position(L2,u->(u==l2)); 
			pf=position(f1,q->(q==l2)); 
			N=mutableMatrix (M1||HS^{pf}); 
			for c from 0 to (n-d+i-1) do rowAdd(N,n-d+i,-HS_(pf,O1_p1_c),c); 
			nm=flatten entries((matrix N)^{n-d+i}); 
			p3=position(nm,g->g!=0); 
			if nm_p3<0 then return -1; 
			if nm_p3>0 then return 1) ) ) ) )

--Outputs the differential relative to the boundary--

relbMap=method();
relbMap(ZZ,PolyhedralComplex):=(i,PC)->(
	n:=PC#"ambient dimension"; 
	d:=dim PC;
	BC:=bComp(PC);
	P1:=polyhedra(d-i,PC);
	P2:=polyhedra(d-i-1,PC);
	SP1:= (if i==0 then set({}) else set(polyhedra(d-i,BC)));
	SP2:=set(polyhedra(d-i-1,BC));
	L1:=select(P1,x->(not member(x,SP1))); 
	L2:=select(P2,y->(not member(y,SP2)));
	O1:=apply(L1,f1->(
		H:=(hyperplanes(f1))#0;
		ref(H)) );
	O2:=apply(L2,f2->(
		H:=(hyperplanes(f2))#0;
		ref(H)) );
	if (n==d and i==0) then return transpose matrix apply(L1,l1->(
		f10=faces(1,l1);
		I1=select(toList(0..(#f10-1)),x->(not member(f10_x,SP2)));
		f1=apply(I1,x->(f10_x));
		HS1:=halfsort(halfspaces(l1),f10) ** QQ;
		HS=HS1^I1;
		apply(L2,l2->(
			if not isFace(l2,l1) then return 0;
			p2=position(L2,u->(u==l2));
			pf=position(f1,q->(q==l2));
			nm=flatten entries(HS^{pf}); 
			p3=position(nm,g->g!=0); 
			H2=(hyperplanes(l2))#0 ** QQ; 
			M2=inverse(H2_(O2_p2))*H2; 
			S1=sort(toList(set(0..n-1)-set({p3}) )); 
			S2=sort(toList(set(0..n-1)-set(O2_p2))); 
			D=det(matrix table(length S1,length S2,(i,j)->(
				if S1_i==S2_j then return 1/1; 
				if member(S1_i,set(S2)) and S1_i != S2_j then return 0/1; 
				if member(S1_i,O2_p2) then return -(M2_(S2))_(i,j) ) ) ); 
			if (-1)^(p3)*nm_p3*D<0 then return -1; 
			if (-1)^(p3)*nm_p3*D>0 then return 1) ) ) ); 
	if (0<i and i<d-1) or (n != d and i==0) then return transpose matrix apply(L1,l1->(
		f10=faces(1,l1);
		I1=select(toList(0..(#f10-1)),x->(not member(f10_x,SP2)));
		f1=apply(I1,x->(f10_x));
		H1=(hyperplanes(l1))#0 ** QQ; 
		HS1=halfsort(halfspaces(l1),f10) ** QQ;
		HS=HS1^I1;
		p1=position(L1,v->(v==l1)); 
		M1=inverse(H1_(O1_p1))*H1; 
		apply(L2,l2->( 
			if not isFace(l2,l1) then return 0; 
			p2=position(L2,u->(u==l2)); 
			pf=position(f1,q->(q==l2)); 
			N=mutableMatrix (M1||HS^{pf}); 
			for c from 0 to (n-d+i-1) do rowAdd(N,n-d+i,-HS_(pf,O1_p1_c),c); 
			nm=flatten entries((matrix N)^{n-d+i}); 
			p3=position(nm,g->g!=0); 
			H2=(hyperplanes(l2))#0 ** QQ; 
			M2=inverse(H2_(O2_p2))*H2; 
			S=sort(toList(set(0..n-1)-set(O1_p1))); 
			S1=sort(toList(set(0..n-1)-(set(O1_p1)+set({p3})) )); 
			S2=sort(toList(set(0..n-1)-set(O2_p2))); 
			D=det(matrix table(length S1,length S2,(i,j)->(
				if S1_i==S2_j then return 1/1; 
				if member(S1_i,set(S2)) and S1_i != S2_j then return 0/1; 
				if member(S1_i,O2_p2) then return -(M2_S2)_(i,j) ) ) ); 
			if (-1)^(position(S,t->(t==p3)))*nm_p3*D<0 then return -1; 
			if (-1)^(position(S,t->(t==p3)))*nm_p3*D>0 then return 1) ) ) );
	if i==d-1 then return transpose matrix apply(L1,l1->(
		f10=faces(1,l1);
		I1=select(toList(0..(#f10-1)),x->(not member(f10_x,SP2)));
		f1=apply(I1,x->(f10_x));
		H1=(hyperplanes(l1))#0 ** QQ; 
		HS1=halfsort(halfspaces(l1),f10) ** QQ; 
		HS=HS1^I1;
		p1=position(L1,v->(v==l1)); 
		M1=inverse(H1_(O1_p1))*H1; 
		apply(L2,l2->( 
			if not isFace(l2,l1) then return 0; 
			p2=position(L2,u->(u==l2)); 
			pf=position(f1,q->(q==l2)); 
			N=mutableMatrix (M1||HS^{pf}); 
			for c from 0 to (n-d+i-1) do rowAdd(N,n-d+i,-HS_(pf,O1_p1_c),c); 
			nm=flatten entries((matrix N)^{n-d+i}); 
			p3=position(nm,g->g!=0); 
			if nm_p3<0 then return -1; 
			if nm_p3>0 then return 1) ) ) ) )




--Constructs Spline complex--
--Inputs: A polyhedral complex P--
--Outputs: the Spline Complex--


splineComplex=method();
splineComplex(ZZ,PolyhedralComplex):=(r,PC)->(
	n:=PC#"ambient dimension"; 
	d:=dim PC;
	S:=QQ[t_0..t_n];
	BC:=bComp(PC);

--Sets up list of interior faces in order of dim--	
	IntF:=for i from 0 to d list(
		Faces=polyhedra(i,PC);
		if i==d then continue Faces;
		BCP=polyhedra(i,BC);
		select(Faces,f->(not member(f,set(BCP) ) ) ) );

--Sets local coordinates on interior faces as the variables corresponding to indices of leading ones in a row echelon form for the matrix of hyperplanes defining the face--
	lCoords:=apply(IntF, j-> apply(j, x->(
		H:=(hyperplanes(x))#0;
		ref(H) ) ) );

--Sets up topological boundary maps d_i from interior dim i faces to interior dim (i-1) faces--
	DD:=for i from 1 to d list(
		L1:=IntF_i;
		L2:=IntF_(i-1);		
		O1:=lCoords_i;
		O2:=lCoords_(i-1);
		if L2=={} then continue (if L1=={} then (transpose matrix{{}})*matrix{{}} else transpose matrix apply(#L1,x->{}) );			
		if i==1 then continue transpose matrix apply(L1,l1->(
			f10=faces(1,l1);
			I1=select(toList(0..(#f10-1)),x->(member(f10_x,set(L2) ) ));
			f1=apply(I1,k->f10_k);
			H1=(hyperplanes(l1))#0 ** QQ; 
			HS1=halfsort(halfspaces(l1),f10) ** QQ; 
			HS=HS1^I1;
			p1=position(L1,v->(v==l1)); 
			M1=inverse(H1_(O1_p1))*H1; 
			apply(L2,l2->( 
				if not isFace(l2,l1) then return 0; 
				p2=position(L2,u->(u==l2)); 
				pf=position(f1,q->(q==l2)); 
				N=mutableMatrix (M1||HS^{pf}); 
				for c from 0 to (n-i-1) do rowAdd(N,n-i,-HS_(pf,O1_p1_c),c); 
				nm=flatten entries((matrix N)^{n-i}); 
				p3=position(nm,g->g!=0); 
				if nm_p3<0 then return -1; 
				if nm_p3>0 then return 1) ) ) );		
		if (1<i and i<d) or (n != d and i==d) then continue transpose matrix apply(L1,l1->(
			f10=faces(1,l1);
			I1=select(toList(0..(#f10-1)),x->(member(f10_x,set(L2)) ));
			f1=apply(I1,k->f10_k);
			H1=(hyperplanes(l1))#0 ** QQ; 
			HS1=halfsort(halfspaces(l1),f10) ** QQ;
			HS=HS1^I1;
			p1=position(L1,v->(v==l1)); 
			M1=inverse(H1_(O1_p1))*H1; 
			apply(L2,l2->( 
				if not isFace(l2,l1) then return 0; 
				p2=position(L2,u->(u==l2)); 
				pf=position(f1,q->(q==l2)); 
				N=mutableMatrix (M1||HS^{pf}); 
				for c from 0 to (n-i-1) do rowAdd(N,n-i,-HS_(pf,O1_p1_c),c); 
				nm=flatten entries((matrix N)^{n-i}); 
				p3=position(nm,g->g!=0); 
				H2=(hyperplanes(l2))#0 ** QQ; 
				M2=inverse(H2_(O2_p2))*H2; 
				S0=sort(toList(set(0..n-1)-set(O1_p1))); 
				S1=sort(toList(set(0..n-1)-(set(O1_p1)+set({p3})) )); 
				S2=sort(toList(set(0..n-1)-set(O2_p2))); 
				D=det(matrix table(length S1,length S2,(i,j)->(
					if S1_i==S2_j then return 1/1; 
					if member(S1_i,set(S2)) and S1_i != S2_j then return 0/1; 
					if member(S1_i,O2_p2) then return -(M2_S2)_(i,j) ) ) ); 
				if (-1)^(position(S0,t->(t==p3)))*nm_p3*D<0 then return -1; 
				if (-1)^(position(S0,t->(t==p3)))*nm_p3*D>0 then return 1) ) ) );
		if (n==d and i==d) then continue transpose matrix apply(L1,l1->(
			f10=faces(1,l1);
			I1=select(toList(0..(#f10-1)),x->(member(f10_x,set(L2) ) ));
			f1=apply(I1,k->f10_k);
			HS1:=halfsort(halfspaces(l1),f10) ** QQ;
			HS=HS1^I1;
			apply(L2,l2->(
				if not isFace(l2,l1) then return 0;
				p2=position(L2,u->(u==l2));
				pf=position(f1,q->(q==l2));
				nm=flatten entries(HS^{pf}); 
				p3=position(nm,g->g!=0); 
				H2=(hyperplanes(l2))#0 ** QQ; 
				M2=inverse(H2_(O2_p2))*H2; 
				S1=sort(toList(set(0..n-1)-set({p3}) )); 
				S2=sort(toList(set(0..n-1)-set(O2_p2))); 
				D=det(matrix table(length S1,length S2,(i,j)->(
					if S1_i==S2_j then return 1/1; 
					if member(S1_i,set(S2)) and S1_i != S2_j then return 0/1; 
					if member(S1_i,O2_p2) then return -(M2_(S2))_(i,j) ) ) ); 
				if (-1)^(p3)*nm_p3*D<0 then return -1; 
				if (-1)^(p3)*nm_p3*D>0 then return 1) ) ) ) ); 

--Initializing Spline Chain Complex--
	C:=new ChainComplex;
	C.ring=S;

--Sets up the maps whose images are the ideals J(f) for each interior face f--
	IFacets=IntF_(d-1);
	E:=apply(d+1,i->(
		Faces=polyhedra(i,PC);
		IFaces=IntF_i;
		if IFaces=={} then return IFaces;
		if i !=d then return(		
		Flist=apply(IFaces,f->(
			select(IFacets,F->contains(F,f))));
		apply(Flist,x->(
			matrix{apply(x,F->(
				(((((hyperplanes(F))#0|(hyperplanes(F))#1)**S)*(transpose vars S))_(0,0))^(r+1) ) )} ) ));
		if i==d then return apply(IFaces,F->(matrix{{}}**S));
		));

--Initializing Modules--
	MList:=for j from 0 to d list(
		if E_j=={} then continue source matrix{{}} **S;
		directSum toSequence(apply(E_j,e->(coker e))));

--Installing Maps--
	for i from 1 to d do(
		C.dd#i=map(MList#(i-1),MList#i,(DD_(i-1))**S) );


--Installing Modules--
	C#0=target C.dd#1;
	for k from 1 to d do(
		C#k=source C.dd#k);

--Initializing and Installing Cellular Chain Comlex--
	T:=chainComplex toSequence DD;

--Returning the Spline and Topological Complexes--	

	{C,T})



--Constructs Spline Complex for Simplicial Complexes--
--Inputs:  A simplicial complex P
--Outputs:  Spline complex for P

splineComplexS=method()
splineComplexS(ZZ,PolyhedralComplex):=(r,SC)->(
	n:=SC#"ambient dimension"; 
	d:=dim SC;
	S:=QQ[t_0..t_n];
	BC:=bComp(SC);

--Sets up list of interior faces in order of dim--	
	IPos:=for i from 0 to d list(
		Faces=polyhedra(i,SC);
		if i==d then continue toList(0..(length(Faces)-1));
		BCP=polyhedra(i,BC);
		positions(Faces,f->(not member(f,set(BCP) ) ) ) );

--Sets up topological boundary maps d_i from interior dim i faces to interior dim (i-1) faces--
	DD:=for i from 1 to d list(
		M=boundaryMap(i,SC);
		L1=IPos_i;
		L2=IPos_(i-1);
		if L2=={} then continue(
			 if L1=={} then (transpose matrix{{}})*matrix{{}} else transpose matrix apply(#L1,x->{}) );
		M_L1^L2
		);
		
--Sets up the maps whose images are the ideals J(f) for each interior face f--
	IFacets=(polyhedra(d-1,SC))_(IPos_(d-1));
	E:=apply(d+1,i->(
		Faces=polyhedra(i,SC);
		IFaces=Faces_(IPos_i);
		if IFaces=={} then return IFaces;
		if i !=d then return(		
		Flist=apply(IFaces,f->(
			select(IFacets,F->contains(F,f))));
		apply(Flist,x->(
			matrix{apply(x,F->(
				(((((hyperplanes(F))#0|(hyperplanes(F))#1)**S)*(transpose vars S))_(0,0))^(r+1) ) )} ) ));
		if i==d then return apply(IFaces,F->(matrix{{}}**S));
		));

--Initializing Spline Chain Complex--
	C:=new ChainComplex;
	C.ring=S;

--Initializing Modules--
	MList:=for j from 0 to d list(
		if E_j=={} then continue source matrix{{}} **S;
		directSum toSequence(apply(E_j,e->(coker e))));

--Installing Maps--
	for i from 1 to d do(
		C.dd#i=map(MList#(i-1),MList#i,(DD_(i-1))**S) );


--Installing Modules--
	C#0=target C.dd#1;
	for k from 1 to d do(
		C#k=source C.dd#k);

--Initializing and Installing Cellular Chain Comlex--
	T:=chainComplex toSequence DD;

--Returning the Spline and Topological Complexes--	

	{C,T})	

--Saves the spline complex to a file--
saveComplex=method();
saveComplex(ChainComplex,String):=(C,N)->(
	S=C.ring;
	d:=rank source vars S;
	L:=toSequence for i from 1 to d-1 list(C.dd#(i));
	H:=openOut N;
	H<<"S="<<toExternalString(S)<<endl;
	H<<"Maps="<<toExternalString(L)<<close )


