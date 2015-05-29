path=append(path,"GitHub/Workshop-2015-Boise/Splines/")
loadPackage("Polyhedra");

pComp=method();
pComp(List, List):=(V,F)->(
	PList := apply(F,x->(convexHull(transpose matrix V_x)));
	polyhedralComplex PList
)

--Computes an echelon form of a matrix, 
--returning a list of columns with leading 1s 
--(code mostly copied from Polyhedra).  
--Note: The rows of the matrix can be viewed as 
--coefficients of linear forms whose common 
--vanishing locus is a hyperplane.  
--Interpreted thus, the list of columns with 
--leading 1s is a list of variables which can 
--be considered as dependent variables on the 
--hyperplane.    
ref=method();
ref(Matrix):=M->(
	n := numcols M;
	s := numrows M;
	m := min(n,s);
	N := mutableMatrix(M**QQ);
	i := 0;
	stopper := 0;
	leading := {};
	while i<m and stopper < n do(
		j := select(1,toList(i..s-1),k->N_(k,stopper) != 0);
		if j != {} then (
			j = j#0;
			leading = append(leading,stopper);
			scan((j+1)..(s-1),k->(
				if N_(k,i) != 0 then(
					a := N_(j,i);
					b := N_(k,i);
					N = rowAdd(N,k,-b/a,j) ) ) );
			if i != j then (N=rowSwap(N,i,j) );
			i = i+1);
		stopper = stopper+1);
	leading
)


--shuffles the rows of halfspaces matrix 
--so that they appear in same order as codim 1 faces of polyhedron --
halfsort = method(); 
halfsort(Sequence,List):=(x,y)->(
	r := numrows(x#0); 
	L := apply(toList(0..r-1), i-> (
		position(toList(0..r-1), j->
		    {(x#1)_(i,0)} == unique flatten entries ((x#0)^{i} * vertices(y_j))
		    )
		)
	    );
	matrix apply(L,i-> flatten entries (x#0)^{i})
	)

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


--The Spline Complex--
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
--Sets local coordinates on interior faces as the 
--variables corresponding to indices of leading ones 
--in a row echelon form for the matrix of hyperplanes defining the face--
	lCoords:=apply(IntF, j-> apply(j, x->(
		H:=(hyperplanes(x))#0;
		ref(H) ) ) );
--Sets up topological boundary maps d_i from interior 
--dim i faces to interior dim (i-1) faces--
	DD:=for i from 1 to d list(
		L1:=IntF_i;
		L2:=IntF_(i-1);		
		O1:=lCoords_i;
		O2:=lCoords_(i-1);
		if L2=={} then continue (if L1=={} then (transpose matrix{{}})*matrix{{}} else transpose matrix apply(#L1,x->{}) );			
		if i==1 then continue transpose matrix apply(L1,l1->(
			f10=faces(1,l1);
			I1=select(#f10,x->(member(f10_x,set(L2) ) ));
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
    
V={{-1, -1, -1}, {3, -1, -1}, {-1, 3, -1}, {-1, -1, 3}, {-4, -4, -4}, {12, -4, -4}, {-4, 12, -4}, {-4, -4, 12}};
F={{0, 1, 2, 3}, {1, 2, 3, 5, 6, 7}, {0, 2, 3, 4, 6, 7}, {0, 1, 3, 4, 5, 7}, {0, 1, 2, 4, 5, 6}};
time P=pComp(V,F);
time CP=splineComplex(1,P);
    