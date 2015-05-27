genericTensor=method()

genericTensor(Ring,ZZ,List):=(R,i,dims) ->(
     d:=product dims;
     ents:=flatten entries genericMatrix(R,R_i,1,d);
     M:=tensorModule(R,dims);
     new M from vector ents
     )

genericTensor(Ring,List):=(R,dims) -> genericTensor(R,0,dims)

getIndex := (R,x) -> (
     M := try monoid R else error "expected a polynomial ring or quotient of one";
     if class x =!= R then error "expected an element of the ring";
     x = try baseName x else error "expected a variable of the ring";
     M.index#x)

genericTensor(Ring,RingElement,List):=(R,x,dims) -> (
     genericTensor(R,getIndex(R,x),dims))

TEST ///
R=QQ[a..z]
genericTensor(R,{3,4})
///

randomTensor=method(TypicalValue=>Tensor)
randomTensor (Ring,List) := (R,dims) -> (
     try random(R) else error "no method for random(Ring)";
     n:=product dims;
     ents := for i in 0..<n list random(R);
     M:=tensorModule(R,dims);
     new M from vector ents
     )
randomTensor (Ring,Thing,List) := (R,deg,dims) -> (
     if not(instance(deg,ZZ) or instance(deg,List)) then error "randomTensor (Ring,Thing,List) expected Thing to be a ZZ or List";
     n:=product dims;
     ents := for i in 0..<n list random(deg,R);
     M:=tensorModule(R,dims);
     new M from vector ents
     )

TEST ///
randomTensor(ZZ,{3,3})
S=QQ[x,y]
randomTensor(S,1,{3,3})
///

end




