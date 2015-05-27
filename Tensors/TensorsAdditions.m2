newPackage(
	"TensorsAdditions",
    	Version => "0.1", 
    	Date => "May 27, 2015",
    	Authors => {
	     {Name => "Hirotachi Abo"},
	     {Name => "Roberto Barrera"},
	     {Name => "Robert Krone"},
	     {Name => "Benjamin Reames"},
	     {Name => "Zach Teitler"}
	     },
    	Headline => "tensors",
	AuxiliaryFiles => true,
	PackageExports => {
	  "Tensors"
	}
)

export {
  symmetrize
}

symmetrize = method()
symmetrize (Tensor) := (T) -> (
  error "not implemented yet";
);


end
