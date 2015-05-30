{*
undocumented {(net,NCGroebnerBasis),
              (net,NCIdeal),
	      Basis}
*}
beginDocumentation()

-------------------------
----- Types
-------------------------
    
doc ///
  Key
    NCAlgebra
  Description
    Text
      This package is used to define and manipulate noncommutative algebras.  Many of the
      commands contain calls to the existing noncommutative algebra package Bergman.
  Subnodes
    "Basic operations on noncommutative algebras"
    "General setup information"
    "Using the Bergman interface"
///

doc ///
   Key
      NCRing
   Headline
      Type of a noncommutative ring
   Description
      Text
         All noncommutative rings have this as an ancestor type.  It is the parent of the
	 types @ TO NCPolynomialRing @ and @ TO NCQuotientRing @. 
      Text
         In addition to defining a ring as a quotient of a @ TO NCPolynomialRing @, some common ways to create
	 NCRings include @ TO skewPolynomialRing @, @ TO endomorphismRing @, and @ TO oreExtension @.      
      
         Let's consider a three dimensional Sklyanin algebra.  We first define the tensor algebra:
      Example
         A = QQ{x,y,z}
      Text
         Then input the defining relations, and put them in an ideal:
      Example
	 f = y*z + z*y - x^2
	 g = x*z + z*x - y^2
	 h = z^2 - x*y - y*x
     	 I=ncIdeal{f,g,h}
      Text
         Next, define the quotient ring (as well as try a few functions on this new ring).  Note that
	 when the quotient ring is defined, a call is made to Bergman to compute the Groebner basis
	 of I (out to a certain degree, should the Groebner basis be infinite).
      Example
	 B=A/I
	 generators B
	 numgens B
	 isCommutative B
	 coefficientRing B
      Text
	 As we can see, x is an element of B.
      Example
         x
      Text
         If we define a new ring containing x, x is now part of that new ring
      Example
      	 C = skewPolynomialRing(QQ,(-1)_QQ,{x,y,z,w}) 
         x
      Text
         We can 'go back' to B using the command @ TO (use, NCRing) @.
      Example
	 use B
	 x
	 use C
      Text
         We can also create an Ore extension.  First define a @ TO NCRingMap @ with @ TO ncMap @.
      Example
	 sigma = ncMap(C,C,{y,z,w,x})
      Text
         Then call the command @ TO oreExtension @.
      Example
	 D = oreExtension(C,sigma,a)
	 generators D
	 numgens D
   SeeAlso
      "Basic operations on noncommutative algebras"
///

doc ///
  Key
    (generators, NCRing)
  Headline
    The list of algebra generators of an NCRing
  Usage
    gensA = generators A
  Inputs
    A : NCRing
  Outputs
    gensA : List
  Description
    Text
       This function returns the generators of an NCRing as a list.  As usual,
       gens is a synonym for generators.
    Example
       A = QQ{x,y,z}
       generators A
       gens A
///

doc ///
  Key
    (numgens, NCRing)
  Headline
    The number of algebra generators of an NCRing
  Usage
    numgensA = numgens A
  Inputs
    A : NCRing
  Outputs
    numgensA : ZZ
  Description
    Text
       This function returns the number of generators of an NCRing.
    Example
       A = QQ{x,y,z}
       numgens A
///
