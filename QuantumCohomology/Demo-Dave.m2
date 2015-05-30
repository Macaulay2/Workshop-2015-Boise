restart
loadPackage("QuantumCohomology");
loadPackage("SchurRings");

--Classical example 
R=qcRing(2,3,"s","q")
time s_{1}^6

--Compare to SchurRings
R=schurRing(QQ,t,3);
time t_{1}^6

--Quantum example
R=qcRing(5,6,"s","q")
time s_{3,3,1}^6*s_{6}

--Compare to Anders
--do what you need to read qcalc
with(qcalc);
with(CodeTools);
Usage(qtoS(S[3,3,1]^6*S[6]));
