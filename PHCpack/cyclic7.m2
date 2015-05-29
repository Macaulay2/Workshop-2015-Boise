-- run this script on the command line as time M2 cyclic7.m2
-- changing the numThreads option will change the execution time
R = CC[z0,z1,z2,z3,z4,z5,z6]; 
c7 = { z0 + z1 + z2 + z3 + z4 + z5 + z6,
 z0*z1 + z1*z2 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z0,
 z0*z1*z2 + z1*z2*z3 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z0 + z6*z0*z1,
 z0*z1*z2*z3 + z1*z2*z3*z4 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z0
+ z5*z6*z0*z1 + z6*z0*z1*z2,
 z0*z1*z2*z3*z4 + z1*z2*z3*z4*z5 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z0 
+ z4*z5*z6*z0*z1 + z5*z6*z0*z1*z2 + z6*z0*z1*z2*z3,
 z0*z1*z2*z3*z4*z5 + z1*z2*z3*z4*z5*z6 + z2*z3*z4*z5*z6*z0 + z3*z4*z5*z6*z0*z1
+ z4*z5*z6*z0*z1*z2 + z5*z6*z0*z1*z2*z3 + z6*z0*z1*z2*z3*z4,
 z0*z1*z2*z3*z4*z5*z6 - 1
}
loadPackage("PHCpack");
(mv,c7q,c7qsols) = mixedVolume(c7,StartSystem=>true)
c7sols = trackPaths(c7,c7q,c7qsols,numThreads=>4)
stdio << "tracked " << mv << " solution paths" << endl;
exit()
