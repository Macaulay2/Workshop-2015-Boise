R = CC[a11,a12,a13,a21,a32];
m = {a11 - a12 - a13 - a32,
  a11*a12 + a11*a13 - a12*a13 + a12*a21 + a11*a32 - a13*a32,
  a11*a12*a13 + a12*a13*a21 + a11*a13*a32 + a13*a21*a32,
  a12 + a13 + a32, a12*a13 + a13*a32}
v11 = random(CC)
v12 = random(CC)
v13 = random(CC)
v21 = random(CC)
v32 = random(CC)
p0 = m#0
p1 = m#1
p2 = m#2
p3 = m#3
p4 = m#4
c0 = p0 - sub(p0, {a11 => v11, a12 => v12, a13 => v13, a21 => v21, a32 => v32} )
c1 = p1 - sub(p1, {a11 => v11, a12 => v12, a13 => v13, a21 => v21, a32 => v32} )
c2 = p2 - sub(p2, {a11 => v11, a12 => v12, a13 => v13, a21 => v21, a32 => v32} )
c3 = p3 - sub(p3, {a11 => v11, a12 => v12, a13 => v13, a21 => v21, a32 => v32} )
c4 = p4 - sub(p4, {a11 => v11, a12 => v12, a13 => v13, a21 => v21, a32 => v32} )
f = {c0, c1, c2, c3, c4}
print f
loadPackage("PHCpack")
sols = solveSystem(f)
#sols
sols/print
