sl3weights=(l) -> (
    v:={};
    for i from 0 to l do (
    for j from i to l do (
    for k from j to l do(
    lambda:= {k,j,i};
    v=append(v, lambda)
    ); ););
    return v
)
