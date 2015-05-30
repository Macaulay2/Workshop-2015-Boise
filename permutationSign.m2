permutationSign = method()
permutationSign(List) := L -> (
    if (sum apply(#L-1, i-> number(drop(L,i+1),j-> j<=i)) % 2) === 0 then 1 else -1
    )

end
--Short snippet of code to compute the sign of a permutation of a list L
--Note that this does NOT check if you actually entered a permutation.
L = {3, 1, 2, 5, 0, 7, 6}

permutationSign(L)
permutationSign(L2)
L = {0,1,2,3}
L2 = {1,0,2,3}

