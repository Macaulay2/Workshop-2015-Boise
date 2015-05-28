listener = openListener "$:8088"
verbose = true

hexdigits := "0123456789ABCDEF"
hext := new HashTable from for i from 0 to 15 list hexdigits#i => i
hex1 := c -> if hext#?c then hext#c else 0
hex2 = (c,d) -> 16 * hex1 c + hex1 d
toHex1 := asc -> ("%",hexdigits#(asc>>4),hexdigits#(asc&15))
toHex := str -> concatenate apply(ascii str, toHex1)

server = () -> (
    stderr << "listening:" << endl;
    while true do (
        local fun;
        wait {listener};
        g := openInOut listener;				    -- this should be interruptable!
        r := read g;
        if verbose then stderr << "request: " << stack lines r << endl;
        r = lines r;
        if #r == 0 then (close g; continue);
        r = first r;
        if match("^GET /fcn1/",r) then (
            s := first select("^GET /fcn1/(.*) ", "\\1", r);
            fun = fcn1;
            )
	  else if match("^GET /fcn2/(.*) ",r) then (
	       s = first select("^GET /fcn2/(.*) ", "\\1", r);
	       fun = fcn2;
	       )
	  else if match("^HEAD /(.*) ",r) then (
	       s = first select("^HEAD /(.*) ", "\\1", r);
	       fun = identity;
	       )
	  else (
	       s = "";
	       fun = identity;
	       );
	  t := select(".|%[0-9A-F]{2,2}", s);
	  u := apply(t, x -> if #x == 1 then x else ascii hex2(x#1, x#2));
	  u = concatenate u;
      g << httpHeaders fun u << close;
	  )
     )

fcn1 = x -> "called fcn1 on " | x
fcn2 = x -> "called fcn2 on " | x
