


runMaple = (runMaplePrefix,script) -> (
    filename := temporaryFileName()|currentTime()|".mpl";
    filename << script << endl << close;
    toString(get(concatenate("!",runMaplePrefix," ",filename)))
);

{*
parseCoefficient = (f,t,Y) -> (
    if t=="" then return value (first select(///([0-9]*q?)S\[///|t|///\]///,"\\1",f));
    t=first select(///([0-9]*q?)S\[///|t|///\]///,"\\1",f);
    if t=="" then return 1_Y;
    if t=="q" then return Y_0;
    if last(t)=="q" then return value( concatenate(substring(0,#t-1,t),"*q") );
    return (1_Y)*(value t);    
);
*}

stringToTerms = (f) -> (
    answer:={};    
    i:=0;
    while i < #f do (
        if f_i == "+" then (
	    answer=append(answer,substring(0,i,f));
	    f=substring(i+1,#f,f);
	    i=0;
	) else i=i+1;
    if i==#f then answer=append(answer,f)
    );
    return answer
);

f="q^2*S[1,1]+2*S[2]+S[]"

parseTerm = (t,Y) -> (
    i:=0;
    while i < #f do (
        if t_i == "S" then break else i=i+1   
    );   
    coeff:=substring(0,i,t);
    if coeff=="" then coeff="1";
    yt:=concatenate("{",substring(i+2,#t-i-3,t),"}");
    if last(coeff)=="*" then coeff=substring(0,#coeff-1,coeff);
    {value yt,(1_Y)*(value coeff)}    
)



{*
getTerms = (s,Y) -> (
    f := replace(" ","",first select(";\n *(.*)","\\1",s));
    T := select(///S\[([0-9,]*)\]///,"\\1",f);
    h := for t in T list (
	if t=="" then {} => parseCoefficient(f,t,Y) else value("{"|t|"}") => parseCoefficient(f,t,Y)    
    );
    new HashTable from h
);
*}

getTerms = (s,Y) -> (
    f := replace(" ","",first select(";\n *(.*)","\\1",s));
    T := stringToTerms f;
    new HashTable from apply(T, t -> parseTerm(t,Y))
)


maple="/Library/Frameworks/Maple.framework/Versions/2015/bin/maple";
qcalclib="/Users/davids/Desktop/Dropbox/Boise/qcalc.mpl";
--mapleScript = ///read("///|qcalclib|///"): with(qcalc): Gr(///|r|","|n|"): qtoS("|expr|");"

--Note: changed on May 28 to match pieriFirstSum and pieriSecondSum
callqcalc = (r,l,expr,Y) -> (
    mapleScript = ///read("///|qcalclib|///"): with(qcalc): Gr(///|toString(r+1)|","|toString(l+r+1)|"): lprint(qtoS("|expr|"));";
    getTerms(runMaple(maple,mapleScript),Y)   
)

callqcalcNoParse = (r,l,expr) -> (
    mapleScript = ///read("///|qcalclib|///"): with(qcalc): Gr(///|toString(r+1)|","|toString(l+r+1)|"): lprint(qtoS("|expr|"));";
    runMaple(maple,mapleScript)   
)

qcalcMonomialMultiplication = (r,l,yt1,yt2,Y) -> (
    syt1:=toString(yt1);
    syt1=substring(1,#syt1-2,syt1);
    syt2:=toString(yt2); 
    syt2=substring(1,#syt2-2,syt2);
    s:=concatenate("S[",syt1,"]*S[",syt2,"]");
    callqcalc(r,l,s,Y)
)




end









<<<<<<< HEAD
-- try something like : callqcalc(3,7,"S[2,1]^4")
Using interface:
o7 = HashTable{{2, 2, 1} => {17q}}
               {3, 1, 1} => {21q}
               {3, 2} => {31q}
               {4, 1} => {17q}
               {4, 4, 4} => {8}

o7 : HashTable

From Maple
> qtoS(S[2,1]^4);
31 q S[3, 2] + 17 q S[4, 1] + 17 q S[2, 2, 1] + 21 q S[3, 1, 1] + 8 S[4, 4, 4]


Test pieri:
r=2, l=4, S[3]*S[3,2,1]

load "OurMapleInterface-Dave.m2";
callqcalc(2,4,"S[3]*S[3,2,1]")
o2 = HashTable{{1, 1} => {q}  }
               {2} => {q}
               {4, 3, 2} => {1}
pieriFirstSum(3,2,4,{3,2,1})
pieriSecondSum(3,2,4,{3,2,1})


callqcalc(5,8,"S[4]*S[6,4,2,2]")
S[8, 6, 2, 2] + S[8, 5, 3, 2] + S[8, 4, 4, 2] + S[7, 6, 3, 2] + S[7, 5, 4, 2]

     + S[6, 6, 4, 2] + S[8, 5, 2, 2, 1] + S[8, 4, 3, 2, 1] + S[7, 6, 2, 2, 1]

     + S[7, 5, 3, 2, 1] + S[7, 4, 4, 2, 1] + S[6, 6, 3, 2, 1]

     + S[6, 5, 4, 2, 1] + S[8, 4, 2, 2, 2] + S[7, 5, 2, 2, 2]

     + S[7, 4, 3, 2, 2] + S[6, 6, 2, 2, 2] + S[6, 5, 3, 2, 2]

     + S[6, 4, 4, 2, 2]

callqcalc(5,8,"S[4]*S[6,4,2,2,1,1]")
q S[3, 2, 1] + q S[4, 1, 1] + q S[3, 1, 1, 1] + S[6, 5, 4, 2, 2, 1]

     + S[6, 6, 3, 2, 2, 1] + S[6, 6, 4, 2, 1, 1] + S[7, 4, 4, 2, 2, 1]

     + S[7, 5, 3, 2, 2, 1] + S[7, 5, 4, 2, 1, 1] + S[7, 6, 2, 2, 2, 1]

     + S[7, 6, 3, 2, 1, 1] + S[8, 4, 3, 2, 2, 1] + S[8, 4, 4, 2, 1, 1]

     + S[8, 5, 2, 2, 2, 1] + S[8, 5, 3, 2, 1, 1] + S[8, 6, 2, 2, 1, 1]



=======
-- try something like : callqcalc(3,7,"S[2,1]^4")
>>>>>>> 6d5b33610119565fc8187456392fb5ba0a13cb63
