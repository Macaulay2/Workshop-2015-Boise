runMaple = (runMaplePrefix,script) -> (
     filename := temporaryFileName()|currentTime()|".mpl";
     filename << script << endl << close;
     toString(get(concatenate("!",runMaplePrefix," ",filename)))
)

getTerms = s -> (
     f := replace(" ","",first select(";\n *(.*)","\\1",s));
     h := for t in select(///S\[([0-9,]*)\]///,"\\1",f) list (
     	 value("{"|t|"}") => select(///([0-9]*q?)S\[///|t|///\]///,"\\1",f)
     );
     new HashTable from h
)

maple="/usr/bin/maple";
qcalclib="/Users/corey/Dropbox/Maple/qcalc.mpl";
mapleScript = ///read("///|qcalclib|///"): with(qcalc): Gr(///|r|","|n|"): qtoS("|expr|");"

callqcalc = (r,n,expr) -> (
    mapleScript = ///read("///|qcalclib|///"): with(qcalc): Gr(///|r|","|n|"): "|expr|";";
    getTerms(runMaple(maple,mapleScript))   
)

-- try something like : callqcalc(3,7,"qtoS(S[2,1]^4")