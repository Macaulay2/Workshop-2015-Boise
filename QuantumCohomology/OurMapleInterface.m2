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
r=3
n=7
expr="S[2,1]^4"
mapleScript = ///read("///|qcalclib|///"): with(qcalc): Gr(///|r|","|n|"): qtoS("|expr|");"
{*
mapleScript = ///read("/Users/corey/Dropbox/Maple/qcalc.mpl"):
with(qcalc):
Gr(3,7):
qtoS(S[2,1]^4);
///
*}

getTerms(runMaple(maple,mapleScript))
