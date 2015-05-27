
runMaple = (script) -> (
     filename := temporaryFileName()|currentTime()|".mpl";
     filename << script << endl << close;
     s:=toString(get("!"|maple|" "|filename))
)

mapleScript = ///read("/Users/davids/Desktop/Dropbox/Boise/qcalc.mpl"):
with(qcalc):
Gr(3,7):
qtoS(S[2,1]^3);

///