
runMaple = (runMaplePrefix,script) -> (
     filename := temporaryFileName()|currentTime()|".mpl";
     filename << script << endl << close;
     s:=toString(get(concatenate("!",runMaplePrefix," ",filename)))
)

runMaplePrefix="/Library/Frameworks/Maple.framework/Versions/2015/bin/maple";

mapleScript = ///read("/Users/davids/Desktop/Dropbox/Boise/qcalc.mpl"):
with(qcalc):
Gr(3,7):
qtoS(S[2,1]^3);
///
