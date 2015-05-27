
runMaple = (script) -> (
     filename := temporaryFileName()|currentTime()|".mpl";
     filename << script << endl << close;
     s:=toString(get("!"|maple|" "|filename))
)

