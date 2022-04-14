function s = str2symbolic( str )
 %takes:    string array/cell array of strings
 %returns:  symbolic array

 %replace newlines to suppress warnings
 str = regexprep(str, "\n", " ");
 
 try 
     % call str2sym (introduced R2017b) if possible
     s = str2sym(str);
 catch
     % try older function if failed
     disp("STR2SYM CRASHED for input: " + str);
     s = sym(str);
 end
end
