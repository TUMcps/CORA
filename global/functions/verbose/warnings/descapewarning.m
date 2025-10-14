function str = descapewarning(str)
% descapewarning - de-escapes the string and returns it as warning string
%
% Syntax:
%    descapewarning(str)
%
% Inputs:
%    str - string
%
% Outputs:
%    str - de-escaped string
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CORAwarning, escapewarning

% Authors:       Tobias Ladner
% Written:       17-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% de-escape
str = strrep(str,'&quot;','"');
str = strrep(str,"&apos;","'");
str = strrep(str,"&newline;",newline);
str = sprintf('[\b%s]\b',str);

end

% ------------------------------ END OF CODE ------------------------------
