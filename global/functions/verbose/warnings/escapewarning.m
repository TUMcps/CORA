function str = escapewarning(str)
% escapewarning - escapes the string to be shown in warning
%
% Syntax:
%    escapewarning(str)
%
% Inputs:
%    str - string
%
% Outputs:
%    str - escaped string
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CORAwarning, descapewarning

% Authors:       Tobias Ladner
% Written:       17-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% de-escape
str = strrep(str,'"','&quot;');
str = strrep(str,"'","&apos;");
str = strrep(str,newline,"&newline;");
str = strrep(str,filesep,'/');

end

% ------------------------------ END OF CODE ------------------------------
