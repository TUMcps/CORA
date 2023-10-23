function s = str2symbolic(str)
% str2symbolic - Convert strings of equations into symbolic equations
%
% Syntax:
%    s = str2symbolic(str)
%
% Inputs:
%    str - string array/cell array of strings
%
% Outputs:
%    s - symbolic array
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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

% ------------------------------ END OF CODE ------------------------------
