function res = isEmptyLine(line)
% isEmptyLine - checks if the given line is empty
%
% Syntax:
%    res = isEmptyLine(line)
%
% Inputs:
%    line - char
%
% Outputs:
%    res = true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: test_codingConventions

% Authors:       Tobias Ladner
% Written:       17-August-2023
% Last update:   06-February-2025 (speed up)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = all(line == ' ');

% ------------------------------ END OF CODE ------------------------------
