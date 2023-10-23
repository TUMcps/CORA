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
% See also: test_docstring

% Authors:       Tobias Ladner
% Written:       17-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get unique chars
chars = unique(line);

% check each unique character
res = true;
for c = 1:length(chars)
    if ~strcmp(chars(c),' ')
        res = false;
        break
    end
end

% ------------------------------ END OF CODE ------------------------------
