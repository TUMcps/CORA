function res = isEmptyCommentLine(line)
% isEmptyCommentLine - checks if the given line is an empty comment line
%
% Syntax:
%    res = isEmptyCommentLine(line)
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

% contains '%'
res = length(chars) >= 1; 

% check each unique character
for c = 1:length(chars)
    if ~strcmp(chars(c),'%') && ~strcmp(chars(c),' ')
        res = false;
        break
    end
end

% ------------------------------ END OF CODE ------------------------------
