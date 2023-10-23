function res = validateAuthorDateLine(line)
% validateAuthorDateLine - checks if a line in the written/last update/revision 
%    block is valid
%
% Syntax:
%    res = validateAuthorDateLine(line)
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
% Written:       18-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if length(line) < 17 || line(17) ~= ' '
    res = false;
    return
end

regex = ' (([0-9][0-9]-|)(January|February|March|April|May|June|July|August|September|October|November|December)-20[0-9][0-9]( \([A-Z][A-Z](, .*|)\)| \([a-z].*\)|)|---) *$';

startIndex = regexp(line,regex,'once');
res = ~isempty(startIndex);

% ------------------------------ END OF CODE ------------------------------
