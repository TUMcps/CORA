function [str,len] = centerString(str,width)
% centerString - centers a string
%
% Syntax:
%    [str,len] = centerString(str,width)
%
% Inputs:
%    str - char, string
%    width - max column width
%
% Outputs:
%    str - centered string
%    len - determined length of input string
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CORAtable

% Authors:       Tobias Ladner
% Written:       05-March-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% compute length of string
if startsWith(str,'<a')
    % only determine length for shown text of hyperlinks
    displayedStr = regexprep(str, '<a.*?>(.*?)</a>', '$1');
else
    displayedStr = str;
end
len = numel(displayedStr);

% centers a string by padding it with white spaces
str = [ ...
    blanks(ceil((width - len)/2)), ...
    str,...
    blanks(floor((width - len)/2)) ...
];

% ------------------------------ END OF CODE ------------------------------
