function str = bracketSubs(str)
% bracketSubs - substitute 'L' and 'R' by opening/closing parenthesis
%
% Syntax:
%    str = bracketSubs(str)
%
% Inputs:
%    str - string
%
% Outputs:
%    str - string
%
% Example:
%    str = 'xL1R';
%    str = bracketSubs(str);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ???
% Written:       ---
% Last update:   01-May-2020 (MW, added header)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate left and right brackets
str = strrep(str,'L','(');
str = strrep(str,'R',')');

% ------------------------------ END OF CODE ------------------------------
