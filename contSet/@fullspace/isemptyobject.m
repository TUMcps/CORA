function res = isemptyobject(fs)
% isemptyobject - checks whether a fullspace object contains any
%    information at all; consequently, the set is interpreted as the empty
%    set 
%
% Syntax:
%    res = isemptyobject(fs)
%
% Inputs:
%    fs - fullspace object
%
% Outputs:
%    res - true/false
%
% Example: 
%    fs = fullspace(2);
%    isemptyobject(fs); % false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       24-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = fs.dimension == 0;

% ------------------------------ END OF CODE ------------------------------
