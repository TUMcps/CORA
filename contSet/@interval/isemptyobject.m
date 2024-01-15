function res = isemptyobject(I)
% isemptyobject - checks whether an interval contains any information at
%    all; consequently, the set is interpreted as the empty set 
%
% Syntax:
%    res = isemptyobject(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - true/false
%
% Example: 
%    I = interval(-1,0);
%    isemptyobject(I); % false
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

res = isnumeric(I.inf) && isempty(I.inf) ...
        && isnumeric(I.sup) && isempty(I.sup);

% ------------------------------ END OF CODE ------------------------------
