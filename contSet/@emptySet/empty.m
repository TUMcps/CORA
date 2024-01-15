function O = empty(n)
% empty - instantiates an empty emptySet
%
% Syntax:
%    O = emptySet.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    O - empty emptySet object
%
% Example: 
%    O = emptySet.empty(2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       09-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call constructor
O = emptySet(n);

% ------------------------------ END OF CODE ------------------------------
