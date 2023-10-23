function c = center(O)
% center - returns the center of an empty set
%
% Syntax:
%    c = center(O)
%
% Inputs:
%    O - emptySet object
%
% Outputs:
%    c - center
%
% Example: 
%    O = emptySet(2);
%    c = center(O);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

c = double.empty(O.dimension,0);

% ------------------------------ END OF CODE ------------------------------
