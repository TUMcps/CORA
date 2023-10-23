function fs = not(O)
% not - overloads '~' operator to compute the complement of an empty set,
%    resulting in a full-dimensional space
%
% Syntax:
%    fs = ~O;
%    fs = not(O);
%
% Inputs:
%    O - emptySet object
%
% Outputs:
%    fs - fullspace object
%
% Example: 
%    O = emptySet(2);
%    fs = ~O;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       07-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% complement is an empty set of same dimension
fs = fullspace(O.dimension);

% ------------------------------ END OF CODE ------------------------------
