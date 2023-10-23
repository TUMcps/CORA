function O = not(fs)
% not - overloads '~' operator to compute the complement of a
%    full-dimensional space, resulting in an empty set
%
% Syntax:
%    O = ~fs;
%    O = not(fs);
%
% Inputs:
%    fs - fullspace object
%
% Outputs:
%    O - emptySet object
%
% Example: 
%    fs = fullspace(2);
%    O = ~fs;
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
O = emptySet(fs.dimension);

% ------------------------------ END OF CODE ------------------------------
