function spec = add(spec1,spec2)
% add - joins two specification objects
%
% Syntax:
%    spec = add(spec1,spec2)
%
% Inputs:
%    spec1 - specification object
%    spec2 - specification object
%
% Outputs:
%    spec - resulting specification object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: specification

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       29-May-2020             
% Last update:   30-April-2023 (MW, massive simplification)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

spec = [spec1;spec2];

% ------------------------------ END OF CODE ------------------------------
