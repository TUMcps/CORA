function cPZ = empty(n)
% empty - instantiates an empty constrained polynomial zonotope
%
% Syntax:
%    cPZ = conPolyZono.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    cPZ - empty constrained polynomial zonotope
%
% Example: 
%    cPZ = conPolyZono.empty(2);
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

if n <= 0
    throw(CORAerror('CORA:wrongValue','first','positive'));
end

cPZ = conPolyZono(zeros(n,0));

% ------------------------------ END OF CODE ------------------------------
