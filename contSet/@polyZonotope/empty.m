function pZ = empty(n)
% empty - instantiates an empty polynomial zonotope
%
% Syntax:
%    pZ = polyZonotope.empty(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    pZ - empty polynomial zonotope
%
% Example: 
%    pZ = polyZonotope.empty(2);
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

pZ = polyZonotope(zeros(n,0));

% ------------------------------ END OF CODE ------------------------------
