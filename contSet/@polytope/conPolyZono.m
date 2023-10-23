function cPZ = conPolyZono(P)
% conPolyZono - convert a polytope to a constrained polynomial zonotope
%
% Syntax:
%    cPZ = conPolyZono(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    A = [1 2; -2 3; -3 0; -1 -4; 1 -3; 3 -1];
%    b = ones(6,1);
%    P = polytope(A,b);
%    cPZ = conPolyZono(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polytope, conZonotope

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call function from conZonotope class
cPZ = conPolyZono(conZonotope(P));
    
% ------------------------------ END OF CODE ------------------------------
