function res = test_ellipsoid_intersectStrip
% test_ellipsoid_intersectStrip - unit test function of intersetStrip
%
% Syntax:
%    res = test_ellipsoid_intersectStrip
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       29-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% simple example
C = [1 0; 0 1; 1 1];
phi = [5; 3; 3];
y = [-2; 2; 2];

Z = [1 2 2 2 6 2 8; 1 2 2 0 5 0 6];
A = [1 0 1 0 0 0];
b = 1;
cZ = conZonotope(Z,A,b);
res_zono = intersectStrip(cZ,C,phi,y);
    
res = true;

end

% ------------------------------ END OF CODE ------------------------------
