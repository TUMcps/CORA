function res = testLong_mptPolytope_contains
% testLong_mptPolytope_contains - unit test function for
%    containment check
%
% Syntax:  
%    res = testLong_mptPolytope_contains()
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

% Author:       Victor Gassmann
% Written:      26-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

for i = 2:2:4
   
    % create random polytope
    P = mptPolytope.generateRandom('Dimension',i);
    
    % check if randomly generated points are inside
    Y = randPoint(P,2*i);
    if ~contains(P,Y,'exact',1e-12)
        res = false;
        break;
    end
    
    % generate zonotope
    Z = zonotope.generateRandom('Dimension',i);
    Ei = ellipsoid(Z,'inner:norm');
    Pz = mptPolytope(Z);
    % check if Ei is in Pz
    if ~contains(Pz,Ei,'exact',Ei.TOL)
        res = false;
        break;
    end
end

%------------- END OF CODE --------------