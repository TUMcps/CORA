function cPZ = conPolyZono(poly)
% conPolyZono - convert a polytope to a constrained polynomial zonotope
%
% Syntax:  
%    cPZ = conPolyZono(poly)
%
% Inputs:
%    poly - mptPolytope object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    P = mptPolytope([1 1; -2 1; 0 -1],[1;1;1]);
%    cPZ = conPolyZono(P)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mptPolytope, conZonotope

% Author:       Niklas Kochdumper
% Written:      21-January-2020
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    cPZ = conPolyZono(conZonotope(poly));
    
end
    
%------------- END OF CODE --------------