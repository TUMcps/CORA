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
%    poly = mptPolytope.generateRandom(2);
%    cPZ = conPolyZono(poly)
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