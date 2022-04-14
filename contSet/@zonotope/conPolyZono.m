function cPZ = conPolyZono(Z)
% conPolyZono - convert a zonotope to a constrained polynomial zonotope
%
% Syntax:  
%    cPZ = conPolyZono(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    Z = zonotope.generateRandom(2);
%    cPZ = conPolyZono(Z)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope, polyZonotope

% Author:       Niklas Kochdumper
% Written:      21-January-2020
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    cPZ = conPolyZono(polyZonotope(Z));
    
end
    
%------------- END OF CODE --------------