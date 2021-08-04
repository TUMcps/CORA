function cPZ = conPolyZono(pZ)
% conPolyZono - convert a polynomial zonotope to a constrained polynomial
%               zonotope
%
% Syntax:  
%    cPZ = conPolyZono(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    pZ = polyZonotope.generateRandom(2);
%    cPZ = conPolyZono(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope

% Author:       Niklas Kochdumper
% Written:      21-January-2020
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    cPZ = conPolyZono(pZ.c,pZ.G,pZ.expMat,pZ.Grest,pZ.id);
    
end
    
%------------- END OF CODE --------------