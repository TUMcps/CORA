function cPZ = conPolyZono(I)
% conPolyZono - convert an interval to a constrained polynomial zonotope
%
% Syntax:
%    cPZ = conPolyZono(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    I = interval([-2;1],[2;4]);
%    cPZ = conPolyZono(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, zonotope, polyZonotope

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

cPZ = conPolyZono(polyZonotope(I));
    
% ------------------------------ END OF CODE ------------------------------
