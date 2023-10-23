function cPZ = conPolyZono(T)
% conPolyZono - convert a Taylor model to a constrained polynomial zonotope
%
% Syntax:
%    cPZ = conPolyZono(T)
%
% Inputs:
%    T - taylm object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    tay = taylm(interval([-2;1],[1;2]));
%    f = @(x) [sin(x(1))*x(2); cos(x(2))];
%    T = f(tay);
%
%    cPZ = conPolyZono(T)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonoBundle, conZonotope

% Authors:       Niklas Kochdumper
% Written:       21-January-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    cPZ = conPolyZono(polyZonotope(T));
    
end
    
% ------------------------------ END OF CODE ------------------------------
