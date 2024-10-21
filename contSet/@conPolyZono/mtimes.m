function cPZ = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix (or
%    interval matrix) with a constrained polynomial zonotope
%
% Syntax:
%    cPZ = factor1 * factor2
%    cPZ = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - conPolyZono object, numeric matrix or scalar
%    factor2 - conPolyZono object, numeric scalar 
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [2 0 2; 0 2 2];
%    E = [1 0 1; 0 1 1; 0 0 0];
%    A = [2 2 4 -4];
%    b = 0;
%    EC = [1 0 1 0; 0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
% 
%    M = [3 1;2 -4];
%    cPZ_ = M * cPZ;
% 
%    figure; hold on;
%    plot(cPZ,[1,2],'r');
%    plot(cPZ_,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/mtimes, plus

% Authors:       Niklas Kochdumper
% Written:       15-May-2018
% Last update:   ---
% Last revision: 04-October-2024 (MW, remove InferiorClasses from conPolyZono)

% ------------------------------ BEGIN CODE -------------------------------

try

    % matrix/scalar * constrained polynomial zonotope
    if isnumeric(factor1)
        c = factor1*factor2.c;
        G = factor1*factor2.G;
        GI = factor1*factor2.GI;
        cPZ = conPolyZono(c,G,factor2.E,factor2.A,factor2.b, ...
            factor2.EC,GI,factor2.id);
        return
    end

    % constrained polynomial zonotope * scalar
    % (note that constrained polynomial zonotope * matrix is not supported)
    if isnumeric(factor2) && isscalar(factor2)
        c = factor2*factor1.c;
        G = factor2*factor1.G;
        GI = factor2*factor1.GI;
        cPZ = conPolyZono(c,G,factor1.E,factor1.A,factor1.b, ...
            factor1.EC,GI,factor1.id);
        return
    end
    
    % use polynomial zonotope method for other cases
    pZ = polyZonotope(factor2.c,factor2.G,factor2.GI,factor2.E);
    pZ_mapped = mtimes(factor1,pZ);
    cPZ = conPolyZono(pZ_mapped.c,pZ_mapped.G,pZ_mapped.E,factor2.A,factor2.b, ...
                      factor2.EC,pZ_mapped.GI,factor2.id);
    return

catch ME
    % check whether different dimension of ambient space
    equalDimCheck(factor1,factor2);
    rethrow(ME);
end

throw(CORAerror('CORA:noops',factor1,factor2));

% ------------------------------ END OF CODE ------------------------------
