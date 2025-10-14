function pZ = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%    interval matrix with a polynomial zonotope
%
% Syntax:
%    pZ = factor1 * factor2
%    pZ = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - polyZonotope object, numeric matrix or scalar
%    factor2 - polyZonotope object, numeric scalar
%
% Outputs:
%    pZ - polyZonotope after the multiplication
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%    matrix = [1 2;-1 3];
%    intMatrix = interval([0.9 1.9;-1.1 2.9],[1.1 2.1;-0.9 3.1]);
%       
%    pZres = matrix*pZ;
%    pZresInt = intMatrix*pZ;
%   
%    figure; hold on;
%    plot(pZ,[1,2],'r');
%    plot(pZres,[1,2],'b');
%    plot(pZresInt,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus, zonotope/mtimes, polyZonotope/matMap

% Authors:       Niklas Kochdumper
% Written:       25-June-2018 
% Last update:   ---
% Last revision: 04-October-2024 (MW, remove InferiorClasses)

% ------------------------------ BEGIN CODE -------------------------------

try
    % matrix/scalar * polynomial zonotope
    if isnumeric(factor1)
        c = factor1*factor2.c;
        G = factor1*factor2.G;
        GI = factor1*factor2.GI;
        pZ = polyZonotope(c,G,GI,factor2.E,factor2.id);
        return
    end

    % polynomial zonotope * scalar
    % (note that polynomial zonotope * matrix is not supported)
    if isnumeric(factor2) && isscalar(factor2)
        c = factor2*factor1.c;
        G = factor2*factor1.G;
        GI = factor2*factor1.GI;
        pZ = polyZonotope(c,G,GI,factor1.E,factor1.id);
        return
    end

catch ME
    % check whether different dimension of ambient space
    equalDimCheck(factor1,factor2);
    rethrow(ME);
end

throw(CORAerror('CORA:noops',factor1,factor2));

% ------------------------------ END OF CODE ------------------------------
