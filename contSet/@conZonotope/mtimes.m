function cZ = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or
%    an interval matrix with a constrained zonotope
%
% Syntax:
%    cZ = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - conZonotope object, numerical or interval matrix
%    factor2 - conZonotope object, numerical or interval matrix 
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%    res = [3 1;2 4] * cZ;
%
%    figure; hold on;
%    plot(cZ,[1,2],'r');
%    plot(res,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Niklas Kochdumper
% Written:       15-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% find the conZonotope object
[cZ,factor2] = findClassArg(factor1,factor2,'conZonotope');

% Call superclass method
if ~isnumeric(factor2)
    throw(CORAerror('CORA:noops',cZ,factor2));
else
    cZ.c = factor2*cZ.c;
    cZ.G = factor2*cZ.G;
end

% ------------------------------ END OF CODE ------------------------------
