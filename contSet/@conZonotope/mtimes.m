function [cZ] = mtimes(factor1,factor2)
% mtimes - Overloaded '.*' operator for the multiplication of a matrix or an
%          interval matrix with a constrained zonotope
%
% Syntax:  
%    [cZ] = times(matrix,cZ)
%
% Inputs:
%    matrix - numerical or interval matrix
%    cZ - conZonotope object 
%
% Outputs:
%    cZ - constrained zonotpe after multiplication with a matrix
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%    cMul = [3 1;2 4] * cZono;
%
%    hold on
%    plot(cZono,[1,2],'r');
%    plot(cMul,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Niklas Kochdumper
% Written:      15-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% Call superclass method
if ~isnumeric(factor1)
    error('conZontope/mtimes: operation not implemented yet!')
else
    cZ = mtimes@zonotope(factor1,factor2);
end

%------------- END OF CODE --------------