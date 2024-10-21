function Z = mtimes(factor1,factor2)
% mtimes - overloaded '*' operator for the multiplication of a matrix or an
%    interval matrix with a zonotope
%
% Syntax:
%    Z = factor1 * factor2
%    Z = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - zonotope object, numeric matrix or scalar
%    factor2 - zonotope object, numeric scalar
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z = zonotope([1 1 0; 0 0 1]);
%    M = [1 2; 1 0];
%    Zmat = M*Z;
% 
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(Zmat,[1,2],'r');
%    
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       30-September-2006 
% Last update:   07-September-2007
%                05-January-2009
%                06-August-2010
%                01-February-2011
%                08-February-2011
%                18-November-2015
% Last revision: 04-October-2024 (MW, remove InferiorClasses from zonotope)

% ------------------------------ BEGIN CODE -------------------------------

try

    % matrix/scalar * zonotope
    if isnumeric(factor1)
        c = factor1*factor2.c;
        G = factor1*factor2.G;
        Z = zonotope(c,G);
        return
    end

    % zonotope * matrix
    if isnumeric(factor2)
        c = factor2*factor1.c;
        G = factor2*factor1.G;
        Z = zonotope(c,G);
        return
    end

catch ME
    % check whether different dimension of ambient space
    equalDimCheck(factor1,factor2);
    rethrow(ME);
end

throw(CORAerror('CORA:noops',factor1,factor2));

% ------------------------------ END OF CODE ------------------------------
