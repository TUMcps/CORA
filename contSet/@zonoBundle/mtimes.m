function zB = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%    interval matrix with a zonotope bundle (see Prop. 1 in [1])
%
% Syntax:
%    zB = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - zonoBundle object or matrix set or matrix
%    factor2 - zonoBundle object or matrix set or matrix 
%
% Outputs:
%    zB - zonoBundle object
%
% References:
%    [1] M. Althoff. "Zonotope bundles for the efficient computation of 
%        reachable sets", 2011
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/mtimes

% Authors:       Matthias Althoff
% Written:       09-November-2010
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine zonotope bundle object
[zB,factor] = findClassArg(factor1,factor2,'zonoBundle');

try

% calculate multiplication for each zonotope
for i=1:zB.parallelSets
    zB.Z{i} = factor*zB.Z{i};
end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check for empty sets
    if representsa_(zB,'emptySet',eps)
        return
    end

    % check whether different dimension of ambient space
    equalDimCheck(factor1,factor2);

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
