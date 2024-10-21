function zB = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%    interval matrix with a zonotope bundle (see Prop. 1 in [1])
%
% Syntax:
%    zB = factor1 * factor2
%    zB = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - zonoBundle object, numeric matrix or scalar
%    factor2 - zonoBundle object, numeric scalar
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       09-November-2010
% Last update:   ---
% Last revision: 04-October-2024 (MW, remove InferiorClasses from zonootpe)

% ------------------------------ BEGIN CODE -------------------------------

try

    % matrix/scalar * zonoBundle
    if isnumeric(factor1)
        % calculate multiplication for each zonotope, see [1]
        list = cell(factor2.parallelSets,1);
        for i=1:factor2.parallelSets
            list{i} = factor1*factor2.Z{i};
        end
        zB = zonoBundle(list);
        return
    end

    % zonoBundle * scalar
    % (note that zonoBundle * matrix is not supported)
    if isnumeric(factor2) && isscalar(factor2)
        % calculate multiplication for each zonotope, see [1]
        list = cell(factor1.parallelSets,1);
        for i=1:factor1.parallelSets
            list{i} = factor1.Z{i}*factor2;
        end
        zB = zonoBundle(list);
        return
    end

catch ME
    % check whether different dimension of ambient space
    equalDimCheck(factor1,factor2);
    rethrow(ME);
end

throw(CORAerror('CORA:noops',factor1,factor2));

% ------------------------------ END OF CODE ------------------------------
