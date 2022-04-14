function [Zbundle] = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%          interval matrix with a zonotope bundle (see Prop. 1 in [1])
%
% Syntax:  
%    [Zbundle] = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - zonotope bundle or matrix set or matrix
%    factor2 - zonotope bundle or matrix set or matrix 
%
% Outputs:
%    Zbundle - Zonotope bundle after multiplication
%
% References:
%    [1] M. Althoff. "Zonotope bundles for the efficient computation of 
%        reachable sets", 2011
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      09-November-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % determine zonotope bundle object
    if isa(factor1,'zonoBundle')
        Zbundle = factor1;
        factor = factor2;
    elseif isa(factor2,'zonoBundle')
        Zbundle = factor2;
        factor = factor1;  
    end

    % calculate multiplication for each zonotope
    for i=1:Zbundle.parallelSets
        Zbundle.Z{i} = factor*Zbundle.Z{i};
    end
end

%------------- END OF CODE --------------