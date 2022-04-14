function [Zbundle] = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a zonotope
%        bundle with a zonotope or with a vector (see Prop. 2 in [1])
%
% Syntax:  
%    [Zbundle] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - zonotope bundle or zonotope or numerical vector
%    summand2 - zonotope bundle or zonotope or numerical vector
%
% Outputs:
%    Zbundle - Zonotpe bundle after Minkowsi addition
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
% Last update:  05-May-2020 (MW, standarized error message)
% Last revision:---

%------------- BEGIN CODE --------------

    % determine zonotope bundle object
    if isa(summand1,'zonoBundle')
        Zbundle=summand1;
        summand=summand2; 
    elseif isa(summand2,'zonoBundle')
        Zbundle=summand2;
        summand=summand1;  
    end

    % over-approximate the set if the summand is a zonotope bundle
    if isa(summand,'zonoBundle')
       summand = zonotope(interval(summand)); 
    end

    % different cases depending on the class of the second summand
    if isa(summand,'zonotope') || isnumeric(summand) || isa(summand,'interval')

        for i = 1:Zbundle.parallelSets
            Zbundle.Z{i} = Zbundle.Z{i} + summand;
        end

    elseif isa(summand,'mptPolytope') || isa(summand,'conZonotope') 

        for i = 1:Zbundle.parallelSets
            Zbundle.Z{i} = Zbundle.Z{i} + zonoBundle(summand);
        end        

    elseif isa(summand,'polyZonotope') || isa(summand,'conPolyZono')

        Zbundle = summand + Zbundle;

    else
        % throw error for given arguments
        error(noops(summand1,summand2));
    end

%------------- END OF CODE --------------