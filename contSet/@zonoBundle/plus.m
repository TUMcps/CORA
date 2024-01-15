function zB = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a zonotope
%        bundle with a zonotope or with a vector (see Prop. 2 in [1])
%
% Syntax:
%    zB = plus(summand1,summand2)
%
% Inputs:
%    summand1 - zonoBundle object, zonotope object, or numerical vector
%    summand2 - zonoBundle object, zonotope object, or numerical vector
%
% Outputs:
%    zB - zonoBundle object after Minkowski addition
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

% Authors:       Matthias Althoff
% Written:       09-November-2010
% Last update:   05-May-2020 (MW, standarized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine zonotope bundle object
[zB,summand] = findClassArg(summand1,summand2,'zonoBundle');

try

    % over-approximate the set if the summand is a zonotope bundle
    if isa(summand,'zonoBundle')
       summand = zonotope(interval(summand)); 
    end
    
    % different cases depending on the class of the second summand
    if isa(summand,'zonotope') || isnumeric(summand) || isa(summand,'interval')
    
        for i = 1:zB.parallelSets
            zB.Z{i} = zB.Z{i} + summand;
        end

    elseif isa(summand,'polytope') || isa(summand,'conZonotope') 
    
        for i = 1:zB.parallelSets
            zB.Z{i} = zB.Z{i} + zonoBundle(summand);
        end        
    
    elseif isa(summand,'polyZonotope') || isa(summand,'conPolyZono')
    
        zB = summand + zB;
    
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',summand1,summand2));
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check whether different dimension of ambient space
    equalDimCheck(zB,summand);

    % check for empty sets
    if representsa_(zB,'emptySet',eps)
        return
    elseif representsa_(summand,'emptySet',eps)
        zB = zonoBundle.empty(dim(zB)); return
    end

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
