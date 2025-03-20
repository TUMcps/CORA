function S_out = plus(zB,S)
% plus - Overloaded '+' operator for the Minkowski addition of a zonotope
%        bundle with a zonotope or with a vector (see Prop. 2 in [1])
%
% Syntax:
%    S_out = zB + S
%    S_out = plus(zB,S)
%
% Inputs:
%    zB - zonoBundle object, numeric
%    S - contSet object, numeric
%
% Outputs:
%    S_out - zonoBundle object after Minkowski addition
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

% ensure that numeric is second input argument
[S_out,S] = reorderNumeric(zB,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < S_out.precedence
    S_out = S + S_out;
    return
end

try

    % over-approximate the set if the summand is a zonotope bundle
    % S is a zonotope bundle or an interval
    if isa(S,'zonoBundle') || isa(S,'interval')
        S = zonotope(interval(S));
        for i = 1:S_out.parallelSets
            S_out.Z{i} = S_out.Z{i} + S;
        end
        return
    end
     % S is a zonotope
    if isa(S,'zonotope')
        for i = 1:S_out.parallelSets
            S_out.Z{i} = S_out.Z{i} + S;
        end
        return
    end
    % S is a numeric, column vector
    if isnumeric(S) && iscolumn(S)
        for i = 1:S_out.parallelSets
            S_out.Z{i} = S_out.Z{i} + S;
        end
        return
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % check whether different dimension of ambient space
    equalDimCheck(S_out,S);

    % check for empty sets
    if representsa_(S_out,'emptySet',eps) || representsa_(S,'emptySet',eps)
        S_out = zonoBundle.empty(dim(S_out));
        return
    end

    % other error...
    rethrow(ME);

end

throw(CORAerror('CORA:noops',S_out,S));

% ------------------------------ END OF CODE ------------------------------
