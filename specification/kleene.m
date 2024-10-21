classdef kleene < uint32
% kleene - three-valued logic
%
% Syntax:
%    phi = kleene.True & ~ (kleene.False | kleene.Unknown)
%    ff = kleene(0) % kleene.False
%    uu = kleene(1) % kleene.Unknown
%    tt = kleene(2) % kleene.True
%
% Inputs:
%    none
%
% Outputs:
%    none
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Benedikt Seidl, Florian Lercher
% Written:       10-August-2022
% Last update:   09-February-2024 (FL, make operators work with arrays and add fromBool)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

enumeration
    True(2), False(0), Unknown(1)
end

methods

    % conjunction
    function out = and(lhs, rhs)
        if ~isequal(size(lhs), size(rhs))
            throw(CORAerror('CORA:dimensionMismatch', lhs, rhs));
        end

        % true if both are true
        trueMask = (lhs == kleene.True) & (rhs == kleene.True);
        out(trueMask) = kleene.True;
        
        % false if one is false
        falseMask = (lhs == kleene.False) | (rhs == kleene.False);
        out(falseMask) = kleene.False;

        % rest is unknown
        out(~trueMask & ~falseMask) = kleene.Unknown;
    end

    % disjunction
    function out = or(lhs, rhs)
        if ~isequal(size(lhs), size(rhs))
            throw(CORAerror('CORA:dimensionMismatch', lhs, rhs));
        end

        % true if one is true
        trueMask = (lhs == kleene.True) | (rhs == kleene.True);
        out(trueMask) = kleene.True;
        
        % false if both are false
        falseMask = (lhs == kleene.False) & (rhs == kleene.False);
        out(falseMask) = kleene.False;

        % rest is unknown
        out(~trueMask & ~falseMask) = kleene.Unknown;
    end

    % negation
    function out = not(phi)
        out(phi == kleene.True) = kleene.False;
        out(phi == kleene.False) = kleene.True;
        out(phi == kleene.Unknown) = kleene.Unknown;
    end
end

methods (Static)
    function val = fromBool(b)
        if b
            val = kleene.True;
        else
            val = kleene.False;
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
