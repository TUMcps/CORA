classdef fourValued < uint32
% fourValued - four-valued logic
%
% Syntax:
%    phi = fourValued.True & ~ (fourValued.False | fourValued.Unknown | fourValued.Inconclusive)
%    ff = fourValued(0) % false
%    uu = fourValued(1) % unknown
%    ii = fourValued(2) % inconclusive
%    tt = fourValued(3) % true 
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
% See also: none

% Authors:       Florian Lercher
% Written:       09-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
enumeration
    True(3), False(0), Unknown(1), Inconclusive(2)
end

methods
    % conjunction
    function out = and(lhs,rhs)
        if ~isequal(size(lhs),size(rhs))
            throw(CORAerror('CORA:dimensionMismatch',lhs,rhs));
        end

        % true if both are true
        trueMask = lhs == fourValued.True & rhs == fourValued.True;
        out(trueMask) = fourValued.True;

        % false if one is false
        falseMask = lhs == fourValued.False | rhs == fourValued.False;
        out(falseMask) = fourValued.False;

        % unknown if at least one is unknown, but none is inconclusive or false
        unknownMask = ~(trueMask | falseMask | lhs == fourValued.Inconclusive | rhs == fourValued.Inconclusive);
        out(unknownMask) = fourValued.Unknown;

        % rest is inconclusive
        inconclusiveMask = ~(trueMask | falseMask | unknownMask);
        out(inconclusiveMask) = fourValued.Inconclusive;
    end

    % disjunction
    function out = or(lhs,rhs)
        if ~isequal(size(lhs),size(rhs))
            throw(CORAerror('CORA:dimensionMismatch',lhs,rhs));
        end

        % true if one is true
        trueMask = lhs == fourValued.True | rhs == fourValued.True;
        out(trueMask) = fourValued.True;

        % false if both are false
        falseMask = lhs == fourValued.False & rhs == fourValued.False;
        out(falseMask) = fourValued.False;

        % unknown if at least one is unknown, but none is inconclusive or true
        unknownMask = ~(trueMask | falseMask | lhs == fourValued.Inconclusive | rhs == fourValued.Inconclusive);
        out(unknownMask) = fourValued.Unknown;

        % rest is inconclusive
        inconclusiveMask = ~(trueMask | falseMask | unknownMask);
        out(inconclusiveMask) = fourValued.Inconclusive;
    end

    % negation
    function out = not(op)
        out(op == fourValued.True) = fourValued.False;
        out(op == fourValued.False) = fourValued.True;
        out(op == fourValued.Unknown) = fourValued.Unknown;
        out(op == fourValued.Inconclusive) = fourValued.Inconclusive;
    end
end

methods (Static)
    function val = fromBool(b)
        if b
            val = fourValued.True;
        else
            val = fourValued.False;
        end
    end

    function val = fromKleene(k)
        if k == kleene.True
            val = fourValued.True;
        elseif k == kleene.False
            val = fourValued.False;
        else
            val = fourValued.Unknown;
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
