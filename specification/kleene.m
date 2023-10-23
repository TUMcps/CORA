classdef kleene < uint32
% kleene - three-valued logic
%
% Syntax:
%    phi = kleene.True & ~ (kleene.False | kleene.Unknown)
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

% Authors:       Benedikt Seidl
% Written:       10-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

enumeration
    True(2), False(0), Unknown(1)
end

methods

    % conjunction
    function out = and(lhs, rhs)
        if kleene.True == lhs && kleene.True == rhs
            out = kleene.True;
        elseif kleene.False == lhs || kleene.False == rhs
            out = kleene.False;
        else
            out = kleene.Unknown;
        end
    end

    % disjunction
    function out = or(lhs, rhs)
        if kleene.True == lhs || kleene.True == rhs
            out = kleene.True;
        elseif kleene.False == lhs && kleene.False == rhs
            out = kleene.False;
        else
            out = kleene.Unknown;
        end
    end

    % negation
    function out = not(phi)
        if kleene.True == phi
            out = kleene.False;
        elseif kleene.False == phi
            out = kleene.True;
        else
            out = kleene.Unknown;
        end
    end

end
end

% ------------------------------ END OF CODE ------------------------------
