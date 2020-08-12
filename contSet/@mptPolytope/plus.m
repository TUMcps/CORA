function [P] = plus(summand1,summand2)
% plus - overloaded '+' operator for the addition of a vector to a
% mptPolytope; Minkowski addition is also supported
%
% Syntax:  
%    [P] = plus(summand1,summand2)
%
% Inputs:
%    summand1 - mptPolytope object or numerical vector
%    summand2 - mptPolytope object or numerical vector
%
% Outputs:
%    P - mptPolytope object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff
% Written:      20-October-2010
% Last update:  24-August-2016
%               05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

    % find a polytope object
    if isa(summand1,'mptPolytope')
        P=summand1;
        summand=summand2;  
    elseif isa(summand2,'mptPolytope')
        P=summand2;
        summand=summand1;  
    end

    % different cases depending on the class of the summand
    if isnumeric(summand)

        try % MPT V3
            P.P = P.P + summand;
        catch % MPT V2
            n = length(summand);
            P.P = range(P.P, eye(n), summand);
        end

    elseif isa(summand,'mptPolytope')

        P.P = P.P + summand.P;

    elseif isa(summand,'zonotope') || isa(summand,'interval') || ...
           isa(summand,'conZonotope') || isa(summand,'zonoBundle')

        P = P + mptPolytope(summand);

    elseif isa(summand,'polyZonotope')

        P = summand + P;

    else
        % throw error for given arguments
        error(noops(summand1,summand2));
    end

%------------- END OF CODE --------------