function P = plus(summand1,summand2)
% plus - overloaded '+' operator for the addition of a vector to a
% mptPolytope; Minkowski addition is also supported
%
% Syntax:  
%    P = plus(summand1,summand2)
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

% determine zonotope object
[P,summand] = findClassArg(summand1,summand2,'mptPolytope');

try

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
        throw(CORAerror('CORA:noops',summand1,summand2));
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check for empty sets
    if isempty(P)
        return
    elseif isemptyobject(summand)
        P = mptPolytope(); return
    end

    % check whether different dimension of ambient space
    equalDimCheck(P,summand);

    % other error...
    rethrow(ME);

end

%------------- END OF CODE --------------