function res = plus(summand1,summand2)
% plus - Overloaded '+' operator for intervals
%
% Syntax:
%    res = plus(summand1,summand2)
%
% Inputs:
%    summand1 - interval object
%    summand2 - interval object
%
% Outputs:
%    res - interval
%
% Example:
%    summand1 = interval([-2;1],[3;2]);
%    summand2 = interval([0;1],[1;2]);
%    summand1 + summand2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff
% Written:       19-June-2015
% Last update:   23-June-2015
%                10-August-2016
%                24-August-2016
%                05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine the interval object
[res,summand] = findClassArg(summand1,summand2,'interval');

try

    % different cases depending on the class of the summand
    if isa(summand,'interval')
    
        res.inf = res.inf + summand.inf;
        res.sup = res.sup + summand.sup;
    
    elseif isnumeric(summand)
    
        res.inf = res.inf + summand;
        res.sup = res.sup + summand;    
    
    elseif isa(summand,'zonotope') || isa(summand,'conZonotope') || ...
           isa(summand,'zonoBundle') || isa(summand,'polyZonotope') || ...
           isa(summand,'polytope') || isa(summand,'conPolyZono')
    
        res = summand + res;
    
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
    if representsa_(res,'emptySet',eps)
        return
    elseif representsa_(summand,'emptySet',eps)
        res = interval.empty(dim(res)); return
    end

    % check whether different dimension of ambient space
    equalDimCheck(res,summand);

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
