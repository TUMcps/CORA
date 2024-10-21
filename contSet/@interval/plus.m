function S_out = plus(I,S)
% plus - overloaded '+' operator for the Minkowski sum of an interval and
%    another set or point
%
% Syntax:
%    S_out = I + S
%    S_out = plus(I,S)
%
% Inputs:
%    I - interval object, numeric
%    S - contSet object, numeric
%
% Outputs:
%    S_out - Minkowski sum
%
% Example:
%    I = interval([-2;1],[3;2]);
%    S = interval([0;1],[1;2]);
%    I + S
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       19-June-2015
% Last update:   23-June-2015
%                10-August-2016
%                24-August-2016
%                05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[S_out,S] = reorderNumeric(I,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < S_out.precedence
    S_out = S + S_out;
    return
end

try

    % interval-interval case
    if isa(S,'interval')
        S_out.inf = S_out.inf + S.inf;
        S_out.sup = S_out.sup + S.sup;
        return
    end
    
    % numeric vector/matrix
    if isnumeric(S)
        S_out.inf = S_out.inf + S;
        S_out.sup = S_out.sup + S;
        return
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % check whether different dimension of ambient space
    equalDimCheck(S_out,S);

    % check for empty sets
    if representsa_(S_out,'emptySet',eps) || representsa_(S,'emptySet',eps)
        S_out = interval.empty(dim(S_out));
        return
    end

    % other error...
    rethrow(ME);

end

throw(CORAerror('CORA:noops',S_out,S));

% ------------------------------ END OF CODE ------------------------------
