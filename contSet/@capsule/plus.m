function S_out = plus(C,S)
% plus - Overloaded '+' operator for the over-approximative Minkowski 
%    addition of two capsules or the exact translation of a capsule by a 
%    vector
%
% Syntax:
%    S_out = C + S
%    S_out = plus(C,S)
%
% Inputs:
%    C - capsule object, numeric
%    S - contSet object, numeric
%
% Outputs:
%    C_out - Minkowski sum
%
% Example: 
%    C = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    summand1 = C;
%    summand2 = [2; 2; 1];
%    C1 = C + summand1;
%    C2 = C + summand2;
%    figure; hold on;
%    plot(C); plot(C1); plot(C2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       04-March-2019 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[S_out,S] = reorderNumeric(C,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < C.precedence
    S_out = S + C;
    return
end

try

    if isa(S,'capsule')
        % add centers
        S_out.c = S_out.c + S.c;
        % only replace generator if the one of the summand is longer
        length_C = norm(S_out.g);
        length_S = norm(S.g);
        if length_C < length_S
            S_out.g = S.g;
            radiusOfGenerator = length_C;
        else
            radiusOfGenerator = length_S; 
        end
        % obtain new radius
        S_out.r = S_out.r + S.r + radiusOfGenerator;
        return
    end
    
    % numeric: add to center
    if isnumeric(S) && iscolumn(S)
        S_out.c = S_out.c + S;
        return
    end

    if isa(S,'interval')
        % shift center
        S_out.c = S_out.c + center(S);
        % enlarge radius by radius of enclosing hyperball of interval
        S_out.r = S_out.r + radius(S);
        return
    end

    if isa(S,'zonotope')
        % enclose zonotope by interval, same method as for interval
        S = interval(S);
        S_out.c = S_out.c + center(S);
        S_out.r = S_out.r + radius(S);
        return
    end
    
    if representsa_(S,'origin',eps)
        % addition of only the origin (note: for some set representations,
        % this is a bit slow, so we put it at the end of the list...)
        return
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check for empty sets
    if representsa_(S_out,'emptySet',eps) || representsa_(S,'emptySet',eps)
        S_out = capsule.empty(dim(S_out));
        return
    end

    % check whether different dimension of ambient space
    equalDimCheck(S_out,S);

    % other error...
    rethrow(ME);

end

throw(CORAerror('CORA:noops',C,S));

% ------------------------------ END OF CODE ------------------------------
