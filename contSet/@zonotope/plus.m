function S_out = plus(Z,S)
% plus - overloaded '+' operator for the Minkowski addition of a zonotope
%    with another set or vector
%
% Syntax:
%    S_out = Z + S
%    S_out = plus(Z,S)
%
% Inputs:
%    Z - zonotope object, numeric
%    S - contSet object, numeric
%
% Outputs:
%    S_out - zonotope after Minkowski addition
%
% Example: 
%    Z = zonotope([1;0],[1 0; 0 1]);
%    Z1 = Z + Z;
%    Z2 = Z + [2;-1];
%
%    figure; hold on;
%    plot(Z,[1,2],'b');
%    plot(Z1,[1,2],'r');
%    plot(Z2,[1,2],'g');
%
% References:
%    [1] M. Althoff. "Reachability analysis and its application to the 
%        safety assessment of autonomous cars", 2010
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       30-September-2006 
% Last update:   23-March-2007
%                14-August-2016
%                04-March-2019
%                13-August-2019
%                14-February-2024 (MW, prevent sum with row vectors/matrices)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[S_out,S] = reorderNumeric(Z,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < S_out.precedence
    S_out = S + S_out;
    return
end

try

    % different cases depending on the class of the second summand
    if isa(S,'zonotope')
        % see Equation 2.1 in [1]
        S_out.c = S_out.c + S.c;
        if isempty(S_out.c)
            S_out = zonotope.empty(dim(S_out));
            return
        end
        S_out.G(:,(end+1):(end+size(S.G,2))) = S.G;
        return
    end
    
    % numeric has to be a scalar or a column vector of correct size
    if isnumeric(S) && iscolumn(S)
        S_out.c = S_out.c + S;
        return
    end
    
    if isa(S,'interval')
        S = zonotope(S);
        S_out.c = S_out.c + S.c;
        S_out.G(:,(end+1):(end+size(S.G,2))) = S.G;
        return
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % check whether different dimension of ambient space
    equalDimCheck(S_out,S);

    % check for empty sets
    if representsa_(S_out,'emptySet',eps) || representsa_(S,'emptySet',1e-10)
        S_out = zonotope.empty(dim(S_out));
        return
    end

    % other error...
    rethrow(ME);

end

throw(CORAerror('CORA:noops',S_out,S));

% ------------------------------ END OF CODE ------------------------------
