function S_out = plus(cZ,S)
% plus - Overloaded '+' operator for the Minkowski addition of a
%    constrained zonotope with another set or vector
%
% Syntax:
%    S_out = cZ + S
%    S_out = plus(cZ,S)
%
% Inputs:
%    cZ - conZonotope object, numeric
%    S - contSet object, numeric
%
% Outputs:
%    S_out - conZonotope object
%
% Example: 
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZono = conZonotope(Z,A,b);
%    cPlus = cZono + [5;4];
%
%    figure; hold on;
%    plot(cZono,[1,2],'FaceColor','r');
%    plot(cPlus,[1,2],'FaceColor','b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:       05-December-2017 
% Last update:   15-May-2018
%                05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[S_out,S] = reorderNumeric(cZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < S_out.precedence
    S_out = S + S_out;
    return
end

try

    % sum with constrained zonotope
    if isa(S,'conZonotope')
        S_out = aux_plus_conZonotope(S_out,S);
        return
    end
        
    % numeric: add to center
    if isnumeric(S) && iscolumn(S)
        S_out.c = S_out.c + S;
        return
    end
        
    % other set representations: convert to conZonotope
    if isa(S,'zonotope') || isa(S,'interval') || ...
           isa(S,'polytope') || isa(S,'zonoBundle')
        S_out = aux_plus_conZonotope(S_out,conZonotope(S));
        return
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % check whether different dimension of ambient space
    equalDimCheck(S_out,S);

    % check for empty sets
    if representsa_(S_out,'emptySet',eps) || representsa_(S,'emptySet',eps)
        S_out = conZonotope.empty(dim(S_out));
        return
    end

    % other error...
    rethrow(ME);

end

throw(CORAerror('CORA:noops',S_out,S));

end


% Auxiliary functions -----------------------------------------------------

function S_out = aux_plus_conZonotope(S_out,S)
% Equation (12) in reference paper [1]

S_out.c = S_out.c + S.c;
S_out.G(:,(end+1):(end+size(S.G,2))) = S.G;
S_out.A = blkdiag(S_out.A, S.A);
S_out.b = [S_out.b; S.b];

S_out.ksi = [];
S_out.R = [];

end

% ------------------------------ END OF CODE ------------------------------
