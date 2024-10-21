function S_out = plus(cPZ,S)
% plus - Overloaded '+' operator for the Minkowski addition of a
%    constrained polynomial zonotope with other sets or points
%
% Syntax:
%    S_out = cPZ + S
%    S_out = plus(cPZ,S)
%
% Inputs:
%    cPZ - conPolyZono object, numeric
%    S - contSet object, numeric
%
% Outputs:
%    S_out - conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [2 2;2 -1];
%    E = [1 0; 0 1];
%    A = [1 1];
%    b = 0;
%    EC = [2 0; 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
% 
%    E = ellipsoid(0.1*[2 1;1 2]);
% 
%    res = cPZ + E; 
%
%    figure; hold on;
%    plot(res,[1,2],'FaceColor','b','Splits',12);
%    plot(E,[1,2],'FaceColor','g');
%    plot(cPZ,[1,2],'FaceColor','r','Splits',12);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, exactPlus

% Authors:       Niklas Kochdumper
% Written:       14-August-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[S_out,S] = reorderNumeric(cPZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < S_out.precedence
    S_out = S + S_out;
    return
end

try
    % different cases for different set represnetations
    if isa(S,'conPolyZono')
        % update states
        S_out.c = S_out.c + S.c;
        S_out.G = [S_out.G, S.G];
        S_out.E = blkdiag(S_out.E,S.E);
        S_out.GI = [S_out.GI, S.GI];
        
        % update constraints
        S_out = priv_updateConstraints(S_out,S_out,S);
        return
    end

    if isa(S,'zonotope')
        S_out.c = S_out.c + S.c;
        S_out.GI = [S_out.GI, S.G];
        return
    end

    % numeric and interval: convert to zonotope
    if (isnumeric(S) && iscolumn(S)) || isa(S,'interval')
        S = zonotope(S);
        S_out.c = S_out.c + S.c;
        S_out.GI = [S_out.GI, S.G];
        return
    end
       
    if isa(S,'ellipsoid') || isa(S,'capsule') || ...
            isa(S,'polyZonotope') || isa(S,'polytope') || ...
            isa(S,'conZonotope') || isa(S,'zonoBundle') || isa(S,'taylm')
        
         % convert to conPolyZono object
         S = conPolyZono(S);
        
         % update states
         S_out.c = S_out.c + S.c;
         S_out.G = [S_out.G, S.G];
         S_out.E = blkdiag(S_out.E,S.E);
         
         S_out.GI = [S_out.GI, S.GI];
         
         % update constraints
         S_out = priv_updateConstraints(S_out,S_out,S);
         return
         
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check whether different dimension of ambient space
    equalDimCheck(S_out,S);

    % check for empty sets
    if representsa_(S_out,'emptySet',eps) || representsa_(S,'emptySet',eps)
        S_out = conPolyZono.empty(dim(S_out)); return
    end

    % other error...
    rethrow(ME);

end

throw(CORAerror('CORA:noops',S_out,S));

% ------------------------------ END OF CODE ------------------------------
