function S_out = plus(pZ,S)
% plus - Overloaded '+' operator for the Minkowski addition of a polynomial
%    zonotope with another set representation or a point
%
% Syntax:
%    pZ = pZ + S
%    pZ = plus(pZ,S)
%
% Inputs:
%    pZ - polyZonotope object, numeric
%    S - contSet object, numeric
%
% Outputs:
%    pZ - polyZonotope object after Minkowski addition
%
% Example: 
%    pZ = polyZonotope([0;0],[2 1 2;0 2 2],[],[1 0 3;0 1 1]);
%    zono = zonotope([6 0.2 0 0.2;0 0 0.2 -0.2]);
%   
%    pZsum = pZ + zono;
%
%    figure; hold on
%    plot(pZ,[1,2],'FaceColor','r');
%    plot(zono,[1,2],'FaceColor','b');
%
%    figure
%    plot(pZsum,[1,2],'FaceColor',[0 0.5 0]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, zonotope/plus

% Authors:       Niklas Kochdumper
% Written:       26-March-2018 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[S_out,S] = reorderNumeric(pZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < S_out.precedence
    S_out = S + S_out;
    return
end

try

    if isa(S,'polyZonotope')
        % compute Minkowski sum
        S_out.c = S_out.c + S.c;
        if isempty(S_out.c)
            S_out = polyZonotope.empty(dim(S_out));
            return
        end
        S_out.G = [S_out.G,S.G];
        S_out.E = blkdiag(S_out.E,S.E);
        S_out.GI = [S_out.GI,S.GI];
        max_id = max(S_out.id);
        if isempty(max_id)
            max_id = 0;
        end
        S_out.id = [S_out.id; max_id + S.id];
        return
    end

    % summand is a numeric vector
    if isnumeric(S) && iscolumn(S)
        S_out.c = S_out.c + S;
        return;
    end
    
    % different cases for the different set representations
    if isa(S,'zonotope')
        S_out.c = S_out.c + center(S);
        S_out.GI = [S_out.GI, generators(S)];
        return
    end
       
    if isa(S,'interval')
        S = zonotope(S);
        S_out.c = S_out.c + center(S);
        S_out.GI = [S_out.GI, generators(S)];
        return
    end
        
    % convert other set representations to polynomial zonotopes
    if isa(S,'polytope') || isa(S,'zonoBundle') || isa(S,'conZonotope')
        S = polyZonotope(S);
        S_out = S_out + S;
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
        S_out = polyZonotope.empty(dim(S_out));
        return
    end

    % other error...
    rethrow(ME);

end

throw(CORAerror('CORA:noops',S_out,S));

% ------------------------------ END OF CODE ------------------------------
