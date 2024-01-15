function cPZ = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a
%    constrained polynomial zonotope with other sets or points
%
% Syntax:
%    cPZ = plus(summand1,summand2)
%
% Inputs:
%    summand1 - conPolyZono object, contSet object, or numerical vector
%    summand2 - conPolyZono object, contSet object, or numerical vector
%
% Outputs:
%    cPZ - conPolyZono object
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

% determine which summand is the conPolyZono object
[cPZ,S] = findClassArg(summand1,summand2,'conPolyZono');

try

    % different cases for different set represnetations
    if isnumeric(S) || isa(S,'interval') || ...
       isa(S,'zonotope')
    
       Z = zonotope(S);
       cPZ.c = cPZ.c + Z.c;
       cPZ.GI = [cPZ.GI, Z.G];
       
    elseif isa(S,'conPolyZono') || isa(S,'polytope') || ...
           isa(S,'zonoBundle') || isa(S,'conZonotope') || ...
           isa(S,'ellipsoid') || isa(S,'capsule') || ...
           isa(S,'polyZonotope') || isa(S,'taylm')
        
         % convert to conPolyZono object
         S = conPolyZono(S);
       
         % update states
         cPZ.c = cPZ.c + S.c;
         cPZ.G = [cPZ.G, S.G];
         cPZ.E = blkdiag(cPZ.E,S.E);
         
         cPZ.GI = [cPZ.GI, S.GI];
         
         % update constraints
         cPZ = updateConstraints(cPZ,cPZ,S);
         
    else
         % throw error for given arguments
         throw(CORAerror('CORA:noops',cPZ,S));
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check whether different dimension of ambient space
    equalDimCheck(cPZ,summand);

    % check for empty sets
    if representsa_(cPZ,'emptySet',eps)
        return
    elseif representsa_(summand,'emptySet',eps)
        cPZ = conPolyZono.empty(dim(cPZ)); return
    end

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
