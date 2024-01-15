function pZ = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a polynomial
%    zonotope with another set representation or a point
%
% Syntax:
%    pZ = plus(summand1,summand2)
%
% Inputs:
%    summand1 - polyZonotope object
%    summand2 - contSet object or numerical vector
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

% determine which summand is the polyZonotope object
[pZ,summand] = findClassArg(summand1,summand2,'polyZonotope');

try

    % summand is a numeric vector
    if isnumeric(summand)
        pZ.c = pZ.c + summand;
        return;
    end
    
    % different cases for the different set representations
    if isa(summand,'zonotope')
        
        pZ.c = pZ.c + center(summand);
        pZ.GI = [pZ.GI, generators(summand)];
       
    elseif isa(summand,'interval')
        
        pZ = pZ + zonotope(summand);
        
    elseif isa(summand,'conPolyZono')
        
        pZ = summand + pZ;
        
    else
        
        % convert other set representations to polynomial zonotopes
        if ~isa(summand,'polyZonotope') 
            if isa(summand,'polytope') || isa(summand,'zonoBundle') || ...
               isa(summand,'conZonotope')
    
                summand = polyZonotope(summand);
            else
                % throw error for given arguments
                throw(CORAerror('CORA:noops',summand1,summand2));
            end
        end
    
        % compute Minkowski sum
        pZ.c = pZ.c + summand.c;
        if isempty(pZ.c)
            pZ = polyZonotope.empty(dim(pZ));
            return
        end
        pZ.G = [pZ.G,summand.G];
        pZ.E = blkdiag(pZ.E,summand.E);
        pZ.GI = [pZ.GI,summand.GI];
        pZ.id = [pZ.id;max(pZ.id) + summand.id];
        
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check whether different dimension of ambient space
    equalDimCheck(pZ,summand);

    % check for empty sets
    if representsa_(pZ,'emptySet',eps)
        return
    elseif representsa_(summand,'emptySet',eps)
        pZ = polyZonotope.empty(dim(pZ)); return
    end

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
