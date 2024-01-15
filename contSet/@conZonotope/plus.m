function cZ = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a
%    constrained zonotope with other set representations
%
% Syntax:
%    cZ = plus(summand1,summand2)
%
% Inputs:
%    summand1 - conZonotope object or numerical vector
%    summand2 - conZonotope object or numerical vector
%
% Outputs:
%    cZ - conZonotope object
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

% find a conZonotope object
[cZ,summand] = findClassArg(summand1,summand2,'conZonotope');

try

    % handle different classes of the second summand
    if isa(summand,'conZonotope')
        
        % Calculate minkowski sum (Equation (12) in reference paper [1])
        cZ.c=cZ.c+summand.c;
        cZ.G(:,(end+1):(end+size(summand.G,2))) = summand.G;
        
        if isempty(cZ.A)
            if ~isempty(summand.A)
                cZ.A = [zeros(size(summand.A,1), ...
                       size(cZ.G,2)-size(summand.A,2)),summand.A]; 
                cZ.b = summand.b;
            end
        else
            if isempty(summand.A)
                cZ.A = [cZ.A,zeros(size(cZ.A,1),size(summand.G,2))];
            else
                cZ.A = blkdiag(cZ.A, summand.A);
                cZ.b = [cZ.b; summand.b];
            end
        end
        
        cZ.ksi = [];
        cZ.R = [];
        
    elseif isnumeric(summand) 
        
        cZ.c = cZ.c + summand;
        
    elseif isa(summand,'zonotope') || isa(summand,'interval') || ...
           isa(summand,'polytope') || isa(summand,'zonoBundle')
        
        cZ = cZ + conZonotope(summand);
        
    elseif isa(summand,'polyZonotope') || isa(summand,'conPolyZono')
        
        cZ = polyZonotope(cZ) + summand;
     
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

    % check whether different dimension of ambient space
    equalDimCheck(cZ,summand);

    % check for empty sets
    if representsa_(cZ,'emptySet',eps)
        return
    elseif representsa_(summand,'emptySet',eps)
        cZ = conZonotope.empty(dim(cZ)); return
    end    

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
