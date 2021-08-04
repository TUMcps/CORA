function Z = cartProd(Z1,Z2)
% cartProd - returns the cartesian product of two zonotopes
%
% Syntax:  
%    Z = cartProd(Z1,Z2)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    zono1 = zonotope.generateRandom(2);
%    zono2 = zonotope.generateRandom(3);
%
%    zono = cartProd(zono1,zono2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      18-May-2011
% Last update:  27-Aug-2019
%               05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

% first or second set is zonotope
if isa(Z1,'zonotope')

    % different cases for different set representations
    if isa(Z2,'zonotope')

        c = [center(Z1);center(Z2)];
        G = blkdiag(generators(Z1),generators(Z2));

        Z = zonotope([c,G]);

    elseif isnumeric(Z2)

        c = [center(Z1);Z2];
        G = [generators(Z1);zeros(size(Z2,1),size(Z1.Z,2)-1)];

        Z = zonotope([c,G]);

    elseif isa(Z2,'interval') 
        Z = cartProd(Z1,zonotope(Z2));
    elseif isa(Z2,'conZonotope')
        Z = cartProd(conZonotope(Z1),Z2);
    elseif isa(Z2,'zonoBundle')
        Z = cartProd(zonoBundle(Z1),Z2);
    elseif isa(Z2,'mptPolytope')
        Z = cartProd(mptPolytope(Z1),Z2);
    elseif isa(Z2,'polyZonotope')
        Z = cartProd(polyZonotope(Z1),Z2);
    elseif isa(Z2,'conPolyZono')
        Z = cartProd(conPolyZono(Z1),Z2);
    else
        % throw error for given arguments
        error(noops(Z1,Z2));
    end

elseif isa(Z2,'zonotope')

    % first argument is a vector
    if isnumeric(Z1)

        c = [Z1;center(Z2)];
        G = [zeros(size(Z1,1),size(Z2.Z,2)-1);generators(Z2)];
        Z = zonotope([c,G]);

    else
        % throw error for given arguments
        error(noops(Z1,Z2));
    end  
    
else
    
    % throw error for given arguments
    error(noops(Z1,Z2));
end  
    
    
    
end

%------------- END OF CODE --------------