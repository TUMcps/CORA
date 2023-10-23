function cZ = cartProd_(cZ,S,varargin)
% cartProd_ - Returns the Cartesian product of a constrained zonotope and
%    other set representations or points
% 
% Syntax:
%    cZ = cartProd_(cZ,S)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    Z = [0 1 2];
%    A = [1 1]; b = 1;
%    cZ = conZonotope(Z,A,b);
%    Z = zonotope([0 1]);
%    cZcart = cartProd(cZ,Z);
%
%    figure; hold on; xlim([0.5 2.5]); ylim([-1.5 1.5]);
%    plot(cZcart,[1,2],'FaceColor','r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/cartProd, zonotope/cartProd_

% Authors:       Niklas Kochdumper
% Written:       10-August-2018
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: 27-March-2023 (MW, rename cartProd_)

% ------------------------------ BEGIN CODE -------------------------------

% first or second set is constrained zonotope
if isa(cZ,'conZonotope')

    % different cases for different set representations
    if isa(S,'conZonotope')

        % new center vector
        c = [cZ.c; S.c];

        % new generator matrix
        G = blkdiag(cZ.G,S.G);

        % new constraint matrix
        h1 = size(cZ.A,1);
        h2 = size(S.A,1);

        m1 = size(cZ.G,2);
        m2 = size(S.G,2);

        if isempty(cZ.A)
           if isempty(S.A)
              A = []; 
           else
              A = [zeros(h2,m1),S.A];
           end
        else
           if isempty(S.A)
              A = [cZ.A,zeros(h1,m2)];
           else
              A = [[cZ.A,zeros(h1,m2)];[zeros(h2,m1),S.A]];
           end
        end

        % new constraint offset
        b = [cZ.b;S.b];

        % generate resulting constrained zonotope
        cZ = conZonotope([c,G],A,b);

    elseif isnumeric(S) || isa(S,'zonotope') || ...
           isa(S,'interval') || isa(S,'polytope') || ...
           isa(S,'zonoBundle')

        cZ = cartProd_(cZ,conZonotope(S),'exact');
        
    elseif isa(S,'polyZonotope') || isa(S,'conPolyZono')
        
        cZ = cartProd_(polyZonotope(cZ),S,'exact');
        
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',cZ,S));
    end

else

    % different cases for different set representations
    if isnumeric(cZ)

        cZ = cartProd_(conZonotope(cZ),S,'exact');

    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',cZ,S));
    end
    
end

% ------------------------------ END OF CODE ------------------------------
