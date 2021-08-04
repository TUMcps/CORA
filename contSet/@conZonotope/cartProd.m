function cZ = cartProd(cZ1,cZ2)
% cartProd - Returns the cartesian product of two conZonotope objects
% 
% Syntax:  
%    cZ = cartProd(cZ1,cZ2)
%
% Inputs:
%    cZ1 - conZonotope object
%    cZ2 - conZonotope object
%
% Outputs:
%    cZ - resulting conZonotope object
%
% Example: 
%    Z = [0 1 2];
%    A = [1 1];
%    b = 1;
%    cZ = conZonotope(Z,A,b);
%    zono = zonotope([0 1]);
%
%    cZcart = cartProd(cZ,zono);
%
%    plot(cZcart,[1,2],'r','Filled',true);
%    xlim([0.5 2.5]);
%    ylim([-1.5 1.5]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/cartProd

% Author:       Niklas Kochdumper
% Written:      10-August-2018
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

    % first or second set is constrained zonotope
    if isa(cZ1,'conZonotope')

        % different cases for different set representations
        if isa(cZ2,'conZonotope')

            % new center vector
            c = [cZ1.Z(:,1); cZ2.Z(:,1)];

            % new generator matrix
            G = blkdiag(cZ1.Z(:,2:end),cZ2.Z(:,2:end));

            % new constraint matrix
            h1 = size(cZ1.A,1);
            h2 = size(cZ2.A,1);

            m1 = size(cZ1.Z,2)-1;
            m2 = size(cZ2.Z,2)-1;

            if isempty(cZ1.A)
               if isempty(cZ2.A)
                  A = []; 
               else
                  A = [zeros(h2,m1),cZ2.A];
               end
            else
               if isempty(cZ2.A)
                  A = [cZ1.A,zeros(h1,m2)];
               else
                  A = [[cZ1.A,zeros(h1,m2)];[zeros(h2,m1),cZ2.A]];
               end
            end

            % new constraint offset
            b = [cZ1.b;cZ2.b];

            % generate resulting constrained zonotope
            cZ = conZonotope([c,G],A,b);

        elseif isnumeric(cZ2) || isa(cZ2,'zonotope') || ...
               isa(cZ2,'interval') || isa(cZ2,'mptPolytope') || ...
               isa(cZ2,'zonoBundle')

            cZ = cartProd(cZ1,conZonotope(cZ2));
            
        elseif isa(cZ2,'polyZonotope') || isa(cZ2,'conPolyZono')
            
            cZ = cartProd(polyZonotope(cZ1),cZ2);
            
        else
            % throw error for given arguments
            error(noops(cZ1,cZ2));
        end

    else

        % different cases for different set representations
        if isnumeric(cZ1)

            cZ = cartProd(conZonotope(cZ1),cZ2);

        else
            % throw error for given arguments
            error(noops(cZ1,cZ2));
        end  
    end
end

%------------- END OF CODE --------------