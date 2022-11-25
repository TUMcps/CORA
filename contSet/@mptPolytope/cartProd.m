function P = cartProd(P1,P2)
% cartProd - Returns the cartesian product of two mptPolytope objects
% 
% Syntax:  
%    P = cartProd(P1,P2)
%
% Inputs:
%    P1 - mptPolytope object
%    P2 - mptPolytope object
%
% Outputs:
%    P - resulting mptPolytope object
%
% Example: 
%    poly1 = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    poly2 = mptPolytope([1;-1],[3;2]);
%
%    poly = cartProd(poly1,poly2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/cartProd

% Author:       Niklas Kochdumper
% Written:      26-November-2019
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

    % first or second set is polytope
    if isa(P1,'mptPolytope')

        % different cases for different set representations
        if isa(P2,'mptPolytope')

            A = blkdiag(P1.P.A,P2.P.A);
            b = [P1.P.b;P2.P.b];
            
            P = mptPolytope(A,b);
            
        elseif isnumeric(P2)
            
            n = size(P2,1);
            
            A = [eye(n);-eye(n)];
            b = [P2;-P2];
            
            P = cartProd(P1,mptPolytope(A,b));

        elseif isa(P2,'zonoBundle') || isa(P2,'zonotope') || ...
               isa(P2,'interval') || isa(P2,'conZonotope')              

            P = cartProd(P1,mptPolytope(P2));
            
        elseif isa(P2,'polyZonotope')
            
            P = cartProd(polyZonotope(P1),P2);
            
        else
            % throw error for given arguments
            error(noops(P1,P2));
        end

    else

        % different cases for different set representations
        if isnumeric(P1)

            n = size(P1,1);
            
            A = [eye(n);-eye(n)];
            b = [P1;-P1];
            
            P = cartProd(mptPolytope(A,b),P2);

        else
            % throw error for given arguments
            error(noops(P1,P2));
        end  
    end
end

%------------- END OF CODE --------------