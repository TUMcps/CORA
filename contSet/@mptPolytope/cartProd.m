function P = cartProd(P1,S)
% cartProd - Returns the cartesian product of two mptPolytope objects
% 
% Syntax:  
%    P = cartProd(P1,S)
%
% Inputs:
%    P1 - mptPolytope object
%    S - contSet object
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

    % pre-processing
    [res,vars] = pre_cartProd('mptPolytope',P1,S);
    
    % check premature exit
    if res
        % if result has been found, it is stored in the first entry of var
        P = vars{1}; return
    else
        % assign variables
        P1 = vars{1}; S = vars{2};
    end


    % first or second set is polytope
    if isa(P1,'mptPolytope')

        % different cases for different set representations
        if isa(S,'mptPolytope')

            A = blkdiag(P1.P.A,S.P.A);
            b = [P1.P.b;S.P.b];
            
            P = mptPolytope(A,b);
            
        elseif isnumeric(S)
            
            n = size(S,1);
            
            A = [eye(n);-eye(n)];
            b = [S;-S];
            
            P = cartProd(P1,mptPolytope(A,b));

        elseif isa(S,'zonoBundle') || isa(S,'zonotope') || ...
               isa(S,'interval') || isa(S,'conZonotope')              

            P = cartProd(P1,mptPolytope(S));
            
        elseif isa(S,'polyZonotope')
            
            P = cartProd(polyZonotope(P1),S);
            
        else
            % throw error for given arguments
            throw(CORAerror('CORA:noops',P1,S));
        end

    else

        % different cases for different set representations
        if isnumeric(P1)

            n = size(P1,1);
            
            A = [eye(n);-eye(n)];
            b = [P1;-P1];
            
            P = cartProd(mptPolytope(A,b),S);

        else
            % throw error for given arguments
            throw(CORAerror('CORA:noops',P1,S));
        end  
    end
end

%------------- END OF CODE --------------