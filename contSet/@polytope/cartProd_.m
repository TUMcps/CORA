function P_out = cartProd_(P,S,type,varargin)
% cartProd_ - Returns the cartesian product of two polytope objects
% 
% Syntax:
%    P_out = cartProd_(P,S,type,varargin)
%
% Inputs:
%    P - polytope object, numerical vector
%    S - contSet object, numerical vector
%    type - 'exact', 'inner', 'outer'
%
% Outputs:
%    P_out - resulting polytope object
%
% Example: 
%    P1 = polytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    P2 = polytope([1;-1],[3;2]);
%
%    P = cartProd(P1,P2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/cartProd, zonotope/cartProd_

% Authors:       Niklas Kochdumper
% Written:       26-November-2019
% Last update:   05-May-2020 (MW, standardized error message)
%                04-April-2022 (VK, adapted to polytope class)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% note: no re-ordering!

% first or second set is polytope
if isa(P,'polytope')

    % different cases for different set representations
    if isa(S,'polytope')

        % block-concatenation of constraint matrices
        A = blkdiag(P.A,S.A);
        Ae = blkdiag(P.Ae,S.Ae);
        % vertical concatenation of offset vectors
        b = [P.b;S.b];
        be = [P.be;S.be];
        
        P_out = polytope(A,b,Ae,be);
        
    elseif isnumeric(S)

        % ensure that not a Cartesian product with a point cloud
        if size(S,2) > 1
            throw(CORAerror('CORA:notSupported',...
                'Cartesian product only supported for single points.'));
        end
        
        % read out dimension
        n = size(S,1);
        % init equality constraints
        Ae = eye(n);
        be = S;
        
        % instantiate polytope
        P_out = cartProd_(P,polytope(zeros(0,n),[],Ae,be),type);
        % resultng polytope is degenerate
        P_out.fullDim.val = false;

    elseif isa(S,'zonoBundle') || isa(S,'zonotope') || ...
           isa(S,'interval') || isa(S,'conZonotope')              

        P_out = cartProd_(P,polytope(S),type);
        
    elseif isa(S,'polyZonotope')
        
        P_out = cartProd_(polyZonotope(P),S,type);
        
    else
        % throw error for given arguments
        throw(CORAerror('CORAerror:noops',P,S));
    end

else

    % different cases for different set representations
    if isnumeric(P)
        % ... second set has to be a polytope, otherwise this function
        % would not have been called

        % ensure that not a Cartesian product with a point cloud
        if size(P,2) > 1
            throw(CORAerror('CORA:notSupported',...
                'Cartesian product only supported for single points.'));
        end

        % read out dimension
        n = size(P,1);
        % init equality constraints
        Ae = eye(n);
        be = P;
        
        P_out = cartProd(polytope(zeros(0,n),[],Ae,be),S);
        % set is degenerate
        P_out.fullDim.val = false;

    else
        % throw error for given arguments
        throw(CORAerror('CORAerror:noops',P,S));
    end  
end

% ------------------------------ END OF CODE ------------------------------
