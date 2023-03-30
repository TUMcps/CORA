function Z = cartProd_(Z,S,varargin)
% cartProd_ - computes the Cartesian product of two zonotopes
%
% Syntax:  
%    Z = cartProd_(Z,S)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z1 = zonotope([-1;1],[1 3 2; -3 0 1]);
%    Z2 = zonotope([0;2;-3],[1 4 -2; 2 0 -1; 0 2 2]);
%    Z = cartProd(Z1,Z2);
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
% Last revision:27-March-2023 (MW, rename cartProd_)

%------------- BEGIN CODE --------------

% first or second set is zonotope
if isa(Z,'zonotope')

    % different cases for different set representations
    if isa(S,'zonotope')
        Z.Z = [[center(Z);center(S)],blkdiag(generators(Z),generators(S))];
    elseif isnumeric(S)
        Z.Z = [[center(Z);S],[generators(Z);zeros(size(S,1),size(Z.Z,2)-1)]];
    elseif isa(S,'interval') 
        Z = cartProd_(Z,zonotope(S),'exact');
    elseif isa(S,'conZonotope')
        Z = cartProd_(conZonotope(Z),S,'exact');
    elseif isa(S,'zonoBundle')
        Z = cartProd_(zonoBundle(Z),S,'exact');
    elseif isa(S,'mptPolytope')
        Z = cartProd_(mptPolytope(Z),S,'exact');
    elseif isa(S,'polyZonotope')
        Z = cartProd_(polyZonotope(Z),S,'exact');
    elseif isa(S,'conPolyZono')
        Z = cartProd_(conPolyZono(Z),S,'exact');
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',Z,S));
    end

elseif isa(S,'zonotope')

    % first argument is a vector
    if isnumeric(Z)
        Z = zonotope([Z;center(S)],[zeros(size(Z,1),size(S.Z,2)-1);generators(S)]);
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',Z,S));
    end  
    
else
    
    % throw error for given arguments
    throw(CORAerror('CORA:noops',Z,S));
    
end

%------------- END OF CODE --------------