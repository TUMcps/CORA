function Z = cartProd(Z,S)
% cartProd - returns the cartesian product of two zonotopes
%
% Syntax:  
%    Z = cartProd(Z,S)
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
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_cartProd('zonotope',Z,S);

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    Z = zonotope(vars{1}); return
else
    % potential re-ordering
    Z = vars{1}; S = vars{2};
end


% first or second set is zonotope
if isa(Z,'zonotope')

    % different cases for different set representations
    if isa(S,'zonotope')
        Z.Z = [[center(Z);center(S)],blkdiag(generators(Z),generators(S))];
    elseif isnumeric(S)
        Z.Z = [[center(Z);S],[generators(Z);zeros(size(S,1),size(Z.Z,2)-1)]];
    elseif isa(S,'interval') 
        Z = cartProd(Z,zonotope(S));
    elseif isa(S,'conZonotope')
        Z = cartProd(conZonotope(Z),S);
    elseif isa(S,'zonoBundle')
        Z = cartProd(zonoBundle(Z),S);
    elseif isa(S,'mptPolytope')
        Z = cartProd(mptPolytope(Z),S);
    elseif isa(S,'polyZonotope')
        Z = cartProd(polyZonotope(Z),S);
    elseif isa(S,'conPolyZono')
        Z = cartProd(conPolyZono(Z),S);
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',Z,S));
    end

elseif isa(S,'zonotope')

    % first argument is a vector
    if isnumeric(Z)
        Z.Z = [[Z;center(S)],[zeros(size(Z,1),size(S.Z,2)-1);generators(S)]];
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',Z,S));
    end  
    
else
    
    % throw error for given arguments
    throw(CORAerror('CORA:noops',Z,S));
    
end  
    
    
    
end

%------------- END OF CODE --------------