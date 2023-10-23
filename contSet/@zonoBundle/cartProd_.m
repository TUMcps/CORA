function zB = cartProd_(zB,S,varargin)
% cartProd_ - Cartesian product of a zonotope bundle and a set
%
% Syntax:
%    zB = cartProd_(zB,S)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object
%
% Outputs:
%    zB - zonoBundle object
%
% Example: 
%    Z1 = zonotope([1;-1],[2 -1 3; 0 1 -1]);
%    Z2 = Z1 + [1;0];
%    zB1 = zonoBundle({Z1,Z2});
%    Z1 = zonotope([0;-2],[1 -3 2 0 1; -2 3 -1 -1 0]);
%    Z2 = Z1 + [-1;0];
%    zB2 = zonoBundle({Z1,Z2});
%    zB = cartProd(zB1,zB2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/cartProd, zonotope/cartProd_

% Authors:       Niklas Kochdumper
% Written:       13-June-2018
% Last update:   24-September-2019 (renaming)
%                05-May-2020 (MW, standardized error message)
% Last revision: 27-March-2023 (MW, rename cartProd_)

% ------------------------------ BEGIN CODE -------------------------------

% first or second set is zonotope
if isa(zB,'zonoBundle')

    % different cases for different set representations
    if isa(S,'zonoBundle')

        % compute cartesian product
        list = cell(zB.parallelSets*S.parallelSets,1);
        counter = 1;

        for i = 1:zB.parallelSets
            for j = 1:S.parallelSets
                list{counter,1} = cartProd_(zB.Z{i},S.Z{j},'exact');
                counter = counter + 1;
            end
        end

        zB = zonoBundle(list);

    elseif isnumeric(S) || isa(S,'zonotope') || isa(S,'interval')
        
        for i = 1:zB.parallelSets
            zB.Z{i} = cartProd_(zB.Z{i},S,'exact'); 
        end

    elseif isa(S,'conZonotope') || isa(S,'polytope')
        
        zB = cartProd_(zB,zonoBundle(S),'exact');
        
    elseif isa(S,'polyZonotope')
        
        zB = cartProd_(polyZonotope(zB),S,'exact');
        
    elseif isa(S,'conPolyZono')
        
        zB = cartProd_(conPolyZono(zB),S,'exact');
        
    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',zB,S));
    end

else

    % different cases for different set representations
    if isnumeric(zB)

        % convert vector to zonotope with center and no generators, then
        % compute Cartesian product with all parallel sets in zonoBundle
        c = zB;
        zB = S;
        for i = 1:zB.parallelSets
            temp = zonotope(c);
            zB.Z{i} = cartProd_(temp,zB.Z{i},'exact');
        end

    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',zB,S));
    end  
end

% ------------------------------ END OF CODE ------------------------------
