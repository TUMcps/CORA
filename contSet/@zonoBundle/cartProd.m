function zB = cartProd(zB1,zB2)
% cartProd - Returns the cartesian product of a zonoBundle with
%                    a zonotope
%
% Syntax:  
%    zB = cartProd(zB1,zB2)
%
% Inputs:
%    zB1 - zonoBundle object
%    zB2 - zonotope object
%
% Outputs:
%    zB - zonoBundle object
%
% Example: 
%    zB1 = zonoBundle.generateRandom(2,3);
%    zB2 = zonoBundle.generateRandom(3,2);
%
%    zB = cartProd(zB1,zB2);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/cartProd

% Author:       Niklas Kochdumper
% Written:      13-June-2018
% Last update:  24-Sep-2019 (renaming)
%               05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------

    % first or second set is zonotope
    if isa(zB1,'zonoBundle')

        % different cases for different set representations
        if isa(zB2,'zonoBundle')

            % compute cartesian product
            list = cell(zB1.parallelSets*zB2.parallelSets,1);
            counter = 1;

            for i = 1:zB1.parallelSets
                for j = 1:zB2.parallelSets
                    list{counter,1} = cartProd(zB1.Z{i},zB2.Z{j});
                    counter = counter + 1;
                end
            end

            zB = zonoBundle(list);

        elseif isnumeric(zB2) || isa(zB2,'zonotope') || isa(zB2,'interval')

            zB = zB1;
            
            for i = 1:zB.parallelSets
               zB.Z{i} = cartProd(zB.Z{i},zB2); 
            end

        elseif isa(zB2,'conZonotope') || isa(zB2,'mptPolytope')
            
            zB = cartProd(zB1,zonoBundle(zB2));
            
        elseif isa(zB2,'polyZonotope')
            
            zB = cartProd(polyZonotope(zB1),zB2);
            
        elseif isa(zB2,'conPolyZono')
            
            zB = cartProd(conPolyZono(zB1),zB2);
            
        else
            % throw error for given arguments
            error(noops(zB1,zB2));
        end

    else

        % different cases for different set representations
        if isnumeric(zB1)

            zB = zB2;
            
            for i = 1:zB2.parallelSets
               zB.Z{i} = cartProd(zB1,zB.Z{i}); 
            end

        else
            % throw error for given arguments
            error(noops(zB1,zB2));
        end  
    end
end

%------------- END OF CODE --------------