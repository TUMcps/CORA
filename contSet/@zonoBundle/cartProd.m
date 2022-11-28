function zB = cartProd(zB,S)
% cartProd - Cartesian product of a zonotope bundle and a set
%
% Syntax:  
%    zB = cartProd(zB,S)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object
%
% Outputs:
%    zB - zonoBundle object
%
% Example: 
%    zB1 = zonoBundle.generateRandom('Dimension',2,'NrZonotopes',3);
%    zB2 = zonoBundle.generateRandom('Dimension',3,'NrZonotopes',2);
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

    % pre-processing
    [res,vars] = pre_cartProd('zonoBundle',zB,S);
    
    % check premature exit
    if res
        % if result has been found, it is stored in the first entry of var
        zB = vars{1}; return
    else
        % potential re-ordering
        zB = vars{1}; S = vars{2};
    end


    % first or second set is zonotope
    if isa(zB,'zonoBundle')

        % different cases for different set representations
        if isa(S,'zonoBundle')

            % compute cartesian product
            list = cell(zB.parallelSets*S.parallelSets,1);
            counter = 1;

            for i = 1:zB.parallelSets
                for j = 1:S.parallelSets
                    list{counter,1} = cartProd(zB.Z{i},S.Z{j});
                    counter = counter + 1;
                end
            end

            zB = zonoBundle(list);

        elseif isnumeric(S) || isa(S,'zonotope') || isa(S,'interval')
            
            for i = 1:zB.parallelSets
                zB.Z{i} = cartProd(zB.Z{i},S); 
            end

        elseif isa(S,'conZonotope') || isa(S,'mptPolytope')
            
            zB = cartProd(zB,zonoBundle(S));
            
        elseif isa(S,'polyZonotope')
            
            zB = cartProd(polyZonotope(zB),S);
            
        elseif isa(S,'conPolyZono')
            
            zB = cartProd(conPolyZono(zB),S);
            
        else
            % throw error for given arguments
            throw(CORAerror('CORA:noops',zB,S));
        end

    else

        % different cases for different set representations
        if isnumeric(zB)

            zB = S;
            for i = 1:S.parallelSets
                zB.Z{i} = cartProd(zB,zB.Z{i}); 
            end

        else
            % throw error for given arguments
            throw(CORAerror('CORA:noops',zB,S));
        end  
    end
end

%------------- END OF CODE --------------