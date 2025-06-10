function res = testLong_conPolyZono_linComb
% testLong_conPolyZono_linComb - unit test function for the 
%    linear combination of constrained polynomial zonotopes
%
% Syntax:
%    res = testLong_conPolyZono_linComb()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conPolyZono/linComb

% Authors:       Niklas Kochdumper
% Written:       03-February-2021
% Last update:   06-June-2025 (TL, fixed unit test)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
splits = 4;

% Random Tests --------------------------------------------------------

% define set representations that are tested
sets = {@conPolyZono.generateRandom,...
        @polyZonotope.generateRandom, ...
        @zonotope.generateRandom, ...
        @conZonotope.generateRandom, ...
        @ellipsoid.generateRandom, ...
        @capsule.generateRandom};

% loop over all test cases
for i = 1:2
    
    % generate random constrained polynomial zonotope
    cPZ1 = conPolyZono.generateRandom('Dimension',2);
    
    % loop over all other set representations
    for j = 1:length(sets)
        
        % generate random object of the current set representation
        S2 = sets{j}('Dimension',2);

        % compute linear combination
        cPZ3 = linComb(cPZ1,S2);

        % get random points inside the two conPolyZono objects
        N12 = 10;
        points1 = randPoint(cPZ1,N12,'extreme');
        points2 = randPoint(S2,N12,'extreme');
        
        % compute random combinations of the points
        N3 = 100;
        points3 = zeros(dim(cPZ3),N3);
        
        for k = 1:N3
            ind1 = randi([1,N12]);
            ind2 = randi([1,N12]);
            % compute random point on the line between point 1 and 2
            points3(:,k) = points1(:,ind1) + ...
                           rand()*(points2(:,ind2) - points1(:,ind1));
        end
        
        % concatenate points
        points = [points1,points2,points3];
        
        % check if all points are inside polygon enclosures
        pgon = polygon(cPZ3,splits);
        assertLoop(contains(pgon,points),i,j)
    end
end

% ------------------------------ END OF CODE ------------------------------
