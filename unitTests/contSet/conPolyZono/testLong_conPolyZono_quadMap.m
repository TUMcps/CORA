function res = testLong_conPolyZono_quadMap
% testLong_conPolyZono_quadMap - unit test function for 
%    quadratic map of constrained polynomial zonotopes
%
% Syntax:
%    res = testLong_conPolyZono_quadMap()
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
% See also: conPolyZono/or

% Authors:       Niklas Kochdumper
% Written:       03-February-2021
% Last update:   ---
% Last revision: 28-March-2025 (TL, fixed test and randPoint can return [])

% ------------------------------ BEGIN CODE -------------------------------

res = true;
splits = 4;

% Test 1: test quadratic map for single set

% loop over all test cases
for i = 1:10
    
    % generate random constrained polynomial zonotopes
    cPZ1 = conPolyZono.generateRandom('Dimension',2);
    
    % generate random matrices
    Q{1} = -2 + 4*rand(2);
    Q{2} = -2 + 4*rand(2);

    % compute quadratic map
    cPZ = quadMap(cPZ1,Q);

    % compute quadractic map for random points inside the set
    N = 10;
    points = zeros(dim(cPZ),N);

    % try to sample N points
    cnt = 0;
    for k = 1:2*N
        % sample point
        p = randPoint(cPZ1,1,'extreme');

        % randPoint might return [] if no point was found
        if isempty(p)
            continue
        end

        % store
        cnt = cnt + 1;
        points(:,cnt) = quadMapPoint(p,p,Q);

        if cnt == N
            break
        end
    end
    
    % check if all points are inside polygon enclosures
    pgon = polygon(cPZ,splits);

    assertLoop(contains(pgon,points),i)
end


% Test 2: test mixed quadratic map involving two conPolyZono objects

% loop over all test cases
for i = 1:10
    
    % generate random constrained polynomial zonotopes
    cPZ1 = conPolyZono.generateRandom('Dimension',2);
    cPZ2 = conPolyZono.generateRandom('Dimension',2);
    
    % generate random matrices
    Q{1} = -2 + 4*rand(2);
    Q{2} = -2 + 4*rand(2);

    % compute quadratic map
    cPZ = quadMap(cPZ1,cPZ2,Q);

    % get random points inside the two conPolyZono objects
    N = 10;
    points1 = randPoint(cPZ1,N,'extreme');
    points2 = randPoint(cPZ2,N,'extreme');

    % add the points
    points = zeros(dim(cPZ),N^2);
    cnt = 1;

    for k = 1:N
        for l = 1:N
            points(:,cnt) = quadMapPoint(points1(:,k),points2(:,l),Q);
            cnt = cnt + 1;
        end
    end
    
    % check if all points are inside polygon enclosures
    pgon = polygon(cPZ,splits);

    assertLoop(contains(pgon,points),i)
end

% ------------------------------ END OF CODE ------------------------------
