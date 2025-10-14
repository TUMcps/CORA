function res = test_polygon_hausdorffDist
% test_polygon_hausdorffDist - unit test function of hausdorffDist between
%       polygons
%
% Syntax:
%    res = test_polygon_hausdorffDist
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
% See also: polygon/hausdorffDist

% Authors:       Niklas Kochdumper
% Written:       29-July-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    tol = 1e-2;

    % polygon-to-polygon Hausdorff-distance
    pgon1 = polygon([0 0 2 2],[0 2 2 0]);
    pgon2 = polygon([1 1 3 3],[1 3 3 1]);

    dist = hausdorffDist(pgon1,pgon2);

    assert(abs(dist - sqrt(2)) < tol);

    % point-to-polygon Hausdorff-distance
    pgon = polygon([0 0 2 2],[0 2 2 0]);
    
    dist = hausdorffDist(pgon,[3;1]);

    assert(abs(dist - 1) < tol)

    % compare with Hausdorff-distance between polytopes
    P1 = polytope.generateRandom('Dimension',2);
    P2 = polytope.generateRandom('Dimension',2);

    distPoly = hausdorffDist(P1,P2);
    distPgon = hausdorffDist(polygon(P1),polygon(P2));

    assert(abs(distPoly - distPgon) < tol | ...
                                max(distPgon,distPoly) < 100*tol);

    % compare with Hausdorff-distance between poylgon and point
    P = polytope.generateRandom('Dimension',2);
    pgon = polygon(P);

    points = randPoint(interval(pgon),100);

    for i = 1:size(points,2)
        distPoly = hausdorffDist(P,points(:,i));
        distPgon = hausdorffDist(pgon,points(:,i));

        assert(abs(distPoly - distPgon) < tol | ...
                                max(distPgon,distPoly) < 100*tol);
    end

    distPoly = hausdorffDist(P,points);
    distPgon = hausdorffDist(pgon,points);

    assert(abs(distPoly - distPgon) < tol);

    % gather results
    res = true;
end

% ------------------------------ END OF CODE ------------------------------
