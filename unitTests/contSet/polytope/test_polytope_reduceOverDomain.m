function res = test_polytope_reduceOverDomain()
% test_polytope_reduceOverDomain - unit test function for reduction of
%    inequality constraints over a given domain
%
% Syntax:
%    res = test_polytope_reduceOverDomain()
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
% See also: none

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       24-November-2022
% Last update:   23-September-2024 (MW, move & refactor)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% loop over different dimensions
for j = 2:5

    % generate list of random halfspaces that are similar
    P_base = polytope.generateRandom('Dimension',j,'NrConstraints',1);
    P_base = normalizeConstraints(P_base);
    
    A = interval(P_base.A' - 0.1*abs(P_base.A'), ...
                 P_base.A' + 0.1*abs(P_base.A'));
    b = interval(P_base.b  - 0.1*abs(P_base.b), ...
                 P_base.b  + 0.1*abs(P_base.b));
    
    nrCon = 10;
    P_A = zeros(nrCon,j); P_b = zeros(nrCon,1);
    for i = 1:nrCon
        P_A(i,:) = randPoint(A)';
        P_b(i) = randPoint(b);
    end
    P = polytope(P_A,P_b);
    
    % generate random domain that intersects the halfspaces
    dom = interval.generateRandom('Dimension',j);
    
    P_base_eq = polytope([],[],P_base.A,P_base.b);
    if ~isIntersecting(P_base_eq,polytope(dom))
        offset = interval(P_base.A * zonotope(dom));
        dom = dom + (-center(dom)) + (P_base.A') * randPoint(offset); 
    end
    
    % group the halfspaces together
    P_inner = reduceOverDomain(P,dom,'inner');
    P_outer = reduceOverDomain(P,dom,'outer');
    
    % check the result for the inner-approximation
    P_intersection = P_inner & polytope(dom);
    
    if ~representsa(P_intersection,'emptySet')
        points = randPoint(P_intersection,10);
        for i = 1:nrCon
            assertLoop(contains(polytope(P_A(i,:),P_b(i)),points),j,i)
        end
    end
    
    % check the result for the outer-approximation
    points = [];
    
    for i = 1:nrCon
        P_intersection = polytope(P_A(i,:),P_b(i)) & polytope(dom);
        if ~representsa(P_intersection,'emptySet')
            points = [points, randPoint(P_intersection,10)];
        end
    end
    
    assertLoop(isempty(points) || all(contains(P_outer,points)),j)
end

% ------------------------------ END OF CODE ------------------------------
