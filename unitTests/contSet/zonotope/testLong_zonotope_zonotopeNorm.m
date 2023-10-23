function res = testLong_zonotope_zonotopeNorm
% testLong_zonotope_zonotopeNorm - unit test function of zonotopeNorm
%
% Syntax:
%    res = testLong_zonotope_zonotopeNorm
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
% See also: -

% Authors:       Adrian Kulmburg
% Written:       06-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

dims = [2 5 10]; % Dimensions to be tested
Ntests = 10; % Number of tests in each case
Npoints = 10; % Number of points to be tested

%% Points inside
for n = dims
    for i = 1:Ntests
        Z = zonotope.generateRandom('Dimension',n);
        P = randPoint(Z,Npoints);
        for i_p = 1:size(P,2)
            zN = zonotopeNorm(Z,P(:,i_p)-Z.center);
            if zN > 1 && ~withinTol(zN,1)
                throw(CORAerror('CORA:testFailed'));
            end
        end
    end
end

%% Points outside
for n = dims
    for i = 1:Ntests
        Z = zonotope.generateRandom('Dimension',n);
        % Lift Z to a space with one more dimension, that way we're sure
        % that the points cannot be contained if we choose them right
        Z_lift = zonotope([Z.center;0], [Z.generators;zeros(1,size(Z.generators,2))]);
        
        P = randPoint(Z,Npoints);
        % Same lift for P, but we displace all points by 10, so that
        % containment has now way of being able to work
        P_lift = [P;zeros(1,size(P,2))] + [zeros(n,1);10];
        
        for i_p = 1:size(P_lift,2)
            zN = zonotopeNorm(Z_lift, P_lift(:,i_p)-Z_lift.center);
            if zN < 1 && ~withinTol(zN,1)
                throw(CORAerror('CORA:testFailed'));
            end
        end
    end
end

%% Check norm-properties
for n = dims
    for i=1:Ntests
        Z = zonotope.generateRandom('Dimension',n);
        
        % sample random points
        p1 = randPoint(Z);
        p2 = randPoint(Z);
        
        % Check triangle inequality
        zN12 = zonotopeNorm(Z,p1+p2);
        zN1 = zonotopeNorm(Z,p1);
        zN2 = zonotopeNorm(Z,p2);
        if ~( zN12 < zN1+zN2 || withinTol(zN12,zN1+zN2) )
            throw(CORAerror('CORA:testFailed'));
        end
        
        % Check symmetry
        zN1 = zonotopeNorm(Z,p1);
        zN1_ = zonotopeNorm(Z,-p1);
        if ~withinTol(zN1,zN1_)
            throw(CORAerror('CORA:testFailed'));
        end
        
        % Check part of the positive definiteness
        zN = zonotopeNorm(Z, zeros(n,1));
        if ~( zN < 0 || withinTol(zN,0) )
            throw(CORAerror('CORA:testFailed'));
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
