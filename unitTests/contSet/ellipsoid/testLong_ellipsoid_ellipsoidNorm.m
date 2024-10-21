function res = testLong_ellipsoid_ellipsoidNorm
% testLong_ellipsoid_ellipsoidNorm - unit test function of
%    ellipsoidNorm
%
% Syntax:
%    res = testLong_ellipsoid_ellipsoidNorm
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

tol = 100*eps; % Tolerance to make sure that all comparisons run smoothly

dims = [2 5 10]; % Dimensions to be tested
Ntests = 10; % Number of tests in each case
Npoints = 10; % Number of points to be tested

%% Points inside
for dim = dims
    for i = 1:Ntests
        E = ellipsoid.generateRandom('Dimension',dim,'IsDegenerate',false);
        P = randPoint(E,Npoints);
        for i_p = 1:size(P,2)
            eN = ellipsoidNorm(E,P(:,i_p)-E.center);
            assert(~(eN > 1) || withinTol(eN,1,tol))
        end
    end
end

%% Points outside
for dim = dims
    for i = 1:Ntests
        E = ellipsoid.generateRandom('Dimension',dim,'IsDegenerate',false);
        
        max_d = norm(E-E.q); % Maximal Euclidean norm of any point in E; we
        % thus just have to dispalce any random point by 10*max_d, to make
        % sure that it won't be inside E.
        
        P = randPoint(E,Npoints) + 10*max_d;
        
        for i_p = 1:size(P,2)
            eN = ellipsoidNorm(E, P(:,i_p)-E.q);
            assert(~(eN < 1) || withinTol(eN,1,tol))
        end
    end
end

%% Check norm-properties
for dim = dims
    for i=1:Ntests
        E = ellipsoid.generateRandom('Dimension',dim,'IsDegenerate',false);
        
        p1 = randPoint(E);
        p2 = randPoint(E);
        
        % Check triangle inequality
        if ~(ellipsoidNorm(E,p1+p2) <= ellipsoidNorm(E,p1)+ellipsoidNorm(E,p2)+tol)
            return;
        end
        
        % Check symmetry
        assert(ellipsoidNorm(E,p1) <= ellipsoidNorm(E,-p1) + tol && ellipsoidNorm(E,p1) >= ellipsoidNorm(E,-p1) - tol)
        
        % Check part of the positive definiteness
        assert(ellipsoidNorm(E, zeros(dim,1)) <= tol)
    end
end

% all ok
res = true;

% ------------------------------ END OF CODE ------------------------------
