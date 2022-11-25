function res = testLongDuration_ellipsoid_ellipsoidNorm
% testLongDuration_ellipsoid_ellipsoidNorm - unit test function of
%    ellipsoidNorm
%
% Syntax:  
%    res = testLongDuration_ellipsoid_ellipsoidNorm
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Adrian Kulmburg
% Written:      06-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE ---------------
tol = 100*eps; % Tolerance to make sure that all comparisons run smoothly

dims = [2 5 10]; % Dimensions to be tested
Ntests = 10; % Number of tests in each case
Npoints = 10; % Number of points to be tested

%% Points inside
for dim = dims
    for i = 1:Ntests
        E = ellipsoid.generateRandom(dim,false);
        P = randPoint(E,Npoints);
        for i_p = 1:size(P,2)
            if ellipsoidNorm(E,P(:,i_p)-E.center) > 1+tol
                res = false; return;
            end
        end
    end
end

%% Points outside
for dim = dims
    for i = 1:Ntests
        E = ellipsoid.generateRandom(dim,false);
        
        max_d = norm(E-E.q); % Maximal Euclidean norm of any point in E; we
        % thus just have to dispalce any random point by 10*max_d, to make
        % sure that it won't be inside E.
        
        P = randPoint(E,Npoints) + 10*max_d;
        
        for i_p = 1:size(P,2)
            if ellipsoidNorm(E, P(:,i_p)-E.q) <= 1-tol
                res = false; return;
            end
        end
    end
end

%% Check norm-properties
for dim = dims
    for i=1:Ntests
        E = ellipsoid.generateRandom(dim,false);
        
        p1 = randPoint(E);
        p2 = randPoint(E);
        
        % Check triangle inequality
        if ~(ellipsoidNorm(E,p1+p2) <= ellipsoidNorm(E,p1)+ellipsoidNorm(E,p2)+tol)
            res = false; return;
        end
        
        % Check symmetry
        if ~(ellipsoidNorm(E,p1) <= ellipsoidNorm(E,-p1) + tol && ellipsoidNorm(E,p1) >= ellipsoidNorm(E,-p1) - tol)
            res = false; return;
        end
        
        % Check part of the positive definiteness
        if ~(ellipsoidNorm(E, zeros(dim,1)) <= tol)
            res = false; return;
        end
    end
end

res = true;

%------------- END OF CODE --------------