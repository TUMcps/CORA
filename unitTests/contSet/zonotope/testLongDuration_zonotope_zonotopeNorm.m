function res = testLongDuration_zonotope_zonotopeNorm
% testLongDuration_zonotope_zonotopeNorm - unit test function of zonotopeNorm
%
% Syntax:  
%    res = testLongDuration_zonotope_zonotopeNorm
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
tol = 1e-10; % Tolerance to make sure that all comparisons run smoothly

dims = [2 5 10]; % Dimensions to be tested
Ntests = 10; % Number of tests in each case
Npoints = 10; % Number of points to be tested

%% Points inside
for dim = dims
    for i = 1:Ntests
        Z = zonotope.generateRandom(dim);
        P = randPoint(Z,Npoints);
        for i_p = 1:size(P,2)
            if zonotopeNorm(Z,P(:,i_p)-Z.center) > 1+tol
                error("Test for point-cloud inside zonotope failed!")
            end
        end
    end
end

%% Points outside
for dim = dims
    for i = 1:Ntests
        Z = zonotope.generateRandom(dim);
        % Lift Z to a space with one more dimension, that way we're sure
        % that the points cannot be contained if we choose them right
        Z_lift = zonotope([Z.center;0], [Z.generators;zeros(1,size(Z.generators,2))]);
        
        P = randPoint(Z,Npoints);
        % Same lift for P, but we displace all points by 10, so that
        % containment has now way of being able to work
        P_lift = [P;zeros(1,size(P,2))] + [zeros(dim,1);10];
        
        for i_p = 1:size(P_lift,2)
            if zonotopeNorm(Z_lift, P_lift(:,i_p)-Z_lift.center) <= 1-tol
                error("Test for point-cloud outside zonotope failed!")
            end
        end
    end
end

%% Check norm-properties
for dim = dims
    for i=1:Ntests
        Z = zonotope.generateRandom(dim);
        
        p1 = randPoint(Z);
        p2 = randPoint(Z);
        
        % Check triangle inequality
        if ~(zonotopeNorm(Z,p1+p2) <= zonotopeNorm(Z,p1)+zonotopeNorm(Z,p2)+tol)
            error("Test for triangle inequality failed!")
        end
        
        % Check symmetry
        if ~(zonotopeNorm(Z,p1) <= zonotopeNorm(Z,-p1) + tol && zonotopeNorm(Z,p1) >= zonotopeNorm(Z,-p1) - tol)
            error("Test for symmetry failed!")
        end
        
        % Check part of the positive definiteness
        if ~(zonotopeNorm(Z, zeros(dim,1)) <= tol)
            error("Test for positive definitness failed!")
        end
    end
end

res = true;

%------------- END OF CODE --------------