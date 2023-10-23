function res = test_ellipsoid_ellipsoidNorm
% test_ellipsoid_ellipsoidNorm - unit test
%
% Syntax:
%    res = test_ellipsoid_ellipsoidNorm
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

% Authors:       Adrian Kulmburg, Victor Gassmann
% Written:       26-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
load cases.mat E_c
tol = 100*eps; % Tolerance to make sure that all comparisons run smoothly

for i=1:length(E_c)
    E = E_c{i}.E1;
    n = length(E.q);
    dim = n;
    %% Points inside
    P = randPoint(E,2*n);
    for i_p = 1:size(P,2)
        eN = ellipsoidNorm(E,P(:,i_p)-E.center);
        if eN > 1 && ~withinTol(eN,1,tol)
            res = false; break;
        end
    end

    %% Points outside

    max_d = norm(E-E.q); % Maximal Euclidean norm of any point in E; we
    % thus just have to displace any random point by 10*max_d, to make
    % sure that it won't be inside E.

    P = randPoint(E,2*n) + 10*max_d;

    for i_p = 1:size(P,2)
        eN = ellipsoidNorm(E, P(:,i_p)-E.q);
        if eN < 1 && ~withinTol(eN,1,tol)
            res = false; break;
        end
    end

    %% Check norm-properties

    p1 = randPoint(E);
    p2 = randPoint(E);

    % Check triangle inequality
    if ~(ellipsoidNorm(E,p1+p2) <= ellipsoidNorm(E,p1)+ellipsoidNorm(E,p2)+tol)
        res = false; break;
    end

    % Check symmetry
    if ~(ellipsoidNorm(E,p1) <= ellipsoidNorm(E,-p1) + tol && ellipsoidNorm(E,p1) >= ellipsoidNorm(E,-p1) - tol)
        res = false; break;
    end

    % Check part of the positive definiteness
    if ~(ellipsoidNorm(E, zeros(dim,1)) <= tol)
        res = false; break;
    end

end

% ------------------------------ END OF CODE ------------------------------
