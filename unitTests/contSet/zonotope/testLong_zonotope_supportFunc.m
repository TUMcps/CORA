function res = testLong_zonotope_supportFunc
% testLong_zonotope_supportFunc - unit test function of support function
%
% Syntax:
%    res = testLong_zonotope_supportFunc
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       11-October-2019
% Last update:   21-April-2023 (VG, more understandable)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

dims = 2:5;
dGen = 5;
steps = 3;
res = true;
TOL = 1e-8;
for i=dims
    n = i;
    for j=1:steps
        m = n + j*dGen;
        %randomly generate zonotope 
        Z = zonotope([randn(n,1),randn(n,m)]);
        %check if [~,x]=minnorm(Z)^2==supportFunc(Z,x)
        %check if [~,x]=norm(Z,2,'exact')^2==supportFunc(Z,x)
        [~,x_min] = minnorm(Z);
        [~,x_max] = norm(Z,2,'exact');
        c = center(Z);
        % subtract center from Zonotope, and check if all 2n random
        % directions (normalized to 1) are >= min-ball-radius and <=
        % max-ball-radius
        r_min = norm(x_min-c,2);
        r_max = norm(x_max-c,2);
        Z0 = -c+Z;
        % generate random directions
        L = randn(n,2*n);
        L = L./vecnorm(L);
        for k=1:size(L,2)
            if supportFunc(Z0,L(:,k))<r_min-TOL || ...
                    supportFunc(Z0,L(:,k))>r_max+TOL
                res = false;
                return;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
