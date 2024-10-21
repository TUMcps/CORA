function res = testLong_ellipsoid_supportFunc
% testLong_ellipsoid_supportFunc - unit test function of supportFunc
%
% Syntax:
%    res = testLong_ellipsoid_supportFunc
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

% Authors:       Victor Gassmann
% Written:       18-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nRuns = 2;
bools = [false,true];
for i=10:5:15
    for j=1:nRuns
        for k=1:2
            %%% generate all variables necessary to replicate results
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            % generate a few random directions
            V = randn(dim(E),2*dim(E));
            %%%
            N = 2*dim(E);
            % normalize
            V = V./sqrt(sum(V.^2,1));
            for m=1:N
                try
                    v = V(:,i);
                    [val_u,x_u] = supportFunc(E,v);
                    % check if "v" points in direction of x
                    assertLoop(v'*(x_u-E.q) >= 0,i,j,k,m)

                    % check if x is a boundary point of E
                    % check not supported for degenerate ellipsoids
                    if isFullDim(E)
                        d_rel = ellipsoidNorm(E,x_u-E.center);
                        assertLoop(d_rel <= 1+E.TOL,i,j,k,m)
                    end
                    % check if x attains val
                    assertLoop(all(withinTol(val_u,v'*x_u,E.TOL)),i,j,k,m)

                    [val_l,x_l] = supportFunc(E,v,'lower');
                    x_l_true = x_u - 2*(x_u-E.q);
                    assertLoop(all(withinTol(x_l,x_l_true,E.TOL)),i,j,k,m)

                    % check if x_l attains val_l
                    assertLoop(all(withinTol(val_l,v'*x_l,E.TOL)),i,j,k,m)

                catch ME
                    if strcmp(ME.identifier,'CORA:solverIssue')
                        disp('Randomly generated ellipsoids caused solver issues! Ignoring...');
                        continue;
                    end
                    rethrow(ME);
                end
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
