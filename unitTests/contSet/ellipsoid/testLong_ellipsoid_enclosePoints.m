function res = testLong_ellipsoid_enclosePoints
% testLong_ellipsoid_enclosePoints - unit test function of
%    enclosePoints
%
% Syntax:
%    res = testLong_ellipsoid_enclosePoints
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

runs = 5;
res = true;
bools = [false,true];
for i=3:5:10
    for j=1:runs
        N = 10*i;
        for k=1:2
            try
                %%% generate all variables necessary to replicate results
                % generate point cloud
                V = randn(i,N);
                r = ceil(i*rand);
                %%%
                % generate non-degenerate and degenerate point clouds
                for m=1:2
                    [Q,~] = qr(V);
                    V = Q'*V;
                    % for bools(m)=false, generate degenerate point cloud
                    V(r,:) = bools(m)*V(r,:);
                    V = Q*V;
                    % test both methods
                    E_cov = ellipsoid.enclosePoints(V,'cov');
                    E_mv = ellipsoid.enclosePoints(V,'min-vol');
                    if ~all(contains(E_cov,V)) || ~all(contains(E_mv,V))
                        res = false;
                        return;
                    end
                end
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

% ------------------------------ END OF CODE ------------------------------
