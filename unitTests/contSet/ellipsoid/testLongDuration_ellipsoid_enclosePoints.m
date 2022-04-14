function res = testLongDuration_ellipsoid_enclosePoints
% testLongDuration_ellipsoid_enclosePoints - unit test function of
%    enclosePoints
%
% Syntax:  
%    res = testLongDuration_ellipsoid_enclosePoints
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

% Author:       Victor Gassmann
% Written:      18-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

runs = 5;
res = true;
bools = [false,true];
for i=5:5:10
    for j=1:runs
        N = 10*i;
        for k=1:2
            try
                % generate point cloud
                V = randn(i,N);
                % generate non-degenerate and degenerate point clouds
                for m=1:2
                    [Q,~] = qr(V);
                    V = Q'*V;
                    r = ceil(i*rand);
                    % for bools(m)=false, generate degenerate point cloud
                    V(r,:) = bools(m)*V(r,:);
                    V = Q*V;
                    % test both methods
                    E_cov = ellipsoid.enclosePoints(V,'cov');
                    E_mv = ellipsoid.enclosePoints(V,'min-vol');
                    if ~in(E_cov,V) || ~in(E_mv,V)
                        res = false;
                        break;
                    end
                end
            catch ME
                if strcmp(ME.identifier,'CORA:solverIssue')
                    disp('Randomly generated ellipsoids caused solver issues! Ignoring...');
                    continue;
                end
                rethrow(ME);
            end
            if ~res
                break;
            end
        end
        if ~res
            break;
        end
    end
    if ~res
        break;
    end
end
if ~res
    disp('testLongDuration_ellipsoid_enclosePoints failed');
else
    disp('testLongDuration_ellipsoid_enclosePoints successful');
end

%------------- END OF CODE --------------