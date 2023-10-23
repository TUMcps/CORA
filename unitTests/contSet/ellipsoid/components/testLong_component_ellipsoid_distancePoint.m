function res = testLong_component_ellipsoid_distancePoint
% testLong_component_ellipsoid_distancePoint - unit test function
%    distancePoint
%
% Syntax:
%    res = testLong_component_ellipsoid_distancePoint
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
% smaller dimensions since vertices and halfspaces are involved
for i=5:5:10
    for j=1:nRuns
        for k=1:2 
            %%% generate all variables necessary to replicate results
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            N = 2*i;
            % generate random point cloud center around origin
            V = randn(i,N);
            r = ceil(i*rand);
            %%%
            for m=1:2
                [Q,~] = qr(V);
                V = Q'*V;
                % for bools(m)=false, generate degenerate point cloud
                V(r,:) = bools(m)*V(r,:);
                V = Q*V;
                D = distance(E,V);
                for i_d = 1:size(D,2)
                    % check if masks match
                    if contains(E,V(:,i_d)) ~= (D(i_d)<=E.TOL)
                        res = false;
                        return;
                    end
                end
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
