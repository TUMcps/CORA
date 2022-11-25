function res = testLongDuration_ellipsoid_distanceDouble
% testLongDuration_ellipsoid_distanceDouble - unit test function of testLongDuration_ellipsoid_distanceDouble
%
% Syntax:  
%    res = testLongDuration_ellipsoid_distanceDouble
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
res = true;
nRuns = 2;
bools = [false,true];
% smaller dimensions since vertices and halfspaces are involved
for i=5:5:10
    for j=1:nRuns
        for k=1:2 
            E = ellipsoid.generateRandom(i,bools(k));
            N = 2*i;
            % generate random point cloud center around origin
            V = randn(i,N);
            for m=1:2
                [Q,~] = qr(V);
                V = Q'*V;
                r = ceil(i*rand);
                % for bools(m)=false, generate degenerate point cloud
                V(r,:) = bools(m)*V(r,:);
                V = Q*V;
                V_c = mat2cell(V,E.dim,ones(1,N));
                D = distance(E,V_c);
                for i_d = 1:size(D,2)
                    % check if masks match
                    if in(E,V(:,i_d)) ~= withinTol(D(i_d),0,E.TOL)
                        res = false;
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
            break;
        end
    end
    if ~res
        break;
    end
end
%------------- END OF CODE --------------