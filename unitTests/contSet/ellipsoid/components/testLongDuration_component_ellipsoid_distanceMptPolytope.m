function res = testLongDuration_ellipsoid_distanceMptPolytope
% testLongDuration_ellipsoid_distanceMptPolytope - unit test function of testLongDuration_ellipsoid_distanceMptPolytope
%
% Syntax:  
%    res = testLongDuration_ellipsoid_distanceMptPolytope
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
            N = 3*i;
            % generate random point cloud center around origin
            V = randn(i,N);
            V = V - mean(V,2);
            B = randPoint(E,N,'extreme');
            m = ceil(N*rand);
            s = E.q + rand*(B(:,m)-E.q);
            r = randi([1,i]);
            IntE = interval(E);
            IntV = interval(mptPolytope(V'));
            for m=1:2
                [U,S,V] = svd(V);
                % for bools(m)=false, generate degenerate point cloud
                S(r,r) = bools(m)*S(r,r);
                V = U*S*V';
                
                V1 = V + s;
                % intersects with E
                P1 = mptPolytope(V1');
                d1 = distance(E,P1);
                if ~withinTol(d1,0,E.TOL)
                    res = false;
                    break;
                end
                
                % does not intersect E
                V2 = V+2*rad(IntE)+2*rad(IntV);
                P2 = mptPolytope(V2');
                d2 = distance(E,P2);
                if withinTol(d2,0,E.TOL)
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
%------------- END OF CODE --------------