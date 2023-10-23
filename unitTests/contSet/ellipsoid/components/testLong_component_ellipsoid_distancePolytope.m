function res = testLong_component_ellipsoid_distancePolytope
% testLong_component_ellipsoid_distancePolytope - unit test
%    function of distance
%
% Syntax:
%    res = testLong_component_ellipsoid_distancePolytope
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
for i=[5,6,7]
    for j=1:nRuns
        for k=1:2 
            %%% generate all variables necessary to replicate results
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            N = 3*i;
            % generate random point cloud center around origin
            V = randn(i,N);
            B = randPoint(E,N,'extreme');
            m = ceil(N*rand);
            fac_s = rand;
            r = randi([1,i]);
            %%%
            V = V - mean(V,2);
            s = E.q + fac_s*(B(:,m)-E.q);
            IntE = interval(E);
            IntV = interval(polytope(V));
            for m=1:2
                [U,S,V] = svd(V);
                % for bools(m)=false, generate degenerate point cloud
                S(r,r) = bools(m)*S(r,r);
                V = U*S*V';
                
                V1 = V + s;
                % intersects with E
                P1 = polytope(V1);
                d1 = distance(E,P1);
                if d1>E.TOL
                    res = false;
                    return;
                end
                
                % does not intersect E
                V2 = V+2*rad(IntE)+2*rad(IntV);
                P2 = polytope(V2);
                d2 = distance(E,P2);
                if d2<=E.TOL
                    res = false;
                    return;
                end
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
