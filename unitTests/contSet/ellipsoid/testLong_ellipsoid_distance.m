function res = testLong_ellipsoid_distance
% testLong_ellipsoid_distance - unit test function of distance
%
% Syntax:
%    res = testLong_ellipsoid_distance
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
% See also: none

% Authors:       Victor Gassmann
% Written:       18-March-2021
% Last update:   ---
% Last revision: 22-September-2024 (MW, integrate components)

% ------------------------------ BEGIN CODE -------------------------------

% ellipsoid
assert(aux_distanceEllipsoid());

% hyperplane
assert(aux_distanceHyperplane());

% polytope (halfspace is implicitly tested)
assert(aux_distancePolytope());

% double
assert(aux_distancePoint());

% combine results
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_distanceEllipsoid()

res = true;
nRuns = 2;
bools = [false,true];
for i=5:5:10
    for j=1:nRuns
        for k=1:2 
            try
                %%% generate all variables necessary to replicate results
                E1 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',false);
                E2 = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
                N = 2*i;
                m = ceil(N*rand);
                B = randPoint(E1,N,'extreme');
                fac_s = rand;
                r = ceil(i*rand);
                %%%
                b = B(:,m);
                s = E1.q + fac_s*(b-E1.q);
                % E1 and E2 are intersecting
                E2 = ellipsoid(E2.Q,s);
                assertLoop(distance(E1,E2)<=E1.TOL,i,j,k)

                IntE1 = interval(E1);
                q3 = E1.q+2.1*rad(IntE1);
                [U,S,~] = svd(E1.Q);
                S(r,r) = bools(k)*S(r,r);
                % guaranteed to not intersect with E1
                E3 = ellipsoid(U*S*U',q3);
                assertLoop(distance(E1,E3) > E1.TOL,i,j,k)

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

function res = aux_distanceHyperplane()

res = true;
nRuns = 2;
bools = [false,true];
for i=5:5:20
    for j=1:nRuns
        for k=1:2 
            %%% generate all variables necessary to replicate results
            E = ellipsoid.generateRandom('Dimension',i,'IsDegenerate',bools(k));
            N = 2*i;
            B = randPoint(E,N,'extreme');
            m = ceil(N*rand);
            fac_s = rand;
            v = randn(i,1);
            %%%
            s = E.q + fac_s*(B(:,m)-E.q);
            v = v/norm(v);
            % guaranteed to intersect
            hyp1 = polytope([],[],v',v'*s);
            assertLoop(distance(E,hyp1)<=E.TOL,i,j,k)

            val = supportFunc(E,v);
            % guaranteed to touch
            hyp2 = polytope([],[],v',val);
            assertLoop(withinTol(distance(E,hyp2),0,E.TOL),i,j,k)

            % guaranteed to not touch or intersect
            hyp3 = polytope([],[],v',val+0.1);
            assertLoop(distance(E,hyp3) > E.TOL,i,j,k)
        end
    end
end

end

function res = aux_distancePolytope()

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
                assertLoop(d1<=E.TOL,i,j,k,m)
                
                % does not intersect E
                V2 = V+2*rad(IntE)+2*rad(IntV);
                P2 = polytope(V2);
                d2 = distance(E,P2);
                assertLoop(d2 > E.TOL,i,j,k,m)
            end
        end
    end
end

end

function res = aux_distancePoint()

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
                    assertLoop(contains(E,V(:,i_d)) == (D(i_d)<=E.TOL),i,j,k,m,i_d)
                end
            end
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
