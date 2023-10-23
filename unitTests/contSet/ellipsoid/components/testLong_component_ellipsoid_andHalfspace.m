function res = testLong_component_ellipsoid_andHalfspace
% testLong_component_ellipsoid_andHalfspace - unit test function of
%    andHalfspace
%
% Syntax:
%    res = testLong_component_ellipsoid_andHalfspace
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
% Written:       17-March-2021
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
            c = randn(E.dim,1);
            B = randPoint(E,2*dim(E),'extreme');
            m = ceil(2*dim(E)*rand);
            fac_s = rand;
            v = randn(dim(E),1);
            %%%
            % normalize tangent
            c = c/norm(c);
            [val,x] = supportFunc(E,c);
            % completely contains E
            h1 = halfspace(c,val);
            if ~((E&h1)==E)
                res = false;
                return;
            end
            % completely outside (ignoring touching)
            h2 = halfspace(-c,-val);
            res2 = E&h2;
   
            if 1/cond(res2.Q)>E.TOL || ~all(withinTol(res2.q,x,E.TOL))
                res = false;
                return;
            end
            % choose random point within E
            s = E.q + fac_s*(B(:,m)-E.q);
            h3 = halfspace(v,v'*s);
            if representsa_(E&h3,'emptySet',eps)
                res = false;
                return;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
