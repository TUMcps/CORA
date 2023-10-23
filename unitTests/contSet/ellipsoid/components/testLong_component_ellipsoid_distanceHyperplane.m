function res = testLong_component_ellipsoid_distanceHyperplane
% testLong_component_ellipsoid_distanceHyperplane - unit test
%    function of distanceHyperplane
%
% Syntax:
%    res = testLong_component_ellipsoid_distanceHyperplane
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
            hyp1 = conHyperplane(v',v'*s);
            if distance(E,hyp1)>E.TOL
                res = false;
                return;
            end
            val = supportFunc(E,v);
            % guaranteed to touch
            hyp2 = conHyperplane(v',val);
            if ~withinTol(distance(E,hyp2),0,E.TOL)
                res = false;
                return;
            end
            % guaranteed to not touch or intersect
            hyp3 = conHyperplane(v',val+0.1);
            if distance(E,hyp3)<=E.TOL
                res = false;
                return;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
