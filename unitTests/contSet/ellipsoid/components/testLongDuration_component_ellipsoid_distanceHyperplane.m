function res = testLongDuration_ellipsoid_distanceHyperplane
% testLongDuration_ellipsoid_distanceHyperplane - unit test function of testLongDuration_ellipsoid_distanceHyperplane
%
% Syntax:  
%    res = testLongDuration_ellipsoid_distanceHyperplane
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
for i=5:5:20
    for j=1:nRuns
        for k=1:2 
            E = ellipsoid.generateRandom(i,bools(k));
            N = 2*i;
            B = randPoint(E,N,'extreme');
            m = ceil(N*rand);
            s = E.q + rand*(B(:,m)-E.q);
            v = randn(i,1);
            v = v/norm(v);
            % guaranteed to intersect
            hyp1 = conHyperplane(v',v'*s);
            if distance(E,hyp1)>E.TOL
                res = false;
                break;
            end
            val = supportFunc(E,v);
            % guaranteed to touch
            hyp2 = conHyperplane(v',val);
            if ~withinTol(distance(E,hyp2),0,E.TOL)
                res = false;
                break;
            end
            % guaranteed to not touch or intersect
            hyp3 = conHyperplane(v',val+0.1);
            if distance(E,hyp3)<=E.TOL
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
%------------- END OF CODE --------------