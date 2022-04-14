function res = testLongDuration_ellipsoid_andHalfspace
% testLongDuration_ellipsoid_andHalfspace - unit test function of testLongDuration_ellipsoid_andHalfspace
%
% Syntax:  
%    res = testLongDuration_ellipsoid_andHalfspace
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
% Written:      17-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nRuns = 2;
bools = [false,true];
for i=10:5:15
    for j=1:nRuns
        for k=1:2
            E = ellipsoid.generateRandom(bools(k),i);
            % generate tangent
            c = randn(E.dim,1);
            c = c/norm(c);
            [val,x] = supportFunc(E,c);
            % completely contains E
            h1 = halfspace(c,val);
            if ~((E&h1)==E)
                res = false;
                break;
            end
            % completely outside (ignoring touching)
            h2 = halfspace(-c,-val);
            res2 = E&h2;
   
            if 1/cond(res2.Q)>E.TOL || ~all(withinTol(res2.q,x,E.TOL))
                res = false;
                break;
            end
            % choose random point within E
            N = 2*dim(E);
            B = randPoint(E,N,'extreme');
            m = ceil(N*rand);
            s = E.q + rand*(B(:,m)-E.q);
            v = randn(dim(E),1);
            h3 = halfspace(v,v'*s);
            if isempty(E&h3)
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