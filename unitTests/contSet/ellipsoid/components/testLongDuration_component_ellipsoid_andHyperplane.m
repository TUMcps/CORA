function res = testLongDuration_ellipsoid_andHyperplane
% testLongDuration_ellipsoid_andHyperplane - unit test function of testLongDuration_ellipsoid_andHyperplane
%
% Syntax:  
%    res = testLongDuration_ellipsoid_andHyperplane
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
            % construct hyperplane that intersects and check
            % select random boundary point
            N = 2*E.dim;
            m = ceil(N*rand);
            B = randPoint(E,N,'extreme');
            b = B(:,m);
            % select sample somewhere between b and E.q
            s = E.q + rand*(b-E.q);
            % generate random direction
            d = randn(E.dim,1);
            d = d/norm(d);
            % generate hyperplane
            H = conHyperplane(d',d'*s);
            if isempty(E&H)
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