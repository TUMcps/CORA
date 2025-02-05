function [res,cert,scaling] = priv_containsPoint(E,S,tol)
% priv_containsPoint - checks whether an ellipsoid contains a point cloud
%
% Syntax:
%    [res,cert,scaling] = priv_containsPoint(E,S)
%
% Inputs:
%    E - ellipsoid object
%    S - point cloud
%    tol - tolerance
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, S is
%           guaranteed to not be contained in E, whereas if res=false and
%           cert=false, nothing can be deduced (E could still be
%           contained in E).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(E - E.center) + E.center contains S.
%
% References:
%    [1] Boyd et al. Convex Optimization (B.2, ex. B.1)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/contains_

% Authors:       Adrian Kulmburg
% Written:       16-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

N = size(S,2);
    
c = center(E);

res = false(1,N);
cert = true(1,N);
scaling = zeros(1,N);

for i=1:N
    scaling(i) = E.ellipsoidNorm(S(:,i)-c);
    if scaling(i) <= 1+tol
        res(i) = true;
    elseif isnan(scaling(i))
        % This can only happen when Q=0 and p=0, in which case we need
        % to manually set scaling and res
        res(i) = 1;
        scaling(i) = 0;
    end
end

% ------------------------------ END OF CODE ------------------------------
