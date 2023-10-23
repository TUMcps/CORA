function c = center(zB)
% center - returns an estimate for the center of a zonotope bundle
%
% Syntax:
%    c = center(zB)
%
% Inputs:
%    zB - zonoBundle object
%
% Outputs:
%    c - center
%
% Example:
%    Z1 = zonotope(zeros(2,1),[1 0.5; -0.2 1]);
%    Z2 = zonotope(ones(2,1),[1 -0.5; 0.2 1]);
%    zB = zonoBundle({Z1,Z2});
%    c = center(zB)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       03-February-2011
% Last update:   24-April-2023 (MW, return empty if empty, otherwise cheby)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if zB.parallelSets == 0
    % fully-empty zonotope bundle
    c = [];
else
    cZ = conZonotope(zB);
    if representsa_(cZ,'emptySet',eps)
        c = double.empty(dim(cZ),0);
    else
        c = center(cZ); 
    end
end

end

% ------------------------------ END OF CODE ------------------------------
