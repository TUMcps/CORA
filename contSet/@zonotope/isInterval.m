function res = isInterval(Z)
% isInterval - checks if a zonotope can be equivalently represented by an
%    interval object (all generators axis-aligned)
%
% Syntax:  
%    res = isInterval(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    res - true/false
%
% Example:
%    Z1 = zonotope([1 0 1;3 2 0]);
%    Z2 = zonotope([1 2 -1;3 4 2]);
% 
%    isInterval(Z1)
%    isInterval(Z2)
% 
%    figure; hold on;
%    plot(Z1,[1,2],'r');
%    plot(Z2,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isempty

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      02-January-2020 
% Last update:  21-April-2020 (speed up)
%               09-September-2020 (MW, 1D case)
% Last revision:---

%------------- BEGIN CODE --------------

res = true;
% one-dimensional zonotopes are always intervals
if dim(Z) == 1
    return
end

G = generators(Z);
for i=1:size(G,2)
    if nnz(G(:,i)) > 1
        % two entries -> not an axis-aligned generator
        res = false;
        return
    end
end

%------------- END OF CODE --------------