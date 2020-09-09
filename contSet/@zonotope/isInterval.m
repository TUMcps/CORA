function res = isInterval(Z)
% isInterval - check if a zonotope represents an interval
%
% Syntax:  
%    res = isInterval(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    res - 1 if zonotope represents an interval, 0 else
%
% Example:
%    zono1 = zonotope([1 0 1;3 2 0]);
%    zono2 = zonotope([1 2 -1;3 4 2]);
%
%    isInterval(zono1)
%    isInterval(zono2)
%
%    figure
%    plot(zono1,[1,2],'r');
%
%    figure
%    plot(zono2,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: isempty

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      02-January-2020 
% Last update:  21-April-2020 (MW, speed up)
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
    
end

%------------- END OF CODE --------------