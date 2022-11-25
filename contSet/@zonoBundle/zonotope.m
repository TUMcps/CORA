function Z = zonotope(obj)
% zonotope - encloses a zonotope bundle with a zonotope
%
% Syntax:  
%    Z = zonotope(obj)
%
% Inputs:
%    obj - zonoBundle object
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    Z{1} = zonotope([1 3 0; 1 0 2]);
%    Z{2} = zonotope([0 2 2; 0 2 -2]);
%    zB = zonoBundle(Z);
%
%    zono = zonotope(zB);
%
%    figure
%    hold on
%    plot(zB,[1,2],'r');
%    plot(zono,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Author:       Niklas Kochdumper
% Written:      01-June-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    Z = obj.Z{1};

%------------- END OF CODE --------------