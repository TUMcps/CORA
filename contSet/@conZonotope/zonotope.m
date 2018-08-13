function res = zonotope(obj)
% zonotope - over-approximates a constrained zonotope with a zonotope
%
% Syntax:  
%    res = zonotope(obj)
%
% Inputs:
%    obj - c-zonotope object
%
% Outputs:
%    res - zonotope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%    zono = zonotope(cZono);
%
%    hold on
%    plotFilled(cZono,[1,2],'r','EdgeColor','none')
%    plot(zono,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      13-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
% remove all constraints of the constrained zonotope
ng = floor((size(obj.Z,2)-1)/size(obj.Z,1)) + 1;

obj = reduce(obj,'girard',ng,0);

% construct the resulting zonotope object
res = zonotope(obj.Z);


%------------- END OF CODE --------------