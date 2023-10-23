function Z = zonotope(C,varargin)
% zonotope - Over-approximates a capsule by a zonotope
%
% Syntax:
%    Z = zonotope(C)
%    Z = zonotope(C,order)
%
% Inputs:
%    C - capsule object
%    order - zonotope order of the resulting zonotope
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    C = capsule([0;0],[-2;2],2);
%    Z1 = zonotope(C,2);
%    Z2 = zonotope(C,5);
%
%    figure; hold on
%    plot(C,[1,2],'k');
%    plot(Z1,[1,2],'g');
%    plot(Z2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       20-November-2019 
% Last update:   21-July-2023 (MW, exact conversion for radius = 0)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default input arguments
order = setDefaultValues({5},varargin);

% compute zonotope enclosing the hypersphere
n = dim(C);

% exact conversion of center and generator
Z = zonotope(C.c,C.g);

if ~withinTol(C.r,0)

    % enclose capsule by ellipsoid and resulting ellipsoid by capsule
    m = order*n-1;
    Z = Z + zonotope(ellipsoid(C.r^2*eye(n)),m,'outer:norm');
    
end

% ------------------------------ END OF CODE ------------------------------
