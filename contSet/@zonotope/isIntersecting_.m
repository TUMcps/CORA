function res = isIntersecting_(Z,S,varargin)
% isIntersecting_ - determines if zonotope intersects a set
%
% Syntax:
%    res = isIntersecting_(Z,S)
%    res = isIntersecting_(Z,S,type)
%
% Inputs:
%    Z - zonotope object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z1 = zonotope([0 1 1 0;0 1 0 1]);
%    Z2 = zonotope([2 -1 1 0;2 1 0 1]);
%    Z3 = zonotope([3.5 -1 1 0;3 1 0 1]);
% 
%    isIntersecting(Z1,Z2)
%    isIntersecting(Z1,Z3)
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z2,[1,2],'g');
% 
%    figure; hold on;
%    plot(Z1,[1,2],'b');
%    plot(Z3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, conZonotope/isIntersecting_

% Authors:       Niklas Kochdumper
% Written:       21-November-2019
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

% call function for other set representations
if isa(S,'halfspace') || isa(S,'conHyperplane') || ...
   isa(S,'polytope') || isa(S,'ellipsoid')

    res = isIntersecting_(S,Z,varargin{:});

else
    
    res = isIntersecting_(conZonotope(Z),S,varargin{:});
end

% ------------------------------ END OF CODE ------------------------------
