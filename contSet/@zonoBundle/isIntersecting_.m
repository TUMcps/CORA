function res = isIntersecting_(zB,S,varargin)
% isIntersecting_ - determines if zonotope bundle intersects a set
%
% Syntax:
%    res = isIntersecting_(zB,S)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([2;2],[4;4]);
%    I2 = interval([3.5;3],[5;5]);
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%
%    isIntersecting(zB,I1)
%    isIntersecting(zB,I2)
%
%    figure; hold on
%    plot(zB,[1,2],'b');
%    plot(I1,[1,2],'g');
%    
%    figure; hold on
%    plot(zB,[1,2],'b');
%    plot(I2,[1,2],'r');
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

    res = isIntersecting_(S,zB,varargin{:});

else
    
    res = isIntersecting_(conZonotope(zB),S,varargin{:});
    
end     

% ------------------------------ END OF CODE ------------------------------
