function zB = or(zB1, varargin)
% or - Computes an over-approximation for the union of zonoBundle objects
%
% Syntax:
%    zB = or(zB1, zB2)
%    zB = or(zB1, ... , zBm)
%    zB = or(zB1, ... , zBm, alg)
%    zB = or(zB1, ... , zBm, alg, order)
%
% Inputs:
%    zB1,...,zBm - zonoBundle objects
%    alg - algorithm used to compute the union ('linProg' or 'tedrake')
%    order - zonotope order of the enclosing zonotope
%
% Outputs:
%    zB - resulting zonoBundle object enclosing the union
%
% Example: 
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%    I = interval([4;4],[5;6]);
% 
%    res = zB | I;
% 
%    figure; hold on;
%    plot(zB,[1,2],'FaceColor','r');
%    plot(I,[1,2],'FaceColor','b');
%    plot(res,[1,2],'k');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/or, zonotope/or

% Authors:       Niklas Kochdumper
% Written:       26-November-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to conZonotope and call conZonotope method
zB1 = conZonotope(zB1);
cZ = or(zB1, varargin{:});
zB = zonoBundle(cZ);
    
% ------------------------------ END OF CODE ------------------------------
