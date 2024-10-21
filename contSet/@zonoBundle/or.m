function S_out = or(zB,S,varargin)
% or - Computes an over-approximation for the union of zonoBundle objects
%
% Syntax:
%    zB = zB | S
%    zB = or(zB,S)
%    zB = or(zB,S,alg,order)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object, numeric
%    alg - algorithm used to compute the union ('linprog' or 'tedrake')
%    order - zonotope order of the enclosing zonotope
%
% Outputs:
%    S_out - resulting zonoBundle object enclosing the union
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
% See also: conZonotope/or

% Authors:       Niklas Kochdumper
% Written:       26-November-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% convert to constrained zonotope, call conZonotope method
S_out = or(conZonotope(zB),S,varargin{:});
% convert back to zonotope bundle
S_out = zonoBundle(S_out);
    
% ------------------------------ END OF CODE ------------------------------
