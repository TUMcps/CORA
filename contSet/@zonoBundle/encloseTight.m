function zB = encloseTight(zB1,zB2,W)
% encloseTight - Generates a zonotope bundle that encloses two zonotopes
%    bundles in a tighter way than standard enclose operation
%
% Syntax:
%    zB = encloseTight(zB1,zB2)
%
% Inputs:
%    zB1 - zonoBundle object
%    zB2 - zonoBundle object
%    W - ???
%
% Outputs:
%    zB - zonotope bundle that encloses both bundles
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       07-December-2010 
% Last update:   25-July-2016 (intervalhull replaced by interval)
%                30-July-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%compute vertices
V1 = vertices(polytope(zB1));
V2 = vertices(polytope(zB2));

%unify vertices
V = vertices([V1, V2]);

for i=1:length(W)
    Z{i} = W{i}*zonotope(interval.enclosePoints(pinv(W{i})*V));
end

Z{end+1} = zonotope.enclosePoints(V,'stursberg');

%instantiate zonotope bundle
zB = zonoBundle(Z);

% ------------------------------ END OF CODE ------------------------------
