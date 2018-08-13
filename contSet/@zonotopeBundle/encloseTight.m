function [Zbundle] = encloseTight(Zbundle1,Zbundle2,W)
% encloseTight - Generates a zonotope bundle that encloses two zonotopes bundles
% in a tighter way than 'enclose'
%
% Syntax:  
%    [Zbundle1] = enclose(Zbundle1,Zbundle2)
%
% Inputs:
%    Zbundle1 - first zonotope bundle
%    Zbundle2 - second zonotope bundle
%
% Outputs:
%    Zbundle1 - zonotope bundle, that encloses both bundles
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      07-December-2010 
% Last update:  25-July-2016 (intervalhull replaced by interval)
%               30-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

%compute vertices
V1 = vertices(polytope(Zbundle1));
V2 = vertices(polytope(Zbundle2));

%unify vertices
V = vertices([get(V1,'V'), get(V2,'V')]);

for i=1:length(W)
    Z{i} = W{i}*zonotope(interval(pinv(W{i})*V));
end

Z{end+1} = zonotope(V);

%instantiate zonotope bundle
Zbundle=zonotopeBundle(Z);



%------------- END OF CODE --------------