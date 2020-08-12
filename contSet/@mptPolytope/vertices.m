function V = vertices(obj)
% vertices - computes the vertices of a mptPolytope
%
% Syntax:  
%    V = vertices(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    V - vertices object
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      01-February-2011
% Last update:  12-June-2015
%               12-August-2016
%               09-May-2018 (NK, changed algorithm for vertex computation)
%               29-June-2018 (MA, all polytopes considered)
% Last revision:---

%------------- BEGIN CODE --------------

obj.P.minHRep();

% for each polytope
v = []; %init set of vertices
for i=1:length(obj.P)
    v_new = lcon2vert(obj.P(i).A,obj.P(i).b,[],[],[],0);
    v = [v; v_new];
end
V = v';

%------------- END OF CODE --------------